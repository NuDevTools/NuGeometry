#include "catch2/catch.hpp"

#include "geom/Element.hh"
#include "geom/LineSegment.hh"
#include "geom/Material.hh"
#include "geom/Ray.hh"
#include "geom/Shape.hh"
#include "geom/Transform3D.hh"
#include "geom/Volume.hh"
#include "geom/World.hh"

#include <map>
#include <memory>
#include <string>

namespace {

NuGeom::Material Mat(const std::string &name, double density) {
    NuGeom::Material m(name, density, 1);
    m.AddElement(NuGeom::Element("Hydrogen", 1, 1), 1);
    return m;
}

std::shared_ptr<NuGeom::PhysicalVolume> Place(const std::string &name, const NuGeom::Material &mat,
                                              const NuGeom::Vector3D &size,
                                              const NuGeom::Vector3D &at) {
    auto lv =
        std::make_shared<NuGeom::LogicalVolume>(name, mat, std::make_shared<NuGeom::Box>(size));
    return std::make_shared<NuGeom::PhysicalVolume>(
        name, lv, NuGeom::Translation3D(at.X(), at.Y(), at.Z()), NuGeom::Rotation3D());
}

/// World: 100^3 of Rock, holding a 40^3 Heavy block and a 2^3 Speck.
std::shared_ptr<NuGeom::LogicalVolume> MakeWorldLV() {
    auto world = std::make_shared<NuGeom::LogicalVolume>(
        "World", Mat("Rock", 2.0), std::make_shared<NuGeom::Box>(NuGeom::Vector3D(100, 100, 100)));
    world->AddDaughter(Place("Heavy", Mat("Heavy", 5.0), {40, 40, 40}, {0, 0, 0}));
    world->AddDaughter(Place("Speck", Mat("Speck", 1.0), {2, 2, 2}, {30, 30, 30}));
    return world;
}

size_t CountNodes(const std::shared_ptr<NuGeom::PhysicalVolume> &pv) {
    size_t n = 1;
    for(const auto &d : pv->Daughters()) n += CountNodes(d);
    return n;
}

} // namespace

TEST_CASE("Prune drops volumes below the mass-fraction threshold", "[World][Prune]") {
    NuGeom::World world(MakeWorldLV());
    REQUIRE(world.NDaughters() == 2);

    // Rock: (1e6 - 64000 - 8) * 2 = 1871984 g; Heavy: 64000*5 = 320000; Speck: 8*1 = 8.
    // Total 2191992 g, so Speck is 3.65e-6 of the mass and Heavy is 0.146.
    NuGeom::World::PruneOptions opts;
    opts.min_mass_fraction = 1e-4;
    auto report = world.Prune(opts);

    CHECK(report.removed_subtrees == 1);
    CHECK(report.removed_nodes == 1);
    CHECK(report.entries.size() == 1);
    CHECK(report.entries.front().volume == "Speck");
    CHECK(report.entries.front().reason == "mass fraction");
    CHECK(report.removed_mass == Approx(8.0));
    // Speck's 8 cm^3 is now Rock (density 2), so the geometry GAINS 16 - 8 = 8 g.
    CHECK(report.mass_delta == Approx(8.0));
    CHECK(world.NDaughters() == 1);
}

TEST_CASE("Prune drops volumes by material name", "[World][Prune]") {
    NuGeom::World world(MakeWorldLV());
    NuGeom::World::PruneOptions opts;
    opts.drop_materials = {"Heavy"};
    auto report = world.Prune(opts);

    CHECK(report.removed_subtrees == 1);
    CHECK(report.entries.front().volume == "Heavy");
    CHECK(report.entries.front().reason == "material");
    CHECK(report.entries.front().replaced_by == "Rock");
    CHECK(report.removed_mass == Approx(320000.0));
    // 64000 cm^3 of Rock (2.0) replaces 320000 g of Heavy: delta = 128000 - 320000.
    CHECK(report.mass_delta == Approx(-192000.0));
    CHECK(world.NDaughters() == 1);
}

TEST_CASE("Prune by material takes the whole subtree", "[World][Prune]") {
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>(
        "World", Mat("Rock", 2.0), std::make_shared<NuGeom::Box>(NuGeom::Vector3D(100, 100, 100)));
    auto shell = Place("Shell", Mat("Shell", 1.0), {40, 40, 40}, {0, 0, 0});
    shell->GetLogicalVolume()->AddDaughter(
        Place("Core", Mat("Core", 8.0), {10, 10, 10}, {0, 0, 0}));
    world_lv->AddDaughter(shell);

    NuGeom::World world(world_lv);
    REQUIRE(CountNodes(world.GetRootPV()) == 3); // world + shell + core

    NuGeom::World::PruneOptions opts;
    opts.drop_materials = {"Shell"};
    auto report = world.Prune(opts);

    CHECK(report.removed_subtrees == 1);
    CHECK(report.removed_nodes == 2); // shell AND its core
    CHECK(CountNodes(world.GetRootPV()) == 1);
}

TEST_CASE("Prune keep_volumes protects a volume from both criteria", "[World][Prune]") {
    NuGeom::World world(MakeWorldLV());
    NuGeom::World::PruneOptions opts;
    opts.min_mass_fraction = 0.5; // would drop everything
    opts.drop_materials = {"Heavy"};
    opts.keep_volumes = {"Heavy", "Speck"};
    auto report = world.Prune(opts);

    CHECK(report.removed_subtrees == 0);
    CHECK(world.NDaughters() == 2);
}

TEST_CASE("Prune with no criteria is a no-op", "[World][Prune]") {
    NuGeom::World world(MakeWorldLV());
    auto report = world.Prune({});
    CHECK(report.removed_subtrees == 0);
    CHECK(report.total_mass > 0.0); // still reported, even though nothing was dropped
    CHECK(world.NDaughters() == 2);
}

TEST_CASE("Pruned volumes no longer appear in the traversal", "[World][Prune]") {
    NuGeom::World world(MakeWorldLV());
    // A ray along +z through (30,30) crosses the Speck.
    NuGeom::Ray ray(NuGeom::Vector3D(30, 30, -60), NuGeom::Vector3D(0, 0, 1), 1.0);

    bool saw_speck = false;
    for(const auto &s : world.GetLineSegments(ray))
        if(s.GetMaterial().Name() == "Speck") saw_speck = true;
    REQUIRE(saw_speck);

    NuGeom::World::PruneOptions opts;
    opts.drop_materials = {"Speck"};
    world.Prune(opts);

    for(const auto &s : world.GetLineSegments(ray)) CHECK(s.GetMaterial().Name() != "Speck");
    // and the space is Rock now
    CHECK(world.FindMaterial(NuGeom::Vector3D(30, 30, 30)).Name() == "Rock");
}

// ---------------------------------------------------------------------------
// The column-density-only traversal shares the sweep with GetLineSegments, so
// it must report exactly the same path length per material.
// ---------------------------------------------------------------------------
TEST_CASE("GetColumnLengths matches GetLineSegments per material", "[World][Column]") {
    NuGeom::World world(MakeWorldLV());

    for(const auto &dir : {NuGeom::Vector3D(0, 0, 1), NuGeom::Vector3D(0.3, 0.2, 1),
                           NuGeom::Vector3D(-0.15, 0.4, 1)}) {
        NuGeom::Ray ray(NuGeom::Vector3D(30, 30, -60), dir, 1.0);

        std::map<std::string, double> from_segments;
        for(const auto &s : world.GetLineSegments(ray))
            from_segments[s.GetMaterial().Name()] += s.Length();

        std::map<std::string, double> from_columns;
        for(const auto &[mat, len] : world.GetColumnLengths(ray)) from_columns[mat->Name()] += len;

        REQUIRE(from_columns.size() == from_segments.size());
        for(const auto &[name, len] : from_segments)
            CHECK(from_columns.at(name) == Approx(len).epsilon(1e-9));
    }
}

TEST_CASE("GetColumnLengths sees pruned volumes disappear", "[World][Column][Prune]") {
    NuGeom::World world(MakeWorldLV());
    NuGeom::Ray ray(NuGeom::Vector3D(30, 30, -60), NuGeom::Vector3D(0, 0, 1), 1.0);

    auto has = [&](const std::string &name) {
        for(const auto &[mat, len] : world.GetColumnLengths(ray))
            if(mat->Name() == name && len > 0) return true;
        return false;
    };
    REQUIRE(has("Speck"));
    NuGeom::World::PruneOptions opts;
    opts.drop_materials = {"Speck"};
    world.Prune(opts);
    CHECK(!has("Speck"));
}
