#include "catch2/catch.hpp"

#include "geom/LineSegment.hh"
#include "geom/Ray.hh"
#include "geom/Volume.hh"
#include "geom/World.hh"

using NuGeom::LogicalVolume;
using NuGeom::PhysicalVolume;

TEST_CASE("Single LogicalVolume", "[Volume]") {
    NuGeom::Material mat("Water", 1.0, 2);
    mat.AddElement(NuGeom::Element("Hydrogen", 1, 1), 2);
    mat.AddElement(NuGeom::Element("Oxygen", 8, 16), 1);
    auto box = std::make_shared<NuGeom::Box>();
    LogicalVolume vol(mat, box);

    SECTION("Volume is correct") {
        CHECK(vol.Volume() == Approx(1));
    }

    SECTION("Mass is correct") {
        CHECK(vol.Mass() == Approx(1));
    }
}

TEST_CASE("Intersect", "[Volume]") {
    NuGeom::Material mat("Water", 1.0, 2);
    mat.AddElement(NuGeom::Element("Hydrogen", 1, 1), 2);
    mat.AddElement(NuGeom::Element("Oxygen", 8, 16), 1);
    auto box = std::make_shared<NuGeom::Box>();
    auto vol = std::make_shared<LogicalVolume>(mat, box);
    PhysicalVolume pvol(vol, NuGeom::Transform3D{}, NuGeom::Transform3D{});

    SECTION("Outside Volume") {
        NuGeom::Ray ray({0, 0, -1}, {0, 0, 1}, 1.0);
        static constexpr double expected_time = 0.5;
        CHECK_THAT(pvol.Intersect(ray), Catch::WithinAbs(expected_time, 1e-8));
    }

    SECTION("Inside Volume") {
        // Ray origin is inside the volume.  Intersect() must return 0.0 so the BVH
        // immediately selects this daughter rather than using the (incorrect) exit time.
        static constexpr double eps = 1e-7;
        NuGeom::Ray ray({0, 0, -0.5 + eps}, {0, 0, 1}, 1.0);
        static constexpr double expected_time = 0.0;
        CHECK_THAT(pvol.Intersect(ray), Catch::WithinAbs(expected_time, 1e-8));
    }
}

TEST_CASE("Line Segments", "[Volume]") {
    NuGeom::Material mat_world("Rock", 2.0, 1);
    mat_world.AddElement(NuGeom::Element("Silicon", 14, 28), 1);
    NuGeom::Material mat_outer("Water", 1.0, 2);
    mat_outer.AddElement(NuGeom::Element("Hydrogen", 1, 1), 2);
    mat_outer.AddElement(NuGeom::Element("Oxygen", 8, 16), 1);
    NuGeom::Material mat_inner("Iron", 7.874, 1);
    mat_inner.AddElement(NuGeom::Element("Iron", 26, 56), 1);

    auto inner_box = std::make_shared<NuGeom::Box>();
    auto inner_vol = std::make_shared<LogicalVolume>(mat_inner, inner_box);
    NuGeom::RotationX3D rot(45 * M_PI / 180.0);
    auto inner_pvol = std::make_shared<PhysicalVolume>(inner_vol, NuGeom::Transform3D{}, rot);

    auto outer_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    auto outer_vol = std::make_shared<LogicalVolume>(mat_outer, outer_box);
    outer_vol->AddDaughter(inner_pvol);
    NuGeom::RotationX3D rot2(90 * M_PI / 180.0);
    auto outer_pvol = std::make_shared<PhysicalVolume>(outer_vol, NuGeom::Transform3D{}, rot2);

    auto world_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 4});
    auto world_lv = std::make_shared<LogicalVolume>(mat_world, world_box);
    world_lv->AddDaughter(outer_pvol);

    // Use World to build the unique PV tree with correct physical mothers.
    NuGeom::World world(world_lv);
    NuGeom::Ray ray({0, 0, -2}, {0, 0, 1}, 1.0);
    auto segments = world.GetLineSegments(ray);

    CHECK(segments.size() == 5);
    CHECK_THAT(segments[0].Length(), Catch::WithinAbs(1, 1e-8));
    CHECK_THAT(segments[1].Length(), Catch::WithinAbs(1 - sqrt(2) / 2, 1e-8));
    CHECK_THAT(segments[2].Length(), Catch::WithinAbs(sqrt(2), 1e-8));
    CHECK_THAT(segments[3].Length(), Catch::WithinAbs(1 - sqrt(2) / 2, 1e-8));
    CHECK_THAT(segments[4].Length(), Catch::WithinAbs(1, 1e-8));

    CHECK_THAT(segments[0].Start().Z(), Catch::WithinAbs(-2, 1e-8));
    CHECK_THAT(segments[1].Start().Z(), Catch::WithinAbs(segments[0].End().Z(), 1e-8));
    CHECK_THAT(segments[2].Start().Z(), Catch::WithinAbs(segments[1].End().Z(), 1e-8));
    CHECK_THAT(segments[3].Start().Z(), Catch::WithinAbs(segments[2].End().Z(), 1e-8));
    CHECK_THAT(segments[4].Start().Z(), Catch::WithinAbs(segments[3].End().Z(), 1e-8));
    CHECK_THAT(segments[4].End().Z(), Catch::WithinAbs(2, 1e-8));
    CHECK_THAT(segments[1].Start().Z(), Catch::WithinAbs(-1, 1e-8));
}
