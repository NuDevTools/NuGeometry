#include "catch2/catch.hpp"

#include "geom/DetectorSim.hh" // defines NuGeom::EnergyRay used by TestGen.hh
#include "geom/InteractionViz.hh"
#include "geom/Logging.hh"
#include "geom/Parser.hh"
#include "geom/Random.hh"
#include "geom/Ray.hh"
#include "geom/TestGen.hh"

#include <string>

namespace {
// CreateLogger throws if the logger already exists; other test cases (and files)
// in the same testsuite may have created it already, so make it idempotent.
void ensure_logger() {
    static bool done = false;
    if(done) return;
    try {
        CreateLogger(false, 0, 1);
    } catch(...) {}
    done = true;
}

NuGeom::World load_simple_boxes() {
    std::string gdml = std::string(NUGEOM_SOURCE_DIR) + "/SimpleBoxes.geom.manual.gdml";
    NuGeom::GDMLParser parser(gdml);
    return parser.GetWorld();
}
} // namespace

TEST_CASE("RayInteractionSim records interactions and outgoing tracks", "[interaction_viz]") {
    ensure_logger();
    NuGeom::Random::Instance().Seed(12345);

    auto world = load_simple_boxes();
    NuGeom::RayInteractionSim sim(std::move(world));

    SECTION("A large cross section forces an interaction with traced secondaries") {
        sim.SetConstantCrossSection(1e-20); // huge → mean free path << detector
        sim.SetSecondaryModel(4, 0.3);

        // Fire along +Z through the centre of the detector where the scintillator
        // blocks sit, starting on the upstream face.
        NuGeom::Ray ray({0, 0, -150}, {0, 0, 1}, 1.0);
        auto event = sim.ShootRay(5.0, ray);

        REQUIRE(event.optical_depth > 0.0);
        REQUIRE(event.interacted);
        REQUIRE_FALSE(event.vertex_material.empty());
        REQUIRE(event.secondaries.size() == 4);

        // Every outgoing particle starts at the interaction vertex and is traced.
        for(const auto &track : event.secondaries) {
            REQUIRE((track.start - event.vertex).Norm() < 1e-9);
            CHECK(std::abs(track.direction.Norm() - 1.0) < 1e-9);
        }
    }

    SECTION("A vanishing cross section yields no interaction") {
        sim.SetConstantCrossSection(0.0);
        NuGeom::Ray ray({0, 0, -150}, {0, 0, 1}, 1.0);
        auto event = sim.ShootRay(5.0, ray);

        CHECK(event.optical_depth == Approx(0.0));
        CHECK_FALSE(event.interacted);
        CHECK(event.secondaries.empty());
        CHECK_FALSE(event.incoming_path.empty());
    }
}

TEST_CASE("RayInteractionSim batch summary is consistent", "[interaction_viz]") {
    ensure_logger();
    NuGeom::Random::Instance().Seed(999);

    auto world = load_simple_boxes();
    NuGeom::RayInteractionSim sim(std::move(world));
    sim.SetConstantCrossSection(1e-21);
    sim.SetSecondaryModel(2, 0.2);

    auto box = sim.GetWorld().GetWorldBox();
    NuGeom::Vector3D center{0, 0, box.min.Z()};
    BeamRayGen gen(1.0, 10.0, center, 30.0, 0.005, box.min.Z());

    auto events = sim.ShootRays(50, gen);
    auto stats = NuGeom::RayInteractionSim::Summarize(events);

    REQUIRE(stats.nrays == 50);
    CHECK(stats.ninteracted <= stats.nrays);
    CHECK(stats.interaction_fraction() >= 0.0);
    CHECK(stats.interaction_fraction() <= 1.0);

    // Secondary count in the summary matches the per-event secondaries.
    size_t manual = 0;
    for(const auto &e : events)
        if(e.interacted) manual += e.secondaries.size();
    CHECK(manual == stats.nsecondaries);
}
