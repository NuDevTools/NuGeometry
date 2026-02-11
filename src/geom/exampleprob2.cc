#include "CLI/CLI.hpp"
#include "geom/DetectorSim.hh"
#include "geom/Material.hh"
#include "geom/TestGen.hh"
#include "geom/Vector3D.hh"
#include "geom/Volume.hh"
#include "geom/World.hh"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"
#include <cstdlib>
#include <ctime>
#include <map>
#include <string>

int main(int argc, char **argv) {
    auto console = spdlog::stdout_color_mt("NuGeom");
    spdlog::set_default_logger(console);
    spdlog::set_pattern("[%n] [%^%l%$] %v");
    CLI::App app("Neutrino Geometry Driver");
    argv = app.ensure_utf8(argv);
    int verbosity = 0;
    std::string outfile = "hit_locations1.txt";
    size_t nevents = 1000;
    double pot = 1e23;
    app.add_flag("-v,--verbose", verbosity, "Increase the verbosity level");
    app.add_option("-f,--file", outfile, "File to write out the locations of the hits.");
    app.add_option("-e,--events", nevents, "Number of events to generate. (Default 1000)");
    app.add_option("-p,--pot", pot, "Number of protons on target (POT). (Default 10^23)");

    try {
        app.parse(argc, argv);
    } catch(const CLI::ParseError &e) { return app.exit(e); }

    if(verbosity == 1)
        spdlog::set_level(spdlog::level::debug);
    else if(verbosity == 2)
        spdlog::set_level(spdlog::level::trace);

    // Define materials in the detector
    NuGeom::Material mat("Water", 1.0, 2);
    NuGeom::Material mat1("Argon", 9.0, 1);
    mat.AddElement(NuGeom::Element("Hydrogen", 1, 1), 2);
    mat.AddElement(NuGeom::Element("Oxygen", 8, 16), 1);
    mat1.AddElement(NuGeom::Element("Argon", 18, 40), 1);

    // Define the interaction geometry
    // Define the inner detector
    auto inner_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{1, 1, 1}); // Define a 1x1x1 box

    auto inner_vol = std::make_shared<NuGeom::LogicalVolume>(mat, inner_box);
    NuGeom::RotationX3D rot(45 * M_PI / 180.0);
    auto inner_pvol =
        std::make_shared<NuGeom::PhysicalVolume>(inner_vol, NuGeom::Transform3D{}, rot);

    // Define the outer detector
    auto outer_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    auto outer_vol = std::make_shared<NuGeom::LogicalVolume>(mat1, outer_box);
    outer_vol->AddDaughter(inner_pvol);
    inner_vol->SetMother(outer_vol);
    NuGeom::RotationX3D rot2(30 * M_PI / 180.0);
    auto outer_pvol =
        std::make_shared<NuGeom::PhysicalVolume>(outer_vol, NuGeom::Transform3D{}, rot2);
    inner_pvol->SetMother(outer_pvol);

    // Define the "World"
    auto world_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 4});
    auto world_vol = std::make_shared<NuGeom::LogicalVolume>(mat, world_box);
    outer_vol->SetMother(world_vol);
    world_vol->AddDaughter(outer_pvol);
    NuGeom::World world(world_vol);

    NuGeom::DetectorSim sim;
    sim.Setup(world);

    std::map<size_t, double> xsec_map = {
        {1000010010, 1e-38}, {1000080160, 1e-38}, {1000180400, 1e-38}};

    double max_prob = 0;
    auto raygen = std::make_shared<TestRayGen>(0, 10, NuGeom::Vector3D{-2, -2, -2},
                                               NuGeom::Vector3D{2, 2, -2});
    TestEventGen gen1(xsec_map);

    std::vector<NuGeom::EnergyRay> rays;
    size_t ntest_rays = 1 << 22;
    for(size_t i = 0; i < ntest_rays; ++i) { rays.push_back(raygen->GetRay()); }
    for(const auto &ray : rays) {
        auto segments = sim.GetLineSegments(ray.second);
        auto materials = sim.GetMaterials(segments);
        auto xsecs = gen1.EvaluateCrossSections(ray.first, materials);
        std::vector<double> probs = sim.EvaluateProbs(segments, xsecs);
        double prob_sum = std::accumulate(probs.begin(), probs.end(), 0.0);
        if(prob_sum > max_prob) { max_prob = prob_sum; }
    }
    sim.SetMaxProb(max_prob);
    spdlog::info("Set max probability to {}", max_prob);
    std::ofstream hist;
    hist.open(outfile);
    if(!hist.is_open()) { throw std::runtime_error("Failed to open hit_locations2.txt"); }
    // for rays const
    size_t ntest_hits = 0;

    while(ntest_hits < nevents) {
        auto ray = raygen->GetRay();
        std::vector<NuGeom::LineSegment> segments0 = sim.GetLineSegments(ray.second);
        std::set<NuGeom::Material> materials = sim.GetMaterials(segments0);
        std::map<NuGeom::Element, double> material_xsec_map =
            gen1.EvaluateCrossSections(ray.first, materials);
        auto [interaction_point, interaction_material] =
            sim.Interaction(segments0, material_xsec_map);
        if(interaction_material != NuGeom::Material()) {
            ntest_hits++;
            // Valid interaction point
            hist << interaction_point.X() << " " << interaction_point.Y() << " "
                 << interaction_point.Z() << " " << interaction_material.Name() << " " << ray.first
                 << "\n";
        } else {
            // No interaction occurred
            spdlog::debug("No interaction occurred for this ray.");
        }
    }
    hist.close();

    return 0;
}
