#include "CLI/CLI.hpp"
#include "geom/DetectorSim.hh"
#include "geom/Material.hh"
#include "geom/Parser.hh"
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
    std::string outfile = "hits.out";
    std::string geomfile = "";
    std::string lar_vol = "ArgonCubeDetector";
    std::string lar_pos = "volArgonCubeDetector_pos";
    std::vector<std::string> parent_pos_names; // parent position refs to accumulate to world frame
    size_t nevents = 1000;
    double pot = 1e21;
    double beam_radius = 200.0; // cm — transverse beam radius on the upstream face
    double sigma_theta = 3e-3;  // rad — RMS angular divergence (3 mrad, DUNE-like)
    app.add_flag("-v,--verbose", verbosity, "Increase the verbosity level");
    app.add_option("-f,--file", outfile, "File to write out the locations of the hits.");
    app.add_option("-e,--events", nevents, "Number of events to generate. (Default 1000)");
    app.add_option("-p,--pot", pot, "Number of protons on target (POT). (Default 10^23)");
    app.add_option("-g,--geometry", geomfile, "File to read the geometry information from");
    app.add_option("--lar-solid", lar_vol,
                   "GDML solid name for the LAr detector (default: ArgonCubeDetector)");
    app.add_option(
        "--lar-pos", lar_pos,
        "GDML define name for the LAr detector position (default: volArgonCubeDetector_pos)");
    app.add_option("--beam-radius", beam_radius,
                   "Transverse beam radius on the upstream face [cm] (default: 200)");
    app.add_option("--sigma-theta", sigma_theta,
                   "RMS beam angular divergence [rad] (default: 3e-3)");
    app.add_option("--parent-pos", parent_pos_names,
                   "GDML position-ref names for each parent volume above the LAr detector "
                   "(summed with --lar-pos to give world-frame beam centre). "
                   "Example: --parent-pos rockBox_lv_pos --parent-pos volDetEnclosure_pos");

    try {
        app.parse(argc, argv);
    } catch(const CLI::ParseError &e) { return app.exit(e); }

    if(verbosity == 1)
        spdlog::set_level(spdlog::level::debug);
    else if(verbosity == 2)
        spdlog::set_level(spdlog::level::trace);

    // Read in the geometry
    NuGeom::GDMLParser parse(geomfile);
    auto world = parse.GetWorld();
    auto box = world.GetWorldBox();

    spdlog::info("BoundingBox: min = {}, max = {}", box.min, box.max);

    // Locate the LAr detector to centre the beam on it.
    // volArgonCubeDetector_pos is relative to its immediate parent (volDetEnclosure), not the
    // world. Sum all parent position offsets to arrive at the world-frame centre.
    NuGeom::Vector3D lar_world_center{0, 0, 0};
    try {
        // Accumulate parent offsets (world → rock → hall → …)
        for(const auto &name : parent_pos_names)
            lar_world_center = lar_world_center + parse.GetPosition(name);
        // Add the LAr's own position within its immediate parent
        lar_world_center = lar_world_center + parse.GetPosition(lar_pos);
        spdlog::info("LAr world-frame centre: ({:.1f}, {:.1f}, {:.1f}) cm", lar_world_center.X(),
                     lar_world_center.Y(), lar_world_center.Z());
    } catch(const std::exception &e) {
        spdlog::warn("Could not locate LAr position ({}); beam centre set to world origin",
                     e.what());
    }

    NuGeom::Vector3D beam_center{lar_world_center.X(), lar_world_center.Y(), box.min.Z()};
    spdlog::info("Beam: centre=({:.1f},{:.1f}) cm  radius={:.1f} cm  sigma_theta={:.4f} rad",
                 beam_center.X(), beam_center.Y(), beam_radius, sigma_theta);

    // Setup DetectorSim
    NuGeom::DetectorSim sim;
    sim.Setup(world);
    sim.SetEventFile(outfile);

    // Set callbacks and initialize the interaction placement (i.e. find max probability)
    // auto raygen = std::make_shared<BeamRayGen>(0, 10, beam_center, beam_radius, sigma_theta,
    //                                            box.min.Z());
    auto raygen = std::make_shared<FileRayGen>(0, 10, "rays.log");
    auto callback = [&]() { return raygen->GetRay(); };
    sim.SetRayGenCallback(callback);

    // TODO: Fix issue with elements in GDML to PDG
    std::map<size_t, double> xsec_map = {
        {1000010010, 1e-38}, {1000080160, 1e-38}, {1000180400, 1e-38}, {1000180390, 1e-38}};
    auto event_gen = std::make_shared<TestEventGen>(xsec_map);
    sim.SetGeneratorCallback(
        [&](double energy, size_t pdg) { return event_gen->CrossSection(energy, pdg); });
    size_t ntrials = 1 << 22;
    sim.Init(ntrials);
    // sim.GenerateEvents(nevents);
    sim.GenerateEvents(pot);

    return 0;
}
