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

/*
// Alternative
double energy = 1.0; // GeV //Remove Later
sim.Setup(world);
// sim.SetMaxProb(1e-3);
double max_prob = 1e-3;
TestEventGen gen1(xsec_map);

// Remove Later
std::vector<NuGeom::Ray> rays;
size_t ntest_rays = 1 << 20;
for(size_t i = 0; i < ntest_rays; ++i) { rays.push_back(raygen->GetRay().second); }
for(const auto &ray : rays) {
    std::vector<NuGeom::LineSegment> segments = sim.GetLineSegments(ray);
    std::set<NuGeom::Material> materials = sim.GetMaterials(segments);
    std::map<NuGeom::Material, double> xsecs = gen1.EvaluateCrossSections(energy, materials);
    std::vector<double> probs = sim.EvaluateProbs(segments, xsecs);
    double prob_sum = std::accumulate(probs.begin(), probs.end(), 0.0);
    if(prob_sum > max_prob) { max_prob = prob_sum; }
    return max_prob;
}

std::ofstream hist;
hist.open("hit_locations2.txt", std::ios::app);
if(!hist.is_open()) { throw std::runtime_error("Failed to open hit_locations2.txt"); }
for(const auto &ray : rays) {
    auto ray = raygen->GetRay();
    std::vector<NuGeom::LineSegment> segments0 = sim.GetLineSegments(ray);
    std::set<NuGeom::Material> materials = sim.GetMaterials(segments0);
    std::map<NuGeom::Material, double> material_xsec_map =
        gen1.EvaluateCrossSections(energy, materials);
    NuGeom::Vector3D interaction_point = sim.Interaction(segments0, material_xsec_map);
    if(interaction_point.X() < 9e9) {
        // Valid interaction point
        hist << interaction_point.X() << " " << interaction_point.Y() << " "
             << interaction_point.Z() << "\n";
    } else {
        // No interaction occurred
        spdlog::info("No interaction occurred for this ray.");
    }
}
hist.close();
}

// Calculate interaction location
std::map<std::string,double> meanfreepaths{
{"Water",1e+7},
{"Air",1e+7},
{"Argon",1e+7}
};

// Shoot Rays and get LineSegments
double prob_max = 0;
for(size_t i = 0; i < ntrials; ++i) {
NuGeom::Ray ray = ShootRay({-2, -2, -2}, {2, 2, -2});
std::vector<NuGeom::LineSegment> segments0;
world->GetLineSegments(ray, segments0);

std::vector<double> probs0(segments0.size());
std::vector<double> probs1(segments0.size());
std::vector<double> seglength0(segments0.size());

std::vector<std::string> material0(segments0.size());

// Calculate probability to interact for each line segment
for (size_t j=0; j<segments0.size(); ++j){
    material0[j]=segments0[j].GetMaterial().Name();

    // NOTE: This only works for l/meanfreepath tiny
    probs0[j]=segments0[j].Length()/meanfreepaths[material0[j]];
    probs1[j]=1-exp(-segments0[j].Length()/meanfreepaths[material0[j]]);

    // For testing only:
    seglength0[j] = segments0[j].Length();
}
double normconst0=std::accumulate(probs0.begin(), probs0.end(), 0.0);
if (normconst0>prob_max) {
    std::cout << "Trial: " << i << " Ray: O(" << ray.Origin() << ") D(" << ray.Direction() <<
")\n"; std::cout << "Hit " << segments0.size() << " segments\n"; double
normconst1=std::accumulate(probs1.begin(), probs1.end(), 0.0); std::cout << "Exact = " << normconst1
<< " Approx = " << normconst0 << "\n"; std::cout << "Old prob max = " << prob_max << "\n"; std::cout
<< "New prob max = " << normconst0 << "\n"; prob_max=normconst0;
}
}

constexpr double safety_factor = 1.5;
prob_max *= safety_factor;

auto rand = NuGeom::Random::Instance();
double nrays = 0;
size_t nhits = 0;
std::ofstream hist;
hist.open("hit_locations.txt");

// Send Ray for check to hit
for(size_t i = 0; i < ntrials; ++i) {
nrays += 1.0/prob_max;
if(i % 10000 == 0) {
    std::cout << "Shot " << nrays << " rays\n";
}
NuGeom::Ray ray = ShootRay({-2, -2, -2}, {2, 2, -2});
std::vector<NuGeom::LineSegment> segments0;

world->GetLineSegments(ray, segments0);

std::vector<double> probs0(segments0.size());
std::vector<double> probs1(segments0.size());
std::vector<double> seglength0(segments0.size());

std::vector<std::string> material0(segments0.size());

// Calculate probability to interact for each line segment
for (size_t j=0; j<segments0.size(); ++j){
    material0[j]=segments0[j].GetMaterial().Name();
    // NOTE: This only works for l/meanfreepath tiny
    probs0[j]=segments0[j].Length()/meanfreepaths[material0[j]];
    // std::cout << -segments0[i].Length()/meanfreepaths[material0[i]] << std::endl;
    probs1[j]=1-exp(-segments0[j].Length()/meanfreepaths[material0[j]]);

    // For testing only:
    seglength0[j] = segments0[j].Length();
}
double normconst0=std::accumulate(probs0.begin(), probs0.end(), 0.0);
if(normconst0 / prob_max < rand.Uniform(0.0, 1.0)) {
    continue;
}
nhits++;
std::cout << nhits << " / " << nrays << "\n";

for(size_t j=0; j<segments0.size(); ++j) probs0[j]/=normconst0;
double sum = 0;
size_t idx = 0;
double r2 = rand.Uniform(0.0, 1.0);
NuGeom::Vector3D position;
for(size_t j = 0; j < probs0.size(); ++j) {
    sum += probs0[j];
    if(sum > r2) {
        idx = j;
        position = segments0[j].Start() + (segments0[j].End() - segments0[j].Start()) * (sum -
r2)/probs0[j]; break;
    }
}

std::cout << "Start: " << segments0[idx].Start() << "\n";
std::cout << "End: " << segments0[idx].End() << "\n";
std::cout << "Sum: " << sum << " R2: " << r2 << " Prob: " << probs0[idx] << " Dist: " << (sum -
r2)/probs0[idx] << "\n"; std::cout << "Hit Location: " << position << "\n";

hist << position.X() << "," << position.Y() << "," << position.Z() << "\n";
if(nhits == 10000) break;
}

hist.close();

return 0;
*/
