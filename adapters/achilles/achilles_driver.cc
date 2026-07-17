// Minimal end-to-end driver: run Achilles as a NuGeometry GeneratorInterface.
//
// Preferred invocation is a single YAML run card:
//
//   achilles_geom_driver <runcard.yml> [--verbosity <level>]
//
// with the schema (Driver section drives NuGeometry; the Achilles section is
// forwarded verbatim -- including !include tags -- to Achilles itself):
//
//   Driver:
//     Geometry: detector.gdml        # required: GDML geometry
//     Rays:
//       File: flux.hepmc             # required: NuHepMC (.hepmc/.hepmc3) or CSV
//       Loop: true                   # optional: rewind at EOF (default false:
//                                    #           running out of rays is an error)
//       DetectorOffset: [dx, dy, dz] # optional: beam->detector translation (cm)
//                                    #           override for NuHepMC files
//     Output: out.hepmc              # optional (default achilles_geom.hepmc)
//     Mode: envelope                 # optional: envelope | total
//     AutoEnergyRange: true          # optional (default true): widen the Achilles
//                                    #   Beam/EnergyRange to cover the flux energies;
//                                    #   false keeps the run card's range as an
//                                    #   energy cut (out-of-range rays are rejected)
//     InitRays: 1000                 # optional: rays for max_prob calibration
//     SafetyFactor: 1.5              # optional
//     Verbosity: info                # optional: trace|debug|info|warn|error|off
//     Run:                           # exactly one of:
//       Events: 1000                 #   stop after N accepted events
//       POT: 1.0e17                  #   stop once accumulated POT reached
//   Achilles:
//     <regular Achilles run card contents>
//
// The legacy flag interface is still accepted:
//
//   achilles_geom_driver --geometry <gdml> --config <achilles.yml> --rays <file>
//                        [--output <out.hepmc>] [--events N | --pot X]
//                        [--mode envelope|total] [--init N] [--safety f]
//                        [--det-offset DX DY DZ] [--loop-rays]
//                        [--no-auto-energy-range] [--verbosity <level>]
//
// Rays are streamed from the flux file one at a time (the file stays open;
// nothing is preloaded), in one of two formats:
//
//   * a NuHepMC (HepMC3 ASCII) flux file (extension .hepmc / .hepmc3), e.g.
//     produced by dune_dk2nu/dk2nu_to_nuhepmc.py.  Rays are in the beam frame;
//     the beam->detector translation recorded in the file is applied per ray.
//
//   * a comma-separated flux file (any other extension) with a header row and
//     one neutrino per line:
//         pid,wgt,E,px,py,pz,x,y,z
//     with pid the neutrino PDG code, wgt the per-sample flux weight, energy E
//     in GeV, (px,py,pz) the ray direction (normalized on load), and (x,y,z)
//     the ray origin in the geometry's length units.

#include "adapters/achilles/AchillesAdapter.hh"
#include "geom/DetectorSim.hh"
#include "geom/FluxSource.hh"
#include "geom/FluxStreamer.hh"
#include "geom/Logging.hh"
#include "geom/Ray.hh"
#include "geom/Vector3D.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string Option(int argc, char **argv, const std::string &key, const std::string &fallback) {
    for(int i = 1; i + 1 < argc; ++i)
        if(key == argv[i]) return argv[i + 1];
    return fallback;
}

bool HasOption(int argc, char **argv, const std::string &key) {
    for(int i = 1; i < argc; ++i)
        if(key == argv[i]) return true;
    return false;
}

// Map a verbosity setting (name or 0-5) onto an spdlog level.
spdlog::level::level_enum ParseVerbosity(const std::string &arg) {
    if(arg == "trace" || arg == "0") return spdlog::level::trace;
    if(arg == "debug" || arg == "1") return spdlog::level::debug;
    if(arg == "info" || arg == "2") return spdlog::level::info;
    if(arg == "warn" || arg == "warning" || arg == "3") return spdlog::level::warn;
    if(arg == "error" || arg == "4") return spdlog::level::err;
    if(arg == "off" || arg == "5") return spdlog::level::off;
    throw std::runtime_error("achilles_geom_driver: unknown verbosity level: " + arg +
                             " (expected trace|debug|info|warn|error|off or 0-5)");
}

// Everything the driver needs for one run, from either input style.
struct DriverConfig {
    std::string geometry;
    std::string rays_file;
    bool loop_rays = false;
    bool have_offset = false;
    NuGeom::Vector3D offset;

    std::string output = "achilles_geom.hepmc";
    std::string mode = "envelope";
    bool auto_energy_range = true;
    std::size_t init_rays = 1000;
    double safety = 1.5;
    std::string verbosity = "info";

    // Exactly one of the two run targets is active; events==0 means POT mode.
    std::size_t events = 0;
    double pot = 0.0;

    // The Achilles section (run-card style) or path (legacy flags).
    std::optional<YAML::Node> achilles_node;
    std::string achilles_config_path;
    std::string runcard_name;
};

DriverConfig FromRunCard(const std::string &path) {
    DriverConfig cfg;
    cfg.runcard_name = path;
    YAML::Node card = YAML::LoadFile(path);

    if(!card["Driver"]) throw std::runtime_error("run card is missing the Driver section: " + path);
    if(!card["Achilles"])
        throw std::runtime_error("run card is missing the Achilles section: " + path);
    const YAML::Node driver = card["Driver"];
    cfg.achilles_node = card["Achilles"];

    if(!driver["Geometry"]) throw std::runtime_error("Driver/Geometry is required");
    cfg.geometry = driver["Geometry"].as<std::string>();

    if(!driver["Rays"]) throw std::runtime_error("Driver/Rays is required");
    const YAML::Node rays = driver["Rays"];
    if(rays.IsScalar()) {
        cfg.rays_file = rays.as<std::string>();
    } else {
        if(!rays["File"]) throw std::runtime_error("Driver/Rays/File is required");
        cfg.rays_file = rays["File"].as<std::string>();
        if(rays["Loop"]) cfg.loop_rays = rays["Loop"].as<bool>();
        if(rays["DetectorOffset"]) {
            const auto v = rays["DetectorOffset"].as<std::vector<double>>();
            if(v.size() != 3)
                throw std::runtime_error("Driver/Rays/DetectorOffset needs three values");
            cfg.offset = NuGeom::Vector3D{v[0], v[1], v[2]};
            cfg.have_offset = true;
        }
    }

    if(driver["Output"]) cfg.output = driver["Output"].as<std::string>();
    if(driver["Mode"]) cfg.mode = driver["Mode"].as<std::string>();
    if(driver["AutoEnergyRange"]) cfg.auto_energy_range = driver["AutoEnergyRange"].as<bool>();
    if(driver["InitRays"]) cfg.init_rays = driver["InitRays"].as<std::size_t>();
    if(driver["SafetyFactor"]) cfg.safety = driver["SafetyFactor"].as<double>();
    if(driver["Verbosity"]) cfg.verbosity = driver["Verbosity"].as<std::string>();

    if(!driver["Run"]) throw std::runtime_error("Driver/Run is required (Events or POT)");
    const YAML::Node run = driver["Run"];
    const bool has_events = static_cast<bool>(run["Events"]);
    const bool has_pot = static_cast<bool>(run["POT"]);
    if(has_events == has_pot)
        throw std::runtime_error("Driver/Run must set exactly one of Events or POT");
    if(has_events)
        cfg.events = run["Events"].as<std::size_t>();
    else
        cfg.pot = run["POT"].as<double>();

    return cfg;
}

DriverConfig FromFlags(int argc, char **argv) {
    DriverConfig cfg;
    cfg.geometry = Option(argc, argv, "--geometry", "");
    cfg.achilles_config_path = Option(argc, argv, "--config", "");
    cfg.rays_file = Option(argc, argv, "--rays", "");
    cfg.output = Option(argc, argv, "--output", cfg.output);
    cfg.mode = Option(argc, argv, "--mode", cfg.mode);
    cfg.init_rays = static_cast<std::size_t>(std::stoul(Option(argc, argv, "--init", "1000")));
    cfg.safety = std::stod(Option(argc, argv, "--safety", "1.5"));
    cfg.verbosity = Option(argc, argv, "--verbosity", cfg.verbosity);
    cfg.loop_rays = HasOption(argc, argv, "--loop-rays");
    cfg.auto_energy_range = !HasOption(argc, argv, "--no-auto-energy-range");
    cfg.runcard_name = cfg.achilles_config_path;

    const std::string pot_s = Option(argc, argv, "--pot", "");
    if(!pot_s.empty())
        cfg.pot = std::stod(pot_s);
    else
        cfg.events = static_cast<std::size_t>(std::stoul(Option(argc, argv, "--events", "1000")));

    if(HasOption(argc, argv, "--det-offset")) {
        for(int i = 1; i + 3 < argc; ++i) {
            if(std::string(argv[i]) == "--det-offset") {
                cfg.offset = NuGeom::Vector3D{std::stod(argv[i + 1]), std::stod(argv[i + 2]),
                                              std::stod(argv[i + 3])};
                cfg.have_offset = true;
                break;
            }
        }
        if(!cfg.have_offset)
            throw std::runtime_error("--det-offset requires three values: DX DY DZ");
    }

    return cfg;
}

void Usage() {
    std::cerr << "usage: achilles_geom_driver <runcard.yml> [--verbosity <level>]\n"
                 "   or: achilles_geom_driver --geometry <gdml> --config <achilles.yml> "
                 "--rays <file> [--output <out.hepmc>] [--events N | --pot X] "
                 "[--mode envelope|total] [--init N] [--safety f] [--det-offset DX DY DZ] "
                 "[--loop-rays] [--no-auto-energy-range] "
                 "[--verbosity trace|debug|info|warn|error|off]\n";
}

} // namespace

int main(int argc, char **argv) {
    try {
        // Run-card mode: first argument is a YAML file (not a --flag).
        DriverConfig cfg;
        if(argc > 1 && argv[1][0] != '-') {
            cfg = FromRunCard(argv[1]);
        } else if(HasOption(argc, argv, "--geometry")) {
            cfg = FromFlags(argc, argv);
        } else {
            Usage();
            return 1;
        }

        if(cfg.geometry.empty() || cfg.rays_file.empty() ||
           (!cfg.achilles_node && cfg.achilles_config_path.empty())) {
            Usage();
            return 1;
        }

        // Command-line verbosity overrides the run card.
        const std::string verb = Option(argc, argv, "--verbosity", cfg.verbosity);

        // Scoped logging: NuGeometry logs through the registered "nugeom"
        // logger, Achilles through the default logger (renamed "achilles"
        // here), both sharing one console sink.  The [%n] pattern field shows
        // which side a message came from.  trace shows per-trial detail;
        // debug shows per-event/per-ray decisions.
        const auto level = ParseVerbosity(verb);
        auto sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        auto nugeom_logger = std::make_shared<spdlog::logger>("nugeom", sink);
        auto achilles_logger = std::make_shared<spdlog::logger>("achilles", sink);
        for(auto &logger : {nugeom_logger, achilles_logger}) {
            logger->set_level(level);
            logger->set_pattern("[%n] [%^%l%$] %v");
        }
        spdlog::register_logger(nugeom_logger);
        spdlog::set_default_logger(achilles_logger);
        NuGeom::RefreshLogger(); // point NuGeom::Log()'s cached handle here

        const bool total = (cfg.mode == "total");
        if(!total && cfg.mode != "envelope")
            throw std::runtime_error("mode must be 'envelope' or 'total', got: " + cfg.mode);
        const auto sampling = total ? NuGeom::DetectorSim::SamplingMode::TotalXSecRetry
                                    : NuGeom::DetectorSim::SamplingMode::EnvelopeNoRetry;

        // Stream rays on demand; the constructor's summary pass reports the
        // energy support so a mismatch with the Achilles beam EnergyRange
        // (which must cover it, in MeV) is visible up front.
        auto streamer = NuGeom::OpenFluxStreamer(cfg.rays_file, cfg.loop_rays,
                                                 cfg.have_offset ? &cfg.offset : nullptr);
        std::cout << "Streaming " << streamer->Count() << " rays per pass from " << cfg.rays_file
                  << " (E in [" << streamer->EMin() << ", " << streamer->EMax()
                  << "] GeV; loop=" << (cfg.loop_rays ? "on" : "off") << "); mode=" << cfg.mode
                  << "\n";

        auto adapter =
            cfg.achilles_node
                ? std::make_shared<NuGeom::AchillesAdapter>(*cfg.achilles_node, cfg.runcard_name)
                : std::make_shared<NuGeom::AchillesAdapter>(cfg.achilles_config_path);
        adapter->SetTotalXSecRetry(total);

        // Align the Achilles Beam/EnergyRange (MeV) with the streamed flux
        // support (GeV) so no ray is rejected as out-of-range.  Opt out with
        // AutoEnergyRange: false / --no-auto-energy-range to keep the run
        // card's range as an energy cut.
        if(cfg.auto_energy_range)
            adapter->EnsureBeamEnergyCoverage(streamer->EMin() * 1000.0, streamer->EMax() * 1000.0);

        NuGeom::DetectorSim sim(cfg.safety);
        sim.Setup(cfg.geometry);
        sim.SetSamplingMode(sampling);
        sim.SetGenerator(adapter);
        sim.SetFluxCallback([&streamer]() { return streamer->Next(); });
        sim.SetEventFile(cfg.output);

        std::cout << "Calibrating max_prob over " << cfg.init_rays << " rays...\n";
        sim.Init(cfg.init_rays);

        if(cfg.events > 0) {
            std::cout << "Generating " << cfg.events << " events...\n";
            sim.GenerateEvents(cfg.events);
        } else {
            std::cout << "Generating until POT = " << cfg.pot << "...\n";
            sim.GenerateEvents(cfg.pot);
        }

        std::cout << "Done. Accumulated POT = " << sim.GetAccumulatedPOT() << " over "
                  << streamer->Loops() << " full flux-file passes; wrote " << cfg.output << "\n";
        return 0;
    } catch(const std::exception &e) {
        std::cerr << "achilles_geom_driver error: " << e.what() << "\n";
        return 1;
    }
}
