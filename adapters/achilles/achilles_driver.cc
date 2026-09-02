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
//                                    #   energy cut (out-of-range rays are rejected).
//                                    #   Needs the flux's energy support, which is
//                                    #   read from the file's
//                                    #   NuGeom.Flux.EnergyRange_GeV attribute or
//                                    #   from FluxEnergyRange below -- the flux is
//                                    #   never scanned to derive it
//     FluxEnergyRange: [lo, hi]      # optional: flux energy support in GeV, for
//                                    #   flux files that do not record it
//     InitRays: 500                  # optional: self-generated probe rays used to
//                                    #   calibrate the layer-1 envelope (NOT flux rays)
//     BurnInRays: 1000               # optional: flux rays consumed (untraced) to fix
//                                    #   the envelope's flux-weight scale before the run
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
//                        [--mode envelope|total] [--init N] [--burn-in N] [--safety f]
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
    // Flux energy support in GeV, when the flux file does not record it.  Only
    // used to widen the Achilles beam range; NuGeometry never scans the file
    // for it.
    bool have_flux_energy_range = false;
    double flux_emin = 0.0, flux_emax = 0.0;
    // Probe rays for the envelope scan.  These are chords through the world
    // box, not flux rays, and each is traced once for the whole energy grid, so
    // a few hundred already pin the maximum column density.
    std::size_t init_rays = 500;
    // Flux rays consumed (without tracing) to establish the envelope's
    // flux-weight scale before any exposure is charged.
    std::size_t burn_in_rays = 1000;
    double safety = 1.5;
    std::string verbosity = "info";

    // Exactly one of the two run targets is active; events==0 means POT mode.
    std::size_t events = 0;
    double pot = 0.0;

    // The run card's root document.  It MUST outlive `achilles_node`, which is
    // a sub-node of it: a yaml-cpp Node shares the document's node arena, and
    // letting the root go while a sub-node is still in use leaves that arena
    // half-owned -- the symptom is a segfault inside free() the next time
    // yaml-cpp grows a vector, typically while Achilles parses its own config.
    std::optional<YAML::Node> runcard_root;
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
    cfg.runcard_root = card;
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
    if(driver["FluxEnergyRange"]) {
        const auto v = driver["FluxEnergyRange"].as<std::vector<double>>();
        if(v.size() != 2) throw std::runtime_error("Driver/FluxEnergyRange needs two values (GeV)");
        cfg.flux_emin = v[0];
        cfg.flux_emax = v[1];
        cfg.have_flux_energy_range = true;
    }
    if(driver["InitRays"]) cfg.init_rays = driver["InitRays"].as<std::size_t>();
    if(driver["BurnInRays"]) cfg.burn_in_rays = driver["BurnInRays"].as<std::size_t>();
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
    cfg.init_rays = static_cast<std::size_t>(std::stoul(Option(argc, argv, "--init", "500")));
    cfg.burn_in_rays =
        static_cast<std::size_t>(std::stoul(Option(argc, argv, "--burn-in", "1000")));
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
                 "[--mode envelope|total] [--init N] [--burn-in N] [--safety f] "
                 "[--det-offset DX DY DZ] "
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

        // Stream rays on demand.  Opening a NuHepMC flux reads one event for
        // the run-level metadata and (only if the converter did not record the
        // ray count) makes one raw byte-count pass -- it never parses the rays.
        auto streamer = NuGeom::OpenFluxStreamer(cfg.rays_file, cfg.loop_rays,
                                                 cfg.have_offset ? &cfg.offset : nullptr);
        std::cout << "Streaming " << streamer->Count() << " rays per pass from " << cfg.rays_file
                  << " (";
        if(streamer->HasEnergyRange())
            std::cout << "E in [" << streamer->EMin() << ", " << streamer->EMax() << "] GeV; ";
        else
            std::cout << "energy range not recorded; ";
        std::cout << "loop=" << (cfg.loop_rays ? "on" : "off") << "); mode=" << cfg.mode << "\n";

        auto adapter =
            cfg.achilles_node
                ? std::make_shared<NuGeom::AchillesAdapter>(*cfg.achilles_node, cfg.runcard_name)
                : std::make_shared<NuGeom::AchillesAdapter>(cfg.achilles_config_path);
        adapter->SetTotalXSecRetry(total);

        // Align the Achilles Beam/EnergyRange (MeV) with the streamed flux
        // support (GeV) so no ray is rejected as out-of-range.  Opt out with
        // AutoEnergyRange: false / --no-auto-energy-range to keep the run
        // card's range as an energy cut.  A flux file that does not record its
        // energy range leaves the run card's range as the authority.
        if(cfg.auto_energy_range) {
            const bool from_file = streamer->HasEnergyRange();
            if(from_file || cfg.have_flux_energy_range) {
                const double lo = from_file ? streamer->EMin() : cfg.flux_emin;
                const double hi = from_file ? streamer->EMax() : cfg.flux_emax;
                adapter->EnsureBeamEnergyCoverage(lo * 1000.0, hi * 1000.0);
            } else {
                std::cout << "AutoEnergyRange requested but '" << cfg.rays_file
                          << "' does not record NuGeom.Flux.EnergyRange_GeV and no "
                             "Driver/FluxEnergyRange was given; keeping the run card's "
                             "Beam/EnergyRange. Rays outside it are rejected as zero cross "
                             "section -- set Driver/FluxEnergyRange (GeV) or regenerate the "
                             "flux to record its range.\n";
            }
        }

        // The envelope scan needs the energies the generator actually supports.
        const auto beam_range = adapter->BeamEnergyRange();
        if(!beam_range)
            throw std::runtime_error(
                "no Achilles Beam/EnergyRange to calibrate the envelope against: set it in the "
                "run card (MeV), or regenerate the flux so it records "
                "NuGeom.Flux.EnergyRange_GeV and leave AutoEnergyRange on");

        NuGeom::DetectorSim sim(cfg.safety);
        sim.Setup(cfg.geometry);
        sim.SetSamplingMode(sampling);
        sim.SetGenerator(adapter);
        sim.SetFluxCallback([&streamer]() { return streamer->Next(); });
        sim.SetEventFile(cfg.output);
        sim.SetEnergyRange(beam_range->first / 1000.0, beam_range->second / 1000.0);
        sim.SetEnvelopeBurnIn(cfg.burn_in_rays);

        // Calibrate on rays NuGeometry generates itself -- chords through the
        // world box scanned over the supported energies -- so the flux file is
        // not consumed (or even read) to find the geometry envelope, and the
        // bound does not depend on how much of the flux's tail a warm-up
        // sample happened to contain.  The remaining flux-weight factor is
        // learned from the burn-in rays and kept up to date during the run.
        std::cout << "Calibrating the envelope on " << cfg.init_rays << " probe rays over E in ["
                  << beam_range->first / 1000.0 << ", " << beam_range->second / 1000.0
                  << "] GeV...\n";
        sim.Init(cfg.init_rays);

        if(cfg.events > 0) {
            std::cout << "Generating " << cfg.events << " events...\n";
            sim.GenerateEvents(cfg.events);
        } else {
            std::cout << "Generating until POT = " << cfg.pot << "...\n";
            sim.GenerateEvents(cfg.pot);
        }

        const auto &stats = sim.GetRunStats();
        std::cout << "Done. Accumulated POT = " << sim.GetAccumulatedPOT() << " over "
                  << streamer->Loops() << " full flux-file passes; wrote " << stats.emitted
                  << " events to " << cfg.output << "\n";
        // The two modes only agree on the physical rate through different
        // estimators when the generator's unweighter is partial (Achilles'
        // Percentile unweighter is): events/POT is the rate in `total` mode,
        // sum(weights)/POT in `envelope` mode.  Print both so a mode-to-mode
        // comparison is made on like for like -- comparing raw event counts
        // shows a systematic <weight> discrepancy that is not statistical.
        if(stats.emitted > 0 && sim.GetAccumulatedPOT() > 0) {
            const double n = static_cast<double>(stats.emitted);
            std::cout << "  events/POT       = " << n / sim.GetAccumulatedPOT() << "\n"
                      << "  sum(weights)/POT = " << stats.sum_weight / sim.GetAccumulatedPOT()
                      << "\n"
                      << "  <weight>         = " << stats.sum_weight / n << "  (the factor by "
                      << "which the two modes' raw event counts differ)\n";
        }
        return 0;
    } catch(const std::exception &e) {
        std::cerr << "achilles_geom_driver error: " << e.what() << "\n";
        return 1;
    }
}
