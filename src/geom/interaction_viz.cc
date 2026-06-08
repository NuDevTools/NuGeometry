/// interaction_viz — shoot rays through a geometry, record interactions, and
/// visualise the outgoing-particle paths.
///
/// A headless companion to the interactive BeamViz viewer: it shoots a beam of
/// rays through a GDML geometry, probabilistically decides whether each ray
/// interacts (using a simple cross-section model to derive a mean free path per
/// material), generates outgoing particles at each interaction vertex, traces
/// them back through the geometry, and writes:
///
///   * a console summary (how many rays interacted, in which materials),
///   * a human-readable event log of every primary path, vertex and secondary,
///   * a self-contained SVG showing the geometry, the incoming rays (grey), the
///     interaction vertices (red) and the outgoing-particle tracks (coloured).
///
/// The SVG needs no external dependencies and opens in any browser, which makes
/// it a convenient way to eyeball what a geometry is doing and to debug ray
/// tracing without a GUI.
///
/// Usage:
///   interaction_viz -g geom.gdml -n 200 --svg event.svg -o events.log
///   interaction_viz -g geom.gdml -n 50 --xsec 5e-24 --secondaries 5 --proj zy

#include "CLI/CLI.hpp"
#include "fmt/format.h"
#include "geom/DetectorSim.hh" // defines NuGeom::EnergyRay used by TestGen.hh
#include "geom/InteractionViz.hh"
#include "geom/Parser.hh"
#include "geom/Random.hh"
#include "geom/TestGen.hh"
#include "geom/World.hh"

#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Event log
// ---------------------------------------------------------------------------
void write_event_log(const std::string &path, const std::vector<NuGeom::InteractionEvent> &events) {
    std::ofstream f(path);
    if(!f) {
        spdlog::error("Cannot open event log for writing: {}", path);
        return;
    }
    auto vec = [](const NuGeom::Vector3D &v) {
        return fmt::format("({:+.4g}, {:+.4g}, {:+.4g})", v.X(), v.Y(), v.Z());
    };
    for(size_t i = 0; i < events.size(); ++i) {
        const auto &e = events[i];
        f << fmt::format("EVENT {} energy={:.4g} interacted={} optical_depth={:.4g}\n", i, e.energy,
                         e.interacted ? 1 : 0, e.optical_depth);
        f << fmt::format("  incoming: origin={} dir={} nseg={}\n", vec(e.incoming.Origin()),
                         vec(e.incoming.Direction()), e.incoming_path.size());
        if(e.interacted) {
            f << fmt::format("  vertex: {} material={}\n", vec(e.vertex), e.vertex_material);
            for(size_t k = 0; k < e.secondaries.size(); ++k) {
                const auto &s = e.secondaries[k];
                f << fmt::format("  secondary {}: pdg={} energy={:.4g} dir={} nseg={} exited={}\n",
                                 k, s.pdg, s.energy, vec(s.direction), s.path.size(),
                                 s.exited_world ? 1 : 0);
            }
        }
    }
    spdlog::info("Wrote event log: {}", path);
}

// ---------------------------------------------------------------------------
// SVG projection
// ---------------------------------------------------------------------------
struct Projector {
    size_t ha{2}, va{1}; // data axis index for horizontal / vertical
    double minh, maxh, minv, maxv;
    double scale, ox, oy; // pixel transform
    double margin, height;

    Projector(const NuGeom::BoundingBox &box, size_t h_axis, size_t v_axis, double px_w,
              double px_h, double m)
        : ha{h_axis}, va{v_axis}, margin{m}, height{px_h} {
        std::array<double, 3> lo{box.min.X(), box.min.Y(), box.min.Z()};
        std::array<double, 3> hi{box.max.X(), box.max.Y(), box.max.Z()};
        minh = lo[ha];
        maxh = hi[ha];
        minv = lo[va];
        maxv = hi[va];
        double sx = (px_w - 2 * margin) / std::max(1e-9, maxh - minh);
        double sy = (px_h - 2 * margin) / std::max(1e-9, maxv - minv);
        scale = std::min(sx, sy); // uniform scale keeps aspect ratio
        ox = margin;
        oy = margin;
    }

    double px(const NuGeom::Vector3D &p) const {
        std::array<double, 3> c{p.X(), p.Y(), p.Z()};
        return ox + (c[ha] - minh) * scale;
    }
    double py(const NuGeom::Vector3D &p) const {
        std::array<double, 3> c{p.X(), p.Y(), p.Z()};
        // SVG y axis points down, so flip the vertical data axis.
        return oy + (maxv - c[va]) * scale;
    }
};

// A small distinct palette for secondary tracks / materials.
const std::vector<std::string> kPalette = {"#e6194b", "#3cb44b", "#4363d8", "#f58231",
                                           "#911eb4", "#46f0f0", "#f032e6", "#bcf60c"};

std::string color_for(const std::string &key) {
    static std::map<std::string, std::string> assigned;
    auto it = assigned.find(key);
    if(it != assigned.end()) return it->second;
    std::string c = kPalette[assigned.size() % kPalette.size()];
    assigned.emplace(key, c);
    return c;
}

void write_svg(const std::string &path, NuGeom::World &world,
               const std::vector<NuGeom::InteractionEvent> &events, int h_axis, int v_axis,
               int max_depth) {
    std::ofstream f(path);
    if(!f) {
        spdlog::error("Cannot open SVG for writing: {}", path);
        return;
    }

    constexpr double W = 1200, H = 800, M = 40;
    auto box = world.GetWorldBox();
    Projector proj(box, static_cast<size_t>(h_axis), static_cast<size_t>(v_axis), W, H, M);

    const char *axis_name[] = {"X", "Y", "Z"};

    f << fmt::format("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{}\" height=\"{}\" "
                     "viewBox=\"0 0 {} {}\">\n",
                     W, H, W, H);
    f << fmt::format("<rect width=\"{}\" height=\"{}\" fill=\"#11131a\"/>\n", W, H);

    // Volume boxes (axis-aligned bounds projected to a rectangle).
    auto bounds = world.GetVolumeBounds(max_depth);
    for(const auto &vb : bounds) {
        NuGeom::Vector3D c0{vb.bb.min}, c1{vb.bb.max};
        double x0 = proj.px(c0), x1 = proj.px(c1);
        double y0 = proj.py(c0), y1 = proj.py(c1);
        std::string col = color_for("mat:" + vb.material.Name());
        f << fmt::format("<rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" "
                         "fill=\"{}\" fill-opacity=\"0.08\" stroke=\"{}\" stroke-opacity=\"0.5\" "
                         "stroke-width=\"1\"/>\n",
                         std::min(x0, x1), std::min(y0, y1), std::abs(x1 - x0), std::abs(y1 - y0),
                         col, col);
    }

    // World box outline.
    f << fmt::format("<rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" "
                     "fill=\"none\" stroke=\"#888\" stroke-width=\"1.5\"/>\n",
                     proj.px(box.min), proj.py(box.max),
                     std::abs(proj.px(box.max) - proj.px(box.min)),
                     std::abs(proj.py(box.min) - proj.py(box.max)));

    // Incoming rays (grey) and outgoing-particle tracks (coloured).
    for(const auto &e : events) {
        for(const auto &seg : e.incoming_path) {
            if(!std::isfinite(seg.Length())) continue;
            f << fmt::format("<line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                             "stroke=\"#9aa\" stroke-width=\"0.8\" stroke-opacity=\"0.5\"/>\n",
                             proj.px(seg.Start()), proj.py(seg.Start()), proj.px(seg.End()),
                             proj.py(seg.End()));
        }
        if(!e.interacted) continue;
        for(size_t k = 0; k < e.secondaries.size(); ++k) {
            std::string col = color_for(fmt::format("sec:{}", e.secondaries[k].pdg));
            for(const auto &seg : e.secondaries[k].path) {
                if(!std::isfinite(seg.Length())) continue;
                f << fmt::format("<line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                                 "stroke=\"{}\" stroke-width=\"1.3\" stroke-opacity=\"0.9\"/>\n",
                                 proj.px(seg.Start()), proj.py(seg.Start()), proj.px(seg.End()),
                                 proj.py(seg.End()), col);
            }
        }
        // Interaction vertex marker.
        f << fmt::format("<circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"3\" fill=\"#ff3030\"/>\n",
                         proj.px(e.vertex), proj.py(e.vertex));
    }

    // Axis label.
    f << fmt::format(
        "<text x=\"{}\" y=\"{}\" fill=\"#bbb\" font-family=\"monospace\" "
        "font-size=\"14\">projection: {} (horizontal) vs {} (vertical), units cm</text>\n",
        M, H - 12, axis_name[h_axis], axis_name[v_axis]);

    f << "</svg>\n";
    spdlog::info("Wrote SVG visualisation: {}", path);
}

int axis_from_char(char c) {
    switch(c) {
    case 'x':
    case 'X':
        return 0;
    case 'y':
    case 'Y':
        return 1;
    case 'z':
    case 'Z':
        return 2;
    default:
        return -1;
    }
}

} // namespace

int main(int argc, char **argv) {
    auto console = spdlog::stdout_color_mt("NuGeom");
    spdlog::set_default_logger(console);
    spdlog::set_pattern("[%n] [%^%l%$] %v");
    spdlog::set_level(spdlog::level::info);

    CLI::App app("interaction_viz: shoot rays, record interactions, visualise outgoing particles");
    argv = app.ensure_utf8(argv);

    std::string geomfile;
    std::string svgfile = "interactions.svg";
    std::string logfile = "interactions.log";
    std::string proj = "zy";
    size_t nrays = 100;
    double emin = 1.0, emax = 10.0;
    double beam_radius = -1.0; // <0 → auto from world size
    double sigma_theta = 0.01;
    double xsec = 1e-24;
    int nsecondaries = 3;
    double cone = 0.35;
    int max_depth = 4;
    int seed = -1;
    bool verbose = false;

    app.add_option("-g,--geometry", geomfile, "GDML geometry file")->required();
    app.add_option("-n,--nrays", nrays, "Number of rays to shoot (default 100)");
    app.add_option("-o,--log", logfile, "Event log output file (default interactions.log)");
    app.add_option("--svg", svgfile, "SVG visualisation output file (default interactions.svg)");
    app.add_option("--proj", proj, "Projection plane: xy, zy, zx, ... (default zy)");
    app.add_option("--emin", emin, "Minimum beam energy [GeV] (default 1)");
    app.add_option("--emax", emax, "Maximum beam energy [GeV] (default 10)");
    app.add_option("--radius", beam_radius, "Beam spot radius [cm] (default: 5% of world size)");
    app.add_option("--sigma", sigma_theta, "Beam angular divergence [rad] (default 0.01)");
    app.add_option("--xsec", xsec, "Constant cross section per nucleus [cm^2] (default 1e-24)");
    app.add_option("--secondaries", nsecondaries, "Outgoing particles per interaction (default 3)");
    app.add_option("--cone", cone, "Secondary cone half-angle [rad] (default 0.35)");
    app.add_option("--max-depth", max_depth, "Volume nesting depth drawn in the SVG (default 4)");
    app.add_option("-s,--seed", seed, "Random seed (default: nondeterministic)");
    app.add_flag("-v,--verbose", verbose, "Verbose (debug) logging");

    try {
        app.parse(argc, argv);
    } catch(const CLI::ParseError &e) { return app.exit(e); }

    if(verbose) spdlog::set_level(spdlog::level::debug);
    if(seed >= 0) NuGeom::Random::Instance().Seed(static_cast<unsigned int>(seed));

    if(proj.size() != 2) {
        spdlog::error("--proj must be two axis letters, e.g. zy");
        return 1;
    }
    int h_axis = axis_from_char(proj[0]);
    int v_axis = axis_from_char(proj[1]);
    if(h_axis < 0 || v_axis < 0 || h_axis == v_axis) {
        spdlog::error("Invalid --proj '{}': use two distinct letters from x/y/z", proj);
        return 1;
    }

    NuGeom::World world;
    try {
        NuGeom::GDMLParser parser(geomfile);
        world = parser.GetWorld();
    } catch(const std::exception &e) {
        spdlog::error("Failed to load geometry '{}': {}", geomfile, e.what());
        return 1;
    }
    auto box = world.GetWorldBox();
    spdlog::info("Loaded geometry {}", geomfile);
    spdlog::info("World box: ({:.1f},{:.1f},{:.1f}) -> ({:.1f},{:.1f},{:.1f}) cm", box.min.X(),
                 box.min.Y(), box.min.Z(), box.max.X(), box.max.Y(), box.max.Z());

    double diag = (box.max - box.min).Norm();
    if(beam_radius < 0) beam_radius = 0.05 * diag;

    // Beam fired along +Z from the upstream (min-Z) face, centred transversely.
    NuGeom::Vector3D center{0.5 * (box.min.X() + box.max.X()), 0.5 * (box.min.Y() + box.max.Y()),
                            box.min.Z()};
    BeamRayGen gen(emin, emax, center, beam_radius, sigma_theta, box.min.Z());

    NuGeom::RayInteractionSim sim(world);
    sim.SetConstantCrossSection(xsec);
    sim.SetSecondaryModel(nsecondaries, cone);

    spdlog::info("Shooting {} rays (E=[{:.2g},{:.2g}] GeV, σ={:.4g}, radius={:.1f} cm, xsec={:.2g} "
                 "cm^2)...",
                 nrays, emin, emax, sigma_theta, beam_radius, xsec);

    std::vector<NuGeom::InteractionEvent> events;
    events.reserve(nrays);
    std::map<std::string, size_t> per_material;
    for(size_t i = 0; i < nrays; ++i) {
        auto [energy, ray] = gen.GetRay();
        auto event = sim.ShootRay(energy, ray);
        if(event.interacted) ++per_material[event.vertex_material];
        events.push_back(std::move(event));
    }

    auto stats = NuGeom::RayInteractionSim::Summarize(events);
    spdlog::info("Interactions: {}/{} rays ({:.1f}%), {} outgoing particles tracked",
                 stats.ninteracted, stats.nrays, 100.0 * stats.interaction_fraction(),
                 stats.nsecondaries);
    for(const auto &[mat, count] : per_material)
        spdlog::info("  {:>8} interactions in {}", count, mat);

    write_event_log(logfile, events);
    write_svg(svgfile, world, events, h_axis, v_axis, max_depth);

    return 0;
}
