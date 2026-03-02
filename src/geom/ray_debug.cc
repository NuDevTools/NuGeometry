/// ray_debug — CLI diagnostic tool for invalid GetLineSegments output.
///
/// Reads a ray log file (same format written by BeamViz), replays each ray
/// through a given GDML geometry, and reports per-segment anomalies:
///
///   INF          — segment length is infinite (Shape::Intersect missed)
///   NAN          — segment length is NaN
///   TINY         — length below --tiny threshold (default 1e-6 cm)
///   DUP_MAT      — same material as the immediately preceding segment
///   DISCONTINUOUS — start of this segment does not match end of previous (> 1e-4 cm gap)
///
/// --diagnose adds four targeted root-cause tests for each bad segment:
///
///   [D1-WorldBound]   Signed-distance and Intersect2 of the world shape at
///                     the bad start point.  Determines whether the ray has
///                     already left the world.
///   [D2-AABBCheck]    Ray–AABB slab test against every volume in the world
///                     hierarchy.  Volumes whose AABB the ray crosses but
///                     whose material never appears in the remaining segments
///                     are flagged as possible BVH or Shape::Intersect misses.
///   [D3-IntersectProbe] For each AABB-positive volume, calls
///                     Shape::Intersect2 in the world shape's local frame to
///                     distinguish "AABB conservative miss" from a genuine
///                     Shape::Intersect failure (returns +∞ even though the
///                     ray entered the bounding volume).
///   [D4-EpsAccum]     Checks whether the immediately preceding segment was
///                     near-zero-length (≤ 1e-4 cm), which suggests the eps
///                     nudge in Volume.cc skipped a thin boundary.
///   [D5-DirNorm]      Reports |d| − 1; a non-unit direction makes every
///                     time-to-distance conversion wrong.
///   [D6-Perturb]      Re-runs GetLineSegments with the bad-start-point
///                     origin nudged ±1e-5 cm along and perpendicular to the
///                     ray to check numerical sensitivity.
///
/// Exit code: 0 if all rays are clean, 1 if any anomaly was found.
///
/// Usage:
///   ray_debug -g geom.gdml -r rays.log --bad-only --diagnose
///   ray_debug -g geom.gdml -r rays.log -i 42 --diagnose -v

#include "CLI/CLI.hpp"
#include "fmt/format.h"
#include "geom/BoundingBox.hh"
#include "geom/LineSegment.hh"
#include "geom/Parser.hh"
#include "geom/Ray.hh"
#include "geom/World.hh"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Ray log I/O (same format as BeamViz: ox oy oz dx dy dz pot per line)
// ---------------------------------------------------------------------------
static std::vector<NuGeom::Ray> load_rays(const std::string &path) {
    std::ifstream f(path);
    if(!f) throw std::runtime_error("Cannot open ray log: " + path);
    std::vector<NuGeom::Ray> out;
    double ox, oy, oz, dx, dy, dz, pot;
    size_t line = 0;
    while(f >> ox >> oy >> oz >> dx >> dy >> dz >> pot) {
        ++line;
        out.emplace_back(NuGeom::Vector3D{ox, oy, oz}, NuGeom::Vector3D{dx, dy, dz}, pot,
                         /*normalize=*/false);
    }
    if(!f.eof() && f.fail())
        throw std::runtime_error(fmt::format("Parse error in {} at line {}", path, line + 1));
    return out;
}

// ---------------------------------------------------------------------------
// Ray–AABB slab intersection
// Returns true if the ray hits the box in the forward direction.
// t_in / t_out are the signed entry/exit distances (t_in may be negative if
// the ray origin is already inside).
// ---------------------------------------------------------------------------
static bool ray_aabb(const NuGeom::Ray &ray, const NuGeom::BoundingBox &bb, double &t_in,
                     double &t_out) {
    const auto &o = ray.Origin();
    const auto &d = ray.Direction();
    const double ox[3] = {o.X(), o.Y(), o.Z()};
    const double dx[3] = {d.X(), d.Y(), d.Z()};
    const double mn[3] = {bb.min.X(), bb.min.Y(), bb.min.Z()};
    const double mx[3] = {bb.max.X(), bb.max.Y(), bb.max.Z()};

    double tmin = -std::numeric_limits<double>::infinity();
    double tmax = std::numeric_limits<double>::infinity();
    for(int i = 0; i < 3; ++i) {
        if(std::abs(dx[i]) < 1e-15) {
            if(ox[i] < mn[i] || ox[i] > mx[i]) return false;
        } else {
            double t1 = (mn[i] - ox[i]) / dx[i];
            double t2 = (mx[i] - ox[i]) / dx[i];
            if(t1 > t2) std::swap(t1, t2);
            tmin = std::max(tmin, t1);
            tmax = std::min(tmax, t2);
            if(tmin > tmax) return false;
        }
    }
    t_in = tmin;
    t_out = tmax;
    return tmax >= 0.0; // hit in the forward half-space
}

// ---------------------------------------------------------------------------
// Per-segment anomaly flags
// ---------------------------------------------------------------------------
struct SegDiag {
    bool inf_len{false};
    bool nan_len{false};
    bool tiny_len{false};
    bool dup_mat{false};
    bool discontinuous{false};

    bool any() const { return inf_len || nan_len || tiny_len || dup_mat || discontinuous; }

    std::string flags_str() const {
        std::string s;
        if(inf_len) s += " INF";
        if(nan_len) s += " NAN";
        if(tiny_len) s += " TINY";
        if(dup_mat) s += " DUP_MAT";
        if(discontinuous) s += " DISCONTINUOUS";
        return s;
    }
};

struct RayDiag {
    size_t index{0};
    bool in_world{false};
    double dir_norm{0.0};
    std::vector<NuGeom::LineSegment> segs;
    std::vector<SegDiag> seg_diags;
    bool has_exception{false};
    std::string exception_msg;
    bool early_exit{false}; // ray's last segment ends inside the world

    bool any_bad() const {
        if(has_exception) return true;
        if(early_exit) return true;
        for(const auto &d : seg_diags)
            if(d.any()) return true;
        return false;
    }

    // Returns index of the first bad segment, or npos.
    size_t first_bad_seg() const {
        for(size_t i = 0; i < seg_diags.size(); ++i)
            if(seg_diags[i].any()) return i;
        return std::string::npos;
    }

    std::string summary_flags() const {
        std::string s;
        if(early_exit) s += " EARLY_EXIT";
        for(size_t si = 0; si < seg_diags.size(); ++si)
            if(seg_diags[si].any()) s += fmt::format(" seg[{}]:{}", si, seg_diags[si].flags_str());
        return s;
    }
};

static RayDiag run_ray(size_t index, const NuGeom::Ray &ray, NuGeom::World &world,
                       double tiny_threshold) {
    RayDiag diag;
    diag.index = index;
    diag.in_world = world.InWorld(ray.Origin());
    diag.dir_norm = ray.Direction().Norm();

    try {
        diag.segs = world.GetLineSegments(ray);
    } catch(const std::exception &e) {
        diag.has_exception = true;
        diag.exception_msg = e.what();
        return diag;
    }

    for(size_t si = 0; si < diag.segs.size(); ++si) {
        const auto &seg = diag.segs[si];
        double len = seg.Length();
        SegDiag sd;
        sd.inf_len = std::isinf(len);
        sd.nan_len = std::isnan(len);
        sd.tiny_len = std::isfinite(len) && (len < tiny_threshold);
        sd.dup_mat = si > 0 && seg.GetMaterial().Name() == diag.segs[si - 1].GetMaterial().Name();
        if(si > 0) {
            auto gap = (diag.segs[si - 1].End() - seg.Start()).Norm();
            sd.discontinuous = gap > 1e-4;
        }
        diag.seg_diags.push_back(sd);
    }
    // EARLY_EXIT: ray's last endpoint is still well inside the world.
    // PruneSegments removes tiny (< 1e-4 cm) segments, so the last endpoint
    // may be up to ~1e-4 cm inside the world boundary.  Use a tolerance to
    // avoid false positives from the pruning artefact.
    if(!diag.segs.empty()) {
        const auto &last = diag.segs.back();
        if(std::isfinite(last.Length())) {
            double sd = world.GetShape(0)->SignedDistance(last.End());
            diag.early_exit = sd < -2e-4; // well inside (past prune tolerance)
        }
    }
    return diag;
}

// ---------------------------------------------------------------------------
// Root-cause diagnostic tests
// ---------------------------------------------------------------------------

// Helper: attempt GetLineSegments from a given start point, return whether all
// segments are finite.
static bool probe_segments(NuGeom::World &world, const NuGeom::Vector3D &origin,
                           const NuGeom::Vector3D &dir, double pot, size_t *out_count = nullptr) {
    try {
        NuGeom::Ray probe(origin, dir, pot, /*normalize=*/false);
        auto segs = world.GetLineSegments(probe);
        if(out_count) *out_count = segs.size();
        for(const auto &s : segs)
            if(!std::isfinite(s.Length())) return false;
        return true;
    } catch(...) { return false; }
}

static void run_diagnostics(const RayDiag &diag, const NuGeom::Ray &ray, NuGeom::World &world) {
    size_t bad_idx = diag.first_bad_seg();
    if(bad_idx == std::string::npos) return;

    const NuGeom::Vector3D p = diag.segs[bad_idx].Start();
    const NuGeom::Vector3D dir = ray.Direction();
    const NuGeom::Ray ray_p(p, dir, ray.POT(), /*normalize=*/false);

    fmt::print("  ┌─ Root-cause diagnostics for first bad segment [{}] ─────────────────\n",
               bad_idx);
    fmt::print("  │  Bad-segment start: ({:+.6g}, {:+.6g}, {:+.6g})\n", p.X(), p.Y(), p.Z());

    // ------------------------------------------------------------------
    // D1 — World boundary: signed distance and Intersect2 at bad point
    // ------------------------------------------------------------------
    fmt::print("  ├─[D1-WorldBound]\n");
    {
        NuGeom::Shape *ws = world.GetShape(0);
        double sdf = ws->SignedDistance(p);
        fmt::print("  │    world.SignedDistance(p) = {:+.6g} cm  → {}\n", sdf,
                   sdf > 1e-6    ? "OUTSIDE WORLD — ray has already exited"
                   : sdf < -1e-6 ? "inside world"
                                 : "on world surface (within 1e-6 cm)");

        auto [te, tx] = ws->Intersect2(ray_p);
        fmt::print("  │    Intersect2 from p:  t_enter={:+.6g}  t_exit={:+.6g}\n", te, tx);
        if(tx < 0)
            fmt::print("  │    => Both entry/exit negative: ray is past the world. "
                       "Prior segment should have stopped at world exit.\n");
        else if(tx < 1e-4)
            fmt::print("  │    => t_exit={:.4g} cm — ray is grazing the world boundary.\n", tx);
        else if(te < 0)
            fmt::print("  │    => Ray is inside the world with {:.4g} cm remaining.\n", tx);
    }

    // ------------------------------------------------------------------
    // D2 — AABB cross-check: which volumes does the ray hit from bad point?
    //      Compares against materials present in segments[bad_idx..].
    // ------------------------------------------------------------------
    fmt::print("  ├─[D2-AABBCheck]\n");
    {
        // Collect materials already produced at or after the bad segment.
        std::set<std::string> covered_mats;
        for(size_t i = bad_idx; i < diag.segs.size(); ++i)
            covered_mats.insert(diag.segs[i].GetMaterial().Name());

        // Use a generous depth to capture nested volumes.
        auto bounds = world.GetVolumeBounds(6);
        fmt::print("  │    Checking {} world-frame AABBs (depth ≤ 6):\n", bounds.size());

        struct AabbHit {
            double t_in, t_out;
            std::string mat;
        };
        std::vector<AabbHit> hits;
        for(const auto &vb : bounds) {
            double ti, to;
            if(ray_aabb(ray_p, vb.bb, ti, to) && to > 0)
                hits.push_back({ti, to, vb.material.Name()});
        }
        std::sort(hits.begin(), hits.end(),
                  [](const AabbHit &a, const AabbHit &b) { return a.t_in < b.t_in; });

        if(hits.empty()) {
            fmt::print("  │    No AABBs hit — ray has exited all known volumes.\n");
            fmt::print("  │    Likely cause: prior boundary was not correctly detected,\n");
            fmt::print("  │    so the segment that should have ended there became infinite.\n");
        } else {
            for(const auto &h : hits) {
                bool missing = covered_mats.find(h.mat) == covered_mats.end();
                fmt::print("  │      t=[{:+.4g},{:+.4g}]  {:<24s}{}\n", h.t_in, h.t_out, h.mat,
                           missing ? " *** NOT in remaining segments — possible BVH/Shape miss ***"
                                   : "");
            }
        }
    }

    // ------------------------------------------------------------------
    // D3 — Shape::Intersect probe: for volumes whose AABB the ray crosses,
    //      call the world-level Shape::Intersect2 to see if the world shape
    //      itself reports a finite intersection (determines if the issue is
    //      in a daughter's Intersect, or in the world-boundary detection).
    // ------------------------------------------------------------------
    fmt::print("  ├─[D3-IntersectProbe]\n");
    {
        // World shape is in world-local frame (same as world frame).
        NuGeom::Shape *ws = world.GetShape(0);
        auto [te, tx] = ws->Intersect2(ray_p);
        fmt::print("  │    World shape Intersect2: t_enter={:+.6g}  t_exit={:+.6g}\n", te, tx);
        if(!std::isfinite(tx))
            fmt::print("  │    => World shape itself returns infinite exit — world shape "
                       "Intersect2 bug or ray truly misses world.\n");
        else
            fmt::print("  │    => World shape is finite — issue is in a daughter volume.\n");

        // Test each direct daughter (idx 1..N) — their shapes are in their own
        // local frames, so we can only report the raw result as a relative check.
        size_t ndau = world.NDaughters();
        if(ndau > 0) {
            fmt::print("  │    Direct daughter shape Intersect (in daughter-local frame):\n");
            fmt::print("  │    Note: these use daughter-local coordinates — "
                       "finite result means a miss in local frame.\n");
            for(size_t i = 1; i <= ndau; ++i) {
                NuGeom::Shape *sh = world.GetShape(i);
                // The ray_p is in world frame; the daughter shape is in its own
                // local frame.  We can still call Intersect to check for ±∞.
                double t = sh->Intersect(ray_p);
                auto [ti2, to2] = sh->Intersect2(ray_p);
                std::string mat = world.GetMaterial(i).Name();
                fmt::print("  │      daughter[{}] ({:<20s}): Intersect={:.6g}  "
                           "Intersect2=({:.6g},{:.6g})\n",
                           i, mat, t, ti2, to2);
            }
        }
    }

    // ------------------------------------------------------------------
    // D4 — Eps accumulation: was the previous segment suspiciously short?
    // ------------------------------------------------------------------
    fmt::print("  ├─[D4-EpsAccum]\n");
    {
        static constexpr double EPS_THRESHOLD =
            1e-4; // cm — matches Volume.cc eps=1e-8 * many calls
        if(bad_idx == 0) {
            fmt::print("  │    Bad segment is the first segment — no predecessor to check.\n");
        } else {
            double prev_len = diag.segs[bad_idx - 1].Length();
            if(std::isfinite(prev_len) && prev_len < EPS_THRESHOLD)
                fmt::print(
                    "  │    *** Previous segment length = {:.4g} cm (< {:.0e} cm threshold)\n"
                    "  │        Likely eps-nudge artifact: a thin volume boundary was\n"
                    "  │        grazed and the 1e-8 cm nudge skipped past it, leaving\n"
                    "  │        the ray in an indeterminate state for the next Intersect.\n",
                    prev_len, EPS_THRESHOLD);
            else
                fmt::print("  │    Previous segment length = {:.6g} cm — no eps accumulation.\n",
                           prev_len);
        }
    }

    // ------------------------------------------------------------------
    // D5 — Direction normalization
    // ------------------------------------------------------------------
    fmt::print("  ├─[D5-DirNorm]\n");
    {
        double err = std::abs(diag.dir_norm - 1.0);
        if(err > 1e-6)
            fmt::print("  │    *** |d| - 1 = {:.6g}  Direction is NOT unit-length.\n"
                       "  │        Every t*direction distance will be scaled by {:.6g}×,\n"
                       "  │        corrupting all boundary-crossing times.\n",
                       err, diag.dir_norm);
        else
            fmt::print("  │    |d| - 1 = {:.2e}  (ok)\n", err);
    }

    // ------------------------------------------------------------------
    // D6 — Perturbation: re-run GetLineSegments from the bad start point
    //      with small nudges to gauge numerical sensitivity.
    // ------------------------------------------------------------------
    fmt::print("  └─[D6-Perturb]\n");
    {
        // Build a perpendicular vector for lateral nudges.
        const NuGeom::Vector3D ex{1, 0, 0}, ey{0, 1, 0};
        NuGeom::Vector3D perp = (std::abs(dir.X()) < 0.9 ? ex : ey);
        // perp = perp - (perp·dir)*dir
        double dot = perp.X() * dir.X() + perp.Y() * dir.Y() + perp.Z() * dir.Z();
        perp = NuGeom::Vector3D{perp.X() - dot * dir.X(), perp.Y() - dot * dir.Y(),
                                perp.Z() - dot * dir.Z()};
        double pn = perp.Norm();
        if(pn > 1e-12) perp = NuGeom::Vector3D{perp.X() / pn, perp.Y() / pn, perp.Z() / pn};

        constexpr double dz = 1e-5; // cm along ray
        constexpr double dr = 1e-5; // cm perpendicular
        struct Probe {
            std::string label;
            NuGeom::Vector3D origin;
        };
        std::vector<Probe> probes = {
            {"forward  +1e-5 cm",
             NuGeom::Vector3D{p.X() + dz * dir.X(), p.Y() + dz * dir.Y(), p.Z() + dz * dir.Z()}},
            {"backward -1e-5 cm",
             NuGeom::Vector3D{p.X() - dz * dir.X(), p.Y() - dz * dir.Y(), p.Z() - dz * dir.Z()}},
            {"lateral  +1e-5 cm",
             NuGeom::Vector3D{p.X() + dr * perp.X(), p.Y() + dr * perp.Y(), p.Z() + dr * perp.Z()}},
        };

        for(const auto &probe : probes) {
            size_t nseg = 0;
            bool ok = probe_segments(world, probe.origin, dir, ray.POT(), &nseg);
            fmt::print("       nudge {:s}: {}  ({} segs)\n", probe.label,
                       ok ? "clean" : "STILL BAD", nseg);
        }
        fmt::print("       If 'forward' is clean but 'backward' is bad (or vice versa),\n"
                   "       the ray is sitting exactly on a surface — eps-sensitivity.\n"
                   "       If all probes are bad, the issue is structural (transform or\n"
                   "       Intersect algorithm), not a boundary grazing artifact.\n");
    }
}

// ---------------------------------------------------------------------------
// Output formatting
// ---------------------------------------------------------------------------
static void print_ray(const RayDiag &diag, const NuGeom::Ray &ray, bool do_diagnose,
                      NuGeom::World &world) {
    const auto &o = ray.Origin(), &d = ray.Direction();
    bool bad = diag.any_bad();

    fmt::print("=== Ray {:>5d}  {} ===\n", diag.index, bad ? "[BAD]" : "[ok]");
    fmt::print("  Origin:    ({:+.10g}, {:+.10g}, {:+.10g})\n", o.X(), o.Y(), o.Z());
    fmt::print("  Direction: ({:+.10g}, {:+.10g}, {:+.10g})  |d|={:.10g}\n", d.X(), d.Y(), d.Z(),
               diag.dir_norm);
    fmt::print("  In world:  {}\n", diag.in_world ? "yes" : "NO  <--- origin outside world");

    if(std::abs(diag.dir_norm - 1.0) > 1e-6)
        fmt::print("  WARNING: direction not unit-length (|d|={:.6g})\n", diag.dir_norm);

    if(diag.has_exception) {
        fmt::print("  EXCEPTION: {}\n\n", diag.exception_msg);
        return;
    }

    fmt::print("  Segments:  {}\n", diag.segs.size());
    fmt::print("  {:>4s}  {:<22s}  {:>14s}  {:>34s}  {:>34s}  {}\n", "idx", "material",
               "length (cm)", "start (x,y,z)", "end (x,y,z)", "flags");
    fmt::print("  {}\n", std::string(130, '-'));

    for(size_t si = 0; si < diag.segs.size(); ++si) {
        const auto &seg = diag.segs[si];
        const auto &sd = diag.seg_diags[si];
        double len = seg.Length();
        const auto &s = seg.Start(), &e = seg.End();

        std::string len_str;
        if(sd.inf_len)
            len_str = "           inf";
        else if(sd.nan_len)
            len_str = "           nan";
        else
            len_str = fmt::format("{:>14.6f}", len);

        std::string start_str = fmt::format("({:+.4g},{:+.4g},{:+.4g})", s.X(), s.Y(), s.Z());
        std::string end_str = (sd.inf_len || sd.nan_len)
                                  ? "???"
                                  : fmt::format("({:+.4g},{:+.4g},{:+.4g})", e.X(), e.Y(), e.Z());
        std::string flags = sd.flags_str();

        fmt::print("  [{:>3d}] {:<22s}  {}  {:>34s}  {:>34s}  {}\n", si, seg.GetMaterial().Name(),
                   len_str, start_str, end_str, flags.empty() ? "" : "<<< " + flags + " >>>");

        if(sd.discontinuous && si > 0) {
            const auto &prev_end = diag.segs[si - 1].End();
            auto gap = (prev_end - s).Norm();
            fmt::print("       ^ gap from previous end ({:+.4g},{:+.4g},{:+.4g}): {:.4g} cm\n",
                       prev_end.X(), prev_end.Y(), prev_end.Z(), gap);
        }
    }

    if(diag.early_exit)
        fmt::print("  *** EARLY_EXIT: last segment ends inside the world — traversal "
                   "terminated before reaching the world boundary.\n"
                   "  *** Likely cause: physical mothers not set (parser omission), so\n"
                   "  *** PhysicalVolume::GetLineSegments cannot unwind up the hierarchy.\n");

    if(do_diagnose && bad) run_diagnostics(diag, ray, world);

    fmt::print("\n");
}

static void print_summary_line(const RayDiag &diag) {
    if(!diag.any_bad()) {
        fmt::print("  {:>5d}  ok     {} segs\n", diag.index, diag.segs.size());
        return;
    }
    std::string issues;
    if(diag.has_exception) {
        issues = "EXCEPTION: " + diag.exception_msg;
    } else {
        issues = diag.summary_flags();
    }
    fmt::print("  {:>5d}  BAD    {} segs  — {}\n", diag.index, diag.segs.size(), issues);
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char **argv) {
    auto console = spdlog::stdout_color_mt("NuGeom");
    spdlog::set_default_logger(console);
    spdlog::set_pattern("[%n] [%^%l%$] %v");
    spdlog::set_level(spdlog::level::warn);

    CLI::App app("ray_debug: diagnose invalid ray segments from a BeamViz ray log");
    argv = app.ensure_utf8(argv);

    std::string geomfile, rayfile = "rays.log";
    int verbosity = 0;
    bool bad_only = false;
    bool summary = false;
    bool do_diag = false;
    int ray_index = -1;
    double tiny_thr = 1e-6;

    app.add_option("-g,--geometry", geomfile, "GDML geometry file")->required();
    app.add_option("-r,--rays", rayfile, "Ray log file (default: rays.log)");
    app.add_flag("-v,--verbose", verbosity,
                 "Verbosity: -v=geometry debug, -vv=trace (enables spdlog inside traversal)");
    app.add_flag("--bad-only", bad_only, "Only print rays that have anomalies");
    app.add_flag("--summary", summary, "One-line-per-ray summary table");
    app.add_flag("--diagnose", do_diag,
                 "Run D1-D6 root-cause tests on the first bad segment of each bad ray");
    app.add_option("-i,--index", ray_index, "Only replay one ray (0-based index)");
    app.add_option("--tiny", tiny_thr, "TINY flag threshold in cm (default 1e-6)");

    try {
        app.parse(argc, argv);
    } catch(const CLI::ParseError &e) { return app.exit(e); }

    if(verbosity == 1)
        spdlog::set_level(spdlog::level::debug);
    else if(verbosity >= 2)
        spdlog::set_level(spdlog::level::trace);

    fmt::print("Geometry : {}\n", geomfile);
    NuGeom::GDMLParser parser(geomfile);
    NuGeom::World world = parser.GetWorld();
    auto bbox = world.GetWorldBox();
    fmt::print("World box: ({:.2f},{:.2f},{:.2f}) -> ({:.2f},{:.2f},{:.2f}) cm\n\n", bbox.min.X(),
               bbox.min.Y(), bbox.min.Z(), bbox.max.X(), bbox.max.Y(), bbox.max.Z());

    fmt::print("Ray log  : {}\n", rayfile);
    std::vector<NuGeom::Ray> rays;
    try {
        rays = load_rays(rayfile);
    } catch(const std::exception &e) {
        fmt::print(stderr, "ERROR: {}\n", e.what());
        return 1;
    }
    fmt::print("Rays     : {}\n\n", rays.size());

    size_t n_bad = 0;
    if(summary) {
        fmt::print("  {:>5s}  {:<7s}  {}\n", "index", "status", "detail");
        fmt::print("  {}\n", std::string(70, '-'));
    }

    for(size_t ri = 0; ri < rays.size(); ++ri) {
        if(ray_index >= 0 && static_cast<int>(ri) != ray_index) continue;

        auto diag = run_ray(ri, rays[ri], world, tiny_thr);
        if(diag.any_bad()) ++n_bad;
        if(bad_only && !diag.any_bad()) continue;

        if(summary)
            print_summary_line(diag);
        else
            print_ray(diag, rays[ri], do_diag, world);
    }

    if(summary) fmt::print("\n");
    fmt::print("Result: {}/{} rays with anomalous segments.\n", n_bad, rays.size());
    return n_bad > 0 ? 1 : 0;
}
