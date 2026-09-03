#include "catch2/catch.hpp"

#include "geom/LineSegment.hh"
#include "geom/Parser.hh"
#include "geom/Random.hh"
#include "geom/Ray.hh"
#include "geom/World.hh"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>

namespace {
// Compare two segment lists for equivalence: same material sequence and
// matching lengths (within tolerance).  Returns "" on match, else a message.
std::string Compare(const std::vector<NuGeom::LineSegment> &a,
                    const std::vector<NuGeom::LineSegment> &b) {
    if(a.size() != b.size())
        return "size " + std::to_string(a.size()) + " vs " + std::to_string(b.size());
    for(size_t i = 0; i < a.size(); ++i) {
        if(a[i].GetMaterial().Name() != b[i].GetMaterial().Name())
            return "material[" + std::to_string(i) + "] " + a[i].GetMaterial().Name() + " vs " +
                   b[i].GetMaterial().Name();
        const double tol = std::max(1e-3, 1e-4 * a[i].Length());
        if(std::abs(a[i].Length() - b[i].Length()) > tol)
            return "length[" + std::to_string(i) + "] " + std::to_string(a[i].Length()) + " vs " +
                   std::to_string(b[i].Length());
    }
    return "";
}
} // namespace

TEST_CASE("GetLineSegments sequential vs sweep", "[.][benchmark]") {
    std::string gdml = std::string(NUGEOM_SOURCE_DIR) + "/nd_hall_with_lar_tms_nosand.gdml";
    NuGeom::GDMLParser parse(gdml);
    auto world = parse.GetWorld();
    auto box = world.GetWorldBox();

    auto vbounds = world.GetVolumeBounds(4);
    double xlo = 0, xhi = 0, ylo = 0, yhi = 0;
    for(const auto &vb : vbounds) {
        if(vb.name == "volArgonCubeDetector") {
            xlo = vb.bb.min.X();
            xhi = vb.bb.max.X();
            ylo = vb.bb.min.Y();
            yhi = vb.bb.max.Y();
            break;
        }
    }

    const double z_start = box.min.Z() + 1.0;
    constexpr size_t nrays = 500;
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist_x(xlo, xhi);
    std::uniform_real_distribution<double> dist_y(ylo, yhi);
    std::uniform_real_distribution<double> dist_ang(-0.15, 0.15);

    std::vector<NuGeom::Ray> rays;
    rays.reserve(nrays);
    for(size_t i = 0; i < nrays; ++i) {
        NuGeom::Vector3D origin(dist_x(rng), dist_y(rng), z_start);
        // Mix of axis-aligned and slightly angled rays (the latter exercise the
        // rotation transforms and mimic angled flux rays).
        NuGeom::Vector3D dir = (i % 2 == 0) ? NuGeom::Vector3D(0, 0, 1)
                                            : NuGeom::Vector3D(dist_ang(rng), dist_ang(rng), 1);
        rays.emplace_back(origin, dir, 1.0);
    }

    // 1. Cross-validate the sweep against the sequential oracle on every ray.
    size_t mismatches = 0, mism_axis = 0, mism_angled = 0;
    std::string first_msg;
    for(size_t i = 0; i < rays.size(); ++i) {
        auto seq = world.GetLineSegmentsSequential(rays[i]);
        auto swp = world.GetLineSegmentsSweep(rays[i]);
        std::string msg = Compare(seq, swp);
        if(!msg.empty()) {
            if(mismatches == 0) {
                first_msg = "ray " + std::to_string(i) + ": " + msg + "\n    seq:";
                for(const auto &s : seq) first_msg += " " + s.GetMaterial().Name();
                first_msg += "\n    swp:";
                for(const auto &s : swp) first_msg += " " + s.GetMaterial().Name();
            }
            ++mismatches;
            (i % 2 == 0) ? ++mism_axis : ++mism_angled;
        }
    }
    std::cout << "  Cross-validation: " << (rays.size() - mismatches) << "/" << rays.size()
              << " rays match  (" << mism_axis << " axis-aligned, " << mism_angled
              << " angled mismatches)\n";
    if(mismatches) std::cout << "  first mismatch -> " << first_msg << "\n";

    // 2. Time both paths.
    auto timeit = [&](const char *label, auto fn) {
        constexpr size_t nsamples = 5;
        double best_ms = 1e30;
        size_t total = 0;
        for(size_t s = 0; s < nsamples; ++s) {
            auto t0 = std::chrono::high_resolution_clock::now();
            total = 0;
            for(const auto &ray : rays) total += fn(ray).size();
            auto t1 = std::chrono::high_resolution_clock::now();
            best_ms = std::min(best_ms, std::chrono::duration<double, std::milli>(t1 - t0).count());
        }
        std::cout << "  " << label << ": best " << best_ms << " ms over " << nrays << " rays  ("
                  << static_cast<double>(nrays) / (best_ms / 1000.0) << " rays/s, " << total
                  << " segs)\n";
    };
    timeit("sequential", [&](const NuGeom::Ray &r) { return world.GetLineSegmentsSequential(r); });
    timeit("sweep     ", [&](const NuGeom::Ray &r) { return world.GetLineSegmentsSweep(r); });

    CHECK(mismatches == 0);
}

// Emit a shared set of rays and our SEQUENTIAL (production) segments for each,
// for comparison against ROOT's TGeoNavigator (tools/check_overlaps_root.C).
// Writes tools/compare_rays.txt (ox oy oz dx dy dz, unit dir) and
// tools/nugeom_segments.txt ("rayIndex  len0:mat0 len1:mat1 ...").
TEST_CASE("Navigator comparison: emit rays + sequential segments", "[.][navcompare]") {
    const std::string root = std::string(NUGEOM_SOURCE_DIR);
    NuGeom::GDMLParser parse(root + "/nd_hall_with_lar_tms_nosand.gdml");
    auto world = parse.GetWorld();
    auto box = world.GetWorldBox();
    auto vbounds = world.GetVolumeBounds(4);
    double xlo = 0, xhi = 0, ylo = 0, yhi = 0;
    for(const auto &vb : vbounds)
        if(vb.name == "volArgonCubeDetector") {
            xlo = vb.bb.min.X();
            xhi = vb.bb.max.X();
            ylo = vb.bb.min.Y();
            yhi = vb.bb.max.Y();
            break;
        }
    const double z0 = box.min.Z() + 1.0;

    // Deterministic ray set: a grid of +z rays through the ArgonCube, plus a few
    // angled ones.  Build (origin, unit-dir) pairs.
    std::vector<std::array<double, 6>> rays;
    const double fx[3] = {0.35, 0.5, 0.65}, fy[3] = {0.35, 0.5, 0.65};
    for(double a : fx)
        for(double b : fy)
            rays.push_back({xlo + a * (xhi - xlo), ylo + b * (yhi - ylo), z0, 0, 0, 1});
    rays.push_back({0.5 * (xlo + xhi), 0.5 * (ylo + yhi), z0, 0.1, 0.05, 1});
    rays.push_back({0.45 * (xlo + xhi) + 0.5 * xlo, 0.5 * (ylo + yhi), z0, -0.08, 0.12, 1});
    rays.push_back({xlo + 0.5 * (xhi - xlo), ylo + 0.4 * (yhi - ylo), z0, 0.03, -0.1, 1});

    std::ofstream fr(root + "/tools/compare_rays.txt");
    std::ofstream fs(root + "/tools/nugeom_segments.txt");
    fr << std::setprecision(12);
    fs << std::setprecision(10);
    for(size_t i = 0; i < rays.size(); ++i) {
        auto &r = rays[i];
        NuGeom::Vector3D o(r[0], r[1], r[2]), d(r[3], r[4], r[5]);
        NuGeom::Ray ray(o, d, 1.0); // ctor normalizes the direction
        const auto &u = ray.Direction();
        fr << o.X() << " " << o.Y() << " " << o.Z() << " " << u.X() << " " << u.Y() << " " << u.Z()
           << "\n";
        auto segs = world.GetLineSegmentsSequential(ray);
        fs << i;
        for(const auto &s : segs) fs << "  " << s.Length() << ":" << s.GetMaterial().Name();
        fs << "\n";
    }
    std::cerr << "  Wrote " << rays.size() << " rays + sequential segments to tools/\n";
    REQUIRE(true);
}

// For ray 7 (which the navigator comparison flags), walk our sequential
// segments and, inside each large G10 span, query FindMaterial (our geometry's
// point-containment) -- if it returns LAr, the geometry is right and the
// sequential traversal is mis-assigning material.
TEST_CASE("Navigator comparison: ray 7 traversal-vs-containment probe", "[.][nav7probe]") {
    const std::string root = std::string(NUGEOM_SOURCE_DIR);
    NuGeom::GDMLParser parse(root + "/nd_hall_with_lar_tms_nosand.gdml");
    auto world = parse.GetWorld();
    NuGeom::Vector3D o(121.0117, -173.973, -29999), d(0, 0, 1);
    NuGeom::Ray ray(o, d, 1.0);
    auto segs = world.GetLineSegmentsSequential(ray);
    double z = 0; // cumulative along ray (dir is +z so z param == path length)
    std::cerr << "  ray 7: probing FindMaterial inside large G10 spans (traversal said G10)\n";
    for(const auto &s : segs) {
        if(s.GetMaterial().Name() == "G10" && s.Length() > 5.0) {
            const auto mid = s.Start() + 0.5 * (s.End() - s.Start());
            std::cerr << "    G10 span z=[" << s.Start().Z() << ", " << s.End().Z()
                      << "] len=" << s.Length() << "  -> FindMaterial(mid)='"
                      << world.FindMaterial(mid).Name() << "' FindVolume='" << world.FindVolume(mid)
                      << "'\n";
        }
        (void)z;
    }
    REQUIRE(true);
}

// Scan rays and report the distinct sweep/containment disagreement "signatures"
// (sweep-material vs contained-material), which localize suspect overlapping /
// non-nested volumes in the GDML.  Run explicitly: testsuite "[overlapscan]".
TEST_CASE("Sweep overlap scan", "[.][overlapscan]") {
    std::string gdml = std::string(NUGEOM_SOURCE_DIR) + "/nd_hall_with_lar_tms_nosand.gdml";
    NuGeom::GDMLParser parse(gdml);
    auto world = parse.GetWorld();
    auto box = world.GetWorldBox();
    auto vbounds = world.GetVolumeBounds(4);
    double xlo = 0, xhi = 0, ylo = 0, yhi = 0;
    for(const auto &vb : vbounds)
        if(vb.name == "volArgonCubeDetector") {
            xlo = vb.bb.min.X();
            xhi = vb.bb.max.X();
            ylo = vb.bb.min.Y();
            yhi = vb.bb.max.Y();
            break;
        }
    const double z_start = box.min.Z() + 1.0;
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist_x(xlo, xhi);
    std::uniform_real_distribution<double> dist_y(ylo, yhi);
    std::uniform_real_distribution<double> dist_ang(-0.15, 0.15);

    struct Sig {
        size_t count = 0;
        double ex_z0 = 0, ex_z1 = 0;
    };
    std::map<std::string, Sig> sigs;
    size_t rays_with_overlap = 0;
    constexpr size_t nrays = 2000;
    for(size_t i = 0; i < nrays; ++i) {
        NuGeom::Vector3D origin(dist_x(rng), dist_y(rng), z_start);
        NuGeom::Vector3D dir = (i % 2 == 0) ? NuGeom::Vector3D(0, 0, 1)
                                            : NuGeom::Vector3D(dist_ang(rng), dist_ang(rng), 1);
        NuGeom::Ray ray(origin, dir, 1.0);
        auto regions = world.CheckSweepConsistency(ray, /*warn=*/false);
        if(!regions.empty()) ++rays_with_overlap;
        for(const auto &r : regions) {
            auto &s = sigs[r.sweep_volume + " [" + r.sweep_material + "] (sweep) vs " +
                           r.contained_volume + " [" + r.contained_material + "] (contained)"];
            if(s.count == 0) {
                s.ex_z0 = r.start.Z();
                s.ex_z1 = r.end.Z();
            }
            ++s.count;
        }
    }
    std::cerr << "  Overlap scan: " << rays_with_overlap << "/" << nrays
              << " rays hit a sweep/containment disagreement.  Distinct signatures:\n";
    for(const auto &kv : sigs)
        std::cerr << "    [" << kv.second.count << "x] " << kv.first << "  (e.g. z=["
                  << kv.second.ex_z0 << ", " << kv.second.ex_z1 << "])\n";
    REQUIRE(true);
}

// Emit the rays where the sweep disagrees with the sequential traversal, plus
// BOTH our segment lists, for a 3-way comparison against ROOT's TGeoNavigator
// (confirms the sweep's residual mismatches are the known CSG single-interval
// gap and not a new traversal bug).
TEST_CASE("Sweep mismatch rays for ROOT 3-way compare", "[.][navmismatch]") {
    const std::string root = std::string(NUGEOM_SOURCE_DIR);
    NuGeom::GDMLParser parse(root + "/nd_hall_with_lar_tms_nosand.gdml");
    auto world = parse.GetWorld();
    auto box = world.GetWorldBox();
    auto vb = world.GetVolumeBounds(4);
    double xlo = 0, xhi = 0, ylo = 0, yhi = 0;
    for(const auto &v : vb)
        if(v.name == "volArgonCubeDetector") {
            xlo = v.bb.min.X();
            xhi = v.bb.max.X();
            ylo = v.bb.min.Y();
            yhi = v.bb.max.Y();
            break;
        }
    const double z0 = box.min.Z() + 1.0;
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dx(xlo, xhi), dy(ylo, yhi), da(-0.15, 0.15);

    std::ofstream fr(root + "/tools/mism_rays.txt"), fseq(root + "/tools/mism_seq.txt"),
        fswp(root + "/tools/mism_swp.txt");
    fr << std::setprecision(12);
    fseq << std::setprecision(10);
    fswp << std::setprecision(10);
    int n = 0;
    for(size_t i = 0; i < 500; ++i) {
        NuGeom::Vector3D o(dx(rng), dy(rng), z0);
        NuGeom::Vector3D d =
            (i % 2 == 0) ? NuGeom::Vector3D(0, 0, 1) : NuGeom::Vector3D(da(rng), da(rng), 1);
        NuGeom::Ray ray(o, d, 1.0);
        auto seq = world.GetLineSegmentsSequential(ray);
        auto swp = world.GetLineSegmentsSweep(ray);
        if(Compare(seq, swp).empty()) continue; // only the mismatching rays
        const auto &u = ray.Direction();
        fr << o.X() << " " << o.Y() << " " << o.Z() << " " << u.X() << " " << u.Y() << " " << u.Z()
           << "\n";
        fseq << n;
        for(auto &s : seq) fseq << "  " << s.Length() << ":" << s.GetMaterial().Name();
        fseq << "\n";
        fswp << n;
        for(auto &s : swp) fswp << "  " << s.Length() << ":" << s.GetMaterial().Name();
        fswp << "\n";
        ++n;
    }
    std::cerr << "  Wrote " << n << " mismatching rays + seq/swp segments\n";
    REQUIRE(true);
}

// Emit ALL 500 benchmark rays + sequential segments for a full ROOT agreement
// measurement (broader than the 12 fixed [navcompare] rays).
TEST_CASE("Emit all benchmark rays + sequential for ROOT compare", "[.][navall]") {
    const std::string root = std::string(NUGEOM_SOURCE_DIR);
    NuGeom::GDMLParser parse(root + "/nd_hall_with_lar_tms_nosand.gdml");
    auto world = parse.GetWorld();
    auto box = world.GetWorldBox();
    auto vb = world.GetVolumeBounds(4);
    double xlo = 0, xhi = 0, ylo = 0, yhi = 0;
    for(const auto &v : vb)
        if(v.name == "volArgonCubeDetector") {
            xlo = v.bb.min.X();
            xhi = v.bb.max.X();
            ylo = v.bb.min.Y();
            yhi = v.bb.max.Y();
            break;
        }
    const double z0 = box.min.Z() + 1.0;
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dx(xlo, xhi), dy(ylo, yhi), da(-0.15, 0.15);
    std::ofstream fr(root + "/tools/all_rays.txt"), fs(root + "/tools/all_seq.txt");
    fr << std::setprecision(12);
    fs << std::setprecision(10);
    for(size_t i = 0; i < 500; ++i) {
        NuGeom::Vector3D o(dx(rng), dy(rng), z0);
        NuGeom::Vector3D d =
            (i % 2 == 0) ? NuGeom::Vector3D(0, 0, 1) : NuGeom::Vector3D(da(rng), da(rng), 1);
        NuGeom::Ray ray(o, d, 1.0);
        const auto &u = ray.Direction();
        fr << o.X() << " " << o.Y() << " " << o.Z() << " " << u.X() << " " << u.Y() << " " << u.Z()
           << "\n";
        auto seq = world.GetLineSegmentsSequential(ray);
        fs << i;
        for(auto &s : seq) fs << "  " << s.Length() << ":" << s.GetMaterial().Name();
        fs << "\n";
    }
    std::cerr << "  Wrote 500 rays + sequential segments\n";
}

// Generic ROOT cross-check emitter (used by the CI gate via tools/root_xcheck.py).
// Env: NUGEOM_XCHECK_GDML (geometry), NUGEOM_XCHECK_OUT (output prefix),
//      NUGEOM_XCHECK_VOLUME (optional volume whose transverse bbox to sample over;
//      defaults to the world box), NUGEOM_XCHECK_NRAYS (default 600).
// Writes <prefix>_rays.txt (ox oy oz dx dy dz, unit dir) and <prefix>_seg.txt
// ("rayIndex  len:mat ..."), using the *default* traversal (the boundary-sweep).
// A matching ROOT TGeoNavigator pass + compare_segments.py form the gate.
TEST_CASE("ROOT cross-check emitter", "[.][xcheck]") {
    const char *gdml = std::getenv("NUGEOM_XCHECK_GDML");
    const char *out = std::getenv("NUGEOM_XCHECK_OUT");
    REQUIRE(gdml != nullptr);
    REQUIRE(out != nullptr);
    const char *vol_env = std::getenv("NUGEOM_XCHECK_VOLUME");
    const char *nrays_env = std::getenv("NUGEOM_XCHECK_NRAYS");
    const size_t nrays = nrays_env ? std::strtoul(nrays_env, nullptr, 10) : 600;

    NuGeom::GDMLParser parse(gdml);
    auto world = parse.GetWorld();
    auto box = world.GetWorldBox();

    // Transverse sampling window.  Default: the xy bbox enclosing all placed
    // volumes (the detector "content" envelope) so we sample where the geometry
    // actually is, not the (often vastly larger) world box.  A named volume's
    // bbox can be requested instead.  Angled rays from the world -z face still
    // sweep across the whole geometry, exercising far-off structure too.
    double xlo = box.min.X(), xhi = box.max.X(), ylo = box.min.Y(), yhi = box.max.Y();
    const auto vbounds = world.GetVolumeBounds(8);
    if(vol_env && vol_env[0]) {
        for(const auto &vb : vbounds)
            if(vb.name == vol_env) {
                xlo = vb.bb.min.X();
                xhi = vb.bb.max.X();
                ylo = vb.bb.min.Y();
                yhi = vb.bb.max.Y();
                break;
            }
    } else if(!vbounds.empty()) {
        xlo = ylo = std::numeric_limits<double>::infinity();
        xhi = yhi = -std::numeric_limits<double>::infinity();
        for(const auto &vb : vbounds) {
            xlo = std::min(xlo, vb.bb.min.X());
            xhi = std::max(xhi, vb.bb.max.X());
            ylo = std::min(ylo, vb.bb.min.Y());
            yhi = std::max(yhi, vb.bb.max.Y());
        }
    }
    const double z0 = box.min.Z() + 1.0;

    std::mt19937 rng(1234);
    std::uniform_real_distribution<double> dx(xlo, xhi), dy(ylo, yhi), da(-0.15, 0.15);
    std::ofstream fr(std::string(out) + "_rays.txt"), fs(std::string(out) + "_seg.txt");
    fr << std::setprecision(12);
    fs << std::setprecision(10);
    for(size_t i = 0; i < nrays; ++i) {
        NuGeom::Vector3D o(dx(rng), dy(rng), z0);
        NuGeom::Vector3D d =
            (i % 2 == 0) ? NuGeom::Vector3D(0, 0, 1) : NuGeom::Vector3D(da(rng), da(rng), 1);
        NuGeom::Ray ray(o, d, 1.0);
        const auto &u = ray.Direction();
        fr << o.X() << " " << o.Y() << " " << o.Z() << " " << u.X() << " " << u.Y() << " " << u.Z()
           << "\n";
        auto segs = world.GetLineSegments(ray); // default traversal (sweep)
        fs << i;
        for(const auto &s : segs) fs << "  " << s.Length() << ":" << s.GetMaterial().Name();
        fs << "\n";
    }
    std::cerr << "  [xcheck] wrote " << nrays << " rays for " << gdml << "\n";
}

// Isolated cost of ONE layer-1 traversal on the DUNE ND hall: wall time and
// heap allocations for GetColumnLengths vs GetLineSegmentsSweep.  Run under
// LD_PRELOAD of an operator-new counter to get the allocation figure.
TEST_CASE("Layer-1 traversal cost", "[.][travcost]") {
    std::string gdml = std::string(NUGEOM_SOURCE_DIR) + "/nd_hall_with_lar_tms_nosand.gdml";
    NuGeom::GDMLParser parse(gdml);
    auto world = parse.GetWorld();

    const char *nr = std::getenv("BENCH_NRAYS");
    const size_t nrays = nr ? std::strtoul(nr, nullptr, 10) : 2000;
    std::vector<NuGeom::Ray> rays;
    rays.reserve(nrays);
    NuGeom::Random::Instance().Seed(20260902);
    for(size_t i = 0; i < nrays; ++i) {
        const double x = NuGeom::Random::Instance().Uniform(-350.0, 350.0);
        const double y = NuGeom::Random::Instance().Uniform(-170.0, 170.0);
        rays.emplace_back(NuGeom::Vector3D(x, y, -30000.0), NuGeom::Vector3D(0, 0, 1), 1.0);
    }

    auto timeit = [&](const char *label, auto fn) {
        double best_ms = 1e30;
        size_t total = 0;
        for(size_t s = 0; s < 3; ++s) {
            auto t0 = std::chrono::high_resolution_clock::now();
            total = 0;
            for(const auto &r : rays) total += fn(r).size();
            auto t1 = std::chrono::high_resolution_clock::now();
            best_ms = std::min(best_ms, std::chrono::duration<double, std::milli>(t1 - t0).count());
        }
        std::cout << "  " << label << ": " << best_ms * 1000.0 / static_cast<double>(nrays)
                  << " us/ray (" << total << " entries)\n";
    };
    timeit("GetColumnLengths  ", [&](const NuGeom::Ray &r) { return world.GetColumnLengths(r); });
    timeit("GetLineSegmentsSwp",
           [&](const NuGeom::Ray &r) { return world.GetLineSegmentsSweep(r); });
    CHECK(true);
}
