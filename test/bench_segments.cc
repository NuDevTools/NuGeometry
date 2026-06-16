#include "catch2/catch.hpp"

#include "geom/LineSegment.hh"
#include "geom/Parser.hh"
#include "geom/Ray.hh"
#include "geom/World.hh"

#include <chrono>
#include <cmath>
#include <iostream>
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
