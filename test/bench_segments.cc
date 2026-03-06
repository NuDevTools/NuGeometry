#include "catch2/catch.hpp"

#include "geom/LineSegment.hh"
#include "geom/Parser.hh"
#include "geom/Ray.hh"
#include "geom/World.hh"

#include <chrono>
#include <iostream>
#include <random>

TEST_CASE("GetLineSegments throughput", "[benchmark]") {
    std::string gdml = std::string(NUGEOM_SOURCE_DIR) + "/nd_hall_with_lar_tms_nosand.gdml";
    NuGeom::GDMLParser parse(gdml);
    auto world = parse.GetWorld();
    auto box = world.GetWorldBox();

    // Find the ArgonCubeDetector bounding box for targeted sampling.
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
    std::cout << "  Sampling: ArgonCubeDetector x=[" << xlo << "," << xhi << "] y=[" << ylo << ","
              << yhi << "]\n";

    // Start rays inside the world, upstream of all geometry.
    const double z_start = box.min.Z() + 1.0;

    constexpr size_t nrays = 500;
    constexpr unsigned seed = 42;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist_x(xlo, xhi);
    std::uniform_real_distribution<double> dist_y(ylo, yhi);

    std::vector<NuGeom::Ray> rays;
    rays.reserve(nrays);
    for(size_t i = 0; i < nrays; ++i) {
        NuGeom::Vector3D origin(dist_x(rng), dist_y(rng), z_start);
        NuGeom::Vector3D dir(0, 0, 1);
        rays.emplace_back(origin, dir, 1.0);
    }

    // Warm up (also ensures BVH is built).
    for(const auto &ray : rays) {
        auto segs = world.GetLineSegments(ray);
        (void)segs;
    }

    constexpr size_t nsamples = 5;
    for(size_t s = 0; s < nsamples; ++s) {
        size_t total = 0;
        auto t0 = std::chrono::high_resolution_clock::now();
        for(const auto &ray : rays) {
            auto segs = world.GetLineSegments(ray);
            total += segs.size();
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << "  Sample " << s << ": " << ms << " ms  (" << total << " segments, "
                  << total / nrays << " segs/ray, " << static_cast<double>(nrays) / (ms / 1000.0)
                  << " rays/s)\n";
    }
    REQUIRE(true);
}
