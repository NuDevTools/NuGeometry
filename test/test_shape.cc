#include "catch2/catch.hpp"
#include "geom/Ray.hh"
#include "geom/Shape.hh"

#include <cmath>
#include <limits>

TEST_CASE("Box", "[Shapes]") {
    SECTION("SDF is correct") {
        NuGeom::Box box;
        NuGeom::Vector3D point;
        CHECK(box.SignedDistance(point) == -0.5);
        point = NuGeom::Vector3D(0.5, 0.5, 0.5);
        CHECK(box.SignedDistance(point) == 0);
        point = NuGeom::Vector3D(1.5, 0.5, 0.5);
        CHECK(box.SignedDistance(point) == 1);
    }

    SECTION("Translated box") {
        NuGeom::Box box{{1, 1, 1}, {}, {1, 2, 3}};
        NuGeom::Vector3D point{1, 2, 3};
        CHECK(box.SignedDistance(point) == -0.5);
        point = NuGeom::Vector3D(1.5, 2.5, 3.4);
        CHECK(box.SignedDistance(point) == 0);
        point = NuGeom::Vector3D(2.5, 2.5, 3.5);
        CHECK(box.SignedDistance(point) == 1);
    }

    SECTION("Rotated box") {
        NuGeom::Box box{{1, 1, 1}, {{0, 0, 1}, 45 * M_PI / 180}};
        NuGeom::Vector3D point;
        CHECK(box.SignedDistance(point) == -0.5);
        point = NuGeom::Vector3D(sqrt(2) / 2, 0, 0);
        CHECK(box.SignedDistance(point) == Approx(0).margin(1e-15));
        point = NuGeom::Vector3D(0, 0, 1.5);
        CHECK(box.SignedDistance(point) == 1);
    }

    SECTION("Volume is correct") {
        NuGeom::Box box;
        CHECK(box.Volume() == 1);
    }

    SECTION("Intersect2 from outside") {
        // Default Box(): half-extents 0.5 in all dims, so z in [-0.5, 0.5]
        NuGeom::Box box;
        NuGeom::Ray ray{{0, 0, -2}, {0, 0, 1}, 1.0};
        auto [t1, t2] = box.Intersect2(ray);
        CHECK(t1 == Approx(1.5));
        CHECK(t2 == Approx(2.5));
    }

    SECTION("Intersect2 from inside") {
        NuGeom::Box box;
        NuGeom::Ray ray{{0, 0, 0}, {0, 0, 1}, 1.0};
        auto [t1, t2] = box.Intersect2(ray);
        CHECK(t1 == Approx(-0.5));
        CHECK(t2 == Approx(0.5));
        // Intersect (first forward hit) should return t2 when inside
        CHECK(box.Intersect(ray) == Approx(0.5));
    }

    SECTION("Intersect2 miss") {
        NuGeom::Box box;
        NuGeom::Ray ray{{2, 0, -2}, {0, 0, 1}, 1.0};
        auto [t1, t2] = box.Intersect2(ray);
        CHECK(!std::isfinite(t1));
        CHECK(!std::isfinite(t2));
        CHECK(!std::isfinite(box.Intersect(ray)));
    }

    SECTION("Intersect2 diagonal ray") {
        // 2x2x2 box; ray from corner at (-1,-1,-1) going (1,1,1)/sqrt(3)
        NuGeom::Box box{{2, 2, 2}};
        NuGeom::Ray ray{{-3, -3, -3}, {1, 1, 1}, 1.0};
        auto [t1, t2] = box.Intersect2(ray);
        // Enters at (-1,-1,-1) which is t = 2*sqrt(3), exits at (1,1,1) t = 4*sqrt(3)
        CHECK(t1 == Approx(2 * std::sqrt(3.0)));
        CHECK(t2 == Approx(4 * std::sqrt(3.0)));
    }

    SECTION("GetBoundingBox") {
        NuGeom::Box box;
        auto bb = box.GetBoundingBox();
        CHECK(bb.IsValid());
        CHECK(bb.min.X() == Approx(-0.5));
        CHECK(bb.max.X() == Approx(0.5));
        CHECK(bb.min.Y() == Approx(-0.5));
        CHECK(bb.max.Y() == Approx(0.5));
        CHECK(bb.min.Z() == Approx(-0.5));
        CHECK(bb.max.Z() == Approx(0.5));
    }

    SECTION("GetBoundingBox for sized box") {
        NuGeom::Box box{{4, 6, 8}};
        auto bb = box.GetBoundingBox();
        CHECK(bb.IsValid());
        CHECK(bb.min.X() == Approx(-2.0));
        CHECK(bb.max.X() == Approx(2.0));
        CHECK(bb.min.Y() == Approx(-3.0));
        CHECK(bb.max.Y() == Approx(3.0));
        CHECK(bb.min.Z() == Approx(-4.0));
        CHECK(bb.max.Z() == Approx(4.0));
    }
}

TEST_CASE("Sphere", "[Shapes]") {
    SECTION("SDF is correct") {
        NuGeom::Sphere sphere;
        NuGeom::Vector3D point;
        CHECK(sphere.SignedDistance(point) == -1);
        point = NuGeom::Vector3D(0, 1, 0);
        CHECK(sphere.SignedDistance(point) == 0);
        point = NuGeom::Vector3D(1, 1, 1);
        CHECK(sphere.SignedDistance(point) == sqrt(3) - 1);
    }

    SECTION("Translated sphere") {
        NuGeom::Sphere sphere{1, {}, {1, 2, 3}};
        NuGeom::Vector3D point{1, 2, 3};
        CHECK(sphere.SignedDistance(point) == -1);
        point = NuGeom::Vector3D(1, 3, 3);
        CHECK(sphere.SignedDistance(point) == 0);
        point = NuGeom::Vector3D(2, 3, 4);
        CHECK(sphere.SignedDistance(point) == sqrt(3) - 1);
    }

    SECTION("Rotated sphere") {
        NuGeom::Sphere sphere{1, {{0, 0, 1}, 45 * M_PI / 180}};
        NuGeom::Vector3D point;
        CHECK(sphere.SignedDistance(point) == -1);
        point = NuGeom::Vector3D(0, 1, 0);
        CHECK(sphere.SignedDistance(point) == 0);
        point = NuGeom::Vector3D(1, 1, 1);
        CHECK(sphere.SignedDistance(point) == sqrt(3) - 1);
    }

    SECTION("Volume is correct") {
        NuGeom::Sphere sphere;
        CHECK(sphere.Volume() == Approx(4.0 * M_PI / 3.0));
    }

    SECTION("Intersect2 from outside") {
        // Sphere radius=1, ray approaching along z-axis from z=-2
        NuGeom::Sphere sphere;
        NuGeom::Ray ray{{0, 0, -2}, {0, 0, 1}, 1.0};
        auto [t1, t2] = sphere.Intersect2(ray);
        CHECK(t1 == Approx(1.0)); // enters at z = -1
        CHECK(t2 == Approx(3.0)); // exits at z = +1
    }

    SECTION("Intersect2 from inside") {
        NuGeom::Sphere sphere;
        NuGeom::Ray ray{{0, 0, 0}, {0, 0, 1}, 1.0};
        auto [t1, t2] = sphere.Intersect2(ray);
        CHECK(t1 == Approx(-1.0)); // back surface
        CHECK(t2 == Approx(1.0));  // front surface
        CHECK(sphere.Intersect(ray) == Approx(1.0));
    }

    SECTION("Intersect2 miss") {
        NuGeom::Sphere sphere;
        NuGeom::Ray ray{{2, 0, 0}, {0, 0, 1}, 1.0}; // passes at x=2, outside radius
        auto [t1, t2] = sphere.Intersect2(ray);
        CHECK(!std::isfinite(t1));
        CHECK(!std::isfinite(t2));
    }

    SECTION("Intersect2 tangent ray") {
        // Ray tangent to the sphere: barely touches at one point
        NuGeom::Sphere sphere;
        NuGeom::Ray ray{{1, 0, -2}, {0, 0, 1}, 1.0}; // touches at (1,0,0)
        auto [t1, t2] = sphere.Intersect2(ray);
        CHECK(std::isfinite(t1));
        CHECK(t1 == Approx(t2).margin(1e-10));
    }

    SECTION("GetBoundingBox") {
        NuGeom::Sphere sphere;
        auto bb = sphere.GetBoundingBox();
        CHECK(bb.IsValid());
        CHECK(bb.min.X() == Approx(-1.0));
        CHECK(bb.max.X() == Approx(1.0));
        CHECK(bb.min.Y() == Approx(-1.0));
        CHECK(bb.max.Y() == Approx(1.0));
        CHECK(bb.min.Z() == Approx(-1.0));
        CHECK(bb.max.Z() == Approx(1.0));
    }
}

TEST_CASE("Cylinder", "[Shapes]") {
    SECTION("SDF is correct") {
        NuGeom::Cylinder cylinder;
        NuGeom::Vector3D point;
        CHECK(cylinder.SignedDistance(point) == -1);
        point = NuGeom::Vector3D(1.0 / sqrt(2), 1.0 / sqrt(2), 0);
        CHECK(cylinder.SignedDistance(point) == Approx(0.0).margin(1e-15));
        point = NuGeom::Vector3D(1, 0, 2);
        CHECK(cylinder.SignedDistance(point) == 1);
    }

    SECTION("Translated cylinder") {
        NuGeom::Cylinder cylinder{1, 1, {}, {1, 2, 3}};
        NuGeom::Vector3D point{1, 2, 3};
        CHECK(cylinder.SignedDistance(point) == -1);
        point = NuGeom::Vector3D(1 + 1.0 / sqrt(2.0), 2 + 1.0 / sqrt(2.0), 3);
        CHECK(cylinder.SignedDistance(point) == Approx(0.0).margin(1e-15));
        point = NuGeom::Vector3D(2, 2, 5);
        CHECK(cylinder.SignedDistance(point) == 1);
    }

    SECTION("Rotated cylinder") {
        NuGeom::Cylinder cylinder{1, 1, {{0, 1, 0}, M_PI / 2}};
        NuGeom::Vector3D point;
        CHECK(cylinder.SignedDistance(point) == -1);
        point = NuGeom::Vector3D(0, 1.0 / sqrt(2), 1.0 / sqrt(2));
        CHECK(cylinder.SignedDistance(point) == Approx(0.0).margin(1e-15));
        point = NuGeom::Vector3D(2, 1, 0);
        CHECK(cylinder.SignedDistance(point) == 1);
    }

    SECTION("Volume is correct") {
        NuGeom::Cylinder cylinder;
        CHECK(cylinder.Volume() == Approx(M_PI));
    }

    SECTION("Intersect2 from side at z=0") {
        // Ray along x-axis at z=0 (inside both z-conventions)
        // Default cylinder: radius=1, height=1.
        // SDF says z in [-1,1]; regardless, z=0 is in the interior.
        NuGeom::Cylinder cylinder;
        NuGeom::Ray ray{{3, 0, 0}, {-1, 0, 0}, 1.0};
        auto [t1, t2] = cylinder.Intersect2(ray);
        CHECK(t1 == Approx(2.0)); // hits curved surface at x=1
        CHECK(t2 == Approx(4.0)); // exits at x=-1
    }

    SECTION("Intersect2 miss") {
        NuGeom::Cylinder cylinder;
        NuGeom::Ray ray{{3, 0, 0}, {0, 0, 1}, 1.0}; // parallel to axis, outside
        auto [t1, t2] = cylinder.Intersect2(ray);
        CHECK(!std::isfinite(t1));
        CHECK(!std::isfinite(t2));
    }

    SECTION("GetBoundingBox") {
        NuGeom::Cylinder cylinder;
        auto bb = cylinder.GetBoundingBox();
        CHECK(bb.IsValid());
        CHECK(bb.min.X() == Approx(-1.0));
        CHECK(bb.max.X() == Approx(1.0));
        CHECK(bb.min.Y() == Approx(-1.0));
        CHECK(bb.max.Y() == Approx(1.0));
        // Bounding box uses half-height convention matching SDF
        CHECK(bb.min.Z() == Approx(-1.0));
        CHECK(bb.max.Z() == Approx(1.0));
    }
}

// ---------------------------------------------------------------------------
// Polyhedra
// ---------------------------------------------------------------------------
TEST_CASE("Polyhedra", "[Shapes]") {
    // Regular hexagonal prism: nsides=6, full 2π, rmax=1, z in [-1, 1]
    std::vector<NuGeom::zplane> planes = {{0.0, 1.0, -1.0}, {0.0, 1.0, 1.0}};
    NuGeom::Polyhedra hex(0.0, 2.0 * M_PI, 6, planes);

    SECTION("Center is inside") {
        CHECK(hex.SignedDistance({0, 0, 0}) < 0);
    }

    SECTION("Point outside along z-axis") {
        CHECK(hex.SignedDistance({0, 0, 2}) > 0);
    }

    SECTION("Point outside along x-axis") {
        CHECK(hex.SignedDistance({2, 0, 0}) > 0);
    }

    SECTION("Point inside/outside along x near vertex") {
        // For a hexagon with vertices at rmax=1 and startphi=0, the vertex
        // at phi=0 is at (1, 0, 0).  Along the x-axis the face boundary
        // passes through x=1 (the vertex), so:
        CHECK(hex.SignedDistance({0.99, 0, 0}) < 0);
        CHECK(hex.SignedDistance({1.01, 0, 0}) > 0);
    }

    SECTION("Intersect2 along z-axis from outside") {
        NuGeom::Ray ray{{0, 0, -3}, {0, 0, 1}, 1.0};
        auto [t1, t2] = hex.Intersect2(ray);
        CHECK(t1 == Approx(2.0)); // enters at z=-1
        CHECK(t2 == Approx(4.0)); // exits at z=+1
    }

    SECTION("Intersect2 from inside") {
        NuGeom::Ray ray{{0, 0, 0}, {0, 0, 1}, 1.0};
        auto [t1, t2] = hex.Intersect2(ray);
        CHECK(t1 < 0); // back cap
        CHECK(t2 > 0); // front cap
        CHECK(t1 == Approx(-1.0));
        CHECK(t2 == Approx(1.0));
    }

    SECTION("Intersect2 miss") {
        // Ray outside the hexagon's XY footprint
        NuGeom::Ray ray{{2, 0, -3}, {0, 0, 1}, 1.0};
        auto [t1, t2] = hex.Intersect2(ray);
        CHECK(!std::isfinite(t1));
        CHECK(!std::isfinite(t2));
    }

    SECTION("GetBoundingBox") {
        auto bb = hex.GetBoundingBox();
        CHECK(bb.IsValid());
        CHECK(bb.min.X() == Approx(-1.0));
        CHECK(bb.max.X() == Approx(1.0));
        CHECK(bb.min.Y() == Approx(-1.0));
        CHECK(bb.max.Y() == Approx(1.0));
        CHECK(bb.min.Z() == Approx(-1.0));
        CHECK(bb.max.Z() == Approx(1.0));
    }

    SECTION("Volume of hexagonal prism") {
        // Regular hexagon area = (3√3/2)*r² = 3√3/2 for r=1
        // Height = 2 → Volume = 3√3
        double expected = 3.0 * std::sqrt(3.0);
        CHECK(hex.Volume() == Approx(expected).epsilon(1e-6));
    }

    SECTION("Partial phi wedge") {
        // Half-hexagon (deltaphi = π, startphi = 0) covers y ≥ 0 half
        NuGeom::Polyhedra halfhex(0.0, M_PI, 6, planes);
        // Origin is on the cut plane (y=0) so SDF=0 there; use offset points
        CHECK(halfhex.SignedDistance({0, 0.1, 0}) < 0);  // inside (y > 0)
        CHECK(halfhex.SignedDistance({0, -0.1, 0}) > 0); // outside (y < 0)
        CHECK(halfhex.Volume() == Approx(hex.Volume() / 2.0).epsilon(1e-6));
    }

    SECTION("Requires at least two planes") {
        std::vector<NuGeom::zplane> one{{0.0, 1.0, 0.0}};
        CHECK_THROWS_AS(NuGeom::Polyhedra(0.0, 2.0 * M_PI, 6, one), std::runtime_error);
    }

    SECTION("Requires at least three sides") {
        CHECK_THROWS_AS(NuGeom::Polyhedra(0.0, 2.0 * M_PI, 2, planes), std::runtime_error);
    }
}

// ---------------------------------------------------------------------------
// Trapezoid
// ---------------------------------------------------------------------------
TEST_CASE("Trapezoid", "[Shapes]") {
    SECTION("Symmetric trapezoid acts like a box") {
        // x1=x2=1, y1=y2=1, z=1 → a 2×2×2 box
        NuGeom::Trapezoid trd(1.0, 1.0, 1.0, 1.0, 1.0);

        SECTION("Center is inside") {
            CHECK(trd.SignedDistance({0, 0, 0}) == Approx(-1.0));
        }

        SECTION("Point outside along z") {
            CHECK(trd.SignedDistance({0, 0, 2}) == Approx(1.0));
        }

        SECTION("Point outside along x") {
            CHECK(trd.SignedDistance({2, 0, 0}) == Approx(1.0));
        }

        SECTION("Corner point on surface") {
            // Nearest point to (1,1,1) is the corner; SDF = 0
            CHECK(trd.SignedDistance({1, 1, 1}) == Approx(0.0));
        }

        SECTION("Intersect2 along z-axis from outside") {
            NuGeom::Ray ray{{0, 0, -3}, {0, 0, 1}, 1.0};
            auto [t1, t2] = trd.Intersect2(ray);
            CHECK(t1 == Approx(2.0));
            CHECK(t2 == Approx(4.0));
        }

        SECTION("Intersect2 from inside") {
            NuGeom::Ray ray{{0, 0, 0}, {0, 0, 1}, 1.0};
            auto [t1, t2] = trd.Intersect2(ray);
            CHECK(t1 == Approx(-1.0));
            CHECK(t2 == Approx(1.0));
        }

        SECTION("Intersect2 along x-axis from outside") {
            NuGeom::Ray ray{{3, 0, 0}, {-1, 0, 0}, 1.0};
            auto [t1, t2] = trd.Intersect2(ray);
            CHECK(t1 == Approx(2.0));
            CHECK(t2 == Approx(4.0));
        }

        SECTION("Intersect2 miss") {
            NuGeom::Ray ray{{3, 0, 0}, {0, 0, 1}, 1.0};
            auto [t1, t2] = trd.Intersect2(ray);
            CHECK(!std::isfinite(t1));
        }

        SECTION("Volume equals box volume") {
            // A 2×2×2 box has volume 8
            CHECK(trd.Volume() == Approx(8.0));
        }

        SECTION("GetBoundingBox") {
            auto bb = trd.GetBoundingBox();
            CHECK(bb.IsValid());
            CHECK(bb.min.X() == Approx(-1.0));
            CHECK(bb.max.X() == Approx(1.0));
            CHECK(bb.min.Y() == Approx(-1.0));
            CHECK(bb.max.Y() == Approx(1.0));
            CHECK(bb.min.Z() == Approx(-1.0));
            CHECK(bb.max.Z() == Approx(1.0));
        }
    }

    SECTION("Asymmetric trapezoid (narrower at -z, wider at +z)") {
        // x1=0.5 (at z=-1), x2=1.5 (at z=+1), symmetric in y
        NuGeom::Trapezoid trd(0.5, 1.5, 0.5, 1.5, 1.0);

        SECTION("Center is inside") {
            CHECK(trd.SignedDistance({0, 0, 0}) < 0);
        }

        SECTION("Point outside along z") {
            CHECK(trd.SignedDistance({0, 0, 2}) > 0);
        }

        SECTION("Point inside at +z face midpoint") {
            // At z=+1 the half-width is x2=1.5; origin x=0 is inside
            CHECK(trd.SignedDistance({0, 0, 1}) == Approx(0.0));
        }

        SECTION("Point outside at -z face, x too large") {
            // At z=-1 the half-width is x1=0.5; x=0.6 is outside
            CHECK(trd.SignedDistance({0.6, 0, -1}) > 0);
        }

        SECTION("x1 and x2 are not swapped") {
            // At z=-1 half-width should be x1=0.5, not x2=1.5
            // x=0.6 at z=-1: outside (narrow face)
            CHECK(trd.SignedDistance({0.6, 0, -1}) > 0);
            // x=0.6 at z=+0.9: inside (wide face has half-width ≈ 1.45 at z=0.9)
            CHECK(trd.SignedDistance({0.6, 0, 0.9}) < 0);
        }

        SECTION("Intersect2 along z-axis") {
            NuGeom::Ray ray{{0, 0, -3}, {0, 0, 1}, 1.0};
            auto [t1, t2] = trd.Intersect2(ray);
            CHECK(t1 == Approx(2.0));
            CHECK(t2 == Approx(4.0));
        }

        SECTION("Volume of frustum") {
            // Correct volume formula for linearly-varying trapezoid:
            // V = (4z/3) * (2*x1*y1 + 2*x2*y2 + x1*y2 + x2*y1)
            // x1=0.5, x2=1.5, y1=0.5, y2=1.5, z=1
            double x1 = 0.5, x2 = 1.5, y1 = 0.5, y2 = 1.5, z = 1.0;
            double expected = (4.0 * z / 3.0) * (2 * x1 * y1 + 2 * x2 * y2 + x1 * y2 + x2 * y1);
            CHECK(trd.Volume() == Approx(expected));
        }

        SECTION("GetBoundingBox uses maximum half-extents") {
            auto bb = trd.GetBoundingBox();
            CHECK(bb.IsValid());
            CHECK(bb.max.X() == Approx(1.5)); // max(x1, x2) = 1.5
            CHECK(bb.min.X() == Approx(-1.5));
            CHECK(bb.max.Z() == Approx(1.0));
            CHECK(bb.min.Z() == Approx(-1.0));
        }
    }
}

// ---------------------------------------------------------------------------
// Random point generator (used by CombinedShape SDF tests)
// ---------------------------------------------------------------------------
class PointGenerator : public Catch::Generators::IGenerator<NuGeom::Vector3D> {
    std::minstd_rand m_rand;
    std::uniform_real_distribution<double> m_dist;
    NuGeom::Vector3D current_point;

  public:
    PointGenerator(double low, double high) : m_rand{std::random_device{}()}, m_dist{low, high} {
        static_cast<void>(next());
    }

    NuGeom::Vector3D const &get() const override;
    bool next() override {
        current_point.X() = m_dist(m_rand);
        current_point.Y() = m_dist(m_rand);
        current_point.Z() = m_dist(m_rand);
        return true;
    }
};

NuGeom::Vector3D const &PointGenerator::get() const {
    return current_point;
}

Catch::Generators::GeneratorWrapper<NuGeom::Vector3D> randomPoint(double low, double high) {
    return Catch::Generators::GeneratorWrapper<NuGeom::Vector3D>(
        std::unique_ptr<Catch::Generators::IGenerator<NuGeom::Vector3D>>(
            new PointGenerator(low, high)));
}

// ---------------------------------------------------------------------------
// CombinedShape — SDF correctness (random-point tests, unchanged)
// ---------------------------------------------------------------------------
TEST_CASE("Combined Shape", "[Shapes]") {
    NuGeom::Vector3D size{2, 2, 2};
    NuGeom::Rotation3D rotation;
    NuGeom::Translation3D translation{0, 0, 2};
    auto box = std::make_shared<NuGeom::Box>(size);
    auto sphere = std::make_shared<NuGeom::Sphere>(1, rotation, translation);
    NuGeom::Vector3D point = GENERATE(take(30, randomPoint(-5, 5)));

    SECTION("Single Union") {
        NuGeom::CombinedShape shape(box, sphere, NuGeom::ShapeBinaryOp::kUnion);
        CHECK(shape.SignedDistance(point) ==
              std::min(box->SignedDistance(point), sphere->SignedDistance(point)));
        NuGeom::CombinedShape shape2(sphere, box, NuGeom::ShapeBinaryOp::kUnion);
        CHECK(shape.SignedDistance(point) == shape2.SignedDistance(point));
    }

    SECTION("Single Intersect") {
        NuGeom::CombinedShape shape(box, sphere, NuGeom::ShapeBinaryOp::kIntersect);
        CHECK(shape.SignedDistance(point) ==
              std::max(box->SignedDistance(point), sphere->SignedDistance(point)));
        NuGeom::CombinedShape shape2(sphere, box, NuGeom::ShapeBinaryOp::kIntersect);
        CHECK(shape.SignedDistance(point) == shape2.SignedDistance(point));
    }

    SECTION("Single Subtraction") {
        NuGeom::CombinedShape shape(box, sphere, NuGeom::ShapeBinaryOp::kSubtraction);
        CHECK(shape.SignedDistance(point) ==
              std::max(-box->SignedDistance(point), sphere->SignedDistance(point)));
        NuGeom::CombinedShape shape2(sphere, box, NuGeom::ShapeBinaryOp::kSubtraction);
        CHECK(shape.SignedDistance(point) != shape2.SignedDistance(point));
    }
}

// ---------------------------------------------------------------------------
// CombinedShape — Intersect2 (fixed geometry, deterministic)
// ---------------------------------------------------------------------------
TEST_CASE("Combined Shape Intersect2 disjoint", "[Shapes]") {
    // box1: 2x2x2 at origin  → z in [-1, 1], Intersect2 interval {2, 4}
    // box2: 2x2x2 at (0,0,3) → z in [ 2, 4], Intersect2 interval {5, 7}
    // Ray: origin (0,0,-3), direction +z
    auto box1 = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    auto box2 = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2}, NuGeom::Rotation3D{},
                                              NuGeom::Translation3D{0, 0, 3});
    NuGeom::Ray ray{{0, 0, -3}, {0, 0, 1}, 1.0};

    SECTION("Union returns first segment") {
        NuGeom::CombinedShape shape(box1, box2, NuGeom::ShapeBinaryOp::kUnion);
        auto [t1, t2] = shape.Intersect2(ray);
        CHECK(t1 == Approx(2.0));
        CHECK(t2 == Approx(4.0));
    }

    SECTION("Intersection misses disjoint boxes") {
        NuGeom::CombinedShape shape(box1, box2, NuGeom::ShapeBinaryOp::kIntersect);
        auto [t1, t2] = shape.Intersect2(ray);
        CHECK(!std::isfinite(t1));
    }

    SECTION("Subtraction right-minus-left returns second segment") {
        // kSubtraction = right - left, so box2 - box1
        NuGeom::CombinedShape shape(box1, box2, NuGeom::ShapeBinaryOp::kSubtraction);
        auto [t1, t2] = shape.Intersect2(ray);
        CHECK(t1 == Approx(5.0));
        CHECK(t2 == Approx(7.0));
    }
}

TEST_CASE("Combined Shape Intersect2 overlapping", "[Shapes]") {
    // box1: 2x2x2 at origin        → z in [-1, 1],   interval {2, 4}
    // box2: 2x2x2 at (0, 0, 1.5)  → z in [0.5, 2.5], interval {3.5, 5.5}
    // Ray: origin (0,0,-3), direction +z
    auto box1 = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    auto box2 = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2}, NuGeom::Rotation3D{},
                                              NuGeom::Translation3D{0, 0, 1.5});
    NuGeom::Ray ray{{0, 0, -3}, {0, 0, 1}, 1.0};

    SECTION("Union spans merged interval") {
        NuGeom::CombinedShape shape(box1, box2, NuGeom::ShapeBinaryOp::kUnion);
        auto [t1, t2] = shape.Intersect2(ray);
        CHECK(t1 == Approx(2.0));
        CHECK(t2 == Approx(5.5));
    }

    SECTION("Intersection returns overlap region") {
        NuGeom::CombinedShape shape(box1, box2, NuGeom::ShapeBinaryOp::kIntersect);
        auto [t1, t2] = shape.Intersect2(ray);
        CHECK(t1 == Approx(3.5));
        CHECK(t2 == Approx(4.0));
    }

    SECTION("Subtraction right-minus-left") {
        // box2 - box1: keep region inside box2 but outside box1
        // box2: [3.5, 5.5], box1: [2, 4]  → remaining: [4, 5.5]
        NuGeom::CombinedShape shape(box1, box2, NuGeom::ShapeBinaryOp::kSubtraction);
        auto [t1, t2] = shape.Intersect2(ray);
        CHECK(t1 == Approx(4.0));
        CHECK(t2 == Approx(5.5));
    }
}

TEST_CASE("Combined Shape Intersect2 from inside", "[Shapes]") {
    // Two 2x2x2 boxes fully overlapping (union = one box)
    auto box1 = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    auto box2 = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    NuGeom::Ray ray{{0, 0, 0}, {0, 0, 1}, 1.0}; // origin inside both

    SECTION("Union from inside: exit at t=1") {
        NuGeom::CombinedShape shape(box1, box2, NuGeom::ShapeBinaryOp::kUnion);
        auto [t1, t2] = shape.Intersect2(ray);
        // Origin is inside; t1 should be ≤ 0, t2 = 1 (exit face)
        CHECK(t1 <= 0.0);
        CHECK(t2 == Approx(1.0));
        CHECK(shape.Intersect(ray) == Approx(1.0));
    }

    SECTION("Intersection from inside: exit at t=1") {
        NuGeom::CombinedShape shape(box1, box2, NuGeom::ShapeBinaryOp::kIntersect);
        auto [t1, t2] = shape.Intersect2(ray);
        CHECK(t1 <= 0.0);
        CHECK(t2 == Approx(1.0));
    }
}

// ---------------------------------------------------------------------------
// CombinedShape — GetBoundingBox
// Note: GetBoundingBox() returns bounds in each child shape's local frame;
// it does NOT account for the child shape's own translation transform.
// These tests therefore use untranslated shapes where local == world frame.
// ---------------------------------------------------------------------------
TEST_CASE("Combined Shape GetBoundingBox", "[Shapes]") {
    // box1: 2x2x2 → [-1,1]³
    // box2: 4x4x4 → [-2,2]³  (both centred at origin, no translation)
    auto box1 = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    auto box2 = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 4});

    SECTION("Union bounding box spans both") {
        NuGeom::CombinedShape shape(box1, box2, NuGeom::ShapeBinaryOp::kUnion);
        auto bb = shape.GetBoundingBox();
        CHECK(bb.IsValid());
        // Union = larger box extent
        CHECK(bb.min.X() == Approx(-2.0));
        CHECK(bb.max.X() == Approx(2.0));
        CHECK(bb.min.Y() == Approx(-2.0));
        CHECK(bb.max.Y() == Approx(2.0));
    }

    SECTION("Intersection bounding box is overlap region") {
        NuGeom::CombinedShape shape(box1, box2, NuGeom::ShapeBinaryOp::kIntersect);
        auto bb = shape.GetBoundingBox();
        CHECK(bb.IsValid());
        // Intersection = smaller box extent
        CHECK(bb.min.X() == Approx(-1.0));
        CHECK(bb.max.X() == Approx(1.0));
    }

    SECTION("Subtraction bounding box is right operand (box2) extent") {
        NuGeom::CombinedShape shape(box1, box2, NuGeom::ShapeBinaryOp::kSubtraction);
        auto bb = shape.GetBoundingBox();
        CHECK(bb.IsValid());
        // Result ⊆ box2 = [-2,2]³
        CHECK(bb.min.X() == Approx(-2.0));
        CHECK(bb.max.X() == Approx(2.0));
    }
}

// ---------------------------------------------------------------------------
// CombinedShape — Volume (Monte Carlo; use wide tolerance)
// Note: Volume() samples the bounding box returned by GetBoundingBox(), which
// works in local (untranslated) coordinates.  Use concentric untranslated
// shapes so the bounding box is correct.
// ---------------------------------------------------------------------------
TEST_CASE("Combined Shape Volume", "[Shapes]") {
    // Box(4x4x4) contains Sphere(r=1) entirely.
    // Intersection = sphere;  volume ≈ 4π/3 ≈ 4.189
    // Bounding box of intersection = [-1,1]³ (sphere's local box), volume 8.
    // MC fills sphere portion: fraction ≈ 4π/3/8 ≈ 0.524.
    auto bigbox = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 4});
    auto sphere = std::make_shared<NuGeom::Sphere>(1.0);
    NuGeom::CombinedShape shape(bigbox, sphere, NuGeom::ShapeBinaryOp::kIntersect);
    // MC with ~4M samples; accept 1% relative tolerance
    CHECK(shape.Volume() == Approx(4.0 * M_PI / 3.0).epsilon(0.01));
}
