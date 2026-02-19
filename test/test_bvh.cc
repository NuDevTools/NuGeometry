#include "catch2/catch.hpp"

#include "geom/BVH.hh"
#include "geom/BoundingBox.hh"
#include "geom/Element.hh"
#include "geom/Material.hh"
#include "geom/Ray.hh"
#include "geom/Shape.hh"
#include "geom/Transform3D.hh"
#include "geom/Volume.hh"
#include "geom/World.hh"

#include <cmath>
#include <limits>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static NuGeom::Material make_mat() {
    NuGeom::Material m("Test", 1.0, 1);
    m.AddElement(NuGeom::Element("Hydrogen", 1, 1), 1);
    return m;
}

/// Brute-force minimum intersection time — mirrors the pre-BVH linear scan.
static double brute_force(const std::vector<std::shared_ptr<NuGeom::PhysicalVolume>> &daughters,
                          const NuGeom::Ray &ray) {
    double tmin = std::numeric_limits<double>::infinity();
    for(const auto &d : daughters) tmin = std::min(tmin, d->Intersect(ray));
    return tmin;
}

/// Brute-force nearest-daughter pointer.
static std::shared_ptr<NuGeom::PhysicalVolume>
brute_force_vol(const std::vector<std::shared_ptr<NuGeom::PhysicalVolume>> &daughters,
                const NuGeom::Ray &ray) {
    double tmin = std::numeric_limits<double>::infinity();
    std::shared_ptr<NuGeom::PhysicalVolume> best;
    for(const auto &d : daughters) {
        double t = d->Intersect(ray);
        if(t < tmin) {
            tmin = t;
            best = d;
        }
    }
    return best;
}

// ---------------------------------------------------------------------------
// BoundingBox::Merge
// ---------------------------------------------------------------------------

TEST_CASE("BoundingBox Merge", "[AABB]") {
    SECTION("Disjoint boxes span combined extent") {
        NuGeom::BoundingBox a{{-1, -1, -1}, {0, 0, 0}};
        NuGeom::BoundingBox b{{1, 1, 1}, {2, 2, 2}};
        auto c = NuGeom::BoundingBox::Merge(a, b);
        CHECK(c.min.X() == Approx(-1.0));
        CHECK(c.min.Y() == Approx(-1.0));
        CHECK(c.min.Z() == Approx(-1.0));
        CHECK(c.max.X() == Approx(2.0));
        CHECK(c.max.Y() == Approx(2.0));
        CHECK(c.max.Z() == Approx(2.0));
        CHECK(c.IsValid());
    }

    SECTION("Merge with self is a no-op") {
        NuGeom::BoundingBox a{{-3, -2, -1}, {4, 5, 6}};
        auto c = NuGeom::BoundingBox::Merge(a, a);
        CHECK(c.min.X() == Approx(a.min.X()));
        CHECK(c.min.Y() == Approx(a.min.Y()));
        CHECK(c.min.Z() == Approx(a.min.Z()));
        CHECK(c.max.X() == Approx(a.max.X()));
        CHECK(c.max.Y() == Approx(a.max.Y()));
        CHECK(c.max.Z() == Approx(a.max.Z()));
    }

    SECTION("Smaller box fully inside larger produces the larger") {
        NuGeom::BoundingBox big{{-5, -5, -5}, {5, 5, 5}};
        NuGeom::BoundingBox small_{{-1, -1, -1}, {1, 1, 1}};
        auto c = NuGeom::BoundingBox::Merge(big, small_);
        CHECK(c.min.X() == Approx(-5.0));
        CHECK(c.max.X() == Approx(5.0));
        CHECK(c.min.Y() == Approx(-5.0));
        CHECK(c.max.Y() == Approx(5.0));
    }

    SECTION("Overlapping boxes") {
        NuGeom::BoundingBox a{{0, 0, 0}, {3, 3, 3}};
        NuGeom::BoundingBox b{{1, 1, 1}, {4, 4, 4}};
        auto c = NuGeom::BoundingBox::Merge(a, b);
        CHECK(c.min.X() == Approx(0.0));
        CHECK(c.max.X() == Approx(4.0));
        CHECK(c.IsValid());
    }
}

// ---------------------------------------------------------------------------
// BoundingBox::Intersect  (slab method)
// ---------------------------------------------------------------------------

TEST_CASE("BoundingBox Intersect", "[AABB]") {
    // Unit box centred at origin: [-0.5, 0.5]^3
    NuGeom::BoundingBox box{{-0.5, -0.5, -0.5}, {0.5, 0.5, 0.5}};

    SECTION("Axial ray along +z hits box") {
        NuGeom::Ray ray{{0, 0, -2}, {0, 0, 1}, 1.0};
        CHECK(box.Intersect(ray));
    }

    SECTION("Axial ray along -z hits box") {
        NuGeom::Ray ray{{0, 0, 2}, {0, 0, -1}, 1.0};
        CHECK(box.Intersect(ray));
    }

    SECTION("Ray offset in x misses along z") {
        NuGeom::Ray ray{{2, 0, -2}, {0, 0, 1}, 1.0};
        CHECK_FALSE(box.Intersect(ray));
    }

    SECTION("Ray offset in y misses along z") {
        NuGeom::Ray ray{{0, 2, -2}, {0, 0, 1}, 1.0};
        CHECK_FALSE(box.Intersect(ray));
    }

    SECTION("Ray origin inside box") {
        NuGeom::Ray ray{{0, 0, 0}, {0, 0, 1}, 1.0};
        CHECK(box.Intersect(ray));
    }

    SECTION("t_far cutoff before box entry rejects hit") {
        // Box enters at t = 1.5; cutoff at t = 1.0 should reject
        NuGeom::Ray ray{{0, 0, -2}, {0, 0, 1}, 1.0};
        CHECK_FALSE(box.Intersect(ray, 0.0, 1.0));
    }

    SECTION("t_far just past box entry accepts hit") {
        // Box enters at t = 1.5; cutoff at t = 2.0 should accept
        NuGeom::Ray ray{{0, 0, -2}, {0, 0, 1}, 1.0};
        CHECK(box.Intersect(ray, 0.0, 2.0));
    }

    SECTION("t_near beyond box exit rejects hit") {
        // Box exits at t = 2.5; start search at t = 3.0 should reject
        NuGeom::Ray ray{{0, 0, -2}, {0, 0, 1}, 1.0};
        CHECK_FALSE(box.Intersect(ray, 3.0));
    }

    SECTION("Ray parallel to x-axis inside yz slab") {
        NuGeom::Ray ray{{-2, 0, 0}, {1, 0, 0}, 1.0};
        CHECK(box.Intersect(ray));
    }

    SECTION("Ray parallel to x-axis outside yz slab") {
        NuGeom::Ray ray{{-2, 1, 0}, {1, 0, 0}, 1.0};
        CHECK_FALSE(box.Intersect(ray));
    }

    SECTION("Diagonal ray through box corner") {
        NuGeom::Ray ray{{-2, -2, -2}, {1, 1, 1}, 1.0};
        CHECK(box.Intersect(ray));
    }

    SECTION("Diagonal ray that misses") {
        // Travels at z = 2 (above the box): never enters z slab
        NuGeom::Ray ray{{-2, 0, 2}, {1, 0, 0}, 1.0};
        CHECK_FALSE(box.Intersect(ray));
    }

    SECTION("Ray pointing away from box") {
        // Box is at z in [-0.5, 0.5]; ray starts at z=-2 going -z
        NuGeom::Ray ray{{0, 0, -2}, {0, 0, -1}, 1.0};
        CHECK_FALSE(box.Intersect(ray, 0.0));
    }
}

// ---------------------------------------------------------------------------
// PhysicalVolume::GetParentBoundingBox
// ---------------------------------------------------------------------------

TEST_CASE("GetParentBoundingBox", "[AABB]") {
    auto mat = make_mat();
    // 2x2x2 box → local AABB = [-1, 1]^3
    auto box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    auto vol = std::make_shared<NuGeom::LogicalVolume>(mat, box);

    SECTION("Identity placement equals shape AABB") {
        NuGeom::PhysicalVolume pv(vol, NuGeom::Transform3D{}, NuGeom::Transform3D{});
        auto bb = pv.GetParentBoundingBox();
        CHECK(bb.IsValid());
        CHECK(bb.min.X() == Approx(-1.0));
        CHECK(bb.max.X() == Approx(1.0));
        CHECK(bb.min.Y() == Approx(-1.0));
        CHECK(bb.max.Y() == Approx(1.0));
        CHECK(bb.min.Z() == Approx(-1.0));
        CHECK(bb.max.Z() == Approx(1.0));
    }

    SECTION("Pure translation in x shifts AABB") {
        // Translate by (5, 0, 0): expect x in [4, 6]
        NuGeom::Translation3D trans(5, 0, 0);
        NuGeom::PhysicalVolume pv(vol, trans, NuGeom::Transform3D{});
        auto bb = pv.GetParentBoundingBox();
        CHECK(bb.IsValid());
        CHECK(bb.min.X() == Approx(4.0));
        CHECK(bb.max.X() == Approx(6.0));
        CHECK(bb.min.Y() == Approx(-1.0));
        CHECK(bb.max.Y() == Approx(1.0));
        CHECK(bb.min.Z() == Approx(-1.0));
        CHECK(bb.max.Z() == Approx(1.0));
    }

    SECTION("Pure translation in z shifts AABB") {
        NuGeom::Translation3D trans(0, 0, -10);
        NuGeom::PhysicalVolume pv(vol, trans, NuGeom::Transform3D{});
        auto bb = pv.GetParentBoundingBox();
        CHECK(bb.IsValid());
        CHECK(bb.min.Z() == Approx(-11.0));
        CHECK(bb.max.Z() == Approx(-9.0));
        // x and y unchanged
        CHECK(bb.min.X() == Approx(-1.0));
        CHECK(bb.max.X() == Approx(1.0));
    }

    SECTION("90-degree rotation around Z swaps X and Y half-extents") {
        // 2x4x2 box: local x in [-1,1], y in [-2,2], z in [-1,1].
        // After 90° around Z: RotZ maps (x,y,z) -> (-y, x, z).
        // Parent AABB: x in [-2,2], y in [-1,1], z in [-1,1].
        auto tall = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 4, 2});
        auto tall_vol = std::make_shared<NuGeom::LogicalVolume>(mat, tall);
        NuGeom::RotationZ3D rot(M_PI / 2);
        NuGeom::PhysicalVolume pv(tall_vol, NuGeom::Transform3D{}, rot);
        auto bb = pv.GetParentBoundingBox();
        CHECK(bb.IsValid());
        CHECK(bb.min.X() == Approx(-2.0).margin(1e-10));
        CHECK(bb.max.X() == Approx(2.0).margin(1e-10));
        CHECK(bb.min.Y() == Approx(-1.0).margin(1e-10));
        CHECK(bb.max.Y() == Approx(1.0).margin(1e-10));
        CHECK(bb.min.Z() == Approx(-1.0));
        CHECK(bb.max.Z() == Approx(1.0));
    }

    SECTION("Translation then rotation around Z") {
        // Translate by (5,0,0), then RotZ(90°).
        // Local corners [-1,1]^3:
        //   step1 (trans by (5,0,0)): x in [4,6], y in [-1,1]
        //   step2 (RotZ 90°): x_new = -y_old in [-1,1], y_new = x_old in [4,6]
        // Parent AABB: x in [-1,1], y in [4,6], z in [-1,1].
        NuGeom::Translation3D trans(5, 0, 0);
        NuGeom::RotationZ3D rot(M_PI / 2);
        NuGeom::PhysicalVolume pv(vol, trans, rot);
        auto bb = pv.GetParentBoundingBox();
        CHECK(bb.IsValid());
        CHECK(bb.min.X() == Approx(-1.0).margin(1e-10));
        CHECK(bb.max.X() == Approx(1.0).margin(1e-10));
        CHECK(bb.min.Y() == Approx(4.0).margin(1e-10));
        CHECK(bb.max.Y() == Approx(6.0).margin(1e-10));
        CHECK(bb.min.Z() == Approx(-1.0));
        CHECK(bb.max.Z() == Approx(1.0));
    }

    SECTION("AABB conservatively contains rotated box") {
        // A 45-degree rotation produces an AABB that is larger than the shape.
        // We only check that the AABB is valid and at least as large as the shape.
        NuGeom::RotationZ3D rot(M_PI / 4);
        NuGeom::PhysicalVolume pv(vol, NuGeom::Transform3D{}, rot);
        auto bb = pv.GetParentBoundingBox();
        CHECK(bb.IsValid());
        // The AABB must contain all 8 shape corners after rotation.
        // Shape half-diagonal in XY = sqrt(2) ≈ 1.414, so |x| and |y| ≤ sqrt(2).
        CHECK(bb.min.X() <= Approx(-1.0));
        CHECK(bb.max.X() >= Approx(1.0));
        CHECK(bb.min.Y() <= Approx(-1.0));
        CHECK(bb.max.Y() >= Approx(1.0));
        CHECK(bb.min.Z() == Approx(-1.0));
        CHECK(bb.max.Z() == Approx(1.0));
    }
}

// ---------------------------------------------------------------------------
// BVH edge cases
// ---------------------------------------------------------------------------

TEST_CASE("BVH empty daughter list", "[BVH]") {
    auto mat = make_mat();
    auto world_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{10, 10, 10});
    auto world_vol = std::make_shared<NuGeom::LogicalVolume>(mat, world_box);

    NuGeom::Ray ray{{0, 0, -5}, {0, 0, 1}, 1.0};
    double time = 0;
    std::shared_ptr<NuGeom::PhysicalVolume> hit;
    CHECK_FALSE(world_vol->RayTrace(ray, time, hit));
}

TEST_CASE("BVH single daughter", "[BVH]") {
    auto mat = make_mat();
    // Unit box at origin: z in [-0.5, 0.5]
    auto inner = std::make_shared<NuGeom::Box>();
    auto inner_vol = std::make_shared<NuGeom::LogicalVolume>(mat, inner);
    auto pv = std::make_shared<NuGeom::PhysicalVolume>(inner_vol, NuGeom::Transform3D{},
                                                       NuGeom::Transform3D{});

    auto world_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{10, 10, 10});
    auto world_vol = std::make_shared<NuGeom::LogicalVolume>(mat, world_box);
    world_vol->AddDaughter(pv);

    SECTION("Hit") {
        NuGeom::Ray ray{{0, 0, -1}, {0, 0, 1}, 1.0};
        double time = 0;
        std::shared_ptr<NuGeom::PhysicalVolume> hit;
        REQUIRE(world_vol->RayTrace(ray, time, hit));
        CHECK(time == Approx(0.5));
        CHECK(hit == pv);
    }

    SECTION("Miss") {
        NuGeom::Ray ray{{5, 0, -1}, {0, 0, 1}, 1.0}; // x=5, outside box
        double time = 0;
        std::shared_ptr<NuGeom::PhysicalVolume> hit;
        CHECK_FALSE(world_vol->RayTrace(ray, time, hit));
    }
}

// ---------------------------------------------------------------------------
// BVH vs brute-force regression: three boxes on the z-axis
// ---------------------------------------------------------------------------

TEST_CASE("BVH matches brute-force: three boxes on z-axis", "[BVH]") {
    auto mat = make_mat();

    auto make_pv = [&](double x, double y, double z) {
        auto b = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{1, 1, 1});
        auto v = std::make_shared<NuGeom::LogicalVolume>(mat, b);
        NuGeom::Translation3D trans(x, y, z);
        return std::make_shared<NuGeom::PhysicalVolume>(v, trans, NuGeom::Transform3D{});
    };

    // d1 at z=0 → z in [-0.5, 0.5]
    // d2 at z=3 → z in [2.5, 3.5]
    // d3 at z=6 → z in [5.5, 6.5]
    auto d1 = make_pv(0, 0, 0);
    auto d2 = make_pv(0, 0, 3);
    auto d3 = make_pv(0, 0, 6);

    auto world_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{20, 20, 20});
    auto world_vol = std::make_shared<NuGeom::LogicalVolume>(mat, world_box);
    world_vol->AddDaughter(d1);
    world_vol->AddDaughter(d2);
    world_vol->AddDaughter(d3);
    const auto &daughters = world_vol->Daughters();

    SECTION("Ray along +z hits d1 first at t=4.5") {
        NuGeom::Ray ray{{0, 0, -5}, {0, 0, 1}, 1.0};
        double time = 0;
        std::shared_ptr<NuGeom::PhysicalVolume> vol;
        REQUIRE(world_vol->RayTrace(ray, time, vol));
        CHECK(time == Approx(brute_force(daughters, ray)));
        CHECK(time == Approx(4.5));
        CHECK(vol == d1);
    }

    SECTION("Ray along -z hits d3 first at t=5.5") {
        NuGeom::Ray ray{{0, 0, 12}, {0, 0, -1}, 1.0};
        double time = 0;
        std::shared_ptr<NuGeom::PhysicalVolume> vol;
        REQUIRE(world_vol->RayTrace(ray, time, vol));
        CHECK(time == Approx(brute_force(daughters, ray)));
        CHECK(time == Approx(5.5));
        CHECK(vol == d3);
    }

    SECTION("Lateral ray hits d2 (middle) in x") {
        // Ray at y=0, z=3 going +x; only d2's z range [2.5,3.5] matches
        NuGeom::Ray ray{{-5, 0, 3}, {1, 0, 0}, 1.0};
        double time = 0;
        std::shared_ptr<NuGeom::PhysicalVolume> vol;
        REQUIRE(world_vol->RayTrace(ray, time, vol));
        CHECK(time == Approx(brute_force(daughters, ray)));
        CHECK(vol == d2);
    }

    SECTION("Ray misses all volumes") {
        NuGeom::Ray ray{{5, 0, -5}, {0, 0, 1}, 1.0}; // x=5, outside all boxes
        double time = 0;
        std::shared_ptr<NuGeom::PhysicalVolume> vol;
        CHECK_FALSE(world_vol->RayTrace(ray, time, vol));
        CHECK(!std::isfinite(brute_force(daughters, ray)));
    }

    SECTION("Diagonal ray hits nearest box") {
        // Angled ray that passes through all three z-slabs; hits d1 first
        NuGeom::Ray ray{{0, 0, -5}, {0.01, 0, 1}, 1.0};
        double time = 0;
        std::shared_ptr<NuGeom::PhysicalVolume> vol;
        REQUIRE(world_vol->RayTrace(ray, time, vol));
        CHECK(time == Approx(brute_force(daughters, ray)).margin(1e-10));
        CHECK(vol == brute_force_vol(daughters, ray));
    }

    SECTION("World::RayTrace returns correct 1-based index") {
        NuGeom::World world(world_vol);
        NuGeom::Ray ray{{0, 0, -5}, {0, 0, 1}, 1.0};
        double dist = 0;
        size_t idx = 0;
        REQUIRE(world.RayTrace(ray, dist, idx));
        CHECK(idx == 1); // d1 is daughters[0], World uses 1-based indexing
        CHECK(dist == Approx(brute_force(daughters, ray)));
    }

    SECTION("World::RayTrace index for d3") {
        NuGeom::World world(world_vol);
        NuGeom::Ray ray{{0, 0, 12}, {0, 0, -1}, 1.0};
        double dist = 0;
        size_t idx = 0;
        REQUIRE(world.RayTrace(ray, dist, idx));
        CHECK(idx == 3); // d3 is daughters[2], 1-based → 3
    }
}

// ---------------------------------------------------------------------------
// BVH vs brute-force regression: many spheres (stress test)
// ---------------------------------------------------------------------------

TEST_CASE("BVH matches brute-force: ten spheres in a row", "[BVH]") {
    auto mat = make_mat();
    auto world_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{200, 200, 200});
    auto world_vol = std::make_shared<NuGeom::LogicalVolume>(mat, world_box);

    // 10 unit spheres at x = 0, 10, 20, ..., 90; y = z = 0
    constexpr int N = 10;
    std::vector<std::shared_ptr<NuGeom::PhysicalVolume>> pvols;
    for(int i = 0; i < N; ++i) {
        auto s = std::make_shared<NuGeom::Sphere>(1.0);
        auto v = std::make_shared<NuGeom::LogicalVolume>(mat, s);
        NuGeom::Translation3D trans(i * 10.0, 0, 0);
        auto pv = std::make_shared<NuGeom::PhysicalVolume>(v, trans, NuGeom::Transform3D{});
        world_vol->AddDaughter(pv);
        pvols.push_back(pv);
    }
    const auto &daughters = world_vol->Daughters();

    struct TestRay {
        double ox, oy, oz, dx, dy, dz;
        const char *label;
    };
    std::vector<TestRay> rays = {
        {-5, 0, 0, 1, 0, 0, "+x from left  → hits sphere 0"},
        {95, 0, 0, -1, 0, 0, "-x from right → hits sphere 9"},
        {45, 0, 0, 1, 0, 0, "+x from mid   → hits sphere 5"},
        {45, 0, 0, -1, 0, 0, "-x from mid   → hits sphere 4"},
        {-5, 5, 0, 1, 0, 0, "y=5 offset    → misses all"},
        {-5, -3, 0, 1, 0, 0, "y=-3 offset   → misses all"},
        {-5, 0, 2, 1, 0, 0, "z=2 offset    → misses all"},
    };

    for(const auto &r : rays) {
        NuGeom::Ray ray{{r.ox, r.oy, r.oz}, {r.dx, r.dy, r.dz}, 1.0};
        double bvh_time = 0;
        std::shared_ptr<NuGeom::PhysicalVolume> bvh_vol;
        bool hit = world_vol->RayTrace(ray, bvh_time, bvh_vol);

        double bf_time = brute_force(daughters, ray);

        if(!std::isfinite(bf_time)) {
            CHECK_FALSE(hit);
        } else {
            REQUIRE(hit);
            CHECK(bvh_time == Approx(bf_time).margin(1e-10));
            CHECK(bvh_vol == brute_force_vol(daughters, ray));
        }
    }
}
