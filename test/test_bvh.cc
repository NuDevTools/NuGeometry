#include "catch2/catch.hpp"

#include "geom/BVH.hh"
#include "geom/BoundingBox.hh"
#include "geom/Element.hh"
#include "geom/LineSegment.hh"
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

static NuGeom::Material make_named_mat(const std::string &name) {
    NuGeom::Material m(name, 1.0, 1);
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
        // Position at (5,0,0) with RotZ(90°) orientation.
        // GDML semantics: position and rotation are independent.
        // Local→parent: R^T * p_local + position.
        // RotZ(90°)^T maps (x,y,z)→(y,-x,z), so local [-1,1]^3 rotates to
        // [-1,1]^3 (cube is symmetric), then shifts by (5,0,0).
        // Parent AABB: x in [4,6], y in [-1,1], z in [-1,1].
        NuGeom::Translation3D trans(5, 0, 0);
        NuGeom::RotationZ3D rot(M_PI / 2);
        NuGeom::PhysicalVolume pv(vol, trans, rot);
        auto bb = pv.GetParentBoundingBox();
        CHECK(bb.IsValid());
        CHECK(bb.min.X() == Approx(4.0).margin(1e-10));
        CHECK(bb.max.X() == Approx(6.0).margin(1e-10));
        CHECK(bb.min.Y() == Approx(-1.0).margin(1e-10));
        CHECK(bb.max.Y() == Approx(1.0).margin(1e-10));
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

// ---------------------------------------------------------------------------
// GetLineSegments: multi-level hierarchy traversal tests
//
// These exercise the from_global coordinate-frame fix in
// LogicalVolume::GetLineSegments.  The BVH inside each LogicalVolume stores
// daughter bounding-boxes in that volume's *local* frame; before the fix the
// ray was searched in world-frame coordinates, causing misses for volumes
// placed far from the world origin.
// ---------------------------------------------------------------------------

TEST_CASE("GetLineSegments: translated daughter returns correct material and length", "[BVH]") {
    // World: 10x10x10 Rock box  (world z in [-5, 5])
    // Child: 2x2x2  Water box placed at z=3  (world z in [2, 4])
    // Ray: (0,0,-5) → +z
    // Expected:  Rock(7) | Water(2) | Rock(1)
    auto rock = make_named_mat("Rock");
    auto water = make_named_mat("Water");

    auto child_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    auto child_lv = std::make_shared<NuGeom::LogicalVolume>(water, child_shape);
    NuGeom::Translation3D child_trans(0, 0, 3);
    auto child_pv =
        std::make_shared<NuGeom::PhysicalVolume>(child_lv, child_trans, NuGeom::Transform3D{});

    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{10, 10, 10});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>(rock, world_shape);
    world_lv->AddDaughter(child_pv);

    NuGeom::World world(world_lv);
    NuGeom::Ray ray{{0, 0, -5}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    REQUIRE(segs.size() == 3);
    CHECK(segs[0].GetMaterial().Name() == "Rock");
    CHECK(segs[1].GetMaterial().Name() == "Water");
    CHECK(segs[2].GetMaterial().Name() == "Rock");
    CHECK_THAT(segs[0].Length(), Catch::WithinAbs(7.0, 1e-4));
    CHECK_THAT(segs[1].Length(), Catch::WithinAbs(2.0, 1e-4));
    CHECK_THAT(segs[2].Length(), Catch::WithinAbs(1.0, 1e-4));
    CHECK_THAT(segs[0].Start().Z(), Catch::WithinAbs(-5.0, 1e-4));
    CHECK_THAT(segs[2].End().Z(), Catch::WithinAbs(5.0, 1e-4));
}

TEST_CASE("GetLineSegments: 3-level hierarchy with large translation (from_global regression)",
          "[BVH]") {
    // World:  200x200x200 Rock     (z in [-100, 100])
    // Middle:  20x20x20  Concrete  placed at z=70 in world  (world z in [60, 80])
    // Inner:    4x4x4   Argon     placed at z=2 in Middle   (Middle-local z in [0,4],
    //                                                         world z in [70, 74])
    //
    // Inner has NO physical mother: when it exits, it calls
    //   Middle's LogicalVolume::GetLineSegments with from_global = middle_m_transform.
    // Without the from_global fix the BVH in Middle would be searched with a
    // world-frame ray (z≈74), while the BVH stores boxes in Middle-local frame
    // (z in [0,4]) → wrong exit time → wrong Concrete trailing segment length.
    //
    // Expected segments (world z):
    //   Rock      -100 → 60    len = 160
    //   Concrete   60 → 70    len =  10
    //   Argon      70 → 74    len =   4
    //   Concrete   74 → 80    len =   6   ← fails without fix
    auto rock = make_named_mat("Rock");
    auto concrete = make_named_mat("Concrete");
    auto argon = make_named_mat("Argon");

    auto inner_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 4});
    auto inner_lv = std::make_shared<NuGeom::LogicalVolume>(argon, inner_shape);
    NuGeom::Translation3D inner_trans(0, 0, 2); // in Middle-local frame
    auto inner_pv =
        std::make_shared<NuGeom::PhysicalVolume>(inner_lv, inner_trans, NuGeom::Transform3D{});

    auto middle_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{20, 20, 20});
    auto middle_lv = std::make_shared<NuGeom::LogicalVolume>(concrete, middle_shape);
    middle_lv->AddDaughter(inner_pv);
    // inner_pv's physical mother is assigned by World's BuildUniqueTree

    NuGeom::Translation3D middle_trans(0, 0, 70); // in world frame
    auto middle_pv =
        std::make_shared<NuGeom::PhysicalVolume>(middle_lv, middle_trans, NuGeom::Transform3D{});

    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{200, 200, 200});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>(rock, world_shape);
    world_lv->AddDaughter(middle_pv);

    NuGeom::World world(world_lv);
    NuGeom::Ray ray{{0, 0, -100}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    // After physical-mother fix the traversal correctly unwinds back to world
    // and emits the trailing Rock segment (z 80→100).
    REQUIRE(segs.size() == 5);
    CHECK(segs[0].GetMaterial().Name() == "Rock");
    CHECK(segs[1].GetMaterial().Name() == "Concrete");
    CHECK(segs[2].GetMaterial().Name() == "Argon");
    CHECK(segs[3].GetMaterial().Name() == "Concrete");
    CHECK(segs[4].GetMaterial().Name() == "Rock");
    CHECK_THAT(segs[0].Length(), Catch::WithinAbs(160.0, 1e-4));
    CHECK_THAT(segs[1].Length(), Catch::WithinAbs(10.0, 1e-4));
    CHECK_THAT(segs[2].Length(), Catch::WithinAbs(4.0, 1e-4));
    CHECK_THAT(segs[3].Length(), Catch::WithinAbs(6.0, 1e-4));
    CHECK_THAT(segs[4].Length(), Catch::WithinAbs(20.0, 1e-4));
}

TEST_CASE("GetLineSegments: ray traverses two sibling daughters", "[BVH]") {
    // World: 20x20x20 Rock
    // SibA:  2x2x2 Water  placed at z=-3  (world z in [-4, -2])
    // SibB:  2x2x2 Air    placed at z= 3  (world z in [ 2,  4])
    // Ray: (0,0,-10) → +z
    // Expected: Rock(6) | Water(2) | Rock(4) | Air(2) | Rock(6)
    auto rock = make_named_mat("Rock");
    auto water = make_named_mat("Water");
    auto air = make_named_mat("Air");

    auto make_daughter = [](const NuGeom::Material &mat,
                            double tz) -> std::pair<std::shared_ptr<NuGeom::LogicalVolume>,
                                                    std::shared_ptr<NuGeom::PhysicalVolume>> {
        auto shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
        auto lv = std::make_shared<NuGeom::LogicalVolume>(mat, shape);
        NuGeom::Translation3D trans(0, 0, tz);
        auto pv = std::make_shared<NuGeom::PhysicalVolume>(lv, trans, NuGeom::Transform3D{});
        return {lv, pv};
    };

    auto [sib_a_lv, sib_a_pv] = make_daughter(water, -3);
    auto [sib_b_lv, sib_b_pv] = make_daughter(air, 3);

    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{20, 20, 20});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>(rock, world_shape);
    world_lv->AddDaughter(sib_a_pv);
    world_lv->AddDaughter(sib_b_pv);

    NuGeom::World world(world_lv);
    NuGeom::Ray ray{{0, 0, -10}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    REQUIRE(segs.size() == 5);
    CHECK(segs[0].GetMaterial().Name() == "Rock");
    CHECK(segs[1].GetMaterial().Name() == "Water");
    CHECK(segs[2].GetMaterial().Name() == "Rock");
    CHECK(segs[3].GetMaterial().Name() == "Air");
    CHECK(segs[4].GetMaterial().Name() == "Rock");
    CHECK_THAT(segs[0].Length(), Catch::WithinAbs(6.0, 1e-4));
    CHECK_THAT(segs[1].Length(), Catch::WithinAbs(2.0, 1e-4));
    CHECK_THAT(segs[2].Length(), Catch::WithinAbs(4.0, 1e-4));
    CHECK_THAT(segs[3].Length(), Catch::WithinAbs(2.0, 1e-4));
    CHECK_THAT(segs[4].Length(), Catch::WithinAbs(6.0, 1e-4));
}

// ---------------------------------------------------------------------------
// Regression tests for adjacent-daughter boundary bug
//
// When two daughters touch (exit of A == entry of B), after exiting A the
// next call to LogicalVolume::GetLineSegments launches a shift_ray that
// starts epsilon *inside* B.  Before the fix, PhysicalVolume::Intersect
// returned B's EXIT time (origin inside → t2), which caused:
//   1. A spurious world-material segment stretching the full length of B.
//   2. B->GetLineSegments called from B's exit → Shape::Intersect = inf.
// After the fix, Intersect returns 0.0 when the origin is inside, so the
// BVH immediately enters B and produces correct segment lengths.
// ---------------------------------------------------------------------------

TEST_CASE("GetLineSegments: two adjacent daughters (touching boundary)", "[BVH]") {
    // World: 30 cm Rock box (z in [-15, 15])
    // Water box: z in [0, 5]   (center z=2.5, halfz=2.5)
    // Air   box: z in [5, 10]  (center z=7.5, halfz=2.5) — touches Water exactly
    // Ray: (0,0,-15) → +z
    // After algorithm fix + pruning: Rock(15) | Water(5) | Air(5) | Rock(5) = 4 segments.
    // No tiny Rock between Water and Air.
    auto rock = make_named_mat("Rock");
    auto water = make_named_mat("Water");
    auto air = make_named_mat("Air");

    auto make_pv = [](const NuGeom::Material &mat, double cz,
                      double halfz) -> std::pair<std::shared_ptr<NuGeom::LogicalVolume>,
                                                 std::shared_ptr<NuGeom::PhysicalVolume>> {
        // Box() takes full extents; store halfz as half-extent → full extent = halfz*2.
        auto shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{5, 5, halfz * 2.0});
        auto lv = std::make_shared<NuGeom::LogicalVolume>(mat, shape);
        NuGeom::Translation3D trans(0, 0, cz);
        auto pv = std::make_shared<NuGeom::PhysicalVolume>(lv, trans, NuGeom::Transform3D{});
        return {lv, pv};
    };

    auto [lv_water, pv_water] = make_pv(water, 2.5, 2.5); // z in [0,  5]
    auto [lv_air, pv_air] = make_pv(air, 7.5, 2.5);       // z in [5, 10]

    // Box({5,5,30}) → half-extent=15 → world z ∈ [-15, 15] (30 cm total).
    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{5, 5, 30});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>(rock, world_shape);
    world_lv->AddDaughter(pv_water);
    world_lv->AddDaughter(pv_air);

    NuGeom::World world(world_lv);
    NuGeom::Ray ray{{0, 0, -15}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    // All segments must be finite and non-negative.
    for(const auto &s : segs) {
        INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
        CHECK(std::isfinite(s.Length()));
        CHECK(s.Length() >= 0.0);
    }

    // Total path length must equal world diameter (30 cm).
    double total = 0.0;
    for(const auto &s : segs) total += s.Length();
    CHECK_THAT(total, Catch::WithinAbs(30.0, 1e-3));

    // Water and Air must each appear with the correct length.
    double water_len = 0.0, air_len = 0.0;
    for(const auto &s : segs) {
        if(s.GetMaterial().Name() == "Water") water_len += s.Length();
        if(s.GetMaterial().Name() == "Air") air_len += s.Length();
    }
    CHECK_THAT(water_len, Catch::WithinAbs(5.0, 1e-3));
    CHECK_THAT(air_len, Catch::WithinAbs(5.0, 1e-3));
}

// ---------------------------------------------------------------------------
// Shared-LogicalVolume regression test
//
// When the same LogicalVolume is placed at two positions via two different
// PhysicalVolumes, parent pointers must be unique per-placement.  Before
// the fix, LogicalVolume::m_mother and daughter PVs' m_mother got
// overwritten by the last AssignPhysicalMothers call, causing rays exiting
// daughters of the first placement to unwind to the second placement's PV.
// ---------------------------------------------------------------------------

TEST_CASE("GetLineSegments: shared LogicalVolume at two positions", "[BVH]") {
    // World:  100x100x100  Rock  (x in [-50, 50])
    // LV_Inner: 10x10x10  Water, contains:
    //   LV_Core:  4x4x4   Iron  (at origin of LV_Inner)
    //
    // PV_Left:  LV_Inner placed at x=-20  (world x in [-25, -15])
    // PV_Right: LV_Inner placed at x=+20  (world x in [ 15,  25])
    //
    // Ray along +x from x=-50: crosses Rock -> Left(Water->Iron->Water) -> Rock ->
    // Right(Water->Iron->Water) -> Rock
    auto rock = make_named_mat("Rock");
    auto water = make_named_mat("Water");
    auto iron = make_named_mat("Iron");

    // Core: 4x4x4 Iron at origin of Inner
    auto core_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 4});
    auto core_lv = std::make_shared<NuGeom::LogicalVolume>("Core", iron, core_shape);
    auto core_pv = std::make_shared<NuGeom::PhysicalVolume>("Core", core_lv, NuGeom::Transform3D{},
                                                            NuGeom::Transform3D{});

    // Inner: 10x10x10 Water with Core daughter -- this LV is SHARED
    auto inner_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{10, 10, 10});
    auto inner_lv = std::make_shared<NuGeom::LogicalVolume>("Inner", water, inner_shape);
    inner_lv->AddDaughter(core_pv);

    // Place Inner at x=-20 and x=+20
    auto pv_left = std::make_shared<NuGeom::PhysicalVolume>(
        "Left", inner_lv, NuGeom::Translation3D(-20, 0, 0), NuGeom::Transform3D{});
    auto pv_right = std::make_shared<NuGeom::PhysicalVolume>(
        "Right", inner_lv, NuGeom::Translation3D(+20, 0, 0), NuGeom::Transform3D{});

    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{100, 100, 100});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>("World", rock, world_shape);
    world_lv->AddDaughter(pv_left);
    world_lv->AddDaughter(pv_right);

    NuGeom::World world(world_lv);

    SECTION("Ray through both placements along +x") {
        NuGeom::Ray ray{{-50, 0, 0}, {1, 0, 0}, 1.0};
        auto segs = world.GetLineSegments(ray);

        for(const auto &s : segs) {
            INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
            CHECK(std::isfinite(s.Length()));
            CHECK(s.Length() >= 0.0);
        }

        double total = 0.0;
        for(const auto &s : segs) total += s.Length();
        CHECK_THAT(total, Catch::WithinAbs(100.0, 1e-3));

        // Ray crosses both placements: 2 * (Water=6 + Iron=4)
        double water_len = 0, iron_len = 0;
        for(const auto &s : segs) {
            if(s.GetMaterial().Name() == "Water") water_len += s.Length();
            if(s.GetMaterial().Name() == "Iron") iron_len += s.Length();
        }
        CHECK_THAT(water_len, Catch::WithinAbs(12.0, 1e-3));
        CHECK_THAT(iron_len, Catch::WithinAbs(8.0, 1e-3));
    }

    SECTION("Ray through both placements along -x") {
        NuGeom::Ray ray{{50, 0, 0}, {-1, 0, 0}, 1.0};
        auto segs = world.GetLineSegments(ray);

        for(const auto &s : segs) {
            INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
            CHECK(std::isfinite(s.Length()));
            CHECK(s.Length() >= 0.0);
        }

        double total = 0.0;
        for(const auto &s : segs) total += s.Length();
        CHECK_THAT(total, Catch::WithinAbs(100.0, 1e-3));

        double water_len = 0, iron_len = 0;
        for(const auto &s : segs) {
            if(s.GetMaterial().Name() == "Water") water_len += s.Length();
            if(s.GetMaterial().Name() == "Iron") iron_len += s.Length();
        }
        CHECK_THAT(water_len, Catch::WithinAbs(12.0, 1e-3));
        CHECK_THAT(iron_len, Catch::WithinAbs(8.0, 1e-3));
    }
}

TEST_CASE("GetLineSegments: 10 adjacent daughters all finite", "[BVH]") {
    // 10 boxes each 2 cm thick, stacked along z: [0,2], [2,4], ..., [18,20].
    // All are touching: exit of slab i == entry of slab i+1.
    // World: 100x5x5 box (z in [-50, 50]).  Ray: (0,0,-50) → +z.
    // Every slab and gap segment must be finite; total = 100 cm.
    auto world_mat = make_named_mat("World");
    // Box({5,5,100}) → half-extent=50 → world z ∈ [-50, 50] (100 cm total).
    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{5, 5, 100});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>(world_mat, world_shape);

    constexpr int N = 10;
    for(int i = 0; i < N; ++i) {
        auto mat = make_named_mat("Slab" + std::to_string(i));
        // Box({4,4,2}) → half-extent z=1 → local z ∈ [-1,1]; full slab = 2 cm.
        auto shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 2}); // halfz=1 → [-1,1]
        auto lv = std::make_shared<NuGeom::LogicalVolume>(mat, shape);
        double cz = 1.0 + 2.0 * i; // centres at z = 1, 3, 5, ..., 19
        NuGeom::Translation3D trans(0, 0, cz);
        auto pv = std::make_shared<NuGeom::PhysicalVolume>(lv, trans, NuGeom::Transform3D{});
        world_lv->AddDaughter(pv);
    }

    NuGeom::World world(world_lv);
    NuGeom::Ray ray{{0, 0, -50}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    for(const auto &s : segs) {
        INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
        CHECK(std::isfinite(s.Length()));
        CHECK(s.Length() >= 0.0);
    }

    double total = 0.0;
    for(const auto &s : segs) total += s.Length();
    CHECK_THAT(total, Catch::WithinAbs(100.0, 0.01));

    // Each slab must contribute 2 cm total (may be split into entry + body due to eps).
    for(int i = 0; i < N; ++i) {
        std::string name = "Slab" + std::to_string(i);
        double slab_len = 0.0;
        for(const auto &s : segs)
            if(s.GetMaterial().Name() == name) slab_len += s.Length();
        CHECK_THAT(slab_len, Catch::WithinAbs(2.0, 1e-3));
    }
}

// ---------------------------------------------------------------------------
// Edge-case tests for world-boundary daughters
//
// These expose two related bugs in LogicalVolume::GetLineSegments:
//
//  Bug A — "infinite last segment":
//    When the last daughter ends exactly at the world exit face, the traversal
//    re-enters world::GetLineSegments from a position at or past the world
//    boundary.  m_shape->Intersect() returns +∞ (both t1 and t2 ≤ 0).
//    There is no guard, so an infinite-length segment is emitted.
//
//  Bug B — "leading zero-length segment":
//    When the first daughter starts exactly at the world entry face, the
//    shift_ray (origin + eps) is already inside the daughter.  BVH returns
//    t = 0 for that daughter, so time += eps = eps.  A world-material segment
//    of length ≈ eps is emitted before the first daughter.
//
// Both bugs manifest together when daughters completely fill the world along
// the ray axis (e.g. two daughters [0,5] and [5,10] in a world [-5,5]).
// ---------------------------------------------------------------------------

TEST_CASE("GetLineSegments: last daughter ends at world exit - no infinite segment", "[BVH]") {
    // World: z ∈ [-10, 10]  (Box({5,5,20}) → halfz=10)
    // Water: center=9, halfz=1 → Box({4,4,2}), world z ∈ [8, 10]
    // Ray: (0,0,-10) → +z
    // Expected: Rock(18) + Water(2), total=20, ALL finite.
    auto rock = make_named_mat("Rock");
    auto water = make_named_mat("Water");

    auto shape_w = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 2}); // halfz=1
    auto lv_w = std::make_shared<NuGeom::LogicalVolume>(water, shape_w);
    NuGeom::Translation3D trans_w(0, 0, 9);
    auto pv_w = std::make_shared<NuGeom::PhysicalVolume>(lv_w, trans_w, NuGeom::Transform3D{});

    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{5, 5, 20}); // halfz=10
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>(rock, world_shape);
    world_lv->AddDaughter(pv_w);

    NuGeom::World world(world_lv);
    NuGeom::Ray ray{{0, 0, -10}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    for(const auto &s : segs) {
        INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
        CHECK(std::isfinite(s.Length()));
        CHECK(s.Length() >= 0.0);
    }

    double total = 0.0;
    for(const auto &s : segs) total += s.Length();
    CHECK_THAT(total, Catch::WithinAbs(20.0, 1e-3));

    double water_len = 0.0;
    for(const auto &s : segs)
        if(s.GetMaterial().Name() == "Water") water_len += s.Length();
    CHECK_THAT(water_len, Catch::WithinAbs(2.0, 1e-3));
}

// ---------------------------------------------------------------------------
// Physical-mother / early-exit regression tests
//
// PhysicalVolume::GetLineSegments, when a nested PV exits and has no more
// siblings, it calls its logical mother's GetLineSegments.  That LV emits
// one more segment and hits "if(!pvol) return" — it has no reference to its
// own physical parent, so the traversal terminates there.  Every volume
// *after* the exited nested region at the grandparent level and above is
// silently dropped.
//
// Fix: the World constructor walks the PV tree and assigns physical mothers
// correctly unwind back through the full physical hierarchy.
// ---------------------------------------------------------------------------

TEST_CASE("GetLineSegments: 3-level deep hierarchy exits world completely", "[BVH]") {
    // World:  200x200x200 Rock  (z in [-100, 100])
    // Hall:   100x100x100 Air   placed at world origin (z in [-50, 50])
    // Cryo:    60x60x60  SSteel placed at hall origin  (world z in [-30, 30])
    // LAr:     40x40x40  LAr   placed at cryo origin   (world z in [-20, 20])
    // Ray: (0,0,-100) → +z
    //
    // Expected after fix and pruning:
    //   Rock(50) | Air(20) | SSteel(10) | LAr(40) | SSteel(10) | Air(20) | Rock(50)
    //   total = 200 cm
    //
    // Without fix: traversal stops after LAr exits — logical cryostat LV emits
    // the outbound SSteel(10) and returns.  Outbound Air(20) and Rock(50) are
    // never produced.  Total = 50+20+10+40+10 = 130 cm (missing 70 cm).
    auto rock = make_named_mat("Rock");
    auto air = make_named_mat("Air");
    auto ssteel = make_named_mat("SSteel");
    auto lar = make_named_mat("LAr");

    // Build from innermost outward.
    auto lar_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{40, 40, 40});
    auto lar_lv = std::make_shared<NuGeom::LogicalVolume>(lar, lar_shape);
    auto lar_pv = std::make_shared<NuGeom::PhysicalVolume>(lar_lv, NuGeom::Transform3D{},
                                                           NuGeom::Transform3D{});

    auto cryo_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{60, 60, 60});
    auto cryo_lv = std::make_shared<NuGeom::LogicalVolume>(ssteel, cryo_shape);
    cryo_lv->AddDaughter(lar_pv);
    auto cryo_pv = std::make_shared<NuGeom::PhysicalVolume>(cryo_lv, NuGeom::Transform3D{},
                                                            NuGeom::Transform3D{});

    auto hall_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{100, 100, 100});
    auto hall_lv = std::make_shared<NuGeom::LogicalVolume>(air, hall_shape);
    hall_lv->AddDaughter(cryo_pv);
    auto hall_pv = std::make_shared<NuGeom::PhysicalVolume>(hall_lv, NuGeom::Transform3D{},
                                                            NuGeom::Transform3D{});

    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{200, 200, 200});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>(rock, world_shape);
    world_lv->AddDaughter(hall_pv);

    // World constructor must assign physical mothers for the fix to take effect.
    NuGeom::World world(world_lv);
    NuGeom::Ray ray{{0, 0, -100}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    for(const auto &s : segs) {
        INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
        CHECK(std::isfinite(s.Length()));
        CHECK(s.Length() >= 0.0);
    }

    double total = 0.0;
    for(const auto &s : segs) total += s.Length();
    CHECK_THAT(total, Catch::WithinAbs(200.0, 1e-3));

    double rock_len = 0, air_len = 0, steel_len = 0, lar_len = 0;
    for(const auto &s : segs) {
        if(s.GetMaterial().Name() == "Rock") rock_len += s.Length();
        if(s.GetMaterial().Name() == "Air") air_len += s.Length();
        if(s.GetMaterial().Name() == "SSteel") steel_len += s.Length();
        if(s.GetMaterial().Name() == "LAr") lar_len += s.Length();
    }
    CHECK_THAT(rock_len, Catch::WithinAbs(100.0, 1e-3)); // 50 in + 50 out
    CHECK_THAT(air_len, Catch::WithinAbs(40.0, 1e-3));   // 20 in + 20 out
    CHECK_THAT(steel_len, Catch::WithinAbs(20.0, 1e-3)); // 10 in + 10 out
    CHECK_THAT(lar_len, Catch::WithinAbs(40.0, 1e-3));
}

TEST_CASE("GetLineSegments: two daughters exactly filling world - all finite", "[BVH]") {
    // World: z ∈ [-5, 5]  (Box({5,5,10}) → halfz=5)
    // Water: center=-2.5, halfz=2.5 → Box({4,4,5}), world z ∈ [-5, 0]  (starts at world entry)
    // Air:   center= 2.5, halfz=2.5 → Box({4,4,5}), world z ∈ [ 0, 5]  (ends at world exit)
    // Ray: (0,0,-5) → +z
    // Expected: Water(5) + Air(5) (with possible tiny Rock eps gaps), total=10, ALL finite.
    auto rock = make_named_mat("Rock");
    auto water = make_named_mat("Water");
    auto air = make_named_mat("Air");

    auto mk = [](const NuGeom::Material &mat, double cz) {
        auto sh = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 5}); // halfz=2.5
        auto lv = std::make_shared<NuGeom::LogicalVolume>(mat, sh);
        NuGeom::Translation3D t(0, 0, cz);
        auto pv = std::make_shared<NuGeom::PhysicalVolume>(lv, t, NuGeom::Transform3D{});
        return std::make_pair(lv, pv);
    };

    auto [lv_w, pv_w] = mk(water, -2.5); // z ∈ [-5, 0]
    auto [lv_a, pv_a] = mk(air, 2.5);    // z ∈ [ 0, 5]

    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{5, 5, 10}); // halfz=5
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>(rock, world_shape);
    world_lv->AddDaughter(pv_w);
    world_lv->AddDaughter(pv_a);

    NuGeom::World world(world_lv);
    NuGeom::Ray ray{{0, 0, -5}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    // Critical: no infinite lengths.
    for(const auto &s : segs) {
        INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
        CHECK(std::isfinite(s.Length()));
        CHECK(s.Length() >= 0.0);
    }

    double total = 0.0;
    for(const auto &s : segs) total += s.Length();
    CHECK_THAT(total, Catch::WithinAbs(10.0, 1e-3));

    double water_len = 0.0, air_len = 0.0;
    for(const auto &s : segs) {
        if(s.GetMaterial().Name() == "Water") water_len += s.Length();
        if(s.GetMaterial().Name() == "Air") air_len += s.Length();
    }
    CHECK_THAT(water_len, Catch::WithinAbs(5.0, 1e-3));
    CHECK_THAT(air_len, Catch::WithinAbs(5.0, 1e-3));
}

// ---------------------------------------------------------------------------
// Daughter ending exactly at parent boundary — traversal must not stop early.
// Before the fix, Shape::Intersect could return ∞ after the eps nudge
// landed outside the parent shape, silently terminating traversal.
// ---------------------------------------------------------------------------
TEST_CASE("GetLineSegments: daughter flush with parent boundary - no early exit", "[BVH]") {
    // World: 10x10x10 box, material Rock
    auto world_mat = make_named_mat("Rock");
    auto world_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{10, 10, 10});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>("world", world_mat, world_box);

    // Outer: 4x4x4 box, material Water, placed at z=+3 (extends z=1..5 = world edge)
    auto outer_mat = make_named_mat("Water");
    auto outer_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 4});
    auto outer_lv = std::make_shared<NuGeom::LogicalVolume>("outer", outer_mat, outer_box);
    NuGeom::Translation3D outer_pos(0, 0, 3);
    auto outer_pv = std::make_shared<NuGeom::PhysicalVolume>("outer_pv", outer_lv, outer_pos,
                                                             NuGeom::Transform3D{});

    // Inner: 2x2x2 box, material Iron, placed at z=+1 inside outer (z=4..5 = outer edge = world
    // edge)
    auto inner_mat = make_named_mat("Iron");
    auto inner_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    auto inner_lv = std::make_shared<NuGeom::LogicalVolume>("inner", inner_mat, inner_box);
    NuGeom::Translation3D inner_pos(0, 0, 1);
    auto inner_pv = std::make_shared<NuGeom::PhysicalVolume>("inner_pv", inner_lv, inner_pos,
                                                             NuGeom::Transform3D{});
    outer_lv->AddDaughter(inner_pv);

    world_lv->AddDaughter(outer_pv);

    NuGeom::World world(world_lv);

    // Ray along z from -5 to +5.
    // World: z=-5..5.  Outer at z=3: z=1..5.  Inner at z=1 inside outer: local z=0..2 -> world
    // z=4..5. Expected: Rock z=-5..1, Water z=1..4, Iron z=4..5 (flush with world edge).
    NuGeom::Ray ray{{0, 0, -5}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    // All segments must be finite.
    for(const auto &s : segs) {
        INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
        CHECK(std::isfinite(s.Length()));
        CHECK(s.Length() >= 0.0);
    }

    // Total path must cover the full world (10 cm).
    double total = 0;
    for(const auto &s : segs) total += s.Length();
    CHECK_THAT(total, Catch::WithinAbs(10.0, 1e-3));

    // The last segment's end must be at z=+5 (world boundary), not inside.
    auto last_end = segs.back().End();
    CHECK_THAT(last_end.Z(), Catch::WithinAbs(5.0, 1e-3));
}

// ---------------------------------------------------------------------------
// Nested daughter flush with BOTH parent and grandparent boundaries.
// Tests multi-level unwinding when Intersect returns inf at each level.
// ---------------------------------------------------------------------------
TEST_CASE("GetLineSegments: nested daughters flush at boundary - multi-level unwind", "[BVH]") {
    auto rock_mat = make_named_mat("Rock");
    auto water_mat = make_named_mat("Water");
    auto iron_mat = make_named_mat("Iron");
    auto lead_mat = make_named_mat("Lead");

    // World: z=-5..5
    auto world_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{10, 10, 10});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>("world", rock_mat, world_box);

    // A: 6x6x7 placed at z=+1.5 -> z=-2..5 (flush with world +z boundary)
    auto a_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{6, 6, 7});
    auto a_lv = std::make_shared<NuGeom::LogicalVolume>("A", water_mat, a_box);
    NuGeom::Translation3D a_pos(0, 0, 1.5);
    auto a_pv =
        std::make_shared<NuGeom::PhysicalVolume>("A_pv", a_lv, a_pos, NuGeom::Transform3D{});

    // B: 4x4x4 placed at z=+1.5 inside A -> A-local z=-0.5..3.5 -> world z=3..5
    // (flush with A +z and world +z)
    auto b_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 4});
    auto b_lv = std::make_shared<NuGeom::LogicalVolume>("B", iron_mat, b_box);
    NuGeom::Translation3D b_pos(0, 0, 1.5);
    auto b_pv =
        std::make_shared<NuGeom::PhysicalVolume>("B_pv", b_lv, b_pos, NuGeom::Transform3D{});
    a_lv->AddDaughter(b_pv);

    // C: 2x2x2 placed at z=+1 inside B -> B-local z=0..2 -> world z=4..5
    // (flush with B, A, and world +z)
    auto c_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    auto c_lv = std::make_shared<NuGeom::LogicalVolume>("C", lead_mat, c_box);
    NuGeom::Translation3D c_pos(0, 0, 1);
    auto c_pv =
        std::make_shared<NuGeom::PhysicalVolume>("C_pv", c_lv, c_pos, NuGeom::Transform3D{});
    b_lv->AddDaughter(c_pv);

    world_lv->AddDaughter(a_pv);

    NuGeom::World world(world_lv);

    // Ray along +z from z=-5 to z=+5 (10 cm).
    // Expected: Rock z=-5..-2, Water z=-2..3, Iron z=3..4, Lead z=4..5 (flush).
    NuGeom::Ray ray{{0, 0, -5}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    for(const auto &s : segs) {
        INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
        CHECK(std::isfinite(s.Length()));
    }

    double total = 0;
    for(const auto &s : segs) total += s.Length();
    CHECK_THAT(total, Catch::WithinAbs(10.0, 1e-3));

    // Last segment must reach the world boundary.
    CHECK_THAT(segs.back().End().Z(), Catch::WithinAbs(5.0, 1e-3));
}

// ---------------------------------------------------------------------------
// Daughter flush with parent boundary in the ray direction, FOLLOWED by a
// sibling.  When the ray exits the inner daughter, the eps nudge lands
// outside the outer parent's shape, and Intersect returns infinity.
// Without the fix, the traversal terminates there, losing the sibling.
// ---------------------------------------------------------------------------
TEST_CASE("GetLineSegments: flush daughter followed by sibling - traversal continues", "[BVH]") {
    auto rock_mat = make_named_mat("Rock");
    auto water_mat = make_named_mat("Water");
    auto iron_mat = make_named_mat("Iron");

    // World: z=-10..10  (Box({10,10,20}) -> half-extent 10)
    auto world_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{10, 10, 20});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>("world", rock_mat, world_box);

    // Outer: 4x4x10 placed at z=-5 -> outer-local z=-5..5 -> world z=-10..0
    auto outer_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 10});
    auto outer_lv = std::make_shared<NuGeom::LogicalVolume>("outer", water_mat, outer_box);

    // Inner: 2x2x4 placed at z=+3 inside outer -> outer-local z=1..5 (flush at +z)
    //   -> world z=-4..0 (flush with outer's +z exit boundary in ray direction)
    auto inner_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 4});
    auto inner_lv = std::make_shared<NuGeom::LogicalVolume>("inner", iron_mat, inner_box);
    NuGeom::Translation3D inner_pos(0, 0, 3);
    auto inner_pv = std::make_shared<NuGeom::PhysicalVolume>("inner_pv", inner_lv, inner_pos,
                                                             NuGeom::Transform3D{});
    outer_lv->AddDaughter(inner_pv);

    NuGeom::Translation3D outer_pos(0, 0, -5);
    auto outer_pv = std::make_shared<NuGeom::PhysicalVolume>("outer_pv", outer_lv, outer_pos,
                                                             NuGeom::Transform3D{});
    world_lv->AddDaughter(outer_pv);

    // Sibling: 4x4x4 at z=+5 -> world z=3..7
    auto sib_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 4});
    auto sib_lv = std::make_shared<NuGeom::LogicalVolume>("sibling", iron_mat, sib_box);
    NuGeom::Translation3D sib_pos(0, 0, 5);
    auto sib_pv =
        std::make_shared<NuGeom::PhysicalVolume>("sib_pv", sib_lv, sib_pos, NuGeom::Transform3D{});
    world_lv->AddDaughter(sib_pv);

    NuGeom::World world(world_lv);

    // Ray along +z.  Path:
    //   Rock(-10..-4), Water(-4..-4), wait...
    //   Actually: Rock(-10), outer entry at -10.
    //   outer-local: ray enters at z=-5, hits inner at z=1, exits inner at z=5 (= outer boundary).
    //   world:  Water(-10..-4), Iron(-4..0), exits outer, Rock(0..3), Iron_sib(3..7), Rock(7..10)
    NuGeom::Ray ray{{0, 0, -10}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    for(const auto &s : segs) {
        INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
        CHECK(std::isfinite(s.Length()));
    }

    double total = 0;
    for(const auto &s : segs) total += s.Length();
    CHECK_THAT(total, Catch::WithinAbs(20.0, 1e-3));

    // Must find the sibling at z>2.
    bool found_sibling = false;
    for(const auto &s : segs) {
        if(s.Start().Z() > 2.0 && s.GetMaterial().Name() == "Iron") found_sibling = true;
    }
    CHECK(found_sibling);
}

// ---------------------------------------------------------------------------
// Exact-fill daughter regression tests
//
// When a daughter has identical dimensions and position as its parent
// (daughter fills parent exactly), the iterative traversal previously
// entered an infinite loop:
//   1. Parent finds daughter at daughter_time=0 (ray inside both)
//   2. Parent emits eps segment, enters daughter
//   3. Daughter's own eps shift compounds to 2*eps, pushing past the
//      identical boundary
//   4. Daughter returns 0 segments → parent doesn't advance → infinite loop
//
// This pattern occurs in DUNE GDML: volHalfDetector → volFieldcage are
// identical boxes at position (0,0,0).
//
// Fix: When daughter_time == 0, enter daughter directly without eps advance.
// If daughter returns 0 segments, compute its exit via Intersect2 and skip.
// ---------------------------------------------------------------------------

TEST_CASE("GetLineSegments: daughter exactly fills parent - no infinite loop", "[BVH]") {
    // Mimics the DUNE pattern: HalfDetector → Fieldcage → LAr
    //
    // World:  100x100x100  Rock
    // Outer:   20x20x20   Air    placed at origin (fills nothing special)
    // Inner:   20x20x20   Water  placed at origin INSIDE Outer (EXACTLY fills Outer)
    // Core:    10x10x10   Iron   placed at origin inside Inner (inset, not flush)
    //
    // Ray: (0,0,-50) → +z
    // Expected: Rock(40) | Water(3) | Iron(10) | Water(3) | ???
    //   Wait — Inner fills Outer exactly, so Air never appears along the ray.
    //   Outer is Air but Inner (Water) fills it completely.
    //   Core (Iron) is inset inside Inner.
    //   Path: Rock(40) | Water(5) | Iron(10) | Water(5) | Rock(40) = 100
    auto rock = make_named_mat("Rock");
    auto air = make_named_mat("Air");
    auto water = make_named_mat("Water");
    auto iron = make_named_mat("Iron");

    // Core: 10x10x10 Iron
    auto core_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{10, 10, 10});
    auto core_lv = std::make_shared<NuGeom::LogicalVolume>("Core", iron, core_shape);
    auto core_pv = std::make_shared<NuGeom::PhysicalVolume>(
        "Core_pv", core_lv, NuGeom::Transform3D{}, NuGeom::Transform3D{});

    // Inner: 20x20x20 Water — EXACTLY same size as Outer, at origin
    auto inner_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{20, 20, 20});
    auto inner_lv = std::make_shared<NuGeom::LogicalVolume>("Inner", water, inner_shape);
    inner_lv->AddDaughter(core_pv);
    auto inner_pv = std::make_shared<NuGeom::PhysicalVolume>(
        "Inner_pv", inner_lv, NuGeom::Transform3D{}, NuGeom::Transform3D{});

    // Outer: 20x20x20 Air — daughter Inner fills it exactly
    auto outer_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{20, 20, 20});
    auto outer_lv = std::make_shared<NuGeom::LogicalVolume>("Outer", air, outer_shape);
    outer_lv->AddDaughter(inner_pv);
    auto outer_pv = std::make_shared<NuGeom::PhysicalVolume>(
        "Outer_pv", outer_lv, NuGeom::Transform3D{}, NuGeom::Transform3D{});

    // World: 100x100x100 Rock
    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{100, 100, 100});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>("World", rock, world_shape);
    world_lv->AddDaughter(outer_pv);

    NuGeom::World world(world_lv);
    NuGeom::Ray ray{{0, 0, -50}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    // Must terminate (no infinite loop) and produce finite segments.
    for(const auto &s : segs) {
        INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
        CHECK(std::isfinite(s.Length()));
        CHECK(s.Length() >= 0.0);
    }

    double total = 0.0;
    for(const auto &s : segs) total += s.Length();
    CHECK_THAT(total, Catch::WithinAbs(100.0, 1e-3));

    // Inner fills Outer exactly, so Air should contribute ~0 (it's replaced by Water).
    // Water: 2 * 5 = 10 cm (half-extent 10 minus half-extent 5 on each side)
    // Iron:  10 cm
    double water_len = 0, iron_len = 0;
    for(const auto &s : segs) {
        if(s.GetMaterial().Name() == "Water") water_len += s.Length();
        if(s.GetMaterial().Name() == "Iron") iron_len += s.Length();
    }
    CHECK_THAT(water_len, Catch::WithinAbs(10.0, 1e-3));
    CHECK_THAT(iron_len, Catch::WithinAbs(10.0, 1e-3));
}

TEST_CASE("GetLineSegments: two exact-fill levels deep - no infinite loop", "[BVH]") {
    // Three nested volumes all with identical 20x20x20 boxes at origin.
    // This tests the double-exact-fill pattern (grandparent = parent = child size).
    //
    // World:  100x100x100  Rock
    // Level1:  20x20x20   Mat_A  at origin
    // Level2:  20x20x20   Mat_B  at origin (fills Level1 exactly)
    // Level3:  20x20x20   Mat_C  at origin (fills Level2 exactly)
    //
    // Ray: (0,0,-50) → +z
    // The deepest volume (Mat_C) should claim the full 20 cm.
    // Expected: Rock(40) | Mat_C(20) | Rock(40) = 100
    auto rock = make_named_mat("Rock");
    auto mat_a = make_named_mat("Mat_A");
    auto mat_b = make_named_mat("Mat_B");
    auto mat_c = make_named_mat("Mat_C");

    auto box20 = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{20, 20, 20});

    auto lv3 = std::make_shared<NuGeom::LogicalVolume>("L3", mat_c, box20);
    auto pv3 = std::make_shared<NuGeom::PhysicalVolume>("L3_pv", lv3, NuGeom::Transform3D{},
                                                        NuGeom::Transform3D{});

    auto lv2 = std::make_shared<NuGeom::LogicalVolume>("L2", mat_b, box20);
    lv2->AddDaughter(pv3);
    auto pv2 = std::make_shared<NuGeom::PhysicalVolume>("L2_pv", lv2, NuGeom::Transform3D{},
                                                        NuGeom::Transform3D{});

    auto lv1 = std::make_shared<NuGeom::LogicalVolume>("L1", mat_a, box20);
    lv1->AddDaughter(pv2);
    auto pv1 = std::make_shared<NuGeom::PhysicalVolume>("L1_pv", lv1, NuGeom::Transform3D{},
                                                        NuGeom::Transform3D{});

    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{100, 100, 100});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>("World", rock, world_shape);
    world_lv->AddDaughter(pv1);

    NuGeom::World world(world_lv);
    NuGeom::Ray ray{{0, 0, -50}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    for(const auto &s : segs) {
        INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
        CHECK(std::isfinite(s.Length()));
        CHECK(s.Length() >= 0.0);
    }

    double total = 0.0;
    for(const auto &s : segs) total += s.Length();
    CHECK_THAT(total, Catch::WithinAbs(100.0, 1e-3));

    // The deepest material (Mat_C) should claim the full 20 cm extent.
    double mat_c_len = 0;
    for(const auto &s : segs)
        if(s.GetMaterial().Name() == "Mat_C") mat_c_len += s.Length();
    CHECK_THAT(mat_c_len, Catch::WithinAbs(20.0, 1e-3));
}

TEST_CASE("GetLineSegments: exact-fill with offset ray - no infinite loop", "[BVH]") {
    // Same exact-fill pattern but with a ray that hits the boundary at an angle,
    // testing that the fix works for non-axis-aligned rays too.
    //
    // World:  100x100x100  Rock
    // Outer:   20x20x20   Air    at origin
    // Inner:   20x20x20   Water  at origin (fills Outer exactly)
    //
    // Ray: (3, 4, -50) → +z  (off-center but still inside the 20x20 cross-section)
    auto rock = make_named_mat("Rock");
    auto air = make_named_mat("Air");
    auto water = make_named_mat("Water");

    auto box20 = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{20, 20, 20});

    auto inner_lv = std::make_shared<NuGeom::LogicalVolume>("Inner", water, box20);
    auto inner_pv = std::make_shared<NuGeom::PhysicalVolume>(
        "Inner_pv", inner_lv, NuGeom::Transform3D{}, NuGeom::Transform3D{});

    auto outer_lv = std::make_shared<NuGeom::LogicalVolume>("Outer", air, box20);
    outer_lv->AddDaughter(inner_pv);
    auto outer_pv = std::make_shared<NuGeom::PhysicalVolume>(
        "Outer_pv", outer_lv, NuGeom::Transform3D{}, NuGeom::Transform3D{});

    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{100, 100, 100});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>("World", rock, world_shape);
    world_lv->AddDaughter(outer_pv);

    NuGeom::World world(world_lv);
    NuGeom::Ray ray{{3, 4, -50}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    for(const auto &s : segs) {
        INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
        CHECK(std::isfinite(s.Length()));
        CHECK(s.Length() >= 0.0);
    }

    double total = 0.0;
    for(const auto &s : segs) total += s.Length();
    CHECK_THAT(total, Catch::WithinAbs(100.0, 1e-3));

    // Water should claim the full 20 cm (Inner fills Outer exactly).
    double water_len = 0;
    for(const auto &s : segs)
        if(s.GetMaterial().Name() == "Water") water_len += s.Length();
    CHECK_THAT(water_len, Catch::WithinAbs(20.0, 1e-3));
}

TEST_CASE("GetLineSegments: exact-fill daughter with sibling - both traversed", "[BVH]") {
    // Exact-fill daughter must not prevent traversal of a sibling volume.
    //
    // World:  100x100x100  Rock
    // Outer:   20x20x20   Air    at z=-20  (world z in [-30, -10])
    // Inner:   20x20x20   Water  at origin inside Outer (fills Outer exactly)
    // Sibling: 10x10x10   Iron   at z=+20  (world z in [15, 25])
    //
    // Ray: (0,0,-50) → +z
    // Expected: Rock(20) | Water(20) | Rock(25) | Iron(10) | Rock(25) = 100
    auto rock = make_named_mat("Rock");
    auto air = make_named_mat("Air");
    auto water = make_named_mat("Water");
    auto iron = make_named_mat("Iron");

    auto box20 = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{20, 20, 20});

    // Inner fills Outer exactly
    auto inner_lv = std::make_shared<NuGeom::LogicalVolume>("Inner", water, box20);
    auto inner_pv = std::make_shared<NuGeom::PhysicalVolume>(
        "Inner_pv", inner_lv, NuGeom::Transform3D{}, NuGeom::Transform3D{});

    auto outer_lv = std::make_shared<NuGeom::LogicalVolume>("Outer", air, box20);
    outer_lv->AddDaughter(inner_pv);
    NuGeom::Translation3D outer_pos(0, 0, -20);
    auto outer_pv = std::make_shared<NuGeom::PhysicalVolume>("Outer_pv", outer_lv, outer_pos,
                                                             NuGeom::Transform3D{});

    // Sibling: separate volume
    auto sib_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{10, 10, 10});
    auto sib_lv = std::make_shared<NuGeom::LogicalVolume>("Sibling", iron, sib_shape);
    NuGeom::Translation3D sib_pos(0, 0, 20);
    auto sib_pv =
        std::make_shared<NuGeom::PhysicalVolume>("Sib_pv", sib_lv, sib_pos, NuGeom::Transform3D{});

    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{100, 100, 100});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>("World", rock, world_shape);
    world_lv->AddDaughter(outer_pv);
    world_lv->AddDaughter(sib_pv);

    NuGeom::World world(world_lv);
    NuGeom::Ray ray{{0, 0, -50}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    for(const auto &s : segs) {
        INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
        CHECK(std::isfinite(s.Length()));
        CHECK(s.Length() >= 0.0);
    }

    double total = 0.0;
    for(const auto &s : segs) total += s.Length();
    CHECK_THAT(total, Catch::WithinAbs(100.0, 1e-3));

    // Water from exact-fill: 20 cm
    double water_len = 0;
    for(const auto &s : segs)
        if(s.GetMaterial().Name() == "Water") water_len += s.Length();
    CHECK_THAT(water_len, Catch::WithinAbs(20.0, 1e-3));

    // Sibling Iron must also be traversed: 10 cm
    double iron_len = 0;
    for(const auto &s : segs)
        if(s.GetMaterial().Name() == "Iron") iron_len += s.Length();
    CHECK_THAT(iron_len, Catch::WithinAbs(10.0, 1e-3));
}

// ---------------------------------------------------------------------------
// AABB-hit but shape-miss regression test
//
// When daughter_time == 0 (BVH reports the ray inside the daughter's AABB)
// but the daughter returns 0 segments and Intersect2 returns (inf, inf)
// — meaning the ray grazes the AABB but misses the actual shape — the
// skip-past path must NOT emit an infinite-length segment.
//
// This occurs in practice when a cylinder (or other non-box shape) has an
// AABB that extends to the corners of its bounding box.  A ray passing
// through the corner of the AABB is "inside" the AABB but outside the
// cylinder.
//
// Bug: the guard was `if(dt2 > eps)` without checking `isfinite(dt2)`.
// Fix: changed to `if(std::isfinite(dt2) && dt2 > eps)`.
// ---------------------------------------------------------------------------

TEST_CASE("GetLineSegments: AABB-hit shape-miss daughter - no infinite segment", "[BVH]") {
    // World:  100x100x100  Rock
    // Parent:  20x20x20   Air   at origin
    // Daughter: Cylinder(r=5, h=10) inside Parent at origin
    //   Cylinder AABB = [-5,-5,-10] to [5,5,10] (full height=20)
    //   A ray at x=4.9, y=4.9 is inside the AABB (4.9 < 5) but outside the
    //   cylinder (sqrt(4.9^2 + 4.9^2) ≈ 6.93 > 5).
    //
    // Ray: (4.9, 4.9, -50) → +z
    // The ray is inside the cylinder's AABB when it enters the parent box,
    // so daughter_time == 0.  But the ray misses the actual cylinder shape.
    // The parent (Air) should fill the full 20 cm with no infinite segments.
    auto rock = make_named_mat("Rock");
    auto air = make_named_mat("Air");
    auto water = make_named_mat("Water");

    // Cylinder: r=5, h=10 (half-height) → full height 20, centered at origin
    auto cyl = std::make_shared<NuGeom::Cylinder>(5.0, 10.0);
    auto cyl_lv = std::make_shared<NuGeom::LogicalVolume>("Cyl", water, cyl);
    auto cyl_pv = std::make_shared<NuGeom::PhysicalVolume>("Cyl_pv", cyl_lv, NuGeom::Transform3D{},
                                                           NuGeom::Transform3D{});

    auto parent_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{20, 20, 20});
    auto parent_lv = std::make_shared<NuGeom::LogicalVolume>("Parent", air, parent_shape);
    parent_lv->AddDaughter(cyl_pv);
    auto parent_pv = std::make_shared<NuGeom::PhysicalVolume>(
        "Parent_pv", parent_lv, NuGeom::Transform3D{}, NuGeom::Transform3D{});

    auto world_shape = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{100, 100, 100});
    auto world_lv = std::make_shared<NuGeom::LogicalVolume>("World", rock, world_shape);
    world_lv->AddDaughter(parent_pv);

    NuGeom::World world(world_lv);

    // Ray at (4.9, 4.9) — inside cylinder AABB but outside cylinder shape
    NuGeom::Ray ray{{4.9, 4.9, -50}, {0, 0, 1}, 1.0};
    auto segs = world.GetLineSegments(ray);

    // Critical: no infinite or NaN segments.
    for(const auto &s : segs) {
        INFO("seg " << s.GetMaterial().Name() << " len=" << s.Length());
        CHECK(std::isfinite(s.Length()));
        CHECK(s.Length() >= 0.0);
    }

    double total = 0.0;
    for(const auto &s : segs) total += s.Length();
    CHECK_THAT(total, Catch::WithinAbs(100.0, 1e-3));

    // The ray misses the cylinder entirely — no Water should appear.
    double water_len = 0;
    for(const auto &s : segs)
        if(s.GetMaterial().Name() == "Water") water_len += s.Length();
    CHECK_THAT(water_len, Catch::WithinAbs(0.0, 1e-3));

    // Air (parent material) should fill the full 20 cm of the parent box.
    double air_len = 0;
    for(const auto &s : segs)
        if(s.GetMaterial().Name() == "Air") air_len += s.Length();
    CHECK_THAT(air_len, Catch::WithinAbs(20.0, 1e-3));
}
