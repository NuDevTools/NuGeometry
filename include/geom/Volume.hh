#pragma once

#include "geom/BoundingBox.hh"
#include "geom/Material.hh"
#include "geom/Shape.hh"
#include <set>

namespace NuGeom {

class BVH;
class LineSegment;
class PhysicalVolume;
class Ray;

/// One boundary-crossing event used by the boundary-sweep traversal.
/// delta = +1 when the volume becomes active at t, -1 when it ends.
/// depth = nesting depth (deeper overrides shallower for material).
struct IntervalEvent {
    double t;
    int delta;
    int depth;
    const Material *material;
    const PhysicalVolume *volume; ///< owning volume (for diagnostics)
};

class LogicalVolume {
  public:
    // All constructors and destructor are defined out-of-line in Volume.cc so that
    // the compiler sees the full BVH definition when instantiating unique_ptr<BVH>.
    LogicalVolume();
    LogicalVolume(Material material, std::shared_ptr<Shape> shape);
    LogicalVolume(std::string name, Material material, std::shared_ptr<Shape> shape);
    ~LogicalVolume();

    std::string Name() const { return m_name; }
    const Material &GetMaterial() const { return m_material; }
    void GetMaterials(std::set<Material> &) const;
    Shape *GetShape() const { return m_shape.get(); }
    const std::vector<std::shared_ptr<PhysicalVolume>> &Daughters() const { return m_daughters; }
    void AddDaughter(std::shared_ptr<PhysicalVolume> daughter) { m_daughters.push_back(daughter); }
    double Volume() const;
    double Mass() const;

    bool InWorld(const Vector3D &) const { return true; }
    bool SphereTrace(const Ray &, double &, size_t &, size_t &) const;
    bool RayTrace(const Ray &, double &, std::shared_ptr<PhysicalVolume> &) const;
    void GetLineSegments(const Ray &, std::vector<LineSegment> &,
                         const Transform3D &from_global = {}) const;

  private:
    double DaughterVolumes() const;
    double DaughterMass() const;
    std::pair<double, size_t> GetSDF(const Vector3D &) const { return {0, 0}; }

    std::string m_name;
    Material m_material;
    std::shared_ptr<Shape> m_shape;
    std::vector<std::shared_ptr<PhysicalVolume>> m_daughters;
    static constexpr size_t m_max_steps{512};
    static constexpr double m_epsilon{1e-4};
    mutable std::shared_ptr<BVH> m_bvh;
};

class PhysicalVolume {
  public:
    PhysicalVolume() = default;
    PhysicalVolume(std::shared_ptr<LogicalVolume> volume, Transform3D trans, Transform3D rot)
        : m_volume{std::move(volume)} {
        m_transform = rot * trans.Inverse();
        m_transform.Decompose(m_scale, m_rot, m_trans);
        is_identity = m_transform.IsIdentity();
        is_translation = m_rot.IsIdentity() && !m_trans.IsIdentity();
    }
    PhysicalVolume(std::string name, std::shared_ptr<LogicalVolume> volume, Transform3D trans,
                   Transform3D rot)
        : m_name{std::move(name)}, m_volume{std::move(volume)} {
        m_transform = rot * trans.Inverse();
        m_transform.Decompose(m_scale, m_rot, m_trans);
        is_identity = m_transform.IsIdentity();
        is_translation = m_rot.IsIdentity() && !m_trans.IsIdentity();
    }
    /// Construct from a pre-computed parent→local transform (used by BuildUniqueTree).
    PhysicalVolume(std::string name, std::shared_ptr<LogicalVolume> volume, Transform3D transform,
                   int /*tag*/)
        : m_name{std::move(name)}, m_volume{std::move(volume)}, m_transform{std::move(transform)} {
        m_transform.Decompose(m_scale, m_rot, m_trans);
        is_identity = m_transform.IsIdentity();
        is_translation = m_rot.IsIdentity() && !m_trans.IsIdentity();
    }
    ~PhysicalVolume() = default;

    std::string Name() const { return m_name; }
    const std::shared_ptr<LogicalVolume> &GetLogicalVolume() const { return m_volume; }
    Transform3D GetTransform() const { return m_transform; }
    std::shared_ptr<PhysicalVolume> Mother() const { return m_mother; }
    void SetMother(std::shared_ptr<PhysicalVolume> mother) { m_mother = std::move(mother); }
    const std::vector<std::shared_ptr<PhysicalVolume>> &Daughters() const {
        return m_own_daughters;
    }
    void AddOwnDaughter(std::shared_ptr<PhysicalVolume> d) {
        m_own_daughters.push_back(std::move(d));
        m_bvh.reset();
    }
    /// Replace this volume's daughter list (used by World::Prune).  Drops the
    /// cached BVH, which was built from the old list.
    void SetOwnDaughters(std::vector<std::shared_ptr<PhysicalVolume>> d) {
        m_own_daughters = std::move(d);
        m_bvh.reset();
    }
    double SignedDistance(const Vector3D &in_point) const {
        auto point = TransformPoint(in_point);
        return m_volume->GetShape()->SignedDistance(point);
    }
    double Intersect(const Ray &in_ray) const;
    /// Returns the AABB of this volume in its parent's coordinate frame.
    BoundingBox GetParentBoundingBox() const;
    bool RayTrace(const Ray &ray, double &time, std::shared_ptr<PhysicalVolume> &pvol) const;
    bool RayTrace(const Ray &ray, double &time, PhysicalVolume *&pvol) const;
    void GetLineSegments(const Ray &, std::vector<LineSegment> &, const Transform3D &) const;

    /// Boundary-sweep collection: append this volume's [enter,exit] crossing
    /// events (and those of every hit descendant) for the world-frame ray.
    /// Each volume's interval is clipped to its parent's active window
    /// [parent_lo, parent_hi] so a daughter that protrudes outside its mother
    /// (a real possibility in detector GDMLs) cannot leak its material into the
    /// gap -- this enforces the same hierarchical containment the sequential
    /// traversal applies.
    /// @param from_global  world -> this volume's parent frame transform.
    void CollectIntervals(const Ray &world_ray, const Transform3D &from_global, int depth,
                          double parent_lo, double parent_hi,
                          std::vector<IntervalEvent> &events) const;

  private:
    Vector3D TransformPoint(const Vector3D &point) const { return m_transform.Apply(point); }
    Vector3D TransformPointInverse(const Vector3D &point) const {
        return m_transform.Inverse().Apply(point);
    }
    Ray TransformRay(const Ray &ray) const;
    Ray TransformRayCached(const Ray &ray) const;
    Ray TransformRayInverse(const Ray &ray) const;
    std::string m_name;
    std::shared_ptr<LogicalVolume> m_volume;
    std::shared_ptr<PhysicalVolume> m_mother;
    std::vector<std::shared_ptr<PhysicalVolume>> m_own_daughters;
    mutable std::shared_ptr<BVH> m_bvh;
    Transform3D m_transform;
    Scale3D m_scale;
    Rotation3D m_rot;
    Translation3D m_trans;
    bool is_identity{}, is_translation{};

    // Geant4-style cached parent->local direction transform.  Within a single
    // ray the input (parent-frame) direction is invariant while only the origin
    // advances, so the rotated local direction and its inverse can be reused
    // across the many leaf Intersect() calls that one ray triggers.  Keyed on
    // the input direction; recomputed when it changes (i.e. once per new ray).
    // mutable + single-threaded, matching the m_bvh lazy-build precedent.
    mutable Vector3D m_dir_cache_in{};
    mutable Vector3D m_dir_cache_local{};
    mutable Vector3D m_dir_cache_inv{};
    mutable bool m_dir_cache_valid = false;
};

} // namespace NuGeom
