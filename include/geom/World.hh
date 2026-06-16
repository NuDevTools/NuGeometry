#pragma once

#include "geom/Material.hh"
#include "geom/Shape.hh"
#include "geom/Volume.hh"
#include <memory>
#include <vector>

namespace NuGeom {

class Ray;
class LineSegment;

class World {
  public:
    World() = default;
    World(std::shared_ptr<LogicalVolume> volume);
    World(Vector3D min, Vector3D max, size_t max_steps = 1000, double epsilon = 1e-4)
        : m_min{min}, m_max{max}, m_max_steps{max_steps}, m_epsilon{epsilon} {}

    Shape *GetShape(size_t idx) const;
    Material GetMaterial(size_t idx) const;
    std::vector<Material> GetMaterials() const;
    BoundingBox GetWorldBox() const;

    bool InWorld(const Vector3D &) const;
    bool SphereTrace(const Ray &, double &, size_t &, size_t &) const;
    bool RayTrace(const Ray &, double &, size_t &) const;
    /// Default traversal (boundary-sweep; NUGEOM_TRAVERSAL=sequential overrides).
    std::vector<LineSegment> GetLineSegments(const Ray &) const;
    /// Sequential boundary-to-boundary traversal (fallback for cross-checking).
    std::vector<LineSegment> GetLineSegmentsSequential(const Ray &) const;
    /// Analytic boundary-sweep (the default): collect all volume crossings in
    /// one descent (Shape::IntersectAll handles multi-interval CSG), sort, and
    /// sweep tracking the innermost (deepest) active volume.  ~1.9x faster than
    /// and ray-for-ray equivalent to the sequential traversal (cross-validated
    /// against it and ROOT's TGeoNavigator).
    std::vector<LineSegment> GetLineSegmentsSweep(const Ray &) const;

    /// A region where the experimental sweep disagrees with the hierarchical
    /// point-in-volume containment (FindMaterial) — a likely overlapping /
    /// non-nested-volume spot in the GDML worth inspecting.
    struct SweepOverlap {
        Vector3D start, end;
        std::string sweep_material;     ///< material the sweep assigned
        std::string contained_material; ///< material at the segment midpoint by containment
        std::string sweep_volume;       ///< volume the sweep selected (deepest active)
        std::string contained_volume;   ///< volume containing the midpoint (placement order)
    };
    /// Run the sweep and flag every segment whose midpoint material disagrees
    /// with the containment oracle.  Logs a warning per discrepancy (likely a
    /// GDML overlap) and returns them.  @param warn  emit log warnings.
    std::vector<SweepOverlap> CheckSweepConsistency(const Ray &ray, bool warn = true) const;
    size_t NDaughters() const { return m_root_pv ? m_root_pv->Daughters().size() : 0; }
    const std::shared_ptr<LogicalVolume> &GetLogicalVolume() const { return m_volume; }
    const std::shared_ptr<PhysicalVolume> &GetRootPV() const { return m_root_pv; }

    /// World-frame axis-aligned bounding box plus material for each placed volume,
    /// collected by a recursive walk of the daughter hierarchy.
    /// @param max_depth  Maximum recursion depth (0 = direct world daughters only).
    struct VolumeBounds {
        BoundingBox bb; ///< axis-aligned bounding box in world frame
        Material material;
        int depth;        ///< nesting depth (0 = direct world daughter)
        std::string name; ///< logical volume name
    };
    std::vector<VolumeBounds> GetVolumeBounds(int max_depth = 4) const;

    /// Return the material at a world-frame point by walking the volume hierarchy.
    Material FindMaterial(const Vector3D &point) const;
    /// Return the name of the deepest volume containing a world-frame point
    /// (placement order resolves overlaps).  Empty if the point is outside.
    std::string FindVolume(const Vector3D &point) const;

    /// Voxelized grid of material indices.
    struct VoxelGrid {
        std::vector<int16_t> data;       ///< material index per voxel, flat [x + nx*(y + ny*z)]
        std::vector<Material> materials; ///< index → Material lookup
        int nx{0}, ny{0}, nz{0};
        BoundingBox bounds;
        int16_t at(int x, int y, int z) const {
            return data[static_cast<size_t>(x + nx * (y + ny * z))];
        }
    };

    /// Voxelize the world geometry.  The longest axis gets @p resolution cells.
    VoxelGrid Voxelize(int resolution) const;

  private:
    std::pair<double, size_t> GetSDF(const Vector3D &) const;
    std::pair<double, size_t> GetSDFNonNeg(const Vector3D &) const;

    Vector3D m_min{}, m_max{};
    size_t m_max_steps{512};
    double m_epsilon{1e-4};
    std::shared_ptr<LogicalVolume> m_volume;
    /// Synthetic root PhysicalVolume wrapping the world LogicalVolume.
    /// Top-level daughters have this as their m_mother, so exiting them
    /// unwinds back to the root PV which uses its own BVH for the next daughter.
    std::shared_ptr<PhysicalVolume> m_root_pv;
};

} // namespace NuGeom
