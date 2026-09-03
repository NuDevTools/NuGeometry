#pragma once

#include "geom/Material.hh"
#include "geom/Shape.hh"
#include "geom/Volume.hh"
#include <memory>
#include <set>
#include <string>
#include <utility>
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

    /// Path length per material along the ray, without materializing
    /// LineSegments.  Layer-1 acceptance needs only sum(n_i sigma_i L_i), so
    /// this shares the sweep with GetLineSegmentsSweep but skips segment
    /// construction and the adjacent-segment merge.  Materials are owned by the
    /// volume tree; the pointers must not outlive this World.
    std::vector<std::pair<const Material *, double>> GetColumnLengths(const Ray &) const;

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

    /// Criteria for dropping placed volumes from the traversal tree.
    ///
    /// Pruning is a *speed* knob, not a physics-neutral one: a dropped volume's
    /// space reverts to its mother's material, so the mother's column density
    /// grows by the dropped shape's volume.  `PruneReport::mass_delta` is that
    /// net change, and it is the number to check before trusting a pruned run.
    struct PruneOptions {
        /// Drop any placed volume whose subtree mass is below this fraction of
        /// the unpruned world mass.  0 disables the mass criterion.
        double min_mass_fraction{0.0};
        /// Drop any placed volume whose material name is in this set, whatever
        /// its mass (e.g. materials the generator has no target for).  The whole
        /// subtree goes with it.
        std::set<std::string> drop_materials{};
        /// Never drop a volume matching these names (checked before the criteria
        /// above), so a light-but-essential volume can be protected.
        std::set<std::string> keep_volumes{};
    };

    struct PruneReport {
        struct Entry {
            std::string volume;      ///< placed volume name
            std::string material;    ///< its material
            std::string replaced_by; ///< mother material now filling the space
            double mass{0};          ///< subtree mass removed [g]
            double fraction{0};      ///< mass / total_mass
            double mass_delta{0};    ///< net geometry mass change from this drop [g]
            size_t nodes{0};         ///< placed volumes removed (this + descendants)
            std::string reason;      ///< "material", "mass fraction"
        };
        size_t removed_subtrees{0};
        size_t removed_nodes{0}; ///< total placed volumes removed
        double total_mass{0};    ///< world mass before pruning [g]
        double removed_mass{0};  ///< sum of dropped subtree masses [g]
        double mass_delta{0};    ///< net change in geometry mass [g]
        std::vector<Entry> entries;
        double RemovedMassFraction() const {
            return total_mass > 0 ? removed_mass / total_mass : 0.0;
        }
        /// Net mass change as a fraction of the original world mass -- the
        /// figure of merit for how much the pruning perturbed the physics.
        double MassDeltaFraction() const { return total_mass > 0 ? mass_delta / total_mass : 0.0; }
    };

    /// Drop placed volumes matching @p opts from the traversal tree, in place.
    /// Returns what was removed.  Call before any ray tracing; cached BVHs of
    /// modified nodes are invalidated automatically.
    PruneReport Prune(const PruneOptions &opts);

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
