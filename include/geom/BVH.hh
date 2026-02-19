#pragma once

#include "geom/BoundingBox.hh"
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace NuGeom {

class PhysicalVolume;
class Ray;

/// Bounding Volume Hierarchy built over a LogicalVolume's daughters.
/// Each leaf corresponds to one PhysicalVolume daughter.  Internal nodes
/// hold the merged AABB of their subtree, allowing fast ray culling.
class BVH {
  public:
    BVH() = default;

    /// Build (or rebuild) the BVH over the given daughter list.
    void Build(const std::vector<std::shared_ptr<PhysicalVolume>> &daughters);

    /// Traverse the BVH and find the nearest daughter intersected by ray.
    /// Returns true and sets time / vol if a hit is found.
    bool Traverse(const Ray &ray, double &time, std::shared_ptr<PhysicalVolume> &vol) const;

    bool IsBuilt() const { return !m_nodes.empty(); }

  private:
    static constexpr size_t kInvalid = std::numeric_limits<size_t>::max();

    struct Node {
        BoundingBox bbox;
        size_t left = kInvalid; // index into m_nodes; kInvalid means none
        size_t right = kInvalid;
        size_t pv_idx = kInvalid; // < kInvalid for leaf nodes (index into *m_daughters)
    };

    /// Recursively build nodes for the subrange [start, end) of leaves.
    /// Returns the index of the newly created node.
    size_t BuildNode(std::vector<std::pair<BoundingBox, size_t>> &leaves, size_t start, size_t end);

    void TraverseNode(size_t idx, const Ray &ray, double &best_time, size_t &best_pv_idx) const;

    std::vector<Node> m_nodes;
    const std::vector<std::shared_ptr<PhysicalVolume>> *m_daughters = nullptr;
};

} // namespace NuGeom
