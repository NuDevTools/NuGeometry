#include "geom/BVH.hh"
#include "geom/Ray.hh"
#include "geom/Volume.hh"

#include <algorithm>
#include <limits>

using NuGeom::BoundingBox;
using NuGeom::BVH;

void BVH::Build(const std::vector<std::shared_ptr<NuGeom::PhysicalVolume>> &daughters) {
    m_daughters = &daughters;
    m_nodes.clear();
    if(daughters.empty()) return;

    // Build leaf entries: (parent-space AABB, daughter index).
    std::vector<std::pair<BoundingBox, size_t>> leaves;
    leaves.reserve(daughters.size());
    for(size_t i = 0; i < daughters.size(); ++i)
        leaves.emplace_back(daughters[i]->GetParentBoundingBox(), i);

    m_nodes.reserve(2 * daughters.size());
    BuildNode(leaves, 0, leaves.size());
}

size_t BVH::BuildNode(std::vector<std::pair<BoundingBox, size_t>> &leaves, size_t start,
                      size_t end) {
    const size_t node_idx = m_nodes.size();
    m_nodes.push_back({}); // reserve slot; recursive calls may reallocate, use index only

    // Merge AABBs for this range.
    BoundingBox combined = leaves[start].first;
    for(size_t i = start + 1; i < end; ++i)
        combined = BoundingBox::Merge(combined, leaves[i].first);

    if(end - start == 1) {
        m_nodes[node_idx].bbox = combined;
        m_nodes[node_idx].pv_idx = leaves[start].second;
        return node_idx;
    }

    // Split along the longest axis at the median centroid.
    NuGeom::Vector3D extent = combined.max - combined.min;
    size_t axis = 0;
    if(extent.Y() > extent.X() && extent.Y() >= extent.Z())
        axis = 1;
    else if(extent.Z() > extent.X() && extent.Z() > extent.Y())
        axis = 2;

    size_t mid = (start + end) / 2;
    std::nth_element(
        leaves.begin() + static_cast<ptrdiff_t>(start),
        leaves.begin() + static_cast<ptrdiff_t>(mid), leaves.begin() + static_cast<ptrdiff_t>(end),
        [axis](const std::pair<BoundingBox, size_t> &a, const std::pair<BoundingBox, size_t> &b) {
            double ca = (a.first.min[axis] + a.first.max[axis]) * 0.5;
            double cb = (b.first.min[axis] + b.first.max[axis]) * 0.5;
            return ca < cb;
        });

    const size_t left_idx = BuildNode(leaves, start, mid);
    const size_t right_idx = BuildNode(leaves, mid, end);

    // Safe to use index after recursive calls (vector may have been reallocated).
    m_nodes[node_idx].bbox = combined;
    m_nodes[node_idx].left = left_idx;
    m_nodes[node_idx].right = right_idx;
    return node_idx;
}

bool BVH::Traverse(const NuGeom::Ray &ray, double &time,
                   std::shared_ptr<NuGeom::PhysicalVolume> &vol) const {
    if(m_nodes.empty()) return false;

    double best_time = std::numeric_limits<double>::infinity();
    size_t best_idx = kInvalid;
    TraverseNode(0, ray, best_time, best_idx);

    if(best_idx == kInvalid) return false;
    time = best_time;
    vol = (*m_daughters)[best_idx];
    return true;
}

void BVH::TraverseNode(size_t idx, const NuGeom::Ray &ray, double &best_time,
                       size_t &best_pv_idx) const {
    const Node &node = m_nodes[idx];

    // Cull entire subtree if its AABB is missed or farther than current best.
    if(!node.bbox.Intersect(ray, 0.0, best_time)) return;

    if(node.pv_idx != kInvalid) {
        // Leaf: exact shape intersection.
        const double t = (*m_daughters)[node.pv_idx]->Intersect(ray);
        if(t < best_time) {
            best_time = t;
            best_pv_idx = node.pv_idx;
        }
        return;
    }

    if(node.left != kInvalid) TraverseNode(node.left, ray, best_time, best_pv_idx);
    if(node.right != kInvalid) TraverseNode(node.right, ray, best_time, best_pv_idx);
}
