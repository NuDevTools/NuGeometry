#include "geom/BVH.hh"
#include "geom/Ray.hh"
#include "geom/Volume.hh"

#include <algorithm>
#include <array>
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

    // Cache non-owning pointers parallel to the daughter list.
    m_daughter_ptrs.clear();
    m_daughter_ptrs.reserve(daughters.size());
    for(const auto &d : daughters) m_daughter_ptrs.push_back(d.get());
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

size_t BVH::TraverseIndex(const NuGeom::Ray &ray, double &time) const {
    if(m_nodes.empty()) return kInvalid;

    constexpr double inf = std::numeric_limits<double>::infinity();
    double best_time = inf;
    size_t best_idx = kInvalid;

    // Raw pointers to skip the vector's bounds-checked operator[] in the loop.
    const Node *const nodes = m_nodes.data();
    PhysicalVolume *const *const dptrs = m_daughter_ptrs.data();

    // Iterative depth-first traversal.  Each node's AABB is tested exactly
    // once (when pushed); the entry time orders the children so the nearer
    // subtree is visited first.  The median-split tree is balanced, so the
    // stack depth is bounded by ~log2(n)+1; 64 covers any realistic geometry.
    struct StackEntry {
        size_t node;
        double t_entry;
    };
    std::array<StackEntry, 64> stack;
    size_t sp = 0;

    if(std::isfinite(nodes[0].bbox.IntersectT(ray, 0.0, best_time))) stack[sp++] = {0, 0.0};

    while(sp > 0) {
        const StackEntry entry = stack[--sp];
        // A closer hit may have been found since this subtree was pushed.
        if(entry.t_entry >= best_time) continue;
        const Node &node = nodes[entry.node];

        if(node.pv_idx != kInvalid) {
            // Leaf: exact shape intersection.
            const double t = dptrs[node.pv_idx]->Intersect(ray);
            if(t < best_time) {
                best_time = t;
                best_idx = node.pv_idx;
            }
            continue;
        }

        double t_left = inf, t_right = inf;
        if(node.left != kInvalid) t_left = nodes[node.left].bbox.IntersectT(ray, 0.0, best_time);
        if(node.right != kInvalid) t_right = nodes[node.right].bbox.IntersectT(ray, 0.0, best_time);

        // Push the farther child first so the nearer one is popped first.
        const bool left_first = t_left <= t_right;
        const size_t far_node = left_first ? node.right : node.left;
        const double far_t = left_first ? t_right : t_left;
        const size_t near_node = left_first ? node.left : node.right;
        const double near_t = left_first ? t_left : t_right;
        if(std::isfinite(far_t)) stack[sp++] = {far_node, far_t};
        if(std::isfinite(near_t)) stack[sp++] = {near_node, near_t};
    }

    if(best_idx == kInvalid) return kInvalid;
    time = best_time;
    return best_idx;
}

bool BVH::Traverse(const NuGeom::Ray &ray, double &time,
                   std::shared_ptr<NuGeom::PhysicalVolume> &vol) const {
    const size_t idx = TraverseIndex(ray, time);
    if(idx == kInvalid) return false;
    vol = (*m_daughters)[idx];
    return true;
}

bool BVH::Traverse(const NuGeom::Ray &ray, double &time, NuGeom::PhysicalVolume *&vol) const {
    const size_t idx = TraverseIndex(ray, time);
    if(idx == kInvalid) return false;
    vol = m_daughter_ptrs[idx];
    return true;
}

void BVH::CollectHits(const NuGeom::Ray &ray, std::vector<size_t> &out) const {
    if(m_nodes.empty()) return;
    const Node *const nodes = m_nodes.data();
    std::array<size_t, 64> stack;
    size_t sp = 0;
    if(std::isfinite(nodes[0].bbox.IntersectT(ray))) stack[sp++] = 0;
    while(sp > 0) {
        const Node &node = nodes[stack[--sp]];
        if(node.pv_idx != kInvalid) {
            out.push_back(node.pv_idx);
            continue;
        }
        if(node.left != kInvalid && std::isfinite(nodes[node.left].bbox.IntersectT(ray)))
            stack[sp++] = node.left;
        if(node.right != kInvalid && std::isfinite(nodes[node.right].bbox.IntersectT(ray)))
            stack[sp++] = node.right;
    }
}
