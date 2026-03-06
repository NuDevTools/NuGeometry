#include "geom/Volume.hh"
#include "geom/BVH.hh"
#include "geom/LineSegment.hh"
#include "geom/Ray.hh"
#include "spdlog/spdlog.h"

#include <limits>
#include <numeric>

using NuGeom::LogicalVolume;
using NuGeom::PhysicalVolume;

LogicalVolume::LogicalVolume() = default;
LogicalVolume::LogicalVolume(Material material, std::shared_ptr<Shape> shape)
    : m_material{std::move(material)}, m_shape{std::move(shape)} {}
LogicalVolume::LogicalVolume(std::string name, Material material, std::shared_ptr<Shape> shape)
    : m_name{std::move(name)}, m_material{std::move(material)}, m_shape{std::move(shape)} {}
// Out-of-line destructor so BVH's full type is visible for unique_ptr<BVH>.
LogicalVolume::~LogicalVolume() = default;

void LogicalVolume::GetMaterials(std::set<Material> &mats) const {
    mats.insert(m_material);
    for(const auto &daughter : m_daughters) { daughter->GetLogicalVolume()->GetMaterials(mats); }
}

double LogicalVolume::Mass() const {
    return Volume() * m_material.Density() + DaughterMass();
}

double AddVolume(double a, const std::shared_ptr<PhysicalVolume> &b) {
    return a + b->GetLogicalVolume()->Volume();
}

double LogicalVolume::DaughterVolumes() const {
    return std::accumulate(m_daughters.begin(), m_daughters.end(), 0.0, AddVolume);
}

double LogicalVolume::Volume() const {
    return m_shape->Volume() - DaughterVolumes();
}

double AddMass(double a, const std::shared_ptr<PhysicalVolume> &b) {
    return a + b->GetLogicalVolume()->Mass();
}

double LogicalVolume::DaughterMass() const {
    return std::accumulate(m_daughters.begin(), m_daughters.end(), 0.0, AddMass);
}

bool LogicalVolume::SphereTrace(const Ray &ray, double &time, size_t &step, size_t &idx) const {
    step = 0;
    Vector3D pos = ray.Origin();
    time = 0;
    auto res = GetSDF(pos);
    while(step < m_max_steps && std::abs(res.first) > m_epsilon && InWorld(pos)) {
        pos = ray.Propagate(time);
        res = GetSDF(pos);
        time += res.first;
        step++;
    }
    if(step == m_max_steps || !InWorld(pos) || res.second == 0) return false;
    idx = res.second;
    return true;
}

bool LogicalVolume::RayTrace(const Ray &ray, double &time,
                             std::shared_ptr<PhysicalVolume> &vol) const {
    if(!m_bvh) {
        m_bvh = std::make_shared<BVH>();
        m_bvh->Build(m_daughters);
    }
    return m_bvh->Traverse(ray, time, vol);
}

void LogicalVolume::GetLineSegments(const Ray &ray, std::vector<LineSegment> &segments,
                                    const Transform3D &from_global) const {
    // Iterative traversal: loop over daughters at this level, only recurse DOWN
    // into daughters (bounded by nesting depth).  This avoids O(N) recursion
    // depth for N daughters which causes stack overflow on complex geometries.
    static constexpr double eps = 1e-8;
    static constexpr size_t kMaxIter = 100000;
    auto current_ray = ray;

    for(size_t iter = 0; iter < kMaxIter; ++iter) {
        auto current_local = Transform3D::ApplyRayDirect(current_ray, from_global);
        auto shift_ray = Ray(current_local.Propagate(eps), current_local.Direction(),
                             current_local.POT(), false);

        // Always compute exit time from this volume's shape.
        double exit_time = m_shape->Intersect(shift_ray);

        // Find closest daughter.
        double daughter_time = std::numeric_limits<double>::infinity();
        std::shared_ptr<PhysicalVolume> pvol = nullptr;
        RayTrace(shift_ray, daughter_time, pvol);

        // Take whichever is closer: daughter entry or volume exit.
        double time;
        if(pvol && daughter_time < exit_time) {
            time = daughter_time;
        } else {
            time = exit_time;
            pvol = nullptr;
        }

        // Ray is already inside the daughter (daughter fills parent exactly).
        if(pvol && daughter_time == 0.0) {
            auto seg_before = segments.size();
            pvol->GetLineSegments(current_ray, segments, from_global);
            if(segments.size() > seg_before) {
                current_ray = Ray(segments.back().End(), ray.Direction(), ray.POT(), false);
                continue;
            }
            // Daughter didn't add segments — skip past it.
            // Intersect2 may return (inf, inf) if the ray grazes the AABB but
            // misses the actual shape; guard against non-finite dt2.
            auto daughter_local = Transform3D::ApplyRayDirect(current_local, pvol->GetTransform());
            auto [dt1, dt2] = pvol->GetLogicalVolume()->GetShape()->Intersect2(daughter_local);
            if(std::isfinite(dt2) && dt2 > eps) {
                segments.emplace_back(current_ray.Origin(), current_ray.Propagate(dt2 + eps),
                                      pvol->GetLogicalVolume()->GetMaterial());
                current_ray = Ray(segments.back().End(), ray.Direction(), ray.POT(), false);
            } else {
                current_ray = Ray(current_ray.Propagate(eps), ray.Direction(), ray.POT(), false);
            }
            continue;
        }

        time += eps;

        if(!std::isfinite(time)) {
            // Shifted ray is outside.  Fall back to the un-shifted ray.
            double fallback = m_shape->Intersect(current_local);
            if(std::isfinite(fallback) && fallback > 0) {
                segments.emplace_back(current_ray.Origin(), current_ray.Propagate(fallback + eps),
                                      m_material);
            }
            return;
        }

        segments.emplace_back(current_ray.Origin(), current_ray.Propagate(time), m_material);
        if(!pvol) return;

        // Enter daughter (recurse DOWN — bounded by nesting depth).
        auto origin = current_ray.Propagate(time);
        auto daughter_ray = Ray(origin, ray.Direction(), ray.POT(), false);
        pvol->GetLineSegments(daughter_ray, segments, from_global);

        // Continue from daughter's exit point.
        if(segments.empty()) return;
        current_ray = Ray(segments.back().End(), ray.Direction(), ray.POT(), false);
    }
    spdlog::warn("LogicalVolume::GetLineSegments: iteration limit reached in '{}'", m_name);
}

bool PhysicalVolume::RayTrace(const Ray &ray, double &time,
                              std::shared_ptr<PhysicalVolume> &pvol) const {
    if(m_own_daughters.empty()) return false;
    if(!m_bvh) {
        m_bvh = std::make_shared<BVH>();
        m_bvh->Build(m_own_daughters);
    }
    return m_bvh->Traverse(ray, time, pvol);
}

NuGeom::BoundingBox PhysicalVolume::GetParentBoundingBox() const {
    BoundingBox local_bb = m_volume->GetShape()->GetTransformedBoundingBox();

    // Forward transform: local → parent  =  m_transform.Inverse()
    const Transform3D forward = m_transform.Inverse();

    const Vector3D cx[2] = {local_bb.min, local_bb.max};
    auto c0 = forward.Apply(cx[0]);
    Vector3D mn = c0, mx = c0;
    for(size_t ix = 0; ix < 2; ++ix)
        for(size_t iy = 0; iy < 2; ++iy)
            for(size_t iz = 0; iz < 2; ++iz) {
                if(ix == 0 && iy == 0 && iz == 0) continue;
                Vector3D corner(cx[ix].X(), cx[iy].Y(), cx[iz].Z());
                auto c = forward.Apply(corner);
                mn = {std::min(mn.X(), c.X()), std::min(mn.Y(), c.Y()), std::min(mn.Z(), c.Z())};
                mx = {std::max(mx.X(), c.X()), std::max(mx.Y(), c.Y()), std::max(mx.Z(), c.Z())};
            }
    return {mn, mx};
}

double PhysicalVolume::Intersect(const Ray &in_ray) const {
    auto ray = TransformRay(in_ray);
    const auto [t1, t2] = m_volume->GetShape()->Intersect2(ray);
    if(t1 > 0) return t1;  // ray enters from outside: return entry time
    if(t2 > 0) return 0.0; // ray origin is inside: enter immediately
    return std::numeric_limits<double>::infinity();
}

void PhysicalVolume::GetLineSegments(const Ray &in_ray, std::vector<LineSegment> &segments,
                                     const Transform3D &from_global) const {
    // Iterative traversal: loop over daughters at this level, only recurse DOWN
    // into daughters (bounded by nesting depth).  Eliminates the mutual recursion
    // (child→mother→child) that caused stack overflow on complex geometries.
    static constexpr double eps = 1e-8;
    static constexpr size_t kMaxIter = 100000;
    auto current_ray = in_ray;
    // from_global for daughters: world → this PV's local frame.
    const auto daughter_fg = from_global * m_transform;

    for(size_t iter = 0; iter < kMaxIter; ++iter) {
        auto local_ray = Transform3D::ApplyRayDirect(current_ray, from_global);
        auto ray = TransformRay(local_ray);
        auto shift_ray = Ray(ray.Propagate(eps), ray.Direction(), ray.POT(), false);

        // Always compute exit time from this volume's shape.
        double exit_time = m_volume->GetShape()->Intersect(shift_ray);

        // Find closest daughter.
        double daughter_time = std::numeric_limits<double>::infinity();
        std::shared_ptr<PhysicalVolume> pvol = nullptr;
        RayTrace(shift_ray, daughter_time, pvol);

        // Take whichever is closer: daughter entry or volume exit.
        double time;
        if(pvol && daughter_time < exit_time) {
            time = daughter_time;
        } else {
            time = exit_time;
            pvol = nullptr;
        }

        // If daughter_time == 0 the ray is already inside the daughter (e.g.
        // the daughter fills the parent exactly).  Enter the daughter directly
        // from the current position — no parent-material segment, no eps
        // advance — so the daughter sees the un-shifted ray and can traverse
        // its full extent without being pushed outside its own boundary.
        if(pvol && daughter_time == 0.0) {
            auto seg_before = segments.size();
            pvol->GetLineSegments(current_ray, segments, daughter_fg);
            if(segments.size() > seg_before) {
                current_ray = Ray(segments.back().End(), in_ray.Direction(), in_ray.POT(), false);
                continue;
            }
            // Daughter returned without adding segments (boundary precision
            // issue).  Compute its exit in our local frame and skip past it.
            // Intersect2 may return (inf, inf) if the ray grazes the AABB but
            // misses the actual shape; guard against non-finite dt2.
            auto daughter_local = Transform3D::ApplyRayDirect(ray, pvol->GetTransform());
            auto [dt1, dt2] = pvol->GetLogicalVolume()->GetShape()->Intersect2(daughter_local);
            if(std::isfinite(dt2) && dt2 > eps) {
                // Emit daughter's material for its full extent
                segments.emplace_back(current_ray.Origin(), current_ray.Propagate(dt2 + eps),
                                      pvol->GetLogicalVolume()->GetMaterial());
                current_ray = Ray(segments.back().End(), in_ray.Direction(), in_ray.POT(), false);
            } else {
                // Daughter is negligible at this point; skip it
                current_ray =
                    Ray(current_ray.Propagate(eps), in_ray.Direction(), in_ray.POT(), false);
            }
            continue;
        }

        time += eps;

        if(!std::isfinite(time)) {
            // Shifted ray is outside this volume.  Fall back to the un-shifted
            // ray to get the actual exit time (handles boundary-flush case).
            double fallback = m_volume->GetShape()->Intersect(ray);
            if(std::isfinite(fallback) && fallback > 0) {
                segments.emplace_back(current_ray.Origin(), current_ray.Propagate(fallback + eps),
                                      m_volume->GetMaterial());
            }
            return;
        }

        auto next_origin = current_ray.Propagate(time);
        segments.emplace_back(current_ray.Origin(), next_origin, m_volume->GetMaterial());
        if(!pvol) return; // exited volume, no daughter ahead

        // Enter daughter (recurse DOWN only — bounded by nesting depth).
        auto daughter_ray = Ray(next_origin, in_ray.Direction(), in_ray.POT(), false);
        pvol->GetLineSegments(daughter_ray, segments, daughter_fg);

        // Continue from daughter's exit point.
        if(segments.empty()) return;
        current_ray = Ray(segments.back().End(), in_ray.Direction(), in_ray.POT(), false);
    }
    spdlog::warn("PhysicalVolume::GetLineSegments: iteration limit reached in '{}'", m_name);
}

NuGeom::Ray PhysicalVolume::TransformRay(const Ray &ray) const {
    if(is_identity)
        return ray;
    else if(is_translation)
        return Transform3D::TranslateRay(ray, m_trans);
    return Transform3D::ApplyRay(ray, m_trans, m_rot);
}

NuGeom::Ray PhysicalVolume::TransformRayInverse(const Ray &ray) const {
    if(is_identity)
        return ray;
    else if(is_translation)
        return Transform3D::TranslateRay(ray, m_trans);
    return Transform3D::ApplyRay(ray, m_transform.Inverse());
}
