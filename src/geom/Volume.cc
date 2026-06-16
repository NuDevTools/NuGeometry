#include "geom/Volume.hh"
#include "geom/BVH.hh"
#include "geom/LineSegment.hh"
#include "geom/Logging.hh"
#include "geom/Ray.hh"

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
    NuGeom::Log().warn("LogicalVolume::GetLineSegments: iteration limit reached in '{}'", m_name);
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

bool PhysicalVolume::RayTrace(const Ray &ray, double &time, PhysicalVolume *&pvol) const {
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
    const auto [t1, t2] = m_volume->GetShape()->Intersect2(TransformRayCached(in_ray));
    if(t1 > 0) return t1;  // ray enters from outside: return entry time
    if(t2 > 0) return 0.0; // ray origin is inside: enter immediately
    return std::numeric_limits<double>::infinity();
}

// Like TransformRay, but reuses a cached rotated direction + its inverse when
// the input direction is unchanged (the common case while a single ray walks
// the geometry), avoiding the per-call direction rotation and the three
// divisions the Ray constructor would otherwise do.  Only the origin transform
// runs every call.  Bit-identical to TransformRay for a given input.
NuGeom::Ray PhysicalVolume::TransformRayCached(const Ray &in_ray) const {
    if(is_identity) return in_ray;

    const Vector3D &in_dir = in_ray.Direction();
    if(!m_dir_cache_valid || in_dir != m_dir_cache_in) {
        const Ray td = TransformRay(in_ray); // direction transform ignores translation
        m_dir_cache_in = in_dir;
        m_dir_cache_local = td.Direction();
        m_dir_cache_inv = td.InvDirection();
        m_dir_cache_valid = true;
    }
    const Vector3D local_origin = is_translation ? m_trans.Apply(in_ray.Origin())
                                                 : m_rot.Apply(m_trans.Apply(in_ray.Origin()));
    return Ray(local_origin, m_dir_cache_local, m_dir_cache_inv, in_ray.POT());
}

void PhysicalVolume::GetLineSegments(const Ray &in_ray, std::vector<LineSegment> &segments,
                                     const Transform3D &from_global) const {
    // Iterative traversal: loop over daughters at this level, only recurse DOWN
    // into daughters (bounded by nesting depth).  Eliminates the mutual recursion
    // (child→mother→child) that caused stack overflow on complex geometries.
    static constexpr double eps = 1e-8;
    static constexpr size_t kMaxIter = 100000;
    auto current_ray = in_ray;
    // from_global for daughters: world → this PV's local frame.  Composition
    // order matters: world->local applies from_global (world->parent) first,
    // then m_transform (parent->local), i.e. m_transform ∘ from_global.  With
    // operator* defined as (A*B).Apply(p) == A.Apply(B.Apply(p)), that is
    // m_transform * from_global.  (Translations commute so this was latent until
    // a rotated, deeply-nested placement like volHalfDetector_R exposed it.)
    const auto daughter_fg = m_transform * from_global;

    // The global ray direction is invariant for the whole traversal, so reuse
    // its precomputed inverse instead of redividing on every reconstruction.
    const Vector3D &gdir = in_ray.Direction();
    const Vector3D &ginv = in_ray.InvDirection();

    // If the ray starts OUTSIDE this volume, advance it to the entry face first.
    // The loop below assumes an interior origin (so Shape::Intersect returns the
    // exit time); with an exterior origin Shape::Intersect returns the *entry*
    // time instead, which made the traversal emit a single segment that stops at
    // the box surface and return.  This happens for flux rays whose decay-vertex
    // origin sits far upstream of the world box.  The exterior region carries no
    // material, so we emit nothing for it — we just step the ray onto the volume.
    {
        auto entry_ray = TransformRay(Transform3D::ApplyRayDirect(current_ray, from_global));
        if(m_volume->GetShape()->SignedDistance(entry_ray.Origin()) > 0) {
            const auto [t_in, t_out] = m_volume->GetShape()->Intersect2(entry_ray);
            if(!(std::isfinite(t_in) && t_in > 0 && t_out > t_in))
                return; // ray never enters this volume
            current_ray = Ray(current_ray.Propagate(t_in + eps), gdir, ginv, in_ray.POT());
        }
    }

    // The shape-local ray direction (and its inverse) is also invariant across
    // iterations — only the origin advances.  Compute it once via the normal
    // transform path, then per iteration transform just the advancing origin.
    // The direction transform ignores translation, so this matches exactly what
    // TransformRay(ApplyRayDirect(...)) produced per-iteration before.
    const Ray probe_local = TransformRay(Transform3D::ApplyRayDirect(in_ray, from_global));
    const Vector3D shape_dir = probe_local.Direction();
    const Vector3D shape_inv = probe_local.InvDirection();
    // world origin → this PV's shape-local frame, mirroring the origin handling
    // of ApplyRayDirect(·, from_global) followed by TransformRay(·).
    const auto to_shape_origin = [&](const Vector3D &world_o) -> Vector3D {
        const Vector3D lo = from_global.Apply(world_o);
        if(is_identity) return lo;
        if(is_translation) return m_trans.Apply(lo);
        return m_rot.Apply(m_trans.Apply(lo));
    };

    for(size_t iter = 0; iter < kMaxIter; ++iter) {
        Ray ray(to_shape_origin(current_ray.Origin()), shape_dir, shape_inv, in_ray.POT());
        Ray shift_ray(ray.Propagate(eps), shape_dir, shape_inv, in_ray.POT());

        // Always compute exit time from this volume's shape.
        double exit_time = m_volume->GetShape()->Intersect(shift_ray);

        // Find closest daughter.
        double daughter_time = std::numeric_limits<double>::infinity();
        PhysicalVolume *pvol = nullptr;
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
                current_ray = Ray(segments.back().End(), gdir, ginv, in_ray.POT());
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
                current_ray = Ray(segments.back().End(), gdir, ginv, in_ray.POT());
            } else {
                // Daughter is negligible at this point; skip it
                current_ray = Ray(current_ray.Propagate(eps), gdir, ginv, in_ray.POT());
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
        auto daughter_ray = Ray(next_origin, gdir, ginv, in_ray.POT());
        pvol->GetLineSegments(daughter_ray, segments, daughter_fg);

        // Continue from daughter's exit point.
        if(segments.empty()) return;
        current_ray = Ray(segments.back().End(), gdir, ginv, in_ray.POT());
    }
    NuGeom::Log().warn("PhysicalVolume::GetLineSegments: iteration limit reached in '{}'", m_name);
}

void PhysicalVolume::CollectIntervals(const Ray &world_ray, const Transform3D &from_global,
                                      int depth, double parent_lo, double parent_hi,
                                      std::vector<IntervalEvent> &events) const {
    // Ray in this volume's shape-local frame (rigid transforms preserve the ray
    // parameter t, so the local entry/exit times are world-frame distances).
    const Ray local = TransformRay(Transform3D::ApplyRayDirect(world_ray, from_global));
    // IntersectAll returns every solid span -- a non-convex CSG (e.g. a window
    // frame with two walls along the ray) has more than one, which a single
    // Intersect2 would collapse to the first wall only.
    const auto spans = m_volume->GetShape()->IntersectAll(local);
    if(spans.empty()) return;

    const Material &mat = m_volume->GetMaterial();
    // Emit each span clipped to the parent's active window (a daughter is only
    // "really" inside the geometry where its mother also contains the ray), and
    // track the overall [span_lo, span_hi] to clip descendants to.  Daughters
    // sit inside this volume's solid (including any gap a frame leaves for them,
    // which lies within the overall span), so the overall span is the right
    // window to pass down.
    double span_lo = std::numeric_limits<double>::infinity();
    double span_hi = -std::numeric_limits<double>::infinity();
    for(const auto &iv : spans) {
        const double enter = iv.first > 0.0 ? iv.first : 0.0; // origin inside => from 0
        const double lo = std::max(enter, parent_lo);
        const double hi = std::min(iv.second, parent_hi);
        if(lo >= hi) continue; // this span does not overlap the mother's window
        events.push_back({lo, +1, depth, &mat, this});
        events.push_back({hi, -1, depth, &mat, this});
        span_lo = std::min(span_lo, lo);
        span_hi = std::max(span_hi, hi);
    }
    if(span_hi <= span_lo) return; // nothing landed in the parent window

    if(m_own_daughters.empty()) return;
    if(!m_bvh) {
        m_bvh = std::make_shared<BVH>();
        m_bvh->Build(m_own_daughters);
    }
    std::vector<size_t> hits;
    m_bvh->CollectHits(local, hits); // daughters' parent-frame AABBs live in `local`'s frame
    const Transform3D daughter_fg = m_transform * from_global; // see GetLineSegments for order
    for(size_t idx : hits) // children clipped to THIS volume's overall span
        m_own_daughters[idx]->CollectIntervals(world_ray, daughter_fg, depth + 1, span_lo, span_hi,
                                               events);
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
