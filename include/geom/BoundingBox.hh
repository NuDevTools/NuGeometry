#pragma once

#include "geom/Ray.hh"
#include "geom/Vector3D.hh"
#include <algorithm>
#include <cmath>
#include <limits>

namespace NuGeom {

struct BoundingBox {
    Vector3D min, max;

    double Volume() const {
        return (max.X() - min.X()) * (max.Y() - min.Y()) * (max.Z() - min.Z());
    }

    bool IsValid() const { return min.X() <= max.X() && min.Y() <= max.Y() && min.Z() <= max.Z(); }

    /// Returns the AABB that contains both a and b.
    static BoundingBox Merge(const BoundingBox &a, const BoundingBox &b) {
        return {Vector3D(std::min(a.min.X(), b.min.X()), std::min(a.min.Y(), b.min.Y()),
                         std::min(a.min.Z(), b.min.Z())),
                Vector3D(std::max(a.max.X(), b.max.X()), std::max(a.max.Y(), b.max.Y()),
                         std::max(a.max.Z(), b.max.Z()))};
    }

    /// Slab-method AABB/ray intersection.
    /// Returns the entry time (clamped to t_near) if the ray hits the box
    /// with t in [t_near, t_far), or +infinity on a miss.
    double IntersectT(const Ray &ray, double t_near = 0.0,
                      double t_far = std::numeric_limits<double>::infinity()) const {
        const auto &origin = ray.Origin();
        const auto &inv_dir = ray.InvDirection();
        // Branchless slab test.  For each axis the two slab-crossing times are
        // (bound - origin) * inv_dir; std::min/std::max fold them into the
        // running [t_near, t_far] interval.  No explicit finite check is needed:
        // for an axis-parallel ray (inv_dir == +/-inf) IEEE arithmetic yields
        // +/-inf when the origin is off the slab face (-> the interval collapses
        // -> miss) and a 0*inf NaN exactly on the face, which std::min/std::max
        // discard so the axis stays unconstrained -- matching the old isfinite
        // branch's "origin inside the slab -> no constraint" behavior.
        // Ray::InvDirection() is precomputed, so this costs no divisions.
        for(size_t axis = 0; axis < 3; ++axis) {
            const double orig = origin[axis];
            const double inv_d = inv_dir[axis];
            const double t0 = (min[axis] - orig) * inv_d;
            const double t1 = (max[axis] - orig) * inv_d;
            t_near = std::max(t_near, std::min(t0, t1));
            t_far = std::min(t_far, std::max(t0, t1));
        }
        return t_near <= t_far ? t_near : std::numeric_limits<double>::infinity();
    }

    /// Returns true if the ray hits the box with t in [t_near, t_far).
    bool Intersect(const Ray &ray, double t_near = 0.0,
                   double t_far = std::numeric_limits<double>::infinity()) const {
        return std::isfinite(IntersectT(ray, t_near, t_far));
    }
};

} // namespace NuGeom
