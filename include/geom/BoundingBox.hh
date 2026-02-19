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

    /// Slab-method AABB/ray intersection test.
    /// Returns true if the ray hits the box with t in [t_near, t_far).
    bool Intersect(const Ray &ray, double t_near = 0.0,
                   double t_far = std::numeric_limits<double>::infinity()) const {
        for(size_t axis = 0; axis < 3; ++axis) {
            const double orig = ray.Origin()[axis];
            const double dir = ray.Direction()[axis];
            const double lo = min[axis];
            const double hi = max[axis];
            if(std::abs(dir) < 1e-12) {
                if(orig < lo || orig > hi) return false;
            } else {
                const double inv_d = 1.0 / dir;
                double t0 = (lo - orig) * inv_d;
                double t1 = (hi - orig) * inv_d;
                if(inv_d < 0) std::swap(t0, t1);
                t_near = std::max(t_near, t0);
                t_far = std::min(t_far, t1);
                if(t_near > t_far) return false;
            }
        }
        return true;
    }
};

} // namespace NuGeom
