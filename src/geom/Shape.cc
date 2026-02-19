#include "geom/Shape.hh"
#include "geom/BoundingBox.hh"
#include "geom/Random.hh"
#include "geom/Ray.hh"
#include "geom/Transform3D.hh"
#include "geom/Vector2D.hh"
#include "geom/Vector3D.hh"
#include "pugixml.hpp"
#include "spdlog/spdlog.h"
#include <algorithm>
#include <limits>
#include <stdexcept>

NuGeom::Location NuGeom::Shape::Contains(const Vector3D &point) const {
    double dist = SignedDistance(point);
    if(dist < 0)
        return NuGeom::Location::kInterior;
    else if(dist > 0)
        return NuGeom::Location::kExterior;
    else
        return NuGeom::Location::kSurface;
}

NuGeom::Vector3D NuGeom::Shape::TransformPoint(const Vector3D &point) const {
    return m_rotation.Apply(m_translation.Apply(point));
}

NuGeom::Ray NuGeom::Shape::TransformRay(const Ray &in_ray) const {
    auto origin = m_rotation.Apply(m_translation.Apply(in_ray.Origin()));
    auto direction = m_rotation.Apply(in_ray.Direction());
    return {origin, direction, in_ray.POT()};
}

NuGeom::BoundingBox NuGeom::Shape::GetTransformedBoundingBox() const {
    if(identity_transform) return GetBoundingBox();

    BoundingBox bb = GetBoundingBox();
    // Forward transform: shape-local → parent.
    // TransformPoint goes parent→local via: rot_inv ∘ trans_inv.
    // The inverse goes local→parent via: trans ∘ rot  =  m_translation.Inverse() ∘
    // m_rotation.Inverse().
    const Transform3D rot_fwd = m_rotation.Inverse();
    const Transform3D trans_fwd = m_translation.Inverse();

    const Vector3D cx[2] = {bb.min, bb.max};
    auto c0 = trans_fwd.Apply(rot_fwd.Apply(cx[0]));
    Vector3D mn = c0, mx = c0;
    for(int ix = 0; ix < 2; ++ix)
        for(int iy = 0; iy < 2; ++iy)
            for(int iz = 0; iz < 2; ++iz) {
                if(ix == 0 && iy == 0 && iz == 0) continue; // already handled
                Vector3D corner(cx[ix].X(), cx[iy].Y(), cx[iz].Z());
                auto c = trans_fwd.Apply(rot_fwd.Apply(corner));
                mn = {std::min(mn.X(), c.X()), std::min(mn.Y(), c.Y()), std::min(mn.Z(), c.Z())};
                mx = {std::max(mx.X(), c.X()), std::max(mx.Y(), c.Y()), std::max(mx.Z(), c.Z())};
            }
    return {mn, mx};
}

double NuGeom::Shape::Intersect(const Ray &in_ray) const {
    const auto [t1, t2] = Intersect2(in_ray);
    if(t1 > 0) return t1;
    if(t2 > 0) return t2;
    return std::numeric_limits<double>::infinity();
}

std::pair<double, double> NuGeom::Shape::Intersect2(const Ray &in_ray) const {
    auto ray = identity_transform ? in_ray : TransformRay(in_ray);
    return Intersect2Impl(ray);
}

std::pair<double, double> NuGeom::Shape::SolveQuadratic(double a, double b, double c) const {
    const double det = b * b - 4 * a * c;
    if(det < 0)
        return {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    double t1 = 2 * c / (-b - sqrt(det));
    double t2 = 2 * c / (-b + sqrt(det));
    t1 = t1 > 0 ? t1 : std::numeric_limits<double>::infinity();
    t2 = t2 > 0 ? t2 : std::numeric_limits<double>::infinity();
    return {t1, t2};
}

// TODO: Do this correctly!!!!
std::unique_ptr<NuGeom::Shape> NuGeom::CombinedShape::Construct(const pugi::xml_node &) {
    auto box1 = std::make_shared<NuGeom::Box>();
    auto box2 = std::make_shared<NuGeom::Box>();
    return std::make_unique<NuGeom::CombinedShape>(box1, box2, NuGeom::ShapeBinaryOp::kUnion);
}

double NuGeom::CombinedShape::SignedDistance(const Vector3D &in_point) const {
    auto point = TransformPoint(in_point);
    double sdf1 = m_left->SignedDistance(point);
    double sdf2 = m_right->SignedDistance(point);
    switch(m_op) {
    case ShapeBinaryOp::kUnion:
        return std::min(sdf1, sdf2);
    case ShapeBinaryOp::kIntersect:
        return std::max(sdf1, sdf2);
    case ShapeBinaryOp::kSubtraction:
        return std::max(-sdf1, sdf2);
    }
    return std::min(sdf1, sdf2); // unreachable
}

NuGeom::BoundingBox NuGeom::CombinedShape::GetBoundingBox() const {
    BoundingBox L = m_left->GetBoundingBox();
    BoundingBox R = m_right->GetBoundingBox();
    switch(m_op) {
    case ShapeBinaryOp::kUnion:
        return {Vector3D(std::min(L.min.X(), R.min.X()), std::min(L.min.Y(), R.min.Y()),
                         std::min(L.min.Z(), R.min.Z())),
                Vector3D(std::max(L.max.X(), R.max.X()), std::max(L.max.Y(), R.max.Y()),
                         std::max(L.max.Z(), R.max.Z()))};
    case ShapeBinaryOp::kIntersect:
        return {Vector3D(std::max(L.min.X(), R.min.X()), std::max(L.min.Y(), R.min.Y()),
                         std::max(L.min.Z(), R.min.Z())),
                Vector3D(std::min(L.max.X(), R.max.X()), std::min(L.max.Y(), R.max.Y()),
                         std::min(L.max.Z(), R.max.Z()))};
    case ShapeBinaryOp::kSubtraction:
        // kSubtraction is right - left; result is contained within right
        return R;
    }
    return L; // unreachable
}

// CSG ray intersection via interval arithmetic.
// Each child's Intersect2 returns a signed [t_enter, t_exit] directly,
// so no probe ray or SignedDistance call is needed.
std::pair<double, double> NuGeom::CombinedShape::Intersect2Impl(const Ray &ray) const {
    static constexpr double eps = 1e-9;
    static constexpr double inf = std::numeric_limits<double>::infinity();

    const auto iA = m_left->Intersect2(ray);
    const auto iB = m_right->Intersect2(ray);
    const double a_in = iA.first, a_out = iA.second;
    const double b_in = iB.first, b_out = iB.second;

    // Is the point at ray parameter t inside the combined shape?
    auto inside = [&](double t) -> bool {
        const bool ia = (t >= a_in && t < a_out);
        const bool ib = (t >= b_in && t < b_out);
        switch(m_op) {
        case ShapeBinaryOp::kUnion:
            return ia || ib;
        case ShapeBinaryOp::kIntersect:
            return ia && ib;
        case ShapeBinaryOp::kSubtraction:
            return ib && !ia; // right − left
        }
        return false; // unreachable
    };

    // Walk all four boundary times sorted; find entry and exit transitions.
    // t_enter = -inf when the origin is already inside (starts_inside).
    double t_enter = inside(0.0) ? -inf : inf;
    double t_exit = inf;

    double times[4] = {a_in, a_out, b_in, b_out};
    std::sort(times, times + 4);
    for(const double t : times) {
        if(!std::isfinite(t)) continue;
        const bool was_inside = inside(t - eps);
        const bool now_inside = inside(t + eps);
        if(!was_inside && now_inside && t_enter == inf) t_enter = t;
        if(was_inside && !now_inside && t_enter != inf) {
            t_exit = t;
            break;
        }
    }

    if(t_enter == inf) return {inf, inf};
    return {t_enter, t_exit};
}

// Monte Carlo volume estimate; result is cached in m_volume.
double NuGeom::CombinedShape::Volume() const {
    if(m_volume != 0) return m_volume;

    const BoundingBox bb = GetBoundingBox();
    const double bb_vol = bb.Volume();
    if(bb_vol <= 0) return 0;

    static constexpr size_t nsamples = 1 << 22; // ~4M samples
    auto rng = Random::Instance();
    size_t count = 0;
    for(size_t i = 0; i < nsamples; ++i) {
        const Vector3D p(rng.Uniform(bb.min.X(), bb.max.X()), rng.Uniform(bb.min.Y(), bb.max.Y()),
                         rng.Uniform(bb.min.Z(), bb.max.Z()));
        if(SignedDistance(p) <= 0) ++count;
    }

    m_volume = bb_vol * static_cast<double>(count) / static_cast<double>(nsamples);
    return m_volume;
}

std::unique_ptr<NuGeom::Shape> NuGeom::Box::Construct(const pugi::xml_node &node) {
    // Load the box parameters
    double x = node.attribute("x").as_double();
    double y = node.attribute("y").as_double();
    double z = node.attribute("z").as_double();
    Vector3D params(x, y, z);

    // Convert the units
    std::string unit = node.attribute("unit").value();
    if(unit == "m") {
        params *= 100;
    } else if(unit == "mm") {
        params /= 10;
    }

    return std::make_unique<NuGeom::Box>(params);
}

NuGeom::BoundingBox NuGeom::Box::GetBoundingBox() const {
    return {-m_params, m_params};
}

double NuGeom::Box::SignedDistance(const Vector3D &in_point) const {
    auto point = TransformPoint(in_point);
    Vector3D q = point.Abs() - m_params;
    return q.Max().Norm() + std::min(q.MaxComponent(), 0.0);
}

std::pair<double, double> NuGeom::Box::Intersect2Impl(const Ray &ray) const {
    static constexpr double inf = std::numeric_limits<double>::infinity();
    const double tx1 = (-m_params.X() - ray.Origin().X()) / ray.Direction().X();
    const double tx2 = (m_params.X() - ray.Origin().X()) / ray.Direction().X();
    const double ty1 = (-m_params.Y() - ray.Origin().Y()) / ray.Direction().Y();
    const double ty2 = (m_params.Y() - ray.Origin().Y()) / ray.Direction().Y();
    const double tz1 = (-m_params.Z() - ray.Origin().Z()) / ray.Direction().Z();
    const double tz2 = (m_params.Z() - ray.Origin().Z()) / ray.Direction().Z();
    const double tmin = std::max({std::min(tx1, tx2), std::min(ty1, ty2), std::min(tz1, tz2)});
    const double tmax = std::min({std::max(tx1, tx2), std::max(ty1, ty2), std::max(tz1, tz2)});
    if(tmin > tmax) return {inf, inf};
    return {tmin, tmax};
}

std::unique_ptr<NuGeom::Shape> NuGeom::Sphere::Construct(const pugi::xml_node &node) {
    // Load the box parameters
    double radius = node.attribute("r").as_double();

    // Convert the units
    std::string unit = node.attribute("unit").value();
    if(unit == "m") {
        radius *= 100;
    } else if(unit == "mm") {
        radius /= 10;
    }

    return std::make_unique<NuGeom::Sphere>(radius);
}

NuGeom::BoundingBox NuGeom::Sphere::GetBoundingBox() const {
    return {-Vector3D(m_radius, m_radius, m_radius), Vector3D(m_radius, m_radius, m_radius)};
}

double NuGeom::Sphere::SignedDistance(const Vector3D &in_point) const {
    auto point = TransformPoint(in_point);
    return point.Norm() - m_radius;
}

std::pair<double, double> NuGeom::Sphere::Intersect2Impl(const Ray &ray) const {
    static constexpr double inf = std::numeric_limits<double>::infinity();
    const double a = ray.Direction() * ray.Direction();
    const double b = 2 * ray.Origin() * ray.Direction();
    const double c = ray.Origin() * ray.Origin() - m_radius;
    const double det = b * b - 4 * a * c;
    if(det < 0) return {inf, inf};
    const double sq = std::sqrt(det);
    const double ta = 2 * c / (-b - sq);
    const double tb = 2 * c / (-b + sq);
    return {std::min(ta, tb), std::max(ta, tb)};
}

// TODO: Handle deltaphi and rmin??
std::unique_ptr<NuGeom::Shape> NuGeom::Cylinder::Construct(const pugi::xml_node &node) {
    // Load the box parameters
    double radius = node.attribute("rmax").as_double();
    double height = node.attribute("z").as_double();

    // Convert the units
    std::string unit = node.attribute("unit").value();
    if(unit == "m") {
        height *= 100;
        radius *= 100;
    } else if(unit == "mm") {
        radius /= 10;
        height /= 10;
    }

    return std::make_unique<NuGeom::Cylinder>(radius, height);
}

// TODO: Ensure that Cylinder is centered at {0, 0, 0} and not {0, 0, height/2}
NuGeom::BoundingBox NuGeom::Cylinder::GetBoundingBox() const {
    return {-Vector3D(m_radius, m_radius, m_height), Vector3D(m_radius, m_radius, m_height)};
}

double NuGeom::Cylinder::SignedDistance(const Vector3D &in_point) const {
    auto point = TransformPoint(in_point);
    Vector2D q = Vector2D(Vector2D(point.X(), point.Y()).Norm(), std::abs(point.Z())) -
                 Vector2D(m_radius, m_height);
    return q.Max().Norm() + std::min(q.MaxComponent(), 0.0);
}

std::pair<double, double> NuGeom::Cylinder::Intersect2Impl(const Ray &ray) const {
    static constexpr double inf = std::numeric_limits<double>::infinity();

    // Infinite cylinder interval in XY
    const double a =
        ray.Direction().X() * ray.Direction().X() + ray.Direction().Y() * ray.Direction().Y();
    const double b =
        2 * (ray.Direction().X() * ray.Origin().X() + ray.Direction().Y() * ray.Origin().Y());
    const double c =
        ray.Origin().X() * ray.Origin().X() + ray.Origin().Y() * ray.Origin().Y() - m_radius;

    double tc1, tc2;
    if(a < 1e-12) {
        if(c > 0) return {inf, inf}; // ray parallel to axis and outside
        tc1 = -inf;
        tc2 = inf;
    } else {
        const double det = b * b - 4 * a * c;
        if(det < 0) return {inf, inf};
        const double sq = std::sqrt(det);
        const double r1 = 2 * c / (-b - sq);
        const double r2 = 2 * c / (-b + sq);
        tc1 = std::min(r1, r2);
        tc2 = std::max(r1, r2);
    }

    // Z-slab [0, m_height] (preserving existing axis convention)
    double tz1, tz2;
    const double dz = ray.Direction().Z();
    if(std::abs(dz) < 1e-12) {
        if(ray.Origin().Z() < 0 || ray.Origin().Z() > m_height) return {inf, inf};
        tz1 = -inf;
        tz2 = inf;
    } else {
        const double t0 = -ray.Origin().Z() / dz;
        const double t1 = (m_height - ray.Origin().Z()) / dz;
        tz1 = std::min(t0, t1);
        tz2 = std::max(t0, t1);
    }

    const double tmin = std::max(tc1, tz1);
    const double tmax = std::min(tc2, tz2);
    if(tmin > tmax) return {inf, inf};
    return {tmin, tmax};
}

NuGeom::Polyhedra::Polyhedra(double startphi, double deltaphi, size_t nsides,
                             const std::vector<zplane> &planes, const Rotation3D &rotation,
                             const Translation3D &translation)
    : Shape(rotation, translation), m_startphi{startphi}, m_deltaphi{deltaphi}, m_nsides{nsides},
      m_planes{planes} {
    if(m_planes.size() < 2) throw std::runtime_error("Polyhedra requires at least two zplanes.");

    if(m_nsides < 3) throw std::runtime_error("Polyhedra must have at least 3 sides.");

    // Sort the zplanes
    std::sort(m_planes.begin(), m_planes.end(),
              [](const zplane &a, const zplane &b) { return a.z < b.z; });

    // Build boundaries
    double dphi_side = m_deltaphi / static_cast<double>(m_nsides);
    std::vector<double> phi(m_nsides);
    for(size_t i = 0; i < m_nsides; ++i) phi[i] = m_startphi + static_cast<double>(i) * dphi_side;

    auto vertex = [](double r, double phi_, double z) -> Vector3D {
        return Vector3D(r * std::cos(phi_), r * std::sin(phi_), z);
    };

    double zmin = m_planes.front().z;
    double zmax = m_planes.back().z;

    // Bottom cap
    m_boundaries.push_back({Vector3D(0, 0, -1), zmin});
    // Top cap
    m_boundaries.push_back({Vector3D(0, 0, 1), -zmax});

    // Faces
    for(size_t iplane = 0; iplane + 1 < m_planes.size(); ++iplane) {
        double z0 = m_planes[iplane].z;
        double z1 = m_planes[iplane + 1].z;

        double r0 = m_planes[iplane].rmax;
        double r1 = m_planes[iplane + 1].rmax;

        for(size_t index1 = 0; index1 < m_nsides; ++index1) {
            size_t index2 = (index1 + 1) % m_nsides;

            // Define corners of trapezoid
            Vector3D vertex0 = vertex(r0, phi[index1], z0);
            Vector3D vertex1 = vertex(r0, phi[index2], z0);
            Vector3D vertex3 = vertex(r1, phi[index1], z1);

            // Construct planes
            Vector3D edge1 = vertex1 - vertex0;
            Vector3D edge2 = vertex3 - vertex0;
            Vector3D normal = edge1.Cross(edge2).Unit();

            m_boundaries.push_back({normal, -normal.Dot(vertex0)});
        }
    }

    // Add wedge cut planes for deltaphi < 2pi
    if(m_deltaphi < 2.0 * M_PI - 1e-6) {
        double phi_end = m_startphi + m_deltaphi;
        // Start plane
        Vector3D nstart(std::sin(m_startphi), -std::cos(m_startphi), 0.0);
        m_boundaries.push_back({nstart, 0.0});
        // End place
        Vector3D nend(-std::sin(phi_end), std::cos(phi_end), 0.0);
        m_boundaries.push_back({nend, 0.0});
    }

    // Ensure normals point outward
    Vector3D test_point(0, 0, (zmin + zmax) / 2);
    for(auto &plane : m_boundaries) {
        double side = plane.normal.Dot(test_point) + plane.distance;
        if(side > 0) {
            plane.normal = -plane.normal;
            plane.distance = -plane.distance;
        }
    }
}

// TODO: Extend to handle rmin
std::unique_ptr<NuGeom::Shape> NuGeom::Polyhedra::Construct(const pugi::xml_node &node) {
    // Convert the units
    std::string unit = node.attribute("lunit").value();
    double lconversion = 1;
    if(unit == "m") {
        lconversion = 100;
    } else if(unit == "mm") {
        lconversion = 0.1;
    } else if(unit != "cm") {
        throw std::runtime_error(fmt::format("Invalid lunit {} found", unit));
    }

    // Convert the angles
    std::string angle = node.attribute("aunit").value();
    double aconversion = 1;
    if(angle == "deg" || angle == "degree") {
        aconversion = M_PI / 180;
    } else if(angle == "rad" || angle == "radian") {
        aconversion = 1;
    } else {
        throw std::runtime_error(fmt::format("Invalid aunit {} found", unit));
    }

    // Load the polyhedra parameters
    auto nsides = node.attribute("numsides").as_ullong();
    double startphi = node.attribute("startphi").as_double() * aconversion;
    double deltaphi = node.attribute("deltaphi").as_double() * aconversion;

    std::vector<zplane> zplanes;
    for(auto zp : node.children("zplane")) {
        zplanes.push_back({zp.attribute("rmin").as_double() * lconversion,
                           zp.attribute("rmax").as_double() * lconversion,
                           zp.attribute("z").as_double() * lconversion});
    }

    return std::make_unique<NuGeom::Polyhedra>(startphi, deltaphi, nsides, zplanes);
}

NuGeom::BoundingBox NuGeom::Polyhedra::GetBoundingBox() const {
    double max_rmax = 0.0;
    for(const auto &p : m_planes) max_rmax = std::max(max_rmax, p.rmax);
    double zmin = m_planes.front().z;
    double zmax = m_planes.back().z;
    return {Vector3D(-max_rmax, -max_rmax, zmin), Vector3D(max_rmax, max_rmax, zmax)};
}

double NuGeom::Polyhedra::SignedDistance(const Vector3D &in_point) const {
    auto point = TransformPoint(in_point);
    double max_dist = -std::numeric_limits<double>::infinity();

    for(const auto &plane : m_boundaries) {
        double dist = plane.normal.Dot(point) + plane.distance;
        max_dist = std::max(max_dist, dist);
    }

    return max_dist;
}

double NuGeom::Polyhedra::Volume() const {
    if(m_planes.size() < 2) return 0.0;

    double total = 0.0;
    const double factor =
        0.5 * static_cast<double>(m_nsides) * std::sin(2.0 * M_PI / static_cast<double>(m_nsides));
    for(size_t i = 0; i + 1 < m_planes.size(); ++i) {
        double z0 = m_planes[i].z;
        double z1 = m_planes[i + 1].z;
        double h = z1 - z0;

        double r0 = m_planes[i].rmax;
        double r1 = m_planes[i + 1].rmax;

        double A0 = factor * r0 * r0;
        double A1 = factor * r1 * r1;

        total += h / 3.0 * (A0 + A1 + std::sqrt(A0 * A1));
    }

    if(m_deltaphi < 2.0 * M_PI) total *= m_deltaphi / (2.0 * M_PI);

    return total;
}

std::pair<double, double> NuGeom::Polyhedra::Intersect2Impl(const Ray &ray) const {
    static constexpr double inf = std::numeric_limits<double>::infinity();
    constexpr double EPS = 1e-12;

    double tmin = -inf;
    double tmax = inf;

    for(const auto &plane : m_boundaries) {
        const double denom = plane.normal.Dot(ray.Direction());
        const double dist = plane.normal.Dot(ray.Origin()) + plane.distance;

        if(std::abs(denom) < EPS) {
            if(dist > 0.0) return {inf, inf}; // parallel and outside
            continue;
        }

        const double t = -dist / denom;
        if(denom < 0.0)
            tmin = std::max(tmin, t); // entering half-space
        else
            tmax = std::min(tmax, t); // exiting half-space

        if(tmin > tmax) return {inf, inf};
    }

    return {tmin, tmax};
}

NuGeom::Trapezoid::Trapezoid(double x1, double x2, double y1, double y2, double z,
                             const Rotation3D &rotation, const Translation3D &translation)
    : Shape(rotation, translation), m_x{x1, x2}, m_y{y1, y2}, m_z{z} {
    m_ax = (m_x.second - m_x.first) / (2.0 * m_z);
    m_bx = 0.5 * (m_x.first + m_x.second);

    m_ay = (m_y.second - m_y.first) / (2.0 * m_z);
    m_by = 0.5 * (m_y.first + m_y.second);
}

// TODO: Extend to handle rmin
std::unique_ptr<NuGeom::Shape> NuGeom::Trapezoid::Construct(const pugi::xml_node &node) {
    // Convert the units
    std::string unit = node.attribute("lunit").value();
    double lconversion = 1;
    if(unit == "m") {
        lconversion = 100;
    } else if(unit == "mm") {
        lconversion = 0.1;
    } else if(unit != "cm") {
        throw std::runtime_error(fmt::format("Invalid lunit {} found", unit));
    }

    // Load the trapezoid parameters
    double x1 = node.attribute("x1").as_double() * lconversion * 0.5;
    double x2 = node.attribute("x2").as_double() * lconversion * 0.5;
    double y1 = node.attribute("y1").as_double() * lconversion * 0.5;
    double y2 = node.attribute("y2").as_double() * lconversion * 0.5;
    double z = node.attribute("z").as_double() * lconversion * 0.5;

    return std::make_unique<NuGeom::Trapezoid>(x1, x2, y1, y2, z);
}

NuGeom::BoundingBox NuGeom::Trapezoid::GetBoundingBox() const {
    double xmax = std::max(m_x.first, m_x.second);
    double ymax = std::max(m_y.first, m_y.second);
    return {Vector3D(-xmax, -ymax, -m_z), Vector3D(xmax, ymax, m_z)};
}

double NuGeom::Trapezoid::SignedDistance(const Vector3D &point) const {
    double dz = std::abs(point.Z()) - m_z;

    double dx1 = point.X() - (m_ax * point.Z() + m_bx);
    double dx2 = -point.X() - (m_ax * point.Z() + m_bx);

    double dy1 = point.Y() - (m_ay * point.Z() + m_by);
    double dy2 = -point.Y() - (m_ay * point.Z() + m_by);

    return std::max({dz, dx1, dx2, dy1, dy2});
}

double NuGeom::Trapezoid::Volume() const {
    double x1 = m_x.first;
    double x2 = m_x.second;
    double y1 = m_y.first;
    double y2 = m_y.second;

    return (4.0 * m_z / 3.0) * (2 * x1 * y1 + 2 * x2 * y2 + x1 * y2 + x2 * y1);
}

std::pair<double, double> NuGeom::Trapezoid::Intersect2Impl(const Ray &ray) const {
    static constexpr double inf = std::numeric_limits<double>::infinity();
    double tmin = -inf;
    double tmax = inf;

    auto update = [&](const Vector3D &normal, double dist) -> bool {
        const double denom = normal.Dot(ray.Direction());
        const double numer = -(normal.Dot(ray.Origin()) + dist);
        if(std::abs(denom) < 1e-12) return numer >= 0; // parallel: ok if inside
        const double t = numer / denom;
        if(denom < 0)
            tmin = std::max(tmin, t);
        else
            tmax = std::min(tmax, t);
        return tmin <= tmax;
    };

    if(!update({0, 0, 1}, -m_z)) return {inf, inf};       // Z+
    if(!update({0, 0, -1}, -m_z)) return {inf, inf};      // Z-
    if(!update({1, 0, -m_ax}, -m_bx)) return {inf, inf};  // X+
    if(!update({-1, 0, -m_ax}, -m_bx)) return {inf, inf}; // X-
    if(!update({0, 1, -m_ay}, -m_by)) return {inf, inf};  // Y+
    if(!update({0, -1, -m_ay}, -m_by)) return {inf, inf}; // Y-

    return {tmin, tmax};
}
