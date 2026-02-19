#pragma once

#include "geom/BoundingBox.hh"
#include "geom/Transform3D.hh"
#include "geom/Vector3D.hh"

#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <utility>

namespace pugi {
class xml_node;
}

namespace NuGeom {

class Ray;

enum class Location { kInterior, kSurface, kExterior };

enum class ShapeBinaryOp { kUnion, kIntersect, kSubtraction };

class Shape {
  public:
    Shape(const Rotation3D &rot = Rotation3D(), const Translation3D &trans = Translation3D())
        : m_rotation{rot.Inverse()}, m_translation{trans.Inverse()} {
        identity_transform = m_rotation.IsIdentity() && m_translation.IsIdentity();
    }
    virtual ~Shape() = default;

    /// Calculates if the shape contains the point
    ///@param point: The point to check the distance from the surface
    ///@return Location: Returns if the point is interior, surface, or exterior
    Location Contains(const Vector3D &) const;

    /// Calculates the signed distance a point is from the surface
    /// Negative values mean inside, positive values mean outside, and 0 means on the surface
    ///@param point: The point to check the distance from the surface
    ///@return double: The signed distance from the surface
    virtual double SignedDistance(const Vector3D &) const = 0;

    /// Finds the first forward intersection of the ray with the shape.
    ///@return The time t ≥ 0, or +∞ if no intersection.
    double Intersect(const Ray &in_ray) const;

    /// Returns the signed entry and exit times along the ray.
    /// t_enter < 0 means the ray origin is already inside the shape.
    /// {+∞, +∞} means the ray misses entirely.
    std::pair<double, double> Intersect2(const Ray &in_ray) const;

    virtual BoundingBox GetBoundingBox() const = 0;

    /// Returns GetBoundingBox() transformed by this shape's own rotation/translation.
    /// Use this when building parent-space AABBs in the BVH.
    BoundingBox GetTransformedBoundingBox() const;

    void SetRotation(const Rotation3D &rot) { m_rotation = rot.Inverse(); }
    void SetTranslation(const Translation3D &trans) { m_translation = trans.Inverse(); }
    virtual double Volume() const = 0;

  protected:
    Vector3D TransformPoint(const Vector3D &) const;
    Ray TransformRay(const Ray &) const;
    std::pair<double, double> SolveQuadratic(double, double, double) const;

  private:
    virtual std::pair<double, double> Intersect2Impl(const Ray &) const = 0;
    Transform3D m_rotation;
    Transform3D m_translation;
    bool identity_transform{false};
};

class ShapeFactory {
    using Constructor = std::function<std::unique_ptr<Shape>(const pugi::xml_node &)>;

    static std::map<std::string, Constructor> &Registry() {
        static std::map<std::string, Constructor> registry;
        return registry;
    }

  public:
    static std::unique_ptr<Shape> Initialize(const std::string &name, const pugi::xml_node &node) {
        auto constructor = Registry().at(name);
        return constructor(node);
    }

    // TODO: Switch to using a logger!
    template <class Derived> static void Register(const std::string &name) {
        if(IsRegistered(name)) std::cerr << name << " is already registered!\n";
        Registry()[name] = Derived::Construct;
    }

    static bool IsRegistered(const std::string &name) {
        return Registry().find(name) != Registry().end();
    }

    static void DisplayShapes() {
        std::cout << "Shapes:\n";
        for(const auto &registered : Registry()) std::cout << "  - " << registered.first << "\n";
    }
};

template <class Derived> class RegistrableShape {
  protected:
    RegistrableShape() = default;
    virtual ~RegistrableShape() {
        if(!m_registered) std::cerr << "Error registering shape\n";
    }

    static bool Register() {
        ShapeFactory::template Register<Derived>(Derived::Name());
        return true;
    }

  private:
    static const bool m_registered;
};
template <class Derived>
const bool RegistrableShape<Derived>::m_registered = RegistrableShape<Derived>::Register();

class CombinedShape : public Shape, RegistrableShape<CombinedShape> {
  public:
    CombinedShape(std::shared_ptr<Shape> left, std::shared_ptr<Shape> right, ShapeBinaryOp op,
                  const Rotation3D &rotation = Rotation3D(),
                  const Translation3D &translation = Translation3D())
        : Shape(rotation, translation), m_left{std::move(left)}, m_right{std::move(right)},
          m_op{op} {}

    static std::string Name() { return "CombinedShape"; }
    static std::unique_ptr<Shape> Construct(const pugi::xml_node &node);

    BoundingBox GetBoundingBox() const override;
    double SignedDistance(const Vector3D &) const override;
    double Volume() const override;

  private:
    std::pair<double, double> Intersect2Impl(const Ray &) const override;
    std::shared_ptr<Shape> m_left, m_right;
    ShapeBinaryOp m_op;
    mutable double m_volume{0};
};

class Box : public Shape, RegistrableShape<Box> {
  public:
    /// Initialize a box with one corner at (-x/2,-y/2,-z/2) and the other at (x/2,y/2,z/2)
    /// Then rotates the box, and translates the box
    ///@param dimensions: The width, depth, and height of the box
    ///@param rot: The rotation matrix of the box
    ///@param trans: The translation of the box from the origin
    Box(const Vector3D &size = Vector3D(1, 1, 1), const Rotation3D &rotation = Rotation3D(),
        const Translation3D &translation = Translation3D())
        : Shape(rotation, translation), m_params{size.X() / 2, size.Y() / 2, size.Z() / 2} {}

    static std::string Name() { return "box"; }
    static std::unique_ptr<Shape> Construct(const pugi::xml_node &node);

    BoundingBox GetBoundingBox() const override;
    double SignedDistance(const Vector3D &) const override;
    double Volume() const override { return m_params.X() * m_params.Y() * m_params.Z() * 8; }

  private:
    std::pair<double, double> Intersect2Impl(const Ray &) const override;
    Vector3D m_params;
};

class Sphere : public Shape, RegistrableShape<Sphere> {
  public:
    /// Initialize a sphere center at the origin with the given radius
    /// Then rotates the sphere, and translates the sphere
    ///@param radius: The radius of the sphere
    ///@param rot: The rotation matrix of the sphere
    ///@param trans: The translation of the sphere from the origin
    Sphere(double radius = 1, const Rotation3D &rotation = Rotation3D(),
           const Translation3D &translation = Translation3D())
        : Shape(rotation, translation), m_radius{radius} {}

    static std::string Name() { return "orb"; }
    static std::unique_ptr<Shape> Construct(const pugi::xml_node &node);

    BoundingBox GetBoundingBox() const override;
    double SignedDistance(const Vector3D &) const override;
    double Volume() const override { return m_radius * m_radius * m_radius * 4 * M_PI / 3.0; }

  private:
    std::pair<double, double> Intersect2Impl(const Ray &) const override;
    double m_radius;
};

class Cylinder : public Shape, RegistrableShape<Cylinder> {
  public:
    /// Initialize a cylinder centered at the origin with the given radius and height
    /// Then rotates the cylinder, and translates the cylinder
    ///@param radius: The radius of the cylinder
    ///@param height: The height of the cylinder
    ///@param rot: The rotation matrix of the cylinder
    ///@param trans: The translation of the cylinder from the origin
    Cylinder(double radius = 1, double height = 1, const Rotation3D &rotation = Rotation3D(),
             const Translation3D &translation = Translation3D())
        : Shape(rotation, translation), m_radius{radius}, m_height{height} {}

    static std::string Name() { return "tube"; }
    static std::unique_ptr<Shape> Construct(const pugi::xml_node &node);

    BoundingBox GetBoundingBox() const override;
    double SignedDistance(const Vector3D &) const override;
    double Volume() const override { return m_radius * m_radius * m_height * M_PI; }

  private:
    std::pair<double, double> Intersect2Impl(const Ray &) const override;
    double m_radius;
    double m_height;
};

/// zplanes for defining polyhedra
struct zplane {
    double rmin, rmax, z;
};

class Polyhedra : public Shape, RegistrableShape<Polyhedra> {
  public:
    /// Initialize a polyhedra centered at the origin with a given set of sides and heights
    /// Then rotates the polyhedra, and translates the polyhedra
    ///@param startphi: The starting phi angle
    ///@param deltaphi: The deltaphi of the shape
    ///@param nsides: The number of sides for the polyhedra
    ///@param planes: The list of z-planes defining the polyhedra
    ///@param rot: The rotation matrix of the polyhedra
    ///@param trans: The translation of the polyhedra from the origin
    Polyhedra(double startphi, double deltaphi, size_t nsides, const std::vector<zplane> &planes,
              const Rotation3D &rotation = Rotation3D(),
              const Translation3D &translation = Translation3D());

    static std::string Name() { return "polyhedra"; }
    static std::unique_ptr<Shape> Construct(const pugi::xml_node &node);

    BoundingBox GetBoundingBox() const override;
    double SignedDistance(const Vector3D &) const override;
    double Volume() const override;

  private:
    std::pair<double, double> Intersect2Impl(const Ray &) const override;

    struct Plane {
        Vector3D normal;
        double distance;
    };

    double m_startphi, m_deltaphi;
    size_t m_nsides;
    std::vector<zplane> m_planes;
    std::vector<Plane> m_boundaries;
};

class Trapezoid : public Shape, RegistrableShape<Trapezoid> {
  public:
    /// Initialize a trapezoid centered at the origin
    /// Then rotates the trapezoid, and translates the trapezoid
    ///@param x1: The x half-length at -z
    ///@param x2: The x half-length at z
    ///@param y1: The y half-length at -z
    ///@param y2: The y half-length at z
    ///@param z: The z half-length
    ///@param rot: The rotation matrix of the polyhedra
    ///@param trans: The translation of the polyhedra from the origin
    Trapezoid(double x1, double x2, double y1, double y2, double z,
              const Rotation3D &rotation = Rotation3D(),
              const Translation3D &translation = Translation3D());

    static std::string Name() { return "trd"; }
    static std::unique_ptr<Shape> Construct(const pugi::xml_node &node);

    BoundingBox GetBoundingBox() const override;
    double SignedDistance(const Vector3D &) const override;
    double Volume() const override;

  private:
    std::pair<double, double> Intersect2Impl(const Ray &) const override;
    std::pair<double, double> m_x, m_y;
    double m_z, m_ax, m_bx, m_ay, m_by;
};

// TODO: Implement details
/*
class Cone : public Shape, RegistrableShape<Cone> {
    public:
        /// Initialize a cone centered at the origin with the given normal vector
        /// Then rotates the cone, and translates the cone
        ///@param normal: The surface normal to the cone
        ///@param rot: The rotation matrix of the cone
        ///@param trans: The translation of the cone from the origin
        Cone(Vector3D normal = Vector3D(),
              const Rotation3D &rotation = Rotation3D(),
              const Translation3D &translation = Translation3D())
            : Shape(rotation, translation), m_normal{normal.Unit()} {}

        static std::string Name() { return "cone"; }
        static std::unique_ptr<Shape> Construct(const pugi::xml_node &node);

        double SignedDistance(const Vector3D&) const override;

    private:
        Vector3D m_normal;

};*/

} // namespace NuGeom
