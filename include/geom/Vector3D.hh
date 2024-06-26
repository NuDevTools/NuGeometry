#pragma once

#include <array>
#include <cmath>

namespace NuGeom {

class Vector3D {
    public:
        Vector3D() : m_vec{} {}
        constexpr Vector3D(double x, double y, double z) : m_vec{x, y, z} {}
        Vector3D(std::array<double, 3> vec) : m_vec{vec} {}
        Vector3D(const Vector3D&) = default;
        Vector3D(Vector3D&&) = default;

        Vector3D& operator=(const Vector3D&) = default;
        Vector3D& operator=(Vector3D&&) = default;

        // const access
        const double& X() const { return m_vec[0]; }
        const double& Y() const { return m_vec[1]; }
        const double& Z() const { return m_vec[2]; }
        const double& R() const { return m_vec[0]; }
        const double& G() const { return m_vec[1]; }
        const double& B() const { return m_vec[2]; }

        // non-const access
        double& X() { return m_vec[0]; }
        double& Y() { return m_vec[1]; }
        double& Z() { return m_vec[2]; }
        double& R() { return m_vec[0]; }
        double& G() { return m_vec[1]; }
        double& B() { return m_vec[2]; }

        // Functions
        double Dot(const Vector3D&) const;
        Vector3D Cross(const Vector3D&) const;
        double Norm2() const { return Dot(*this); }
        double Norm() const { return sqrt(Norm2()); }
        Vector3D Unit() const;
        Vector3D Abs() const { return {std::abs(X()), std::abs(Y()), std::abs(Z())}; }
        Vector3D Max(const Vector3D& = Vector3D()) const;
        double MaxComponent() const;

        // Operators
        friend Vector3D operator*(double, const Vector3D&);
        friend Vector3D operator/(const Vector3D&, double);
        friend Vector3D operator/(double, const Vector3D&);
        const double& operator[](size_t i) const { return m_vec[i]; }
        double& operator[](size_t i) { return m_vec[i]; }
   
        bool operator==(const Vector3D &other) const {
            return m_vec == other.m_vec;
        }
        bool operator!=(const Vector3D &other) const {
            return !(*this == other);
        }

        Vector3D& operator*=(double scale) {
            m_vec[0] *= scale;
            m_vec[1] *= scale;
            m_vec[2] *= scale;

            return *this;
        }
        Vector3D& operator/=(double scale) {
            return *this *= 1.0/scale;
        }
        Vector3D& operator+=(const Vector3D &other) {
            m_vec[0] += other.m_vec[0];
            m_vec[1] += other.m_vec[1];
            m_vec[2] += other.m_vec[2];

            return *this;
        }
        Vector3D& operator-=(const Vector3D &other) {
            return *this += -other;
        }
        Vector3D operator*(double scale) const {
            return Vector3D{*this} *= scale;
        }
        double operator*(const Vector3D &other) const {
            return Dot(other);
        }
        Vector3D operator+(const Vector3D &other) const {
            return Vector3D{*this} += other;
        }
        Vector3D operator-(const Vector3D &other) const {
            return Vector3D{*this} -= other;
        }
        Vector3D operator-() const {
            return {-m_vec[0], -m_vec[1], -m_vec[2]};
        }

        template<typename OStream>
        friend OStream& operator<<(OStream &os, const Vector3D &vec) {
            os << "Vector3D(" << vec.X() << ", " << vec.Y() << ", " << vec.Z() << ")";
            return os;
        }

    private:
        std::array<double, 3> m_vec;
};

Vector3D operator*(double, const Vector3D&);
Vector3D operator/(const Vector3D&, double);
Vector3D operator/(double, const Vector3D&);

static constexpr Vector3D UnitX = Vector3D(1, 0, 0);
static constexpr Vector3D UnitY = Vector3D(0, 1, 0);
static constexpr Vector3D UnitZ = Vector3D(0, 0, 1);

}
