#pragma once

#include "geom/Vector3D.hh"

namespace NuGeom {

class Ray {
public:
    Ray() = default;
    Ray(const Vector3D &origin, const Vector3D &direction, double pot, bool normalize=true) 
        : m_origin{origin}, m_direction{direction}, m_pot{pot} {
            if(normalize) m_direction = m_direction.Unit();
        }

    Vector3D Origin() const { return m_origin; }
    Vector3D Direction() const { return m_direction; }
    Vector3D Propagate(double t) const { return m_origin + t*m_direction; }
    double POT() const { return m_pot; }

private:
    Vector3D m_origin, m_direction;
    double m_pot;
};

}
