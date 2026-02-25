#pragma once

#include "geom/Vector3D.hh"

namespace NuGeom {

class FourMomentum {
  public:
    FourMomentum() : m_vec{} {}
    FourMomentum(double energy,double px, double py, double pz) : m_vec{energy, px, py, pz} {}
    FourMomentum(std::array<double, 4> vec) : m_vec{vec} {}
    FourMomentum(const FourMomentum &) = default;
    FourMomentum(FourMomentum &&) = default;

    FourMomentum &operator=(const FourMomentum &) = default;
    FourMomentum &operator=(FourMomentum &&) = default;

    const double &Energy() const { return m_vec[0]; }
    const double &Px() const { return m_vec[1]; }
    const double &Py() const { return m_vec[2]; }
    const double &Pz() const { return m_vec[3]; }
    
    double Mass2() const {
        double p2 = Px() * Px() + Py() * Py() + Pz() * Pz();
        return Energy() * Energy() - p2;
    }
    double Mass() const { 
        double m2 = Mass2();
        return std::sqrt(std::abs(m2));
    }
    Vector3D Momentum() const { return {Px(), Py(), Pz()}; }
    //TODO: Need to return 3 position
    

  private:
    std::array<double, 4> m_vec;

};
}
