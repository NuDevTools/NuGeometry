#pragma once

#include "geom/DetectorSim.hh"
#include "geom/Element.hh"
#include "geom/Material.hh"
#include "geom/Random.hh"
#include "geom/Ray.hh"
#include "geom/Vector3D.hh"
#include "geom/Volume.hh"
#include "geom/World.hh"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <map>
#include <string>


class TestRayGen {
  public:
    TestRayGen(double emin, double emax, NuGeom::Vector3D corner1, NuGeom::Vector3D corner2)
        : m_emin{emin}, m_emax{emax} {
        auto [xmin, xmax] = std::minmax(corner1.X(), corner2.X());
        auto [ymin, ymax] = std::minmax(corner1.Y(), corner2.Y());
        auto [zmin, zmax] = std::minmax(corner1.Z(), corner2.Z());
        m_xmin = xmin;
        m_xmax = xmax;
        m_ymin = ymin;
        m_ymax = ymax;
        m_zmin = zmin;
        m_zmax = zmax;
    }

    NuGeom::EnergyRay GetRay() const {
        auto ray = ShootRay();
        double energy = NuGeom::Random::Instance().Uniform(m_emin, m_emax);
        return {energy, ray};
    }

  private:
    NuGeom::Ray ShootRay() const {
        auto rand = NuGeom::Random::Instance();

        NuGeom::Vector3D position{rand.Uniform(m_xmin, m_xmax), rand.Uniform(m_ymin, m_ymax),
                                  rand.Uniform(m_zmin, m_zmax)};
        double costheta = rand.Uniform(0.0, 1.0);
        double sintheta = std::sqrt(1 - costheta * costheta);
        double phi = rand.Uniform(0.0, 2 * M_PI);
        NuGeom::Vector3D direction{sintheta * cos(phi), sintheta * sin(phi), costheta};

        return NuGeom::Ray(position, direction, 1e2);
    }

    double m_emin, m_emax;
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;
};

class TestEventGen {
  public:
    TestEventGen(std::map<size_t, double> xsec) : m_xsec{xsec} {}

    double CrossSection(double energy, size_t pdg) const {
        double sigma0 = m_xsec.at(pdg);
        double x0 = std::log10(0.04);
        double sigmax = 1;
        return energy * sigma0 * std::exp(-pow((std::log10(energy) - x0) / sigmax, 2));
    }

    std::map<NuGeom::Element, double>
    EvaluateCrossSections(double energy, const std::set<NuGeom::Material> &mats) const {
        std::map<NuGeom::Element, double> result;
        for(const auto &mat : mats) {
            for(size_t i = 0; i < mat.NElements(); ++i) {
                auto elem = mat.Elements()[i];
                double xsec = CrossSection(energy, elem.PDG());
                result[elem] = xsec;
            }
        }
        return result;
    }

    std::map<size_t, double> GetXSecs(double energy, const std::set<NuGeom::Element> &mats) {
        std::map<size_t, double> result;
        for(const auto &elm : mats) {
            if(result.find(elm.PDG()) == result.end()) {
                result[elm.PDG()] = CrossSection(energy, elm.PDG());
            }
        }
        return result;
    }

  private:
    std::map<size_t, double> m_xsec;
};
