#pragma once

#include "geom/Element.hh"
#include "geom/Material.hh"
#include "geom/Random.hh"
#include "geom/Ray.hh"
#include "geom/Vector3D.hh"
#include "spdlog/spdlog.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <map>
#include <set>
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

/// Ray generator modelling a diverging neutrino beam with a circular transverse profile.
///
/// Origins are sampled uniformly in a disk of radius @p beam_radius centered on @p beam_center
/// on the upstream face at z = @p z_plane.  Directions are drawn with independent Gaussian
/// angular deviations (σ = @p sigma_theta rad) around the +z beam axis.
class BeamRayGen {
  public:
    /// @param emin/emax     Neutrino energy range [GeV]
    /// @param beam_center   Beam axis centre (X, Y used; Z overridden by z_plane)
    /// @param beam_radius   Transverse RMS radius of the beam spot on the sampling plane [cm]
    /// @param sigma_theta   RMS angular divergence of the beam [rad]
    /// @param z_plane       z-coordinate of the upstream sampling plane [cm]
    BeamRayGen(double emin, double emax, NuGeom::Vector3D beam_center, double beam_radius,
               double sigma_theta, double z_plane)
        : m_emin{emin}, m_emax{emax}, m_radius{beam_radius}, m_sigma{sigma_theta}, m_z{z_plane},
          m_center{beam_center} {}

    NuGeom::EnergyRay GetRay() const {
        auto rand = NuGeom::Random::Instance();

        // Uniform sampling in a disk (sqrt gives uniform area distribution)
        double r = m_radius * std::sqrt(rand.Uniform(0.0, 1.0));
        double phi = rand.Uniform(0.0, 2 * M_PI);
        NuGeom::Vector3D origin{m_center.X() + r * std::cos(phi), m_center.Y() + r * std::sin(phi),
                                m_z};

        // Gaussian angular spread around the +z beam axis
        double theta_x = rand.Normal(0.0, m_sigma);
        double theta_y = rand.Normal(0.0, m_sigma);
        NuGeom::Vector3D direction{theta_x, theta_y, 1.0}; // normalised by Ray ctor

        double energy = rand.Uniform(m_emin, m_emax);
        return {energy, NuGeom::Ray(origin, direction, 1.0)};
    }

  private:
    double m_emin, m_emax, m_radius, m_sigma, m_z;
    NuGeom::Vector3D m_center;
};

class FileRayGen {
  public:
    FileRayGen(double emin, double emax, const std::string &filename) : m_emin{emin}, m_emax{emax} {
        char err_buf[256] = {};
        rays = load_rays(filename, err_buf, sizeof(err_buf));
    }
    NuGeom::EnergyRay GetRay() const {
        auto rand = NuGeom::Random::Instance();
        auto ray = rand.Sample(1, rays);
        double energy = rand.Uniform(m_emin, m_emax);
        return {energy, ray[0]};
    }

  private:
    static std::vector<NuGeom::Ray> load_rays(const std::string &path, char *err_buf,
                                              size_t err_sz) {
        std::vector<NuGeom::Ray> out;
        std::ifstream f(path);
        if(!f) {
            std::snprintf(err_buf, err_sz, "Cannot open: %s", path.c_str());
            return out;
        }
        double ox, oy, oz, dx, dy, dz, pot;
        while(f >> ox >> oy >> oz >> dx >> dy >> dz >> pot)
            out.emplace_back(NuGeom::Vector3D{ox, oy, oz}, NuGeom::Vector3D{dx, dy, dz}, pot,
                             /*normalize=*/false);
        return out;
    }
    double m_emin, m_emax;
    std::vector<NuGeom::Ray> rays;
};

class TestEventGen {
  public:
    TestEventGen(std::map<size_t, double> xsec) : m_xsec{xsec} {}

    double CrossSection(double energy, size_t pdg) const {
        static std::set<size_t> invalid_elms;
        if(m_xsec.find(pdg) == m_xsec.end()) {
            if(invalid_elms.find(pdg) == invalid_elms.end()) {
                spdlog::warn("Element {} has no entry, returning zero", pdg);
                invalid_elms.emplace(pdg);
            }
            return 0;
        }
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
