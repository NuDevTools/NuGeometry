#include "geom/InteractionViz.hh"
#include "geom/Random.hh"

#include "spdlog/spdlog.h"

#include <cmath>
#include <limits>
#include <numeric>

using namespace NuGeom;

RayInteractionSim::RayInteractionSim(World world) : m_world{std::move(world)} {}

void RayInteractionSim::SetConstantCrossSection(double sigma) {
    m_const_xsec = sigma;
    m_xsec = nullptr; // fall back to the constant model
}

void RayInteractionSim::SetSecondaryModel(int n, double cone_angle) {
    m_nsecondaries = n;
    m_cone_angle = cone_angle;
    m_secondaries = nullptr; // fall back to the built-in model
}

double RayInteractionSim::MeanFreePath(double energy, const Material &mat) const {
    double inv_mfp = 0.0;
    for(const auto &elm : mat.Elements()) {
        double sigma = m_xsec ? m_xsec(energy, elm.PDG()) : m_const_xsec;
        inv_mfp += mat.NumberDensity(elm) * sigma;
    }
    if(inv_mfp <= 0.0) return std::numeric_limits<double>::infinity();
    return 1.0 / inv_mfp;
}

std::vector<SecondaryParticle>
RayInteractionSim::DefaultSecondaries(double energy, const Vector3D &dir, const Vector3D &) const {
    std::vector<SecondaryParticle> out;
    if(m_nsecondaries <= 0) return out;

    Vector3D axis = dir.Unit();
    // Build an orthonormal basis (axis, u, v) for sampling directions in a cone.
    Vector3D helper = std::abs(axis.X()) < 0.9 ? Vector3D{1, 0, 0} : Vector3D{0, 1, 0};
    Vector3D u = axis.Cross(helper).Unit();
    Vector3D v = axis.Cross(u);

    auto rng = Random::Instance();
    for(int i = 0; i < m_nsecondaries; ++i) {
        double theta = std::abs(rng.Normal(0.0, m_cone_angle));
        double phi = rng.Uniform(0.0, 2 * M_PI);
        Vector3D d =
            std::cos(theta) * axis + std::sin(theta) * (std::cos(phi) * u + std::sin(phi) * v);
        // Share the incoming energy among the secondaries (rough, for display).
        double frac = rng.Uniform(0.1, 1.0);
        out.push_back({i, frac * energy / m_nsecondaries, d.Unit()});
    }
    return out;
}

InteractionEvent RayInteractionSim::ShootRay(double energy, const Ray &ray) const {
    InteractionEvent event;
    event.incoming = ray;
    event.energy = energy;
    event.incoming_path = m_world.GetLineSegments(ray);

    // Per-segment optical depth τᵢ = Lᵢ / λᵢ.
    std::vector<double> taus;
    taus.reserve(event.incoming_path.size());
    for(const auto &seg : event.incoming_path) {
        double mfp = MeanFreePath(energy, seg.GetMaterial());
        double len = seg.Length();
        double tau = (std::isfinite(len) && std::isfinite(mfp) && mfp > 0) ? len / mfp : 0.0;
        taus.push_back(tau);
    }
    event.optical_depth = std::accumulate(taus.begin(), taus.end(), 0.0);

    // Sample the optical depth of the first interaction; no interaction if it
    // exceeds the total optical depth available along the ray.
    double t_target = -std::log(1.0 - Random::Instance().Uniform(0.0, 1.0));
    if(t_target > event.optical_depth || event.optical_depth <= 0.0) {
        spdlog::debug("No interaction (τ_target={:.4g} > τ_total={:.4g})", t_target,
                      event.optical_depth);
        return event;
    }

    // Locate the interacting segment and the vertex within it.
    double cumulative = 0.0;
    for(size_t i = 0; i < event.incoming_path.size(); ++i) {
        if(taus[i] <= 0.0) continue;
        if(cumulative + taus[i] >= t_target) {
            const auto &seg = event.incoming_path[i];
            double frac = (t_target - cumulative) / taus[i];
            event.vertex = seg.Start() + (seg.End() - seg.Start()) * frac;
            event.vertex_material = seg.GetMaterial().Name();
            event.interacted = true;
            break;
        }
        cumulative += taus[i];
    }
    if(!event.interacted) return event; // numerical safety net

    spdlog::debug("Interaction at ({:.3g},{:.3g},{:.3g}) in {}", event.vertex.X(), event.vertex.Y(),
                  event.vertex.Z(), event.vertex_material);

    // Generate the outgoing particles and trace each through the geometry.
    auto secondaries = m_secondaries ? m_secondaries(energy, ray.Direction(), event.vertex)
                                     : DefaultSecondaries(energy, ray.Direction(), event.vertex);
    for(const auto &sp : secondaries) {
        ParticleTrack track;
        track.pdg = sp.pdg;
        track.energy = sp.energy;
        track.start = event.vertex;
        track.direction = sp.direction.Unit();
        Ray outgoing(event.vertex, track.direction, ray.POT());
        track.path = m_world.GetLineSegments(outgoing);
        track.exited_world = !track.path.empty() && std::isfinite(track.path.back().Length());
        event.secondaries.push_back(std::move(track));
    }
    return event;
}

InteractionStats RayInteractionSim::Summarize(const std::vector<InteractionEvent> &events) {
    InteractionStats stats;
    stats.nrays = events.size();
    for(const auto &e : events) {
        if(e.interacted) {
            ++stats.ninteracted;
            stats.nsecondaries += e.secondaries.size();
        }
    }
    return stats;
}
