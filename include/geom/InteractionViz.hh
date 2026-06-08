#pragma once

#include "geom/LineSegment.hh"
#include "geom/Material.hh"
#include "geom/Ray.hh"
#include "geom/Vector3D.hh"
#include "geom/World.hh"

#include <functional>
#include <string>
#include <vector>

namespace NuGeom {

/// A single outgoing (secondary) particle produced at an interaction vertex,
/// together with the path it traces through the detector geometry.
struct ParticleTrack {
    int pdg{0};                    ///< PDG code (used for labelling / colouring)
    double energy{0.0};            ///< energy carried away [GeV]
    Vector3D start{};              ///< production point (the interaction vertex)
    Vector3D direction{};          ///< initial unit direction
    std::vector<LineSegment> path; ///< material segments traversed until leaving the world
    bool exited_world{false};      ///< whether the particle reached the world boundary
};

/// The kinematics of one secondary emitted at an interaction.  Produced by a
/// user-supplied generator so real physics can be plugged in; a simple built-in
/// generator is provided for geometry debugging / visualisation.
struct SecondaryParticle {
    int pdg{0};
    double energy{0.0};
    Vector3D direction{}; ///< need not be normalised; tracing normalises it
};

/// Generator of secondary particles: given the incoming energy, direction and
/// the interaction vertex, return the outgoing particles to be traced.
using SecondaryGenerator = std::function<std::vector<SecondaryParticle>(
    double energy, const Vector3D &direction, const Vector3D &vertex)>;

/// The complete record of shooting one ray into the detector.
struct InteractionEvent {
    Ray incoming{{0, 0, 0}, {0, 0, 1}, 1.0};
    double energy{0.0};
    std::vector<LineSegment> incoming_path; ///< material segments of the primary ray
    bool interacted{false};
    Vector3D vertex{};           ///< interaction point (valid iff interacted)
    std::string vertex_material; ///< material name at the vertex
    double optical_depth{0.0};   ///< total interaction weight Σ Lᵢ/λᵢ along the ray
    std::vector<ParticleTrack> secondaries;
};

/// Summary statistics over a batch of events.
struct InteractionStats {
    size_t nrays{0};
    size_t ninteracted{0};
    size_t nsecondaries{0};
    double interaction_fraction() const {
        return nrays ? static_cast<double>(ninteracted) / static_cast<double>(nrays) : 0.0;
    }
};

/// Shoots rays through a World, decides probabilistically whether each ray
/// interacts (using a cross-section model to derive a mean free path per
/// material), and traces the resulting outgoing particles back through the
/// geometry.  Intended for visualising and debugging detector geometry.
class RayInteractionSim {
  public:
    /// Cross section [cm^2] per nucleus as a function of (energy [GeV], target PDG).
    using XSecCallback = std::function<double(double, size_t)>;

    explicit RayInteractionSim(World world);

    /// Set the cross-section model.  When unset, a constant cross section is
    /// used (see @ref SetConstantCrossSection); the default is sized so that
    /// interactions are visible in typical detector geometries.
    void SetCrossSection(XSecCallback xsec) { m_xsec = std::move(xsec); }

    /// Use a single energy/PDG independent cross section [cm^2].  This is a
    /// debugging knob, not a physical model: larger values shorten the mean
    /// free path and produce more interactions.
    void SetConstantCrossSection(double sigma);

    /// Set the secondary-particle generator.  When unset, the built-in
    /// forward-cone generator is used (see @ref SetSecondaryModel).
    void SetSecondaryGenerator(SecondaryGenerator gen) { m_secondaries = std::move(gen); }

    /// Configure the built-in secondary generator: emit @p n particles within a
    /// cone of half-angle @p cone_angle [rad] about the incoming direction, each
    /// carrying a random fraction of the incoming energy.
    void SetSecondaryModel(int n, double cone_angle);

    /// Shoot a single ray and return the full interaction record.
    InteractionEvent ShootRay(double energy, const Ray &ray) const;

    /// Shoot @p nrays using @p gen, whose `GetRay()` returns an `EnergyRay`
    /// `{energy, ray}` (e.g. BeamRayGen / TestRayGen).
    template <typename Gen>
    std::vector<InteractionEvent> ShootRays(size_t nrays, const Gen &gen) const {
        std::vector<InteractionEvent> events;
        events.reserve(nrays);
        for(size_t i = 0; i < nrays; ++i) {
            auto er = gen.GetRay();
            events.push_back(ShootRay(er.first, er.second));
        }
        return events;
    }

    /// Aggregate summary statistics for a batch of events.
    static InteractionStats Summarize(const std::vector<InteractionEvent> &events);

    const World &GetWorld() const { return m_world; }

    /// The built-in forward-cone secondary generator (used when none is set).
    std::vector<SecondaryParticle> DefaultSecondaries(double energy, const Vector3D &dir,
                                                      const Vector3D &vertex) const;

  private:
    /// Mean free path [cm] of a particle of @p energy in @p mat; +inf if the
    /// material is inert under the current cross-section model.
    double MeanFreePath(double energy, const Material &mat) const;

    World m_world;
    XSecCallback m_xsec;
    SecondaryGenerator m_secondaries;
    double m_const_xsec{1e-24}; // cm^2, ~1 barn — visible in dense materials
    int m_nsecondaries{3};
    double m_cone_angle{0.35}; // rad
};

} // namespace NuGeom
