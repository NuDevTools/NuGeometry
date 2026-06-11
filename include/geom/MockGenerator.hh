#pragma once

#include "geom/GeneratorInterface.hh"

#include <functional>
#include <utility>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenVertex.h"
#pragma GCC diagnostic pop

namespace NuGeom {

// Minimal GeneratorInterface implementation for tests and for driving
// DetectorSim from code that only knows how to compute a total xsec.
//
// Supplies a caller-provided xsec function and, for each event, attaches a
// single outgoing "neutrino" (status 1) to the primary vertex as a stand-in
// final state.  Real adapters should do channel selection and populate the
// full final state here.
class MockGenerator : public GeneratorInterface {
  public:
    using XSecFn = std::function<double(double energy, int nu_pdg, int target_pdg)>;

    // Single cross-section fn is fine for a mock: no partial unweighter,
    // so TotalXSec and EnvelopeXSec return the same value.  Real adapters
    // would return sigma_tot vs max_w respectively.
    explicit MockGenerator(XSecFn xsec) : m_xsec{std::move(xsec)} {}

    void InitRun(std::shared_ptr<HepMC3::GenRunInfo> run_info) override {
        run_info->set_weight_names({"CV"});
        run_info->tools().push_back({"NuGeomMockGenerator", "0.1", ""});
    }

    double TotalXSec(double energy, int nu_pdg, int target_pdg) override {
        return m_xsec ? m_xsec(energy, nu_pdg, target_pdg) : 0.0;
    }
    double EnvelopeXSec(double energy, int nu_pdg, int target_pdg) override {
        return TotalXSec(energy, nu_pdg, target_pdg);
    }

    bool GenerateEvent(HepMC3::GenEvent &evt) override {
        evt.weights() = {1.0};
        // Primary vertex + incoming nu + target pre-filled by DetectorSim.
        // Attach a placeholder outgoing neutrino to exercise the pipeline.
        auto vtxs = evt.vertices();
        if(vtxs.empty()) return false;
        auto &vtx = vtxs.front();
        auto ins = vtx->particles_in();
        if(ins.empty()) return false;
        auto nu_in = ins.front();
        auto nu_out = std::make_shared<HepMC3::GenParticle>(nu_in->momentum(), nu_in->pid(), 1);
        vtx->add_particle_out(nu_out);
        return true;
    }

  private:
    XSecFn m_xsec;
};

} // namespace NuGeom
