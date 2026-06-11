#include "catch2/catch.hpp"

#include "geom/DetectorSim.hh"
#include "geom/FluxSource.hh"
#include "geom/Material.hh"
#include "geom/MockGenerator.hh"
#include "geom/Shape.hh"
#include "geom/Volume.hh"
#include "geom/World.hh"

#include <cstddef>
#include <map>
#include <memory>

namespace {

// Generator whose partial unweighter rejects every odd call (reject, emit,
// reject, emit, ...) -- an unweighting efficiency of exactly 1/2, realized
// deterministically.  Consistently with that, the envelope (max-weight) cross
// section is twice the total: sigma / max_w = 1/2.  Counting calls lets tests
// pin down the per-mode trial discipline that keeps the cross-section
// normalization sum(w)/N_trials correct: TotalXSecRetry must retry the same
// vertex until emit, while EnvelopeNoRetry must take a single shot per vertex.
class AlternatingGenerator : public NuGeom::GeneratorInterface {
  public:
    static constexpr double kSigma = 1e-38;

    void InitRun(std::shared_ptr<HepMC3::GenRunInfo>) override {}
    double TotalXSec(double, int, int) override { return kSigma; }
    double EnvelopeXSec(double, int, int) override { return 2 * kSigma; }
    bool GenerateEvent(HepMC3::GenEvent &evt) override {
        ++m_calls;
        if(m_calls % 2 != 0) return false;
        evt.weights() = {1.0};
        return true;
    }
    size_t Calls() const { return m_calls; }

  private:
    size_t m_calls = 0;
};

struct DetSimFixture {
    static constexpr double flux_weight = 0.5;
    static constexpr double ray_pot = 2.0;

    static void SetupSim(NuGeom::DetectorSim &sim) {
        NuGeom::Material water("Water", 1.0, 2);
        water.AddElement(NuGeom::Element("Hydrogen", 1, 1), 2);
        water.AddElement(NuGeom::Element("Oxygen", 8, 16), 1);
        auto box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{10, 10, 10});
        auto lvol = std::make_shared<NuGeom::LogicalVolume>(water, box);
        NuGeom::World world(lvol);

        sim.Setup(world);
        sim.SetFluxCallback([]() {
            NuGeom::FluxSample fs;
            fs.energy = 1.0;
            fs.pdg = 14;
            fs.ray = NuGeom::Ray({0, 0, -100}, {0, 0, 1}, ray_pot);
            fs.flux_weight = flux_weight;
            return fs;
        });
        // Far below any ray's total interaction probability, so the layer-1
        // accept ratio exceeds 1 and every thrown ray is accepted
        // deterministically.
        sim.SetMaxProb(1e-20);
    }
};

} // namespace

// POT accounting: for a constant flux weight w, constant Ray::POT() = p, and
// a given max_prob M, each thrown ray contributes w*p/M to m_pot regardless
// of accept/reject outcome.  After N throws we expect m_pot == N*w*p/M.
TEST_CASE("DetectorSim POT accounting is deterministic per thrown ray", "[DetectorSim][POT]") {
    NuGeom::Material water("Water", 1.0, 2);
    water.AddElement(NuGeom::Element("Hydrogen", 1, 1), 2);
    water.AddElement(NuGeom::Element("Oxygen", 8, 16), 1);

    auto box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{10, 10, 10});
    auto lvol = std::make_shared<NuGeom::LogicalVolume>(water, box);
    NuGeom::World world(lvol);

    const double flux_weight = 0.5;
    const double ray_pot = 2.0;

    NuGeom::DetectorSim sim(1.0); // safety factor 1: SetMaxProb stays as-is
    sim.Setup(world);
    sim.SetFluxCallback([&]() {
        NuGeom::FluxSample fs;
        fs.energy = 1.0;
        fs.pdg = 14;
        fs.ray = NuGeom::Ray({0, 0, -100}, {0, 0, 1}, ray_pot);
        fs.flux_weight = flux_weight;
        return fs;
    });
    auto gen = std::make_shared<NuGeom::MockGenerator>([](double, int, int) { return 1e-38; });
    sim.SetGenerator(gen);
    sim.SetMaxProb(1.0);

    const double M = sim.GetMaxProb();
    REQUIRE(M == Approx(1.0));

    const double increment = flux_weight * ray_pot / M;
    const double target_pot = 100 * increment;
    sim.GenerateEvents(target_pot);

    // Every thrown ray charges `increment`; the loop exits on the first
    // throw whose running total crosses target_pot.  Accumulated POT lies in
    // [target_pot, target_pot + increment).
    const double acc = sim.GetAccumulatedPOT();
    CHECK(acc >= target_pot);
    CHECK(acc < target_pot + increment);
}

// EnvelopeNoRetry: one generator trial per accepted vertex.  A rejected trial
// produces no event but the thrown ray still charges POT, so the rejection
// inefficiency stays in the normalization (this is what lets the generator
// count rejected trials toward sum(w)/N_trials for the cross section).
TEST_CASE("DetectorSim EnvelopeNoRetry takes a single trial per vertex",
          "[DetectorSim][SamplingMode]") {
    NuGeom::DetectorSim sim(1.0);
    DetSimFixture::SetupSim(sim);
    sim.SetSamplingMode(NuGeom::DetectorSim::SamplingMode::EnvelopeNoRetry);
    auto gen = std::make_shared<AlternatingGenerator>();
    sim.SetGenerator(gen);

    const double increment = DetSimFixture::flux_weight * DetSimFixture::ray_pot / sim.GetMaxProb();
    sim.GenerateEvents(std::size_t{4});

    // Every other trial is rejected, so emitting 4 events takes 8 trials and
    // 8 thrown rays: each rejection consumed (and charged) its own ray.
    CHECK(gen->Calls() == 8);
    CHECK(sim.GetAccumulatedPOT() == Approx(8 * increment));
}

// TotalXSecRetry: the generator is retried on the same vertex until it emits,
// so exactly one ray is consumed per event and the rejections cost no POT.
TEST_CASE("DetectorSim TotalXSecRetry retries the same vertex until emit",
          "[DetectorSim][SamplingMode]") {
    NuGeom::DetectorSim sim(1.0);
    DetSimFixture::SetupSim(sim);
    sim.SetSamplingMode(NuGeom::DetectorSim::SamplingMode::TotalXSecRetry);
    auto gen = std::make_shared<AlternatingGenerator>();
    sim.SetGenerator(gen);

    const double increment = DetSimFixture::flux_weight * DetSimFixture::ray_pot / sim.GetMaxProb();
    sim.GenerateEvents(std::size_t{4});

    // Two trials per vertex (reject, emit) but only one thrown ray per event.
    CHECK(gen->Calls() == 8);
    CHECK(sim.GetAccumulatedPOT() == Approx(4 * increment));
}

// The two sampling modes bookkeep the unweighter inefficiency in different
// places -- per-trial rejections at fixed POT cost (EnvelopeNoRetry) vs. a
// smaller layer-1 envelope with retries (TotalXSecRetry) -- but the physical
// event rate per POT must come out identical.  With the deterministic
// 1/2-efficiency generator (envelope = 2 sigma) and max_prob calibrated by
// Init() on a single-ray flux (safety factor 1), both modes are fully
// deterministic and the rates must agree exactly.
TEST_CASE("DetectorSim sampling modes produce the same events per POT",
          "[DetectorSim][SamplingMode]") {
    const std::size_t nevents = 50;
    std::map<NuGeom::DetectorSim::SamplingMode, double> rate;

    for(auto mode : {NuGeom::DetectorSim::SamplingMode::EnvelopeNoRetry,
                     NuGeom::DetectorSim::SamplingMode::TotalXSecRetry}) {
        NuGeom::DetectorSim sim(1.0);
        DetSimFixture::SetupSim(sim);
        sim.SetSamplingMode(mode);
        auto gen = std::make_shared<AlternatingGenerator>();
        sim.SetGenerator(gen);
        // Calibrate max_prob from the (constant) flux: with safety factor 1
        // it lands exactly on the per-ray probability for the mode's cross
        // section, so layer-1 accepts every thrown ray.
        sim.Init(10);
        sim.GenerateEvents(nevents);
        rate[mode] = static_cast<double>(nevents) / sim.GetAccumulatedPOT();
    }

    CHECK(rate[NuGeom::DetectorSim::SamplingMode::EnvelopeNoRetry] ==
          Approx(rate[NuGeom::DetectorSim::SamplingMode::TotalXSecRetry]));
}
