#include "catch2/catch.hpp"

#include "geom/DetectorSim.hh"
#include "geom/FluxSource.hh"
#include "geom/Material.hh"
#include "geom/MockGenerator.hh"
#include "geom/Random.hh"
#include "geom/Shape.hh"
#include "geom/Volume.hh"
#include "geom/World.hh"

#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
#include <numeric>
#include <vector>

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

// Init must reject a configuration where the cross section vanishes over the
// whole flux: max_prob would be exactly zero, and ProduceEvent would then
// divide by it (0/0 -> NaN) and spin forever.  This is the failure mode that a
// GeV/MeV energy-unit mismatch produced in TotalXSecRetry mode.
TEST_CASE("DetectorSim::Init throws when every ray has zero interaction probability",
          "[DetectorSim][Init]") {
    NuGeom::DetectorSim sim(1.0);
    DetSimFixture::SetupSim(sim);
    sim.SetMaxProb(0.0); // clear the fixture's seed so Init must find max_prob itself
    auto gen = std::make_shared<NuGeom::MockGenerator>([](double, int, int) { return 0.0; });
    sim.SetGenerator(gen);
    CHECK_THROWS_WITH(sim.InitFromFlux(10), Catch::Contains("max_prob"));
}

// POT accounting: POT is the beam exposure the ray represents, so each thrown
// ray contributes p/M to m_pot regardless of its flux weight or accept/reject
// outcome (the flux weight belongs in the flux/acceptance, not in POT).  After
// N throws we expect m_pot == N*p/M.
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

    const double increment = ray_pot / M; // POT does not carry the flux weight
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

    const double increment = DetSimFixture::ray_pot / sim.GetMaxProb(); // POT excludes flux weight
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

    const double increment = DetSimFixture::ray_pot / sim.GetMaxProb(); // POT excludes flux weight
    sim.GenerateEvents(std::size_t{4});

    // Two trials per vertex (reject, emit) but only one thrown ray per event.
    CHECK(gen->Calls() == 8);
    CHECK(sim.GetAccumulatedPOT() == Approx(4 * increment));
}

// Flux-averaged cross-section consistency: the event rate produced by the
// full ray machinery must imply the same flux-averaged cross section that the
// generator reports.  The flux is Σw / POT (w = per-ray flux weight), so each
// ray's interaction rate is w * N_col * sigma while it represents POT_per_ray
// of exposure.  Hence
//   events/POT = <w * N_col * sigma> / POT_per_ray
//   <sigma>_ray = (events/POT) * POT_per_ray / (N_col * <w>)
// with N_col computed here analytically (Avogadro, density, chord length) —
// independent of Material::NumberDensity.  Dividing out <w> here (not folding
// it into POT) is the correct normalization; a units mismatch anywhere in the
// chain (cm vs mm, cm^2 vs nb, per-cm^3 densities) shows up as orders of
// magnitude between <sigma>_ray and <sigma>_gen.
namespace {
// Water cube from DetSimFixture: rho = 1 g/cm^3, Box{10,10,10} = 10 cm chord.
// n_H = (2/18) N_A, n_O = (1/18) N_A  =>  n_total = (3/18) N_A per cm^3.
constexpr double kAvogadro = 6.02214076e23;
constexpr double kWaterAtomsPerCm3 = 3.0 / 18.0 * kAvogadro;
constexpr double kChordCm = 10.0;
constexpr double kColumnDensity = kWaterAtomsPerCm3 * kChordCm; // atoms / cm^2
} // namespace

TEST_CASE("Flux-averaged xsec from event rate matches the generator (monoenergetic)",
          "[DetectorSim][XSec]") {
    constexpr double sigma = 1e-38; // cm^2, constant
    NuGeom::DetectorSim sim(1.0);
    DetSimFixture::SetupSim(sim);
    auto gen = std::make_shared<NuGeom::MockGenerator>([](double, int, int) { return sigma; });
    sim.SetGenerator(gen);

    // Constant flux: Init lands max_prob exactly on the per-ray probability,
    // so every thrown ray is accepted and the rate is deterministic.
    sim.InitFromFlux(5);

    constexpr std::size_t nevents = 100;
    sim.GenerateEvents(nevents);

    const double events_per_pot = static_cast<double>(nevents) / sim.GetAccumulatedPOT();
    // Constant flux weight, so <w> = flux_weight; divide it out to recover sigma.
    const double sigma_ray =
        events_per_pot * DetSimFixture::ray_pot / (kColumnDensity * DetSimFixture::flux_weight);

    CHECK(sigma_ray == Approx(sigma).epsilon(0.02));
}

TEST_CASE("Flux-averaged xsec from event rate matches the generator (spectrum)",
          "[DetectorSim][XSec]") {
    // sigma(E) = sigma0 * E over a three-energy cycling flux: the generator's
    // flux average is sigma0 * <E>; the ray machinery has to reproduce it
    // through Init + accept/reject + POT accounting.
    constexpr double sigma0 = 1e-38; // cm^2 / GeV
    const std::vector<double> energies = {1.0, 2.0, 3.0};

    NuGeom::DetectorSim sim(1.0);
    DetSimFixture::SetupSim(sim);
    auto counter = std::make_shared<std::size_t>(0);
    sim.SetFluxCallback([counter, energies]() {
        NuGeom::FluxSample fs;
        fs.energy = energies[(*counter)++ % energies.size()];
        fs.pdg = 14;
        fs.ray = NuGeom::Ray({0, 0, -100}, {0, 0, 1}, DetSimFixture::ray_pot);
        fs.flux_weight = DetSimFixture::flux_weight;
        return fs;
    });
    auto gen = std::make_shared<NuGeom::MockGenerator>(
        [](double energy, int, int) { return sigma0 * energy; });
    sim.SetGenerator(gen);

    NuGeom::Random::Instance().Seed(20260613);
    sim.InitFromFlux(30); // covers the full energy cycle; max_prob = P(E_max)

    constexpr std::size_t nevents = 600;
    sim.GenerateEvents(nevents);

    const double events_per_pot = static_cast<double>(nevents) / sim.GetAccumulatedPOT();
    // Constant flux weight, so <w> = flux_weight; divide it out to recover sigma.
    const double sigma_ray =
        events_per_pot * DetSimFixture::ray_pot / (kColumnDensity * DetSimFixture::flux_weight);

    double sigma_gen = 0; // equal-weight flux average of the generator xsec
    for(const double e : energies) sigma_gen += gen->TotalXSec(e, 14, 1000080160);
    sigma_gen /= static_cast<double>(energies.size());

    // The headline check is order-of-magnitude agreement (catches unit bugs);
    // with this seed the stochastic acceptance keeps it within a few percent.
    CHECK(sigma_ray / sigma_gen > 0.1);
    CHECK(sigma_ray / sigma_gen < 10.0);
    CHECK(sigma_ray == Approx(sigma_gen).epsilon(0.10));
}

// Varying flux weights: with constant geometry/xsec but per-ray weights drawn
// from {0.05, 0.15} (mean 0.1), the absolute normalization must satisfy
//   events/POT = <w> * N_col * sigma / POT_per_ray.
// This is the case the constant-weight tests above CANNOT catch: if the flux
// weight is (incorrectly) folded into the POT charge, events/POT comes out a
// factor 1/<w> = 10x too high here, because the per-ray weights no longer
// cancel against the single max_prob envelope.
TEST_CASE("Flux-averaged xsec is correctly normalized for varying flux weights",
          "[DetectorSim][XSec]") {
    constexpr double sigma = 1e-38; // cm^2, constant
    const std::vector<double> weights = {0.05, 0.15};
    const double mean_w = (weights[0] + weights[1]) / 2.0; // 0.1

    NuGeom::DetectorSim sim(1.0);
    DetSimFixture::SetupSim(sim);
    auto counter = std::make_shared<std::size_t>(0);
    sim.SetFluxCallback([counter, weights]() {
        NuGeom::FluxSample fs;
        fs.energy = 1.0;
        fs.pdg = 14;
        fs.ray = NuGeom::Ray({0, 0, -100}, {0, 0, 1}, DetSimFixture::ray_pot);
        fs.flux_weight = weights[(*counter)++ % weights.size()];
        return fs;
    });
    auto gen = std::make_shared<NuGeom::MockGenerator>([](double, int, int) { return sigma; });
    sim.SetGenerator(gen);

    NuGeom::Random::Instance().Seed(20260613);
    sim.InitFromFlux(20); // max_prob = max(w) * N_col * sigma

    constexpr std::size_t nevents = 2000;
    sim.GenerateEvents(nevents);

    const double events_per_pot = static_cast<double>(nevents) / sim.GetAccumulatedPOT();
    // Recover sigma using the mean flux weight; the buggy (weight-in-POT)
    // normalization would land this 1/<w> = 10x high.
    const double sigma_ray = events_per_pot * DetSimFixture::ray_pot / (kColumnDensity * mean_w);

    CHECK(sigma_ray == Approx(sigma).epsilon(0.08));
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
        sim.InitFromFlux(10);
        sim.GenerateEvents(nevents);
        rate[mode] = static_cast<double>(nevents) / sim.GetAccumulatedPOT();
    }

    CHECK(rate[NuGeom::DetectorSim::SamplingMode::EnvelopeNoRetry] ==
          Approx(rate[NuGeom::DetectorSim::SamplingMode::TotalXSecRetry]));
}

// ...and the case that breaks that agreement in a real run.  AlternatingGenerator
// above is a *strict* unweighter: every emitted event has weight exactly 1.  Real
// generators use partial unweighting -- Achilles' default is a Percentile
// unweighter whose max_w is the 99th percentile of the weight distribution, so
// ~1% of trials are overweight, accepted with probability 1, and emitted with
// weight max(1, w/max_w) > 1.
//
// The consequence, derived in DetectorSim::ReportRunSummary:
//   * TotalXSecRetry retries until emit, so it produces exactly one event per
//     layer-1 accept: events/POT is the physical rate.
//   * EnvelopeNoRetry takes a single shot that emits with probability
//     E[min(1, w/max_w)] = sigma/(max_w * <weight>): its event count is LOW by
//     <weight>, and sum(weights)/POT is the physical rate instead.
// So the two modes' RAW EVENT COUNTS differ by exactly <weight>, while their
// weighted rates agree.  This is the mechanism behind the observed
// mode-to-mode event-count discrepancy.
namespace {

// Two-valued weight distribution with an explicit overweight tail.  A fraction
// `tail_frac` of trials draw weight `tail * max_w`, the rest draw `max_w/2`.
// Both TotalXSec and EnvelopeXSec are exact for this distribution, so any
// mode-to-mode difference comes purely from the overweight handling.
class PartialUnweightGenerator : public NuGeom::GeneratorInterface {
  public:
    static constexpr double kMaxW = 2e-38; // envelope (the "percentile")
    static constexpr double kTailFrac = 0.1;
    static constexpr double kTailRatio = 4.0; // overweight trials sit at 4 * max_w

    // E[w] = 0.9 * (max_w/2) + 0.1 * (4 max_w) = 0.85 * max_w
    static constexpr double kSigma = (1 - kTailFrac) * 0.5 * kMaxW + kTailFrac * kTailRatio * kMaxW;
    // E[min(1, w/max_w)] = 0.9 * 0.5 + 0.1 * 1 = 0.55
    static constexpr double kAcceptProb = (1 - kTailFrac) * 0.5 + kTailFrac * 1.0;
    // <weight | emitted> = E[w]/max_w / E[min(1, w/max_w)]
    static constexpr double kMeanWeight = (kSigma / kMaxW) / kAcceptProb;

    void InitRun(std::shared_ptr<HepMC3::GenRunInfo>) override {}
    double TotalXSec(double, int, int) override { return kSigma; }
    double EnvelopeXSec(double, int, int) override { return kMaxW; }

    bool GenerateEvent(HepMC3::GenEvent &evt) override {
        const double ratio = NuGeom::Random::Instance().Uniform(0.0, 1.0) < kTailFrac
                                 ? kTailRatio // overweight: accepted with probability 1
                                 : 0.5;
        if(ratio < NuGeom::Random::Instance().Uniform(0.0, 1.0)) return false;
        evt.weights() = {std::max(1.0, ratio)}; // the partial-unweight convention
        return true;
    }
};

} // namespace

TEST_CASE("A partial unweighter makes the two modes' raw event counts differ by <weight>",
          "[DetectorSim][SamplingMode]") {
    constexpr std::size_t nevents = 20000;
    std::map<NuGeom::DetectorSim::SamplingMode, NuGeom::DetectorSim::RunStats> stats;

    for(auto mode : {NuGeom::DetectorSim::SamplingMode::EnvelopeNoRetry,
                     NuGeom::DetectorSim::SamplingMode::TotalXSecRetry}) {
        NuGeom::Random::Instance().Seed(20260807);
        NuGeom::DetectorSim sim(1.0);
        DetSimFixture::SetupSim(sim);
        sim.SetSamplingMode(mode);
        sim.SetGenerator(std::make_shared<PartialUnweightGenerator>());
        sim.InitFromFlux(10);
        sim.GenerateEvents(nevents);
        stats[mode] = sim.GetRunStats();
    }

    const auto &env = stats[NuGeom::DetectorSim::SamplingMode::EnvelopeNoRetry];
    const auto &tot = stats[NuGeom::DetectorSim::SamplingMode::TotalXSecRetry];

    const double count_rate_env = static_cast<double>(env.emitted) / env.pot;
    const double count_rate_tot = static_cast<double>(tot.emitted) / tot.pot;

    // The counts disagree, by the mean emitted weight.
    CHECK(env.sum_weight / static_cast<double>(env.emitted) ==
          Approx(PartialUnweightGenerator::kMeanWeight).epsilon(0.03));
    CHECK(count_rate_tot / count_rate_env ==
          Approx(PartialUnweightGenerator::kMeanWeight).epsilon(0.03));

    // ...while the *right* estimator per mode does agree: sum(weights)/POT in
    // EnvelopeNoRetry, events/POT in TotalXSecRetry.
    CHECK(env.sum_weight / env.pot == Approx(count_rate_tot).epsilon(0.03));
}

// The custom-ray envelope: Init() must bound every flux ray's acceptance weight
// without ever having seen the flux, and the run must land on the same rate the
// flux-calibrated envelope gives.  A deliberately hostile flux (weights spanning
// three decades, presented in increasing order so the early rays badly
// under-represent the tail) is exactly the case a fixed warm-up envelope clips.
TEST_CASE("Custom-ray Init gives an unbiased rate on a heavy-tailed flux", "[DetectorSim][Init]") {
    constexpr double sigma = 1e-38;
    const std::vector<double> weights = {0.001, 0.01, 0.1, 1.0};
    const double mean_w =
        std::accumulate(weights.begin(), weights.end(), 0.0) / static_cast<double>(weights.size());

    NuGeom::Random::Instance().Seed(20260807);
    NuGeom::DetectorSim sim(1.0);
    DetSimFixture::SetupSim(sim);
    auto counter = std::make_shared<std::size_t>(0);
    sim.SetFluxCallback([counter, weights]() {
        NuGeom::FluxSample fs;
        fs.energy = 1.0;
        fs.pdg = 14;
        fs.ray = NuGeom::Ray({0, 0, -100}, {0, 0, 1}, DetSimFixture::ray_pot);
        fs.flux_weight = weights[(*counter)++ % weights.size()];
        return fs;
    });
    sim.SetGenerator(
        std::make_shared<NuGeom::MockGenerator>([](double, int, int) { return sigma; }));

    // No flux involved: chords through the world box over an energy grid.
    sim.SetEnergyRange(0.5, 5.0);
    sim.SetFluxSpecies({14});
    sim.Init(64);

    constexpr std::size_t nevents = 4000;
    sim.GenerateEvents(nevents);

    const double events_per_pot = static_cast<double>(nevents) / sim.GetAccumulatedPOT();
    const double sigma_ray = events_per_pot * DetSimFixture::ray_pot / (kColumnDensity * mean_w);
    CHECK(sigma_ray == Approx(sigma).epsilon(0.05));

    // The probe rays are chords of the world box, so they bound the flux ray's
    // column density: the envelope only ever grows to track the flux weight.
    CHECK(sim.GetRunStats().envelope_growths == weights.size());
}

// Init() must fail loudly rather than leave a zero envelope for ProduceEvent to
// divide by (0/0 -> NaN, then an infinite loop).
TEST_CASE("Custom-ray Init throws when the cross section vanishes", "[DetectorSim][Init]") {
    NuGeom::DetectorSim sim(1.0);
    DetSimFixture::SetupSim(sim);
    sim.SetMaxProb(0.0);
    sim.SetGenerator(std::make_shared<NuGeom::MockGenerator>([](double, int, int) { return 0.0; }));
    sim.SetEnergyRange(0.5, 5.0);
    CHECK_THROWS_WITH(sim.Init(8), Catch::Contains("cross section is zero"));
}

TEST_CASE("Custom-ray Init requires an energy range", "[DetectorSim][Init]") {
    NuGeom::DetectorSim sim(1.0);
    DetSimFixture::SetupSim(sim);
    sim.SetGenerator(
        std::make_shared<NuGeom::MockGenerator>([](double, int, int) { return 1e-38; }));
    CHECK_THROWS_WITH(sim.Init(8), Catch::Contains("energy range"));
}
