// Cross-check Achilles' sigma_tot SPLINE against the integral of its own
// generated event WEIGHTS, at fixed neutrino energies.
//
// This settles which side of the mode discrepancy documented in
// mode_discrepancy_findings.md is at fault.  NuGeometry's two sampling modes
// source sigma_tot from different places:
//
//   TotalXSecRetry  -> EventGen::TotalXSec(E)      (the spline)
//   EnvelopeNoRetry -> LastEvent().Weight() * VertexEnvelope, averaged over
//                      trials with rejected trials counted as zero
//
// and their event spectra disagree above ~4 GeV.  Both are measured here on the
// same configuration, at fixed energy, with no geometry or flux involved.
//
//   xsec_scan <runcard.yml> [trials_per_energy]
#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventGen.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/ParticleInfo.hh"

#include "yaml-cpp/yaml.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {
    if(argc < 2) {
        std::cerr << "usage: xsec_scan <runcard.yml> [trials_per_energy]\n";
        return 1;
    }
    const std::string card = argv[1];
    const std::size_t ntrials = argc > 2 ? std::strtoul(argv[2], nullptr, 10) : 20000;

    YAML::Node root = YAML::LoadFile(card);
    // The driver run cards nest the Achilles config under "Achilles:".
    YAML::Node ach = root["Achilles"] ? root["Achilles"] : root;

    achilles::EventGen gen(ach, {}, card);
    gen.Initialize();

    const achilles::PID nu{14};
    const achilles::PID target{1000180400}; // Ar40
    const double M = gen.VertexEnvelope(nu, target);
    constexpr double cm2_to_pb = 1.0 / achilles::Constant::TO_CM2;

    std::cout << "\nVertexEnvelope (max_w) = " << M * cm2_to_pb << " pb\n";
    std::cout << "trials per energy      = " << ntrials << "\n\n";
    std::cout
        << "   E [GeV]   sigma_spline [pb]   sigma_weights [pb]   +/-      ratio   emit_frac\n";

    for(double E : {1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0}) {
        const double e_mev = E * 1000.0;
        const double spline = gen.TotalXSec(e_mev, nu, target) * cm2_to_pb;

        // Same estimator the adapter uses in EnvelopeNoRetry: every trial is one
        // sample, a rejected trial contributes zero, an emitted one contributes
        // Weight() * M.  The mean over trials estimates sigma.
        achilles::FourVector p4{e_mev, 0.0, 0.0, e_mev};
        double sum = 0.0, sum2 = 0.0;
        std::size_t emitted = 0;
        std::vector<double> ws;
        ws.reserve(ntrials);
        for(std::size_t i = 0; i < ntrials; ++i) {
            gen.InjectRay(p4, nu, target);
            double w = 0.0;
            if(gen.GenerateSingleEvent()) {
                w = gen.LastEvent().Weight() * M * cm2_to_pb;
                ++emitted;
            }
            sum += w;
            sum2 += w * w;
            ws.push_back(w);
        }
        // Tail heaviness: what fraction of the total does the top 1% of trials
        // carry?  A heavy tail means a small sample systematically under-counts.
        std::sort(ws.begin(), ws.end());
        const std::size_t top = ws.size() / 100;
        double tail = 0.0;
        for(std::size_t i = ws.size() - top; i < ws.size(); ++i) tail += ws[i];
        const double tail_frac = sum > 0 ? tail / sum : 0.0;
        const double wmax = ws.empty() ? 0.0 : ws.back();
        const double n = static_cast<double>(ntrials);
        const double mean = sum / n;
        const double err = std::sqrt(std::max(0.0, sum2 / n - mean * mean) / n);
        std::cout.precision(5);
        std::cout << std::fixed;
        std::cout << "   " << E << "        " << spline << "             " << mean << "        "
                  << err << "   " << (spline > 0 ? mean / spline : 0.0) << "   "
                  << static_cast<double>(emitted) / n << "   top1%=" << tail_frac
                  << "   wmax/mean=" << (mean > 0 ? wmax / mean : 0.0) << "\n";
    }
    gen.Finalize();
    return 0;
}
