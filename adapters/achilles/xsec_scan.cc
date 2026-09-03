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

    const std::vector<double> grid{1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};

    // Interleaved pass: energies are visited in round-robin instead of in
    // blocks, so Achilles' single adaptive max_w sees the same mixture it does
    // in a real flux run.  If the per-energy emit probability drops relative to
    // the blocked scan below, the shared scalar max_w is the culprit.
    if(std::getenv("XSEC_SCAN_INTERLEAVE")) {
        std::vector<std::size_t> acc(grid.size(), 0), emi(grid.size(), 0);
        std::vector<double> sw(grid.size(), 0.0);
        for(std::size_t i = 0; i < ntrials * grid.size(); ++i) {
            const std::size_t k = i % grid.size();
            const double e_mev = grid[k] * 1000.0;
            achilles::FourVector p4{e_mev, 0.0, 0.0, e_mev};
            gen.InjectRay(p4, nu, target);
            ++acc[k];
            if(gen.GenerateSingleEvent()) {
                ++emi[k];
                sw[k] += gen.LastEvent().Weight() * M * cm2_to_pb;
            }
        }
        std::cout << "INTERLEAVED (shared adaptive max_w sees all energies)\n";
        std::cout << "   E [GeV]   sigma_spline   sigma_weights   ratio   emit_frac\n";
        for(std::size_t k = 0; k < grid.size(); ++k) {
            const double spl = gen.TotalXSec(grid[k] * 1000.0, nu, target) * cm2_to_pb;
            const double mean = sw[k] / static_cast<double>(acc[k]);
            std::cout.precision(5);
            std::cout << std::fixed << "   " << grid[k] << "        " << spl << "        " << mean
                      << "     " << (spl > 0 ? mean / spl : 0.0) << "   "
                      << static_cast<double>(emi[k]) / static_cast<double>(acc[k]) << "\n";
        }
        gen.Finalize();
        return 0;
    }

    // Direction scan: fixed energy, injected off-axis.  Both other modes always
    // inject along +z, whereas a flux run injects each ray's actual lab
    // direction through the GeometryBeam rotation.  Flux angles here are tiny
    // (median 0.18 deg, max 3.7 deg) but the scan goes well past that to expose
    // any pathology in the rotation path.
    if(const char *ang_e = std::getenv("XSEC_SCAN_ANGLE")) {
        const double E = std::strtod(ang_e, nullptr);
        const double e_mev = E * 1000.0;
        const double spline = gen.TotalXSec(e_mev, nu, target) * cm2_to_pb;
        std::cout << "ANGLE SCAN at E = " << E << " GeV (sigma_spline = " << spline << " pb)\n";
        std::cout << "   theta [deg]   sigma_weights   ratio   emit_frac\n";
        // XSEC_SCAN_ANGLES overrides the list, so the same angle can be
        // repeated (does the collapse follow the ANGLE?) or the order reversed
        // (does it follow POSITION IN THE SEQUENCE, i.e. generator state?).
        std::vector<double> angles{0.0, 0.2, 1.0, 2.0, 4.0, 10.0, 30.0, 90.0, 180.0};
        if(const char *lst = std::getenv("XSEC_SCAN_ANGLES")) {
            angles.clear();
            std::string t(lst);
            std::size_t pos = 0;
            while(pos < t.size()) {
                std::size_t c = t.find(',', pos);
                if(c == std::string::npos) c = t.size();
                angles.push_back(std::strtod(t.substr(pos, c - pos).c_str(), nullptr));
                pos = c + 1;
            }
        }
        for(double deg : angles) {
            const double th = deg * M_PI / 180.0;
            achilles::FourVector p4{e_mev, e_mev * std::sin(th), 0.0, e_mev * std::cos(th)};
            double sum = 0.0;
            std::size_t emitted = 0;
            for(std::size_t i = 0; i < ntrials; ++i) {
                gen.InjectRay(p4, nu, target);
                if(gen.GenerateSingleEvent()) {
                    sum += gen.LastEvent().Weight() * M * cm2_to_pb;
                    ++emitted;
                }
            }
            const double n = static_cast<double>(ntrials);
            std::cout.precision(5);
            std::cout << std::fixed << "   " << deg << "        " << sum / n << "     "
                      << (spline > 0 ? sum / n / spline : 0.0) << "   "
                      << static_cast<double>(emitted) / n << "\n";
        }
        gen.Finalize();
        return 0;
    }

    for(double E : grid) {
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
