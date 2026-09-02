#include "catch2/catch.hpp"

#include "geom/FluxStreamer.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "HepMC3/Attribute.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/WriterAscii.h"
#pragma GCC diagnostic pop

#include <filesystem>
#include <fstream>
#include <string>

namespace {

// Write a small CSV flux file (header + rows) and return its path.
std::string WriteCSV(const std::string &name, const std::string &body) {
    const auto path = (std::filesystem::temp_directory_path() / name).string();
    std::ofstream out(path);
    out << "pid,wgt,E,px,py,pz,x,y,z\n" << body;
    return path;
}

const std::string kThreeRays = "14,0.5,1.0,0,0,1,0,0,-100\n"
                               "-14,0.25,2.5,0,0,1,1,2,-100\n"
                               "12,1.0,0.5,0,1,0,3,4,-100\n";

// Write a minimal NuHepMC-style flux file with `n` events: parent pion in,
// neutrino (status 1) out of a decay vertex.  POT and the beam->detector
// translation are recorded as run-info attributes, as the dk2nu converter does.
// `survey` adds the NuGeom.Flux.NRays / EnergyRange_GeV attributes the current
// converter writes, which let the streamer open the file without scanning it.
std::string WriteHepMCFlux(const std::string &name, int n, bool survey = false) {
    const auto path = (std::filesystem::temp_directory_path() / name).string();
    auto run = std::make_shared<HepMC3::GenRunInfo>();
    run->add_attribute("NuHepMC.Exposure.POT", std::make_shared<HepMC3::StringAttribute>("100.0"));
    run->add_attribute("NuGeom.BeamToDetector.Translation_cm",
                       std::make_shared<HepMC3::StringAttribute>("10,20,30"));
    run->add_attribute("NuGeom.LengthUnit", std::make_shared<HepMC3::StringAttribute>("cm"));
    if(survey) {
        run->add_attribute("NuGeom.Flux.NRays",
                           std::make_shared<HepMC3::StringAttribute>(std::to_string(n)));
        run->add_attribute("NuGeom.Flux.EnergyRange_GeV",
                           std::make_shared<HepMC3::StringAttribute>("1.0," + std::to_string(n)));
    }

    HepMC3::WriterAscii writer(path, run);
    for(int i = 0; i < n; ++i) {
        HepMC3::GenEvent evt(run, HepMC3::Units::GEV, HepMC3::Units::CM);
        evt.set_event_number(i);
        evt.weights() = {0.5 + i};

        const double e = 1.0 + i;
        auto parent =
            std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, e, e), 211, 4); // pi+ in
        auto nu = std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, e, e), 14,
                                                        1); // numu out
        auto vtx = std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector(i, 0, 0, 0));
        vtx->add_particle_in(parent);
        vtx->add_particle_out(nu);
        evt.add_vertex(vtx);
        writer.write_event(evt);
    }
    writer.close();
    return path;
}

} // namespace

TEST_CASE("CSVFluxStreamer streams rays in order", "[FluxStreamer]") {
    const auto path = WriteCSV("nugeom_streamer_order.csv", kThreeRays);
    NuGeom::CSVFluxStreamer streamer(path, /*loop=*/false);

    CHECK(streamer.Count() == 3);
    CHECK(streamer.EMin() == Approx(0.5));
    CHECK(streamer.EMax() == Approx(2.5));

    auto fs = streamer.Next();
    CHECK(fs.pdg == 14);
    CHECK(fs.energy == Approx(1.0));
    CHECK(fs.flux_weight == Approx(0.5));
    CHECK(fs.ray.POT() == Approx(1.0));

    fs = streamer.Next();
    CHECK(fs.pdg == -14);
    CHECK(fs.energy == Approx(2.5));

    fs = streamer.Next();
    CHECK(fs.pdg == 12);
    CHECK(fs.energy == Approx(0.5));
}

TEST_CASE("CSVFluxStreamer errors at end-of-file when not looping", "[FluxStreamer]") {
    const auto path = WriteCSV("nugeom_streamer_eof.csv", kThreeRays);
    NuGeom::CSVFluxStreamer streamer(path, /*loop=*/false);
    for(int i = 0; i < 3; ++i) streamer.Next();
    CHECK_THROWS_WITH(streamer.Next(), Catch::Contains("end of flux file"));
}

TEST_CASE("CSVFluxStreamer rewinds at end-of-file when looping", "[FluxStreamer]") {
    const auto path = WriteCSV("nugeom_streamer_loop.csv", kThreeRays);
    NuGeom::CSVFluxStreamer streamer(path, /*loop=*/true);

    // Seven draws cover two full passes plus one ray; the wrap restarts at
    // the first ray each time.
    std::vector<int> pdgs;
    for(int i = 0; i < 7; ++i) pdgs.push_back(streamer.Next().pdg);
    CHECK(pdgs == std::vector<int>{14, -14, 12, 14, -14, 12, 14});
    CHECK(streamer.Loops() == 2);
}

TEST_CASE("CSVFluxStreamer skips malformed rows and empty lines", "[FluxStreamer]") {
    const auto path = WriteCSV("nugeom_streamer_bad.csv", "14,0.5,1.0,0,0,1,0,0,-100\n"
                                                          "\n"
                                                          "not,a,valid,row\n"
                                                          "12,1.0,0.5,0,1,0,3,4,-100\n");
    NuGeom::CSVFluxStreamer streamer(path, /*loop=*/false);
    CHECK(streamer.Count() == 2);
    CHECK(streamer.Next().pdg == 14);
    CHECK(streamer.Next().pdg == 12);
}

TEST_CASE("CSVFluxStreamer rejects files without rays", "[FluxStreamer]") {
    const auto path = WriteCSV("nugeom_streamer_empty.csv", "");
    CHECK_THROWS_WITH(NuGeom::CSVFluxStreamer(path, false), Catch::Contains("no rays"));
}

TEST_CASE("HepMCFluxStreamer streams NuHepMC rays with POT split and offset", "[FluxStreamer]") {
    const auto path = WriteHepMCFlux("nugeom_streamer_flux.hepmc", 4);
    NuGeom::HepMCFluxStreamer streamer(path, /*loop=*/false);

    CHECK(streamer.Count() == 4);
    CHECK(streamer.TotalPOT() == Approx(100.0));
    CHECK(streamer.Offset().X() == Approx(10.0));
    // No recorded energy range and the streamer never reads the rays up front,
    // so it must report the range as unknown rather than invent one.
    CHECK_FALSE(streamer.HasEnergyRange());

    for(int i = 0; i < 4; ++i) {
        auto fs = streamer.Next();
        CHECK(fs.pdg == 14);
        CHECK(fs.energy == Approx(1.0 + i));
        CHECK(fs.flux_weight == Approx(0.5 + i));
        // Per-ray POT = file POT / ray count, so a full pass reproduces the
        // recorded exposure.
        CHECK(fs.ray.POT() == Approx(25.0));
        // Vertex position plus the recorded beam->detector translation.
        CHECK(fs.ray.Origin().X() == Approx(10.0 + i));
        CHECK(fs.ray.Origin().Y() == Approx(20.0));
        CHECK(fs.ray.Origin().Z() == Approx(30.0));
    }
    CHECK_THROWS_WITH(streamer.Next(), Catch::Contains("end of flux file"));
}

TEST_CASE("HepMCFluxStreamer loops over the file when requested", "[FluxStreamer]") {
    const auto path = WriteHepMCFlux("nugeom_streamer_flux_loop.hepmc", 2);
    NuGeom::HepMCFluxStreamer streamer(path, /*loop=*/true);

    std::vector<double> energies;
    for(int i = 0; i < 5; ++i) energies.push_back(streamer.Next().energy);
    CHECK(energies == std::vector<double>{1.0, 2.0, 1.0, 2.0, 1.0});
    CHECK(streamer.Loops() == 2);
}

TEST_CASE("OpenFluxStreamer dispatches on extension", "[FluxStreamer]") {
    const auto csv = WriteCSV("nugeom_streamer_dispatch.csv", kThreeRays);
    const auto hepmc = WriteHepMCFlux("nugeom_streamer_dispatch.hepmc", 2);

    auto s1 = NuGeom::OpenFluxStreamer(csv, false);
    CHECK(dynamic_cast<NuGeom::CSVFluxStreamer *>(s1.get()) != nullptr);
    auto s2 = NuGeom::OpenFluxStreamer(hepmc, false);
    CHECK(dynamic_cast<NuGeom::HepMCFluxStreamer *>(s2.get()) != nullptr);
}

// The converter records the ray count and energy range so the streamer can open
// a multi-GB flux file without touching it: the count still has to normalize
// POT per ray exactly as the byte-scan fallback does.
TEST_CASE("HepMCFluxStreamer takes the ray count and energy range from the file",
          "[FluxStreamer]") {
    const auto path = WriteHepMCFlux("nugeom_streamer_survey.hepmc", 4, /*survey=*/true);
    NuGeom::HepMCFluxStreamer streamer(path, /*loop=*/false);

    CHECK(streamer.Count() == 4);
    REQUIRE(streamer.HasEnergyRange());
    CHECK(streamer.EMin() == Approx(1.0));
    CHECK(streamer.EMax() == Approx(4.0));
    CHECK(streamer.Next().ray.POT() == Approx(25.0));
}

// ...and the byte-scan fallback must agree with it ray for ray on a file the
// converter did not survey.
TEST_CASE("HepMCFluxStreamer counts rays without the recorded count", "[FluxStreamer]") {
    const auto surveyed = WriteHepMCFlux("nugeom_streamer_count_a.hepmc", 7, /*survey=*/true);
    const auto bare = WriteHepMCFlux("nugeom_streamer_count_b.hepmc", 7, /*survey=*/false);

    NuGeom::HepMCFluxStreamer a(surveyed, /*loop=*/false);
    NuGeom::HepMCFluxStreamer b(bare, /*loop=*/false);
    CHECK(a.Count() == b.Count());
    CHECK(b.Count() == 7);
    CHECK(a.Next().ray.POT() == Approx(b.Next().ray.POT()));
}
