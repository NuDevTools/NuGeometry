#pragma once

#include "geom/FluxSource.hh"
#include "geom/Vector3D.hh"

#include <fstream>
#include <limits>
#include <memory>
#include <string>

namespace HepMC3 {
class ReaderAscii;
}

namespace NuGeom {

// Streams FluxSamples one at a time from a flux file, keeping the file open
// instead of loading every ray into memory.
//
// TryNext() yields the next ray and returns false at end-of-file.  Next() is
// the convenience used as a DetectorSim flux callback: at end-of-file it
// either rewinds to the first ray (loop = true) or throws std::runtime_error
// (loop = false).
//
// Construction performs one light pass over the file to gather a summary
// (ray count, energy range, and for NuHepMC files the POT and beam->detector
// offset) without retaining the rays; the count keeps the per-ray POT
// normalization of NuHepMC fluxes identical to a full load.
class FluxStreamer {
  public:
    virtual ~FluxStreamer() = default;

    // Next ray in file order; false once the end of the file is reached
    // (subsequent calls keep returning false until Next() rewinds).
    virtual bool TryNext(FluxSample &fs) = 0;

    FluxSample Next();

    bool Looping() const { return m_loop; }
    std::size_t Count() const { return m_count; }
    double EMin() const { return m_emin; }
    double EMax() const { return m_emax; }
    // Completed passes over the file (only grows when looping).
    std::size_t Loops() const { return m_loops; }

  protected:
    FluxStreamer(std::string path, bool loop) : m_path{std::move(path)}, m_loop{loop} {}
    // Rewind to the first ray.
    virtual void Rewind() = 0;

    std::string m_path;
    bool m_loop;
    std::size_t m_count = 0;
    double m_emin = std::numeric_limits<double>::infinity();
    double m_emax = -std::numeric_limits<double>::infinity();
    std::size_t m_loops = 0;
};

// Comma-separated flux file with a header row and one neutrino per line:
//     pid,wgt,E,px,py,pz,x,y,z
// (pid = PDG code, wgt = per-sample flux weight, E in GeV, direction
// normalized on load, origin in the geometry's length units).  Malformed or
// blank rows are skipped.  Each ray carries Ray::POT() = 1.
class CSVFluxStreamer : public FluxStreamer {
  public:
    CSVFluxStreamer(const std::string &path, bool loop);
    bool TryNext(FluxSample &fs) override;

  protected:
    void Rewind() override;

  private:
    std::ifstream m_in;
};

// NuHepMC (HepMC3 ASCII) flux file, e.g. from dune_dk2nu/dk2nu_to_nuhepmc.py.
// Same per-ray semantics as LoadHepMCFlux (see HepMCFluxReader.hh): rays are
// stored in the beam frame and the recorded beam->detector translation (or
// `offset_override`, in cm) is applied to each origin; Ray::POT() =
// total POT / total ray count so a full pass reproduces the file exposure.
class HepMCFluxStreamer : public FluxStreamer {
  public:
    HepMCFluxStreamer(const std::string &path, bool loop,
                      const Vector3D *offset_override = nullptr);
    ~HepMCFluxStreamer() override;
    bool TryNext(FluxSample &fs) override;

    double TotalPOT() const { return m_pot; }
    const Vector3D &Offset() const { return m_offset; }

  protected:
    void Rewind() override;

  private:
    // First pass: metadata (POT, offset) + ray count + energy range.
    void Summarize(const Vector3D *offset_override);

    double m_pot = 0.0;
    double m_pot_per_ray = 1.0;
    // Transverse flux-window area (cm^2) recorded by the converter; folded into
    // each sample's window_area so the per-cm^2 weight integrates over the
    // detector face.  1.0 (no-op) when the file predates the attribute.
    double m_window_area = 1.0;
    Vector3D m_offset;
    std::unique_ptr<HepMC3::ReaderAscii> m_reader;
};

// Open a flux file by extension: NuHepMC for .hepmc/.hepmc3, CSV otherwise.
// For NuHepMC files `offset_override` (if non-null, cm) replaces the recorded
// beam->detector translation.
std::unique_ptr<FluxStreamer> OpenFluxStreamer(const std::string &path, bool loop,
                                               const Vector3D *offset_override = nullptr);

} // namespace NuGeom
