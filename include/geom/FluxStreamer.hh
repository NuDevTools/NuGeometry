#pragma once

#include "geom/FluxSource.hh"
#include "geom/Vector3D.hh"

#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

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
// Construction gathers a summary (ray count, energy range, and for NuHepMC
// files the POT and beam->detector offset) without retaining the rays; the
// count keeps the per-ray POT normalization of NuHepMC fluxes identical to a
// full load.  How much of the file that costs is format-specific -- see
// HepMCFluxStreamer, which never parses the file up front.
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
    // Whether EMin()/EMax() are known.  A streamer that never reads the whole
    // file up front can only report them if the file records them.
    bool HasEnergyRange() const { return m_emin <= m_emax; }
    // Completed passes over the file (only grows when looping).
    std::size_t Loops() const { return m_loops; }

    // Mean / total flux weight over the file, when known without reading it.
    // NuHepMC files record these (NuGeom.Flux.MeanWeight / TotalWeight); they
    // let importance sampling be armed -- and its payoff reported -- before a
    // single ray has been parsed.  HasWeightSummary() is false otherwise.
    bool HasWeightSummary() const { return m_mean_weight > 0.0; }
    double MeanFluxWeight() const { return m_mean_weight; }
    double TotalFluxWeight() const { return m_total_weight; }

  protected:
    FluxStreamer(std::string path, bool loop) : m_path{std::move(path)}, m_loop{loop} {}
    // Rewind to the first ray.
    virtual void Rewind() = 0;
    // CachedFluxStreamer wraps another streamer and must be able to rewind it.
    friend class CachedFluxStreamer;

    std::string m_path;
    bool m_loop;
    std::size_t m_count = 0;
    double m_emin = std::numeric_limits<double>::infinity();
    double m_emax = -std::numeric_limits<double>::infinity();
    std::size_t m_loops = 0;
    double m_mean_weight = 0.0;
    double m_total_weight = 0.0;
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
//
// Nothing is parsed up front.  Construction reads exactly one event -- enough
// for the run-level metadata -- and then needs only the ray count, to turn the
// file's total POT into a per-ray exposure.  That count comes from the
// NuGeom.Flux.NRays run attribute when the converter recorded it, and
// otherwise from a raw byte scan for event-record lines, which is ~20x cheaper
// than building a GenEvent per ray (17 s -> under 1 s on a 585 MB DUNE ND
// flux).  The energy range is likewise taken from NuGeom.Flux.EnergyRange_GeV
// if present and reported as unknown if not (HasEnergyRange()); deriving it
// would require reading every ray, which is exactly what this class avoids.
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
    // Read the run-level metadata (POT, offset, window area, and the optional
    // ray-count / energy-range hints) from the file's first event.
    void ReadMetadata(const Vector3D *offset_override);

    double m_pot = 0.0;
    double m_pot_per_ray = 1.0;
    // Transverse flux-window area (cm^2) recorded by the converter; folded into
    // each sample's window_area so the per-cm^2 weight integrates over the
    // detector face.  1.0 (no-op) when the file predates the attribute.
    double m_window_area = 1.0;
    Vector3D m_offset;
    std::unique_ptr<HepMC3::ReaderAscii> m_reader;
};

// Decorator that keeps parsed rays in memory, and can draw them in proportion
// to their flux weight.  Both behaviours are opt-in and independent of the
// underlying format.
//
// CACHING (Options::cache).  The first pass pulls rays from the wrapped
// streamer and stores a compact record for each; every later pass is served
// from RAM.  Parsing a HepMC3 ASCII event costs ~24 us, which dominates the
// run once layer-1 rejection is cheap, and a long run re-reads the file many
// times (a 1e17-POT envelope run threw 55M rays against a 961k-ray file, ~57
// passes).  `max_bytes` caps the footprint: a file that does not fit stops
// being cached mid-pass and the streamer transparently degrades to pure
// streaming for the rest of the run -- so enabling this can never make a large
// file fail, only stop helping.
//
// IMPORTANCE SAMPLING (Options::importance).  Instead of walking the file in
// order, draw ray i with probability w_i / sum(w) and hand it the MEAN weight
// sum(w)/N in place of its own.  That is the textbook importance-sampling
// reweighting and leaves the rate estimator unbiased:
//
//     E_p[(W/N) * A * C] = sum_i (w_i/W) (W/N) A C_i = (1/N) sum_i w_i A C_i
//
// which is exactly the uniform-draw expectation, while Ray::POT() stays
// filePOT/N so a full N draws still reproduce the file exposure.  The payoff is
// that the acceptance envelope's flux factor becomes mean(w) rather than
// max(w) -- a factor ~459 on the DUNE ND flux -- so the same POT needs ~459x
// fewer rays.  The cost is that rays are no longer visited in file order and
// some are never drawn at all, so leave it off when every ray must be used.
// It needs random access and therefore requires caching.
class CachedFluxStreamer : public FluxStreamer {
  public:
    struct Options {
        bool cache{true};         ///< keep parsed rays in memory
        bool importance{false};   ///< draw rays proportional to flux weight
        std::size_t max_bytes{0}; ///< cache cap; 0 = unlimited
        /// Fill the cache at construction instead of on the first pass.
        /// Importance sampling REQUIRES this: the acceptance envelope only ever
        /// grows, so if the run begins on real weights it locks in a max(w)
        /// -sized envelope and the later mean-weight rays are then accepted
        /// ~max/mean times too rarely -- the speedup is lost entirely.
        /// Defaults to on whenever `importance` is set.
        bool preload{false};
    };

    CachedFluxStreamer(std::unique_ptr<FluxStreamer> inner, Options opts);
    bool TryNext(FluxSample &fs) override;

    bool Caching() const { return m_caching; }
    bool Importance() const { return m_importance; }
    std::size_t CachedRays() const { return m_rays.size(); }
    /// Mean flux weight over the cached rays (the weight handed to every ray
    /// once importance sampling is live).  0 before the first pass completes.
    double MeanWeight() const { return m_mean_w; }
    /// Bytes held by the cache.
    std::size_t CacheBytes() const { return m_rays.size() * sizeof(Cached); }

  protected:
    void Rewind() override;

  private:
    /// Compact per-ray record.  window_area and POT-per-ray are file-level
    /// constants, so they are not stored per ray.
    struct Cached {
        double energy, ox, oy, oz, dx, dy, dz, weight;
        int pdg;
    };

    void Emit(const Cached &c, FluxSample &fs, double weight_override) const;
    /// Build the cumulative weight table used to draw rays proportional to w.
    void BuildSampler();
    /// Drain the wrapped streamer into the cache up front.
    void Preload();

    std::unique_ptr<FluxStreamer> m_inner;
    std::vector<Cached> m_rays;
    std::vector<double> m_cumulative; ///< running sum of weights, for sampling
    bool m_caching{true};             ///< still filling the cache (false once capped)
    bool m_importance{false};
    bool m_filled{false}; ///< first pass complete, cache authoritative
    std::size_t m_max_bytes{0};
    std::size_t m_pos{0};   ///< read cursor when serving sequentially
    std::size_t m_draws{0}; ///< draws made in importance mode
    double m_mean_w{0.0};
    double m_total_w{0.0};
    double m_window_area{1.0};
    double m_pot_per_ray{1.0};
};

// Open a flux file by extension: NuHepMC for .hepmc/.hepmc3, CSV otherwise.
// For NuHepMC files `offset_override` (if non-null, cm) replaces the recorded
// beam->detector translation.
std::unique_ptr<FluxStreamer> OpenFluxStreamer(const std::string &path, bool loop,
                                               const Vector3D *offset_override = nullptr);

// As above, wrapped in a CachedFluxStreamer when `opts` asks for caching or
// importance sampling; the bare streamer otherwise.
std::unique_ptr<FluxStreamer> OpenFluxStreamer(const std::string &path, bool loop,
                                               const Vector3D *offset_override,
                                               const CachedFluxStreamer::Options &opts);

} // namespace NuGeom
