#!/usr/bin/env python3
"""Convert dk2nu ROOT flux files to a NuHepMC (HepMC3 ASCII) flux file.

No ROOT dependency: branches are read with uproot and events are written with
pyhepmc (which bundles HepMC3). The result is consumed by NuGeometry, whose
FluxSource/GeneratorInterface already speak NuHepMC.

Each dk2nu entry stores a vector<NuRay> with three precomputed rays:
    [0] random decay direction   [1] near detector   [2] far detector
We emit one HepMC3 event per entry for a chosen ray location (default ND = 1).
The ray already carries the neutrino energy, direction and importance weight
from g4lbne, so calcEnuWgt / decay_through_point (the pydk2nu/ROOT path) is not
needed.

Weight units (IMPORTANT)
------------------------
dk2nu's nuray.wgt is calcEnuWgt's sangdet*emrat^2: the number of neutrinos
through a sphere of radius kRDET = 100 cm, i.e. through a 1 m-radius disk of
area pi*(100 cm)^2 = pi m^2. The main NuGeometry code (and dk2nu_to_csv.py via
pydk2nu) work in per-cm^2 flux, so we convert:
    w_cm2 = nimpwt * nuray.wgt * 1e-4 / pi
where 1e-4 is m^2 -> cm^2 and 1/pi turns "per pi-m^2 disk" into "per m^2".
Skipping this scaling (writing nimpwt*nuray.wgt) overstates the flux by
pi*1e4 ~= 31416; see compare_flux.py for the end-to-end cross-check.

Coordinate convention (IMPORTANT)
---------------------------------
dk2nu decay vertices live in the BEAM frame: target at the origin, +z along
the beam, units cm. They are far upstream of the detector in the GDML's
hall-local frame. We deliberately write rays in their NATIVE beam frame and
record the frame + (optional) beam->detector translation as run-level metadata,
so NuGeometry can apply the transform once at load time rather than having a
hidden offset baked into the file. See the NuGeom.* run attributes below.

Per-event central weight is the dk2nu flux weight, converted to per-cm^2
    w = decay.nimpwt * nuray.wgt[location] * 1e-4 / pi
(the same 1/pi that DUNE's plottingExamples.py applies; it is the disk-area
normalization of calcEnuWgt, not merely a per-degree plotting factor).

Flux window (IMPORTANT)
-----------------------
nuray[location] is precomputed for a single point (ND = on-axis (0,0,57400) cm),
so emitting it verbatim makes every ray converge to that point. The detector
sim would then trace one chord and recover only the per-cm^2 flux rate -- low by
~the detector face area (~2.4e5 cm^2 for ND-LAr). To fix the normalization we
spread each ray's origin by a uniform random (dx,dy) over a transverse flux
window (--window-width/--window-height, beam frame; default = the ND-LAr face)
and record NuGeom.FluxWindow.Area_cm2. NuGeometry multiplies the per-cm^2 weight
by that area, so summing flux_weight * area * P_interaction over the spread rays
Monte-Carlo integrates the rate over the face -- mirroring how GENIE's GNuMIFlux
throws neutrinos at random positions across its flux window. The nuray weight is
reused at the shifted point (a uniform-flux-over-window approximation, excellent
at the ~574 m ND baseline); the exact per-position weight would need pydk2nu's
calcEnuWgt (see compare_flux.py) and is left as a future refinement.
"""

import argparse
import math
import sys

import numpy as np
import uproot
import pyhepmc
from pyhepmc import GenEvent, GenParticle, GenVertex, GenRunInfo, FourVector, Units

# nuray vector index -> human label, used for the default output naming and the
# NuGeom.FluxLocation run attribute.
LOCATIONS = {0: "random", 1: "ND", 2: "FD"}

# nuray.wgt is calcEnuWgt's flux through a 1 m-radius disk (area pi*(100 cm)^2 =
# pi m^2). Multiply by this to get the per-cm^2 weight the main code expects:
# 1e-4 (m^2 -> cm^2) and 1/pi (per pi-m^2 disk -> per m^2).
WGT_TO_PER_CM2 = 1e-4 / math.pi

# Default beam->detector translation (cm) recorded for NuGeometry to apply AT
# LOAD. dk2nu rays are in the BEAM frame; the GDML
# (nd_hall_with_lar_tms_nosand.gdml) uses a hall-local frame whose ND-LAr active
# argon (volArgonInner) center is at (-50.0, -58.25, 667.0) cm. The ND rays are
# aimed at the documented on-axis ND nominal point (0, 0, 57400) cm in beam
# coordinates, so the translation that maps beam-frame coords onto the detector
# is (detector_center_hall - ND_nominal_beam):
#     (-50.0, -58.25, 667.0) - (0, 0, 57400) = (-50.0, -58.25, -56733.0).
# The z term dominates and is well-determined (574 m on-axis ND); the x/y terms
# are the GDML LAr-center offsets. Override with --det-offset if the geometry or
# the beam-axis<->hall-origin mapping is refined.
DEFAULT_DET_OFFSET = (-50.0, -58.25, -56733.0)

# Transverse flux window (cm), in the BEAM frame, over which ray origins are
# spread (see the "Flux window" note below). The dk2nu ND ray (location 1) is
# precomputed for the single on-axis point (0,0,57400) cm; emitting it verbatim
# makes EVERY ray converge to that one point, so the downstream detector
# simulation samples a single chord and recovers only the per-cm^2 flux rate --
# low by ~the detector face area (~2.4e5 cm^2 for ND-LAr). We instead shift each
# ray transversely by a uniform random (dx,dy) within this window, so the rays
# Monte-Carlo sample the detector face; the recorded area lets NuGeometry turn
# the per-cm^2 weight back into a total rate. Default = the ND-LAr active argon
# (ArgonCubeDetector75) transverse face, 700 cm (x) by 342.424 cm (y), centered
# on the beam axis (the LAr center maps from beam (0,0,57400)).
DEFAULT_WINDOW_WIDTH = 700.0
DEFAULT_WINDOW_HEIGHT = 342.424

# Rest masses (GeV) of the decay parents that appear in DUNE beam files, keyed
# by |PDG|. Used to give the incoming parent particle a valid on-shell energy.
# Unknown parents fall back to massless (E = |p|).
PARENT_MASS = {
    211: 0.13957,   # pi+/-
    321: 0.49368,   # K+/-
    130: 0.49761,   # K0L
    310: 0.49761,   # K0S
    311: 0.49761,   # K0
    13: 0.10566,    # mu+/-
    2212: 0.93827,  # proton
    2112: 0.93957,  # neutron
}

# Branches read per chunk. Kept minimal so iteration stays fast on ~1M entries.
BRANCHES = [
    "dk2nu/decay/decay.ntype",
    "dk2nu/decay/decay.nimpwt",
    "dk2nu/decay/decay.vx",
    "dk2nu/decay/decay.vy",
    "dk2nu/decay/decay.vz",
    "dk2nu/decay/decay.ptype",
    "dk2nu/decay/decay.pdpx",
    "dk2nu/decay/decay.pdpy",
    "dk2nu/decay/decay.pdpz",
    "dk2nu/nuray/nuray.px",
    "dk2nu/nuray/nuray.py",
    "dk2nu/nuray/nuray.pz",
    "dk2nu/nuray/nuray.E",
    "dk2nu/nuray/nuray.wgt",
    "dk2nu/tgtexit/tgtexit.tpx",
    "dk2nu/tgtexit/tgtexit.tpy",
    "dk2nu/tgtexit/tgtexit.tpz",
    "dk2nu/tgtexit/tgtexit.tptype",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Convert dk2nu ROOT files to a NuHepMC (HepMC3 ASCII) flux file."
    )
    p.add_argument("files", nargs="+", help="Input dk2nu ROOT files")
    p.add_argument("-o", "--output", default=None,
                   help="Output HepMC3 file (default: flux_<location>.hepmc3)")
    p.add_argument("--location", type=int, default=1, choices=(0, 1, 2),
                   help="nuray index to emit: 0=random, 1=ND, 2=FD (default: 1)")
    p.add_argument("--flavor", type=int, default=None,
                   help="Keep only this neutrino PDG (e.g. 14 for numu). "
                        "Default: keep all flavors.")
    p.add_argument("--max-events", type=int, default=None,
                   help="Stop after writing this many events (default: no limit)")
    p.add_argument("--det-offset", type=float, nargs=3,
                   metavar=("DX", "DY", "DZ"), default=DEFAULT_DET_OFFSET,
                   help="Beam->detector translation in cm, recorded as metadata "
                        "for NuGeometry to apply AT LOAD. Rays themselves are "
                        "left in the beam frame. Default %(default)s maps the "
                        "beam-frame on-axis ND (0,0,57400) onto the LAr center "
                        "of nd_hall_with_lar_tms_nosand.gdml.")
    p.add_argument("--window-width", type=float, default=DEFAULT_WINDOW_WIDTH,
                   metavar="W",
                   help="Transverse flux-window width in x (cm, beam frame) over "
                        "which ray origins are spread. Default %(default)s (ND-LAr "
                        "face). Set to 0 to disable spreading (rays converge to a "
                        "point and the recorded area is 1 cm^2 -- the per-cm^2 flux "
                        "rate, NOT integrated over the detector face).")
    p.add_argument("--window-height", type=float, default=DEFAULT_WINDOW_HEIGHT,
                   metavar="H",
                   help="Transverse flux-window height in y (cm, beam frame). "
                        "Default %(default)s (ND-LAr face).")
    p.add_argument("--step", type=int, default=50000,
                   help="uproot iteration chunk size (default: 50000)")
    return p.parse_args()


def total_pot(files):
    pot = 0.0
    for f in files:
        with uproot.open(f"{f}:dkmetaTree") as meta:
            pot += float(np.sum(meta["dkmeta/pots"].array(library="np")))
    return pot


def survey_rays(files, loc, flavor, max_events, step):
    """Count the rays that will be written and find their energy range.

    HepMC3 ASCII stores the run info ahead of the first event, so these have to
    be known before the writing loop starts. NuGeometry reads them back as
    NuGeom.Flux.NRays / NuGeom.Flux.EnergyRange_GeV / NuGeom.Flux.TotalWeight,
    which is what lets its flux streamer normalize POT per ray, check the
    generator's energy coverage, and size the importance-sampling payoff
    without ever scanning the converted file. Skipping this survey
    only costs NuGeometry a cheap byte-count pass at open time, so the filters
    below must mirror the writing loop's exactly or the POT normalization
    silently drifts.
    """
    sources = [f"{f}:dk2nuTree" for f in files]
    branches = ["dk2nu/decay/decay.ntype", "dk2nu/decay/decay.nimpwt",
                "dk2nu/nuray/nuray.E", "dk2nu/nuray/nuray.wgt"]
    n = 0
    emin, emax = np.inf, -np.inf
    # Sum of the per-cm^2 flux weights. NuGeometry needs mean(w) BEFORE it
    # throws its first ray: its acceptance envelope only ever grows, so a run
    # that starts on raw weights locks in a max(w)-sized envelope and can never
    # benefit from importance sampling afterwards.
    total_w = 0.0
    for chunk in uproot.iterate(sources, branches, step_size=step, library="np"):
        ntype = chunk["dk2nu/decay/decay.ntype"]
        nimpwt = chunk["dk2nu/decay/decay.nimpwt"]
        E = np.array([r[loc] for r in chunk["dk2nu/nuray/nuray.E"]])
        wgt = np.array([r[loc] for r in chunk["dk2nu/nuray/nuray.wgt"]])

        keep = np.ones(len(ntype), dtype=bool)
        if flavor is not None:
            keep &= (ntype.astype(int) == flavor)
        w = nimpwt * wgt * WGT_TO_PER_CM2
        keep &= np.isfinite(w) & (w > 0.0)

        if max_events is not None and n + int(keep.sum()) >= max_events:
            # Honour the same cut-off as the writing loop: keep only the first
            # (max_events - n) surviving rays of this chunk.
            idx = np.flatnonzero(keep)[: max_events - n]
            keep = np.zeros(len(ntype), dtype=bool)
            keep[idx] = True
            n += int(keep.sum())
            if keep.any():
                emin = min(emin, float(E[keep].min()))
                emax = max(emax, float(E[keep].max()))
                total_w += float(w[keep].sum())
            break

        n += int(keep.sum())
        if keep.any():
            emin = min(emin, float(E[keep].min()))
            emax = max(emax, float(E[keep].max()))
            total_w += float(w[keep].sum())
    return n, emin, emax, total_w


def make_run_info(files, pot, location, offset, window, nrays=0,
                  energy_range=None, total_weight=0.0):
    width, height = window
    area = width * height if (width > 0.0 and height > 0.0) else 1.0
    ri = GenRunInfo()
    ri.tools = [GenRunInfo.ToolInfo(
        name="dk2nu_to_nuhepmc.py", version="1.1",
        description="dk2nu->NuHepMC flux converter (uproot+pyhepmc, no ROOT)")]
    ri.weight_names = ["flux"]
    # String-valued attributes (HepMC3 ASCII stores attributes as strings).
    ri.attributes = {
        # NuHepMC standard exposure (GC.4): POT this flux sample corresponds to.
        "NuHepMC.Exposure.POT": repr(pot),
        # NuGeometry coordinate metadata: the load-time transform contract.
        "NuGeom.CoordinateFrame": "beam",
        "NuGeom.LengthUnit": "cm",
        "NuGeom.FluxLocation": LOCATIONS[location],
        # Beam->detector translation (cm) to be APPLIED AT LOAD. Identity here
        # unless --det-offset was given; NuGeometry may also override it.
        "NuGeom.BeamToDetector.Translation_cm":
            ",".join(repr(float(o)) for o in offset),
        # Transverse flux-window area (cm^2) the ray origins were spread over.
        # NuGeometry folds this into the per-cm^2 flux weight so the event rate
        # is integrated over the detector face (see the converter docstring).
        "NuGeom.FluxWindow.Area_cm2": repr(float(area)),
        "NuGeom.FluxWindow.Size_cm": f"{float(width)!r},{float(height)!r}",
        "NuGeom.SourceFiles": ";".join(files),
    }
    # Lets NuGeometry's flux streamer open the file without reading it: the
    # count normalizes POT per ray, the range tells the generator which beam
    # energies it must cover. See survey_rays().
    if nrays:
        ri.attributes["NuGeom.Flux.NRays"] = str(int(nrays))
    if energy_range is not None:
        emin, emax = energy_range
        ri.attributes["NuGeom.Flux.EnergyRange_GeV"] = f"{float(emin)!r},{float(emax)!r}"
    # Flux-weight summary. NuGeometry uses mean(w) to arm importance sampling
    # from the very first ray (see CachedFluxStreamer): drawing rays with
    # probability w_i/sum(w) and handing each the mean weight is unbiased, and
    # shrinks the acceptance envelope by max(w)/mean(w) -- ~459x on the DUNE ND
    # flux -- so the same POT needs that many fewer rays.
    if total_weight > 0.0:
        ri.attributes["NuGeom.Flux.TotalWeight"] = repr(float(total_weight))
        if nrays:
            ri.attributes["NuGeom.Flux.MeanWeight"] = repr(float(total_weight) / int(nrays))
    return ri


def main():
    args = parse_args()
    out = args.output or f"flux_{LOCATIONS[args.location]}.hepmc3"
    loc = args.location

    window = (args.window_width, args.window_height)
    spread = args.window_width > 0.0 and args.window_height > 0.0
    half_w = 0.5 * args.window_width
    half_h = 0.5 * args.window_height
    rng = np.random.default_rng()

    pot = total_pot(args.files)
    print(f"Total POT: {pot:.6e}")
    print(f"Emitting nuray[{loc}] ({LOCATIONS[loc]}) -> {out}")
    if spread:
        print(f"Spreading ray origins over a {args.window_width} x "
              f"{args.window_height} cm flux window "
              f"(area {args.window_width * args.window_height:.6e} cm^2)")
    else:
        print("Flux window disabled (--window-width/height <= 0): rays converge "
              "to a point and the recorded area is 1 cm^2. The event rate will be "
              "the per-cm^2 flux rate, NOT integrated over the detector face.")
    if any(args.det_offset):
        print(f"Recording beam->detector translation (cm, applied at load): "
              f"{tuple(args.det_offset)}")
    else:
        print("No beam->detector offset given; recording identity. NuGeometry "
              "applies the transform at load.")

    n_survey, e_lo, e_hi, w_tot = survey_rays(args.files, loc, args.flavor,
                                              args.max_events, args.step)
    mean_w = w_tot / n_survey if n_survey else 0.0
    print(f"Survey: {n_survey} rays, E in [{e_lo:.6g}, {e_hi:.6g}] GeV, "
          f"mean weight {mean_w:.6g} /cm^2")

    run_info = make_run_info(args.files, pot, loc, args.det_offset, window,
                             nrays=n_survey,
                             energy_range=(e_lo, e_hi) if n_survey else None,
                             total_weight=w_tot)
    sources = [f"{f}:dk2nuTree" for f in args.files]

    n_written = 0
    n_skipped = 0
    with pyhepmc.open(out, "w") as writer:
        for chunk in uproot.iterate(sources, BRANCHES, step_size=args.step,
                                    library="np"):
            ntype = chunk["dk2nu/decay/decay.ntype"]
            nimpwt = chunk["dk2nu/decay/decay.nimpwt"]
            vx = chunk["dk2nu/decay/decay.vx"]
            vy = chunk["dk2nu/decay/decay.vy"]
            vz = chunk["dk2nu/decay/decay.vz"]
            ptype = chunk["dk2nu/decay/decay.ptype"]
            pdpx = chunk["dk2nu/decay/decay.pdpx"]
            pdpy = chunk["dk2nu/decay/decay.pdpy"]
            pdpz = chunk["dk2nu/decay/decay.pdpz"]
            # nuray.* are jagged (3 entries each); index the chosen location.
            E = np.array([r[loc] for r in chunk["dk2nu/nuray/nuray.E"]])
            px = np.array([r[loc] for r in chunk["dk2nu/nuray/nuray.px"]])
            py = np.array([r[loc] for r in chunk["dk2nu/nuray/nuray.py"]])
            pz = np.array([r[loc] for r in chunk["dk2nu/nuray/nuray.pz"]])
            wgt = np.array([r[loc] for r in chunk["dk2nu/nuray/nuray.wgt"]])
            tpx = chunk["dk2nu/tgtexit/tgtexit.tpx"]
            tpy = chunk["dk2nu/tgtexit/tgtexit.tpy"]
            tpz = chunk["dk2nu/tgtexit/tgtexit.tpz"]
            tptype = chunk["dk2nu/tgtexit/tgtexit.tptype"]

            for i in range(len(ntype)):
                if args.flavor is not None and int(ntype[i]) != args.flavor:
                    continue
                w = float(nimpwt[i] * wgt[i] * WGT_TO_PER_CM2)
                if not np.isfinite(w) or w <= 0.0:
                    n_skipped += 1
                    continue

                evt = GenEvent(Units.GEV, Units.CM)
                evt.event_number = n_written
                evt.run_info = run_info

                # Decay vertex (beam frame, cm): parent hadron in, neutrino out.
                # HepMC3's ASCII reader requires a vertex to have both an
                # incoming and an outgoing particle, and this is also the
                # physically faithful picture of the meson decay.
                pmag2 = float(pdpx[i]**2 + pdpy[i]**2 + pdpz[i]**2)
                pmass = PARENT_MASS.get(abs(int(ptype[i])), 0.0)
                parent_E = np.sqrt(pmag2 + pmass * pmass)
                parent = GenParticle(
                    FourVector(float(pdpx[i]), float(pdpy[i]), float(pdpz[i]),
                               parent_E),
                    int(ptype[i]), 4)  # status 4: incoming beam/decay parent

                # Outgoing neutrino (status 1) carrying the chosen ray's
                # 4-momentum. Its production vertex is the ray origin that
                # NuGeometry consumes; direction = momentum, energy = E.
                nu = GenParticle(
                    FourVector(float(px[i]), float(py[i]), float(pz[i]),
                               float(E[i])),
                    int(ntype[i]), 1)

                # Spread the ray transversely across the flux window: shifting
                # the decay vertex by (dx,dy,0) rigidly translates the whole ray
                # (direction unchanged), so its crossing point at the detector
                # plane moves by exactly (dx,dy). The precomputed nuray weight is
                # kept as-is -- a uniform-flux-over-window approximation that is
                # excellent at the ND (the window is ~m-scale, the baseline
                # ~574 m), and exactly recovers the missing face-area factor.
                dx = rng.uniform(-half_w, half_w) if spread else 0.0
                dy = rng.uniform(-half_h, half_h) if spread else 0.0
                vtx = GenVertex(FourVector(float(vx[i]) + dx, float(vy[i]) + dy,
                                           float(vz[i]), 0.0))
                vtx.add_particle_in(parent)
                vtx.add_particle_out(nu)
                evt.add_vertex(vtx)

                evt.weights = [w]
                evt.attributes = {
                    "NuGeom.flux_weight": repr(w),
                    "dk2nu.nimpwt": repr(float(nimpwt[i])),
                    "dk2nu.nuray_wgt": repr(float(wgt[i])),
                    # Parent hadron provenance at target exit.
                    "dk2nu.tgtexit.tptype": str(int(tptype[i])),
                    "dk2nu.tgtexit.tp": ",".join(
                        repr(float(v)) for v in (tpx[i], tpy[i], tpz[i])),
                }
                writer.write(evt)
                n_written += 1
                if args.max_events is not None and n_written >= args.max_events:
                    print(f"Reached --max-events={args.max_events}, stopping.")
                    print(f"Wrote {n_written} events to {out} "
                          f"({n_skipped} skipped)")
                    return

    print(f"Wrote {n_written} events to {out}"
          + (f" ({n_skipped} skipped: non-finite or non-positive weight)"
             if n_skipped else ""))
    # NuGeom.Flux.NRays is what NuGeometry divides the file POT by, so a survey
    # that disagrees with what was written would misnormalize every run.
    if n_written != n_survey:
        print(f"WARNING: survey predicted {n_survey} rays but {n_written} were "
              f"written; NuGeom.Flux.NRays is wrong and the per-ray POT will be "
              f"off by {n_survey / max(n_written, 1):.6f}. Please report this.")


if __name__ == "__main__":
    main()
