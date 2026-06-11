#!/usr/bin/env python3
"""Convert dk2nu ROOT flux files to a NuHepMC (HepMC3 ASCII) flux file.

No ROOT dependency: branches are read with uproot and events are written with
pyhepmc (which bundles HepMC3). The result is consumed by NuGeometry, whose
FluxSource/GeneratorInterface already speak NuHepMC.

Each dk2nu entry stores a vector<NuRay> with three precomputed rays:
    [0] random decay direction   [1] near detector   [2] far detector
We emit one HepMC3 event per entry for a chosen ray location (default ND = 1).
The ray already carries the neutrino energy, direction and per-cm^2 importance
weight from g4lbne, so calcEnuWgt / decay_through_point (the pydk2nu/ROOT path)
is not needed.

Coordinate convention (IMPORTANT)
---------------------------------
dk2nu decay vertices live in the BEAM frame: target at the origin, +z along
the beam, units cm. They are far upstream of the detector in the GDML's
hall-local frame. We deliberately write rays in their NATIVE beam frame and
record the frame + (optional) beam->detector translation as run-level metadata,
so NuGeometry can apply the transform once at load time rather than having a
hidden offset baked into the file. See the NuGeom.* run attributes below.

Per-event central weight is the dk2nu flux weight
    w = decay.nimpwt * nuray.wgt[location]
(the 1/pi factor in DUNE's plottingExamples.py is a per-degree plotting
normalization, not part of the ray weight, so it is not applied here).
"""

import argparse
import sys

import numpy as np
import uproot
import pyhepmc
from pyhepmc import GenEvent, GenParticle, GenVertex, GenRunInfo, FourVector, Units

# nuray vector index -> human label, used for the default output naming and the
# NuGeom.FluxLocation run attribute.
LOCATIONS = {0: "random", 1: "ND", 2: "FD"}

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
    p.add_argument("--step", type=int, default=50000,
                   help="uproot iteration chunk size (default: 50000)")
    return p.parse_args()


def total_pot(files):
    pot = 0.0
    for f in files:
        with uproot.open(f"{f}:dkmetaTree") as meta:
            pot += float(np.sum(meta["dkmeta/pots"].array(library="np")))
    return pot


def make_run_info(files, pot, location, offset):
    ri = GenRunInfo()
    ri.tools = [GenRunInfo.ToolInfo(
        name="dk2nu_to_nuhepmc.py", version="1.0",
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
        "NuGeom.SourceFiles": ";".join(files),
    }
    return ri


def main():
    args = parse_args()
    out = args.output or f"flux_{LOCATIONS[args.location]}.hepmc3"
    loc = args.location

    pot = total_pot(args.files)
    print(f"Total POT: {pot:.6e}")
    print(f"Emitting nuray[{loc}] ({LOCATIONS[loc]}) -> {out}")
    if any(args.det_offset):
        print(f"Recording beam->detector translation (cm, applied at load): "
              f"{tuple(args.det_offset)}")
    else:
        print("No beam->detector offset given; recording identity. NuGeometry "
              "applies the transform at load.")

    run_info = make_run_info(args.files, pot, loc, args.det_offset)
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
                w = float(nimpwt[i] * wgt[i])
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

                vtx = GenVertex(FourVector(float(vx[i]), float(vy[i]),
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


if __name__ == "__main__":
    main()
