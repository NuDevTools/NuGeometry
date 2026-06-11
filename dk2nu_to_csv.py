#!/usr/bin/env python3
"""Convert dk2nu ROOT files to a CSV flux file for NuGeometry.

For each decay we sample one (or more) target points uniformly inside a
bounding box that encloses the DUNE near detector (LAr + TMS) in *beam-frame*
coordinates (cm), and use dk2nu's calcEnuWgt to compute the neutrino energy
and per-cm^2 importance weight aimed at that point. The starting position is
the sampled point and the direction is (point - decay_vertex), normalized.

Output columns:
    pid, wgt, x, y, z, dx, dy, dz, E

  pid       PDG code of the neutrino
  wgt       per-POT, per-ray flux weight in cm^-2 (already includes nimpwt and
            the 1/N_per_decay Monte-Carlo factor; sum over rows / total_POT
            gives the integrated flux through the sampled volume)
  x,y,z     starting position in beam-frame cm
  dx,dy,dz  unit direction vector
  E         neutrino energy in GeV

Note on coordinates: dk2nu decay vertices are in the beam frame (target at
origin, z along the beam axis). The GDML nd_hall_with_lar_tms_nosand.gdml uses
hall-local coordinates and does not carry a beam->hall transform, so we sample
directly in beam frame. The defaults below match the on-axis ND at z ~ 574 m.
"""

import argparse
import csv
import sys

import numpy as np
import pydk2nu as dk

# DUNE ND on-axis envelope in beam-frame cm. The ND-LAr active volume is
# roughly 7 m (x) x 3 m (y) x 5 m (z), and TMS extends ~5 m further downstream;
# this box loosely encloses both. Override with --xmin/--xmax/... if needed.
DEFAULT_BBOX = {
    "xmin": -350.0, "xmax":  350.0,
    "ymin": -150.0, "ymax":  150.0,
    "zmin": 57000.0, "zmax": 58500.0,
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Convert dk2nu ROOT files to a NuGeometry CSV flux file."
    )
    p.add_argument("files", nargs="+", help="Input dk2nu ROOT files")
    p.add_argument("-o", "--output", default="flux.csv",
                   help="Output CSV file (default: flux.csv)")
    p.add_argument("--rays-per-decay", type=int, default=1,
                   help="Number of target points to sample per decay (default: 1)")
    p.add_argument("--max-rays", type=int, default=None,
                   help="Stop after writing this many rays (default: no limit)")
    p.add_argument("--seed", type=int, default=0, help="RNG seed (default: 0)")
    p.add_argument("--validate", action="store_true",
                   help="After CSV writing, plot the bbox-sampled energy spectrum "
                        "against the on-axis reference (a la example_off_axis.py) "
                        "and report integrated flux/POT for each flavor.")
    p.add_argument("--validate-plot", default="validation.pdf",
                   help="Output path for the --validate plot (default: validation.pdf)")
    for k, v in DEFAULT_BBOX.items():
        p.add_argument(f"--{k}", type=float, default=v,
                       help=f"Sampling bbox {k} in beam-frame cm (default: {v})")
    return p.parse_args()


# DUNE ND on-axis reference point used by dk2nu/python/example_off_axis.py.
ON_AXIS_REF_CM = (0.0, 0.0, 57400.0)

# PDG codes -> labels for flavors that show up in DUNE beam files.
NU_LABELS = {12: "nue", -12: "nuebar", 14: "numu", -14: "numubar"}


def validate_against_on_axis(sampled_pid, sampled_E, sampled_wgt,
                             reader, total_pot, plot_path):
    """Plot bbox-sampled spectra against the on-axis reference.

    The bbox-sampled rays are a Monte-Carlo estimate of the average flux per
    cm^2 inside the user's sampling box (since each row's wgt was scaled by
    1/rays_per_decay, summing over rows reproduces the volume average).
    The on-axis reference is the flux at the single point used by
    dk2nu/python/example_off_axis.py: ``decay_all_through_point((0,0,57400))``.
    For a sampling box that hugs the on-axis ND the two spectra should agree.
    """
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    print(f"[validate] running on-axis reference at {ON_AXIS_REF_CM} cm "
          "(matches dk2nu/python/example_off_axis.py)")
    ref = reader.decay_all_through_point(ON_AXIS_REF_CM)  # Nx3: pid, E, wgt
    ref_pid = ref[:, 0].astype(int)
    ref_E = ref[:, 1]
    ref_wgt = ref[:, 2]

    sampled_pid = np.asarray(sampled_pid, dtype=int)
    sampled_E = np.asarray(sampled_E, dtype=float)
    sampled_wgt = np.asarray(sampled_wgt, dtype=float)

    # Match the binning used by example_off_axis.py.
    bins = np.linspace(0.0, 10.0, 200)
    dE = bins[1] - bins[0]

    print(f"[validate] {len(ref)} decays in on-axis reference, "
          f"{len(sampled_pid)} bbox-sampled rays")
    print(f"[validate] integrated flux/POT (cm^-2/POT):")
    print(f"[validate]   {'flavor':>8s}  {'bbox-sampled':>14s}  "
          f"{'on-axis ref':>14s}  {'ratio':>8s}")

    flavors_to_plot = [(14, "numu"), (-14, "numubar"),
                       (12, "nue"), (-12, "nuebar")]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for ax, (pid, label) in zip(axes.flat, flavors_to_plot):
        ref_mask = ref_pid == pid
        smp_mask = sampled_pid == pid
        ref_int = ref_wgt[ref_mask].sum() / total_pot
        smp_int = sampled_wgt[smp_mask].sum() / total_pot
        ratio = (smp_int / ref_int) if ref_int > 0 else float("nan")
        print(f"[validate]   {label:>8s}  {smp_int:>14.6e}  "
              f"{ref_int:>14.6e}  {ratio:>8.3f}")

        ax.hist(ref_E[ref_mask], bins=bins,
                weights=ref_wgt[ref_mask] / (total_pot * dE),
                histtype="step", lw=2, label="on-axis ref")
        ax.hist(sampled_E[smp_mask], bins=bins,
                weights=sampled_wgt[smp_mask] / (total_pot * dE),
                histtype="step", lw=2, ls="--", label="bbox sampled")
        ax.set_title(label)
        ax.set_yscale("log")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    for ax in axes[1]:
        ax.set_xlabel(r"$E_\nu$ (GeV)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"flux /cm$^2$ /GeV /POT")

    fig.suptitle("dk2nu_to_csv.py validation: bbox sampling vs on-axis "
                 "reference (example_off_axis.py)")
    fig.tight_layout()
    fig.savefig(plot_path, bbox_inches="tight")
    print(f"[validate] wrote validation plot to {plot_path}")


def main():
    args = parse_args()

    if args.xmin >= args.xmax or args.ymin >= args.ymax or args.zmin >= args.zmax:
        sys.exit("Invalid bounding box: min must be < max for all axes.")

    rng = np.random.default_rng(args.seed)
    lo = np.array([args.xmin, args.ymin, args.zmin])
    hi = np.array([args.xmax, args.ymax, args.zmax])

    reader = dk.dk2nuFileReader(args.files)
    total_pot = sum(m.pots for m in reader.metas())
    print(f"Total POT: {total_pot:.6e}")
    print(f"Sampling box (beam frame, cm): "
          f"x[{args.xmin},{args.xmax}] "
          f"y[{args.ymin},{args.ymax}] "
          f"z[{args.zmin},{args.zmax}]")

    n_per_decay = args.rays_per_decay
    inv_n = 1.0 / n_per_decay
    n_written = 0
    n_failed = 0

    # Per-ray accumulators used by --validate. Cheap when off; we still keep
    # them unconditionally so the main loop has a single shape.
    sampled_pid: list[int] = []
    sampled_E: list[float] = []
    sampled_wgt: list[float] = []

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pid", "wgt", "x", "y", "z", "dx", "dy", "dz", "E"])

        for entry in reader.decays():
            if args.max_rays is not None and n_written >= args.max_rays:
                break

            vx, vy, vz = entry.decay.vx, entry.decay.vy, entry.decay.vz
            pid = entry.decay.ntype

            for _ in range(n_per_decay):
                point = lo + rng.random(3) * (hi - lo)
                try:
                    E, wgt = entry.decay_through_point((point[0], point[1], point[2]))
                except RuntimeError:
                    n_failed += 1
                    continue

                d = point - np.array([vx, vy, vz])
                norm = np.linalg.norm(d)
                if norm == 0.0:
                    n_failed += 1
                    continue
                d /= norm

                row_wgt = wgt * inv_n
                writer.writerow([
                    pid, row_wgt,
                    point[0], point[1], point[2],
                    d[0], d[1], d[2],
                    E,
                ])
                sampled_pid.append(pid)
                sampled_E.append(E)
                sampled_wgt.append(row_wgt)
                n_written += 1
                if args.max_rays is not None and n_written >= args.max_rays:
                    break

    print(f"Wrote {n_written} rays to {args.output}"
          + (f" ({n_failed} skipped)" if n_failed else ""))

    if args.validate:
        validate_against_on_axis(sampled_pid, sampled_E, sampled_wgt,
                                 reader, total_pot, args.validate_plot)


if __name__ == "__main__":
    main()
