#!/usr/bin/env python3
"""ROOT cross-check gate: verify NuGeometry's traversal matches ROOT's
TGeoNavigator ray-for-ray on a set of GDML geometries.

For each geometry it:
  1. runs the [xcheck] test to emit rays + our (default boundary-sweep) segments,
  2. steps ROOT's TGeoNavigator over the same rays (check_overlaps_root.C),
  3. compares the (material, length) segment lists after the same prune/merge.

Exits non-zero if any ray on any geometry disagrees, so it can gate CI.

Usage:
  tools/root_xcheck.py [--build DIR] [--nrays N] [--root ROOT_EXE] [gdml ...]

With no gdml args it cross-checks every *.gdml in the repo root.  ROOT is found
via --root, $ROOT_EXE, or `root.exe`/`root` on PATH.
"""
import argparse, os, shutil, subprocess, sys
from collections import defaultdict
from pathlib import Path

THR = 1e-4  # cm, matches PruneSegments


def merge(path):
    rays = {}
    for line in open(path):
        x = line.split()
        if not x:
            continue
        out = []
        for t in x[1:]:
            l, m = t.rsplit(":", 1)
            l = float(l)
            if l < THR:
                continue
            if out and out[-1][1] == m:
                out[-1][0] += l
            else:
                out.append([l, m])
        rays[int(x[0])] = out
    return rays


def compare(ours, root):
    """Return (n_match, n_total, first_mismatch_msg)."""
    n_ok = 0
    first = ""
    for i in sorted(ours):
        a, b = ours[i], root.get(i, [])
        ok = len(a) == len(b) and all(
            am == bm and abs(al - bl) <= max(1e-2, 1e-3 * al)
            for (al, am), (bl, bm) in zip(a, b)
        )
        if ok:
            n_ok += 1
        elif not first:
            # summarize the per-material disagreement
            ma, mb = defaultdict(float), defaultdict(float)
            for l, m in a:
                ma[m] += l
            for l, m in b:
                mb[m] += l
            diffs = {
                m: round(ma.get(m, 0) - mb.get(m, 0), 3)
                for m in set(ma) | set(mb)
                if abs(ma.get(m, 0) - mb.get(m, 0)) > 0.05
            }
            first = f"ray {i}: nseg ours={len(a)} root={len(b)} matdiff(cm)={diffs}"
    return n_ok, len(ours), first


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default="build", help="build dir (has test/nugeometry-testsuite)")
    ap.add_argument("--nrays", type=int, default=400)
    ap.add_argument("--root", default=os.environ.get("ROOT_EXE", ""))
    ap.add_argument("gdml", nargs="*")
    a = ap.parse_args()

    repo = Path(__file__).resolve().parent.parent
    tools = repo / "tools"
    testbin = Path(a.build) / "test" / "nugeometry-testsuite"
    if not testbin.exists():
        sys.exit(f"test binary not found: {testbin} (build with -DENABLE_TESTING=ON)")
    root_exe = a.root or shutil.which("root.exe") or shutil.which("root")
    if not root_exe:
        sys.exit("ROOT not found (set --root or $ROOT_EXE, or put root/root.exe on PATH)")

    gdmls = [Path(g) for g in a.gdml] or sorted(repo.glob("*.gdml"))
    if not gdmls:
        sys.exit("no GDML geometries found")

    work = tools / "xcheck_work"
    work.mkdir(exist_ok=True)
    results, failed = [], False
    for gdml in gdmls:
        gdml = gdml.resolve()
        out = work / gdml.stem
        env = dict(os.environ,
                   NUGEOM_XCHECK_GDML=str(gdml),
                   NUGEOM_XCHECK_OUT=str(out),
                   NUGEOM_XCHECK_NRAYS=str(a.nrays))
        # 1. our segments
        r = subprocess.run([str(testbin), "[xcheck]"], env=env,
                           capture_output=True, text=True)
        if r.returncode != 0 or not (out.with_name(out.name + "_seg.txt")).exists():
            print(f"  {gdml.name:48s} EMIT FAILED\n{r.stdout}\n{r.stderr}")
            failed = True
            continue
        # 2. ROOT navigation
        macro = (f'navigate_rays("{gdml}","{out}_rays.txt","{out}_root.txt")')
        rr = subprocess.run([root_exe, "-l", "-b", "-q",
                             str(tools / "check_overlaps_root.C"), "-e", macro],
                            capture_output=True, text=True)
        if not Path(f"{out}_root.txt").exists():
            print(f"  {gdml.name:48s} ROOT FAILED\n{rr.stdout[-2000:]}\n{rr.stderr[-2000:]}")
            failed = True
            continue
        # 3. compare
        n_ok, n_tot, first = compare(merge(f"{out}_seg.txt"), merge(f"{out}_root.txt"))
        ok = n_ok == n_tot
        failed = failed or not ok
        results.append((gdml.name, n_ok, n_tot, "" if ok else first))

    print("\n==== ROOT cross-check summary ====")
    for name, n_ok, n_tot, first in results:
        tag = "OK " if n_ok == n_tot else "FAIL"
        print(f"  [{tag}] {name:52s} {n_ok}/{n_tot}")
        if first:
            print(f"         {first}")
    print("==================================")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
