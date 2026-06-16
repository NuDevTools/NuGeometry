#!/usr/bin/env python3
"""Compare NuGeometry sequential segments vs ROOT TGeoNavigator segments.

Both files have lines "rayIndex  len:mat len:mat ...".  We apply NuGeometry's
PruneSegments rule (drop <1um segments, merge adjacent same-material) to BOTH
sides so the only remaining differences are genuine geometry/material disagreements.
"""
import sys

THR = 1e-4  # cm, matches PruneSegments kPruneThreshold


def load(path):
    rays = {}
    for line in open(path):
        p = line.split()
        if not p:
            continue
        rays[int(p[0])] = [(float(t.rsplit(":", 1)[0]), t.rsplit(":", 1)[1]) for t in p[1:]]
    return rays


def prune(segs):
    out = []
    for l, m in segs:
        if l < THR:
            continue
        if out and out[-1][1] == m:
            out[-1][0] += l
        else:
            out.append([l, m])
    return out


a = load(sys.argv[1] if len(sys.argv) > 1 else "tools/nugeom_segments.txt")
b = load(sys.argv[2] if len(sys.argv) > 2 else "tools/root_segments.txt")

nmatch = 0
for i in sorted(a):
    A, B = prune(a[i]), prune(b.get(i, []))
    tol = lambda x: max(1e-2, 1e-3 * x)  # 0.1mm or 0.1%
    ok = len(A) == len(B) and all(
        am == bm and abs(al - bl) <= tol(al) for (al, am), (bl, bm) in zip(A, B)
    )
    if ok:
        nmatch += 1
        print(f"ray {i:2d}: MATCH  ({len(A)} segments, total {sum(l for l,_ in A):.2f} cm)")
    else:
        print(f"ray {i:2d}: DIFFER  nugeom={len(A)} segs, root={len(B)} segs")
        for k in range(max(len(A), len(B))):
            sa = f"{A[k][0]:.4f}:{A[k][1]}" if k < len(A) else "-"
            sb = f"{B[k][0]:.4f}:{B[k][1]}" if k < len(B) else "-"
            flag = "" if (k < len(A) and k < len(B) and A[k][1] == B[k][1]
                          and abs(A[k][0] - B[k][0]) <= tol(A[k][0])) else "  <--"
            print(f"        nugeom={sa:32s} root={sb:32s}{flag}")

print(f"\n{nmatch}/{len(a)} rays match after prune/merge")
