#!/usr/bin/env python3
"""Statistical + graphical comparison of DetectorSim's two sampling modes.

Shape convention (see mode_discrepancy_findings.md):
  Both modes' *weight-weighted* distributions are the physical differential
  rate: the partial unweighter emits with prob min(1,r) and weight max(1,r),
  so E[w] ~ f*r ~ the true density in BOTH modes.  What differs is the
  normalization estimator: events/POT in total, sum(w)/POT in envelope.
So: compare w-weighted SHAPES, and compare rates with the per-mode estimator.
"""
import sys, math, numpy as np

MU = 0.1056583745

def parse(fn, limit=None):
    """Stream a HepMC3 ASCII file -> per-event kinematics."""
    ev = None
    out = {k: [] for k in ('w','Enu','Emu','Q2','cos','x','y','z','mat')}
    nu = mu = None
    def flush():
        if ev is None or nu is None or mu is None: return
        (pxn,pyn,pzn,en) = nu
        (pxm,pym,pzm,em) = mu
        pn = math.sqrt(pxn*pxn+pyn*pyn+pzn*pzn)
        pm = math.sqrt(pxm*pxm+pym*pym+pzm*pzm)
        if pn <= 0 or pm <= 0: return
        cos = (pxn*pxm+pyn*pym+pzn*pzm)/(pn*pm)
        q2 = 2*en*em - 2*(pxn*pxm+pyn*pym+pzn*pzm) - MU*MU   # -(p_nu - p_mu)^2
        out['w'].append(ev['w']);   out['Enu'].append(en);  out['Emu'].append(em)
        out['Q2'].append(q2);       out['cos'].append(cos)
        out['x'].append(ev['x']);   out['y'].append(ev['y']); out['z'].append(ev['z'])
        out['mat'].append(ev['mat'])
    with open(fn) as f:
        for line in f:
            t = line[0]
            if t == 'E' and line.startswith('E '):
                flush()
                if limit and len(out['w']) >= limit: break
                ev = {'w':1.0,'x':np.nan,'y':np.nan,'z':np.nan,'mat':'?'}
                nu = mu = None
            elif ev is None:
                continue
            elif t == 'W':
                ev['w'] = float(line.split()[1])
            elif t == 'A':
                p = line.split()
                if len(p) >= 4 and p[2] == 'lab_pos':
                    ev['x'],ev['y'],ev['z'] = float(p[3]),float(p[4]),float(p[5])
                elif len(p) >= 4 and p[2] == 'NuGeom.material':
                    ev['mat'] = p[3]
            elif t == 'P':
                p = line.split()
                pdg, status = int(p[3]), int(p[9])
                vec = (float(p[4]),float(p[5]),float(p[6]),float(p[7]))
                if status == 4 and abs(pdg) in (12,14,16): nu = vec
                elif status == 1 and abs(pdg) == 13:       mu = vec
        flush()
    return {k: np.asarray(v) if k != 'mat' else np.asarray(v) for k,v in out.items()}

def summary_from_log(log):
    import re
    t = open(log).read()
    m = re.search(r'Run summary \((\w+) mode\): (\d+) events / (\d+) thrown rays, POT = ([\d.e+-]+)', t)
    return dict(mode=m.group(1), n=int(m.group(2)), rays=int(m.group(3)), pot=float(m.group(4)))
