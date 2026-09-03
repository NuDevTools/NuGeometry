#!/usr/bin/env python3
import sys, math, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats
from mode_compare import parse, summary_from_log

ENV, TOT = '#2a78d6', '#eb6834'          # categorical slots 1,2
SURF, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#dcdcd8'

tag = sys.argv[1] if len(sys.argv) > 1 else 'ab'
env = parse(f'{tag}_env.hepmc'); tot = parse(f'{tag}_tot.hepmc')
senv = summary_from_log(f'{tag}_env.log'); stot = summary_from_log(f'{tag}_tot.log')

print(f"envelope: {len(env['w'])} events parsed  (log says {senv['n']}), POT={senv['pot']:.4e}")
print(f"total   : {len(tot['w'])} events parsed  (log says {stot['n']}), POT={stot['pot']:.4e}")

# ---- rates, each with its own correct estimator -------------------------
r_env = env['w'].sum()/senv['pot'];  e_env = r_env/math.sqrt(len(env['w']))
r_tot = len(tot['w'])/stot['pot'];   e_tot = r_tot/math.sqrt(len(tot['w']))
R = r_env/r_tot; eR = R*math.hypot(1/math.sqrt(len(env['w'])), 1/math.sqrt(len(tot['w'])))
print(f"\nRATE  envelope sum(w)/POT = {r_env:.4e} +/- {e_env:.2e}")
print(f"RATE  total    events/POT = {r_tot:.4e} +/- {e_tot:.2e}")
print(f"R = {R:.4f} +/- {eR:.4f}   ({abs(R-1)/eR:.1f} sigma from 1)\n")

for d in (env, tot):
    d['theta'] = np.degrees(np.arccos(np.clip(d['cos'], -1, 1)))
VARS = [('Enu', r'$E_\nu$ [GeV]'), ('Emu', r'$E_\mu$ [GeV]'),
        ('Q2',  r'$Q^2$ [GeV$^2$]'), ('theta', r'$\theta_\mu$ [deg]'),
        ('z',   'vertex $z$ [cm]'), ('x',  'vertex $x$ [cm]')]

fig = plt.figure(figsize=(15, 10), facecolor=SURF)
gs = GridSpec(5, 3, height_ratios=[3,1,0.9,3,1], hspace=0.08, wspace=0.28)
results = []
for i,(k,xlabel) in enumerate(VARS):
    row, col = (0 if i < 3 else 3), i % 3
    ax = fig.add_subplot(gs[row, col]); axr = fig.add_subplot(gs[row+1, col], sharex=ax)
    a, b = env[k], tot[k]
    lo, hi = np.percentile(np.concatenate([a,b]), [0.5, 99.5])
    bins = np.linspace(lo, hi, 31); c = 0.5*(bins[1:]+bins[:-1])

    h1,_ = np.histogram(a, bins=bins, weights=env['w'])
    s1,_ = np.histogram(a, bins=bins, weights=env['w']**2)
    # envelope: weighted.  total: UNWEIGHTED (see the module docstring).
    h2,_ = np.histogram(b, bins=bins)
    s2,_ = np.histogram(b, bins=bins)
    W1, W2 = h1.sum(), h2.sum()
    n1, e1 = h1/W1, np.sqrt(s1)/W1
    n2, e2 = h2/W2, np.sqrt(s2)/W2

    ok = (e1**2 + e2**2) > 0
    chi2 = (((n1-n2)**2)/(e1**2+e2**2))[ok].sum(); dof = ok.sum()-1
    p_chi = stats.chi2.sf(chi2, dof)
    p_ks  = stats.ks_2samp(a, b).pvalue
    results.append((k, chi2, dof, p_chi, p_ks))

    ax.step(c, n1, where='mid', color=ENV, lw=2, label='envelope (NoRetry)')
    ax.step(c, n2, where='mid', color=TOT, lw=2, label='total (Retry)')
    ax.set_facecolor(SURF); ax.grid(True, color=GRID, lw=0.6, alpha=.8); ax.set_axisbelow(True)
    ax.tick_params(labelbottom=False, colors=INK2, labelsize=8)
    ax.set_ylabel('fraction / bin', color=INK2, fontsize=9)
    ax.set_title(f'{xlabel}   $\\chi^2$/dof={chi2:.0f}/{dof}, p={p_chi:.2f}',
                 color=INK, fontsize=10, pad=6)
    for s in ax.spines.values(): s.set_color(GRID)
    if i == 0: ax.legend(frameon=False, fontsize=9, labelcolor=INK2)

    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(n2>0, n1/n2, np.nan)
        rerr  = np.where(n2>0, np.sqrt((e1/n2)**2 + (n1*e2/n2**2)**2), np.nan)
    axr.axhline(1, color=INK2, lw=1, ls='--')
    axr.errorbar(c, ratio, yerr=rerr, fmt='o', ms=3, color=ENV, ecolor=ENV, elinewidth=1)
    axr.set_facecolor(SURF); axr.grid(True, color=GRID, lw=0.6, alpha=.8); axr.set_axisbelow(True)
    axr.set_ylim(0.6, 1.4); axr.set_xlabel(xlabel, color=INK2, fontsize=9)
    axr.set_ylabel('env/tot', color=INK2, fontsize=8)
    axr.tick_params(colors=INK2, labelsize=8)
    for s in axr.spines.values(): s.set_color(GRID)

fig.suptitle(f'DetectorSim sampling modes: envelope sum(w) vs total counts  '
             f'(POT = {stot["pot"]:.3g}; {len(env["w"])} envelope vs {len(tot["w"])} total events)',
             color=INK, fontsize=12, y=0.965)
fig.savefig(f'{tag}_mode_compare.png', dpi=140, facecolor=SURF, bbox_inches='tight')
print(f"{'var':6s} {'chi2':>8s} {'dof':>4s} {'p(chi2)':>9s} {'p(KS)':>9s}")
for k,c2,d,pc,pk in results: print(f"{k:6s} {c2:8.1f} {d:4d} {pc:9.3f} {pk:9.3f}")
print(f"\nwrote {tag}_mode_compare.png")

# ---- material composition (categorical, weight-weighted) ---------------
import collections
mats = sorted(set(env['mat']) | set(tot['mat']))
print(f"\n{'material':16s} {'envelope':>12s} {'total':>12s}   (weight-weighted fraction)")
rows=[]
for m in mats:
    f1 = env['w'][env['mat']==m].sum()/env['w'].sum()
    f2 = (tot['mat']==m).sum()/len(tot['mat'])
    n1 = (env['mat']==m).sum(); n2 = (tot['mat']==m).sum()
    if n1+n2 >= 10: rows.append((m,f1,f2,n1,n2))
    print(f"{m:16s} {f1:12.5f} {f2:12.5f}")
# 2-sample chi2 on raw counts of the well-populated categories
obs = np.array([[r[3] for r in rows],[r[4] for r in rows]], dtype=float)
chi2c, pc, dofc, _ = stats.chi2_contingency(obs)
print(f"material composition: chi2 = {chi2c:.1f}/{dofc}, p = {pc:.3f}")
