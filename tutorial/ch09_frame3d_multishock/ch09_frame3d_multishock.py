# -*- coding: utf-8 -*-
"""
Chapter 09 - 3D steel moment frame: column removal + mainshock/aftershock sequence
==================================================================================
Units: kN, m, s  (so stress is kN/m^2 = kPa, mass is t = kN*s^2/m).  g = 9.81 m/s^2

What the script does
  * Builds a 3-story, 2x2-bay 3D steel moment frame with fiber sections
    (Steel02) and a rigid floor diaphragm (rigidDiaphragm).
  * Runs 4 cases, each in a freshly built model:
        Intact  - all columns present
        S1      - corner column (axis X1/Y1) removed in story 1
        S2      - the same corner column removed in story 2
        S3      - the same corner column removed in story 3
  * For every case: gravity -> eigen -> one multi-shock time history
    (mainshock + free-vibration gap + aftershock + short free-vibration tail)
    in ONE Path time series along global X.
  * Compares peak and residual inter-story drifts.

Tag convention (chapter 03, 3D version)
  node     = 10000*floor + 100*axisY + axisX      (floor 0 = base, e.g. 101 .. 30303)
  column   = 1000 + 100*story + 10*axisY + axisX  (1xxx, e.g. 1111 = story 1, Y1, X1)
  beam X   = 2000 + 100*floor + 10*axisY + axisX  (2xxx, axisX = start axis of the span)
  beam Y   = 2050 + 100*floor + 10*axisY + axisX  (2xxx, +50 marks the Y direction,
                                                    axisY = start axis of the span)
"""
import os
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import openseespy.opensees as ops

T_START = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, 'figs')
os.makedirs(FIGS, exist_ok=True)
REC_DIR = os.path.join(HERE, '..', 'records')

# %% 0. Parameters ------------------------------------------------------------
G = 9.81                      # m/s^2
NX, NY, NF = 3, 3, 3          # axes in X, axes in Y, number of floors (= stories)
BAY = 5.0                     # bay width in both directions [m]
H = 3.2                       # story height [m]
Q = 7.0                       # seismic gravity load D + 0.25L incl. self weight [kN/m^2]

# Steel S355
E_S = 2.0e8                   # Young's modulus [kN/m^2]
G_S = 7.7e7                   # shear modulus [kN/m^2]
FY = 355.0e3                  # yield stress [kN/m^2]
B_H = 0.01                    # strain-hardening ratio

# Columns: square hollow section 250x250x12
COL_B, COL_T = 0.250, 0.012
# Beams: IPE270
BM_D, BM_BF, BM_TF, BM_TW = 0.270, 0.135, 0.0102, 0.0066

# Ground motion sequence (one record after the other in ONE time series)
DT = 0.02                     # analysis step = record step (coarse on purpose) [s]
SF_MAIN = 1.5                 # mainshock scale factor  (El Centro NS, PGA 0.32 g)
MAIN_DUR = 15.0               # mainshock duration used [s]
GAP = 5.0                     # free-vibration gap between the shocks [s]
SF_AFTER = 0.5                # aftershock scale factor (Northridge LOS270, PGA 0.47 g)
AFTER_DUR = 10.0              # aftershock duration used [s]
TAIL = 5.0                    # free vibration after the aftershock [s]
ZETA = 0.05                   # Rayleigh damping ratio

MAX_HALVING = 4               # substep halving levels on non-convergence
CASE_BUDGET = 40.0            # wall-clock seconds allowed per time-history case

# Removal cases: (label, removed column as (story, axisX, axisY) or None)
CASES = [('Intact', None), ('S1', (1, 1, 1)), ('S2', (2, 1, 1)), ('S3', (3, 1, 1))]

# Tags of materials, sections, transformations, integrations
MAT_STEEL = 1
SEC_COL, SEC_BM = 1, 2
TR_COL, TR_BX, TR_BY = 1, 2, 3
INT_COL, INT_BM = 1, 2


def node_tag(floor, ix, iy):
    """Node tag = 10000*floor + 100*axisY + axisX."""
    return 10000 * floor + 100 * iy + ix


def col_tag(story, ix, iy):
    """Column tag = 1000 + 100*story + 10*axisY + axisX."""
    return 1000 + 100 * story + 10 * iy + ix


def beam_tag(floor, ix, iy, direction):
    """Beam tag: 2000 + 100*floor + 10*axisY + axisX (+50 for Y-direction beams)."""
    return 2000 + 100 * floor + 10 * iy + ix + (50 if direction == 'Y' else 0)


def xyz(floor, ix, iy):
    return (ix - 1) * BAY, (iy - 1) * BAY, floor * H


# Section properties needed for torsion (GJ) and for the text
A_BOX = COL_B**2 - (COL_B - 2 * COL_T)**2
a_m = COL_B - COL_T                                   # mid-line width of the box
J_BOX = 2.0 * COL_T * a_m**4 / (2.0 * a_m)            # Bredt: 2 t a^2 b^2 / (a+b)
J_IPE = (2.0 * BM_BF * BM_TF**3 + (BM_D - 2 * BM_TF) * BM_TW**3) / 3.0


# %% 1. Model (one function, called once per case) ----------------------------
def build_model(removed=None):
    """Build the 3D frame. removed = (story, axisX, axisY) of the column to omit."""
    ops.wipe()                                        # clear any previous model
    ops.model('basic', '-ndm', 3, '-ndf', 6)          # 3D model, 6 DOF per node

    # --- nodes, supports, masses
    for fl in range(NF + 1):
        for iy in range(1, NY + 1):
            for ix in range(1, NX + 1):
                tag = node_tag(fl, ix, iy)
                ops.node(tag, *xyz(fl, ix, iy))       # node(tag, X, Y, Z)
                if fl == 0:
                    ops.fix(tag, 1, 1, 1, 1, 1, 1)    # fixed base: all 6 DOF restrained
                else:
                    # tributary floor area: corner 1/4, edge 1/2, interior 1 panel
                    ax = 0.5 if ix in (1, NX) else 1.0
                    ay = 0.5 if iy in (1, NY) else 1.0
                    m = Q * ax * ay * BAY**2 / G
                    # mass(tag, mX, mY, mZ, IX, IY, IZ): horizontal mass only;
                    # the diaphragm constraint turns it into floor torsional inertia
                    ops.mass(tag, m, m, 0.0, 0.0, 0.0, 0.0)

    # --- rigid floor diaphragm: centre node (axis X2/Y2) is the retained node
    for fl in range(1, NF + 1):
        master = node_tag(fl, 2, 2)
        slaves = [node_tag(fl, ix, iy) for iy in range(1, NY + 1)
                  for ix in range(1, NX + 1) if (ix, iy) != (2, 2)]
        # rigidDiaphragm(perpDirn, retainedNode, *constrainedNodes):
        # perpDirn=3 -> diaphragm in the X-Y plane; ties UX, UY, RZ of the floor
        ops.rigidDiaphragm(3, master, *slaves)

    # --- material: Steel02 = Giuffre-Menegotto-Pinto steel with kinematic hardening
    # Steel02(tag, Fy, E, b, R0, cR1, cR2)
    ops.uniaxialMaterial('Steel02', MAT_STEEL, FY, E_S, B_H, 18.0, 0.925, 0.15)

    # --- column fiber section: square hollow box, 4 wall patches
    # section('Fiber', tag, '-GJ', GJ): 3D fiber section + elastic torsion stiffness
    ops.section('Fiber', SEC_COL, '-GJ', G_S * J_BOX)
    b, t = COL_B / 2.0, COL_T
    # patch('rect', mat, nFibY, nFibZ, yI, zI, yJ, zJ): rectangle from corner I to J
    ops.patch('rect', MAT_STEEL, 2, 8, b - t, -b, b, b)          # wall at +y
    ops.patch('rect', MAT_STEEL, 2, 8, -b, -b, -b + t, b)        # wall at -y
    ops.patch('rect', MAT_STEEL, 8, 2, -b + t, -b, b - t, -b + t)  # wall at -z
    ops.patch('rect', MAT_STEEL, 8, 2, -b + t, b - t, b - t, b)    # wall at +z

    # --- beam fiber section: IPE270, depth along local y (= global Z, see transf)
    ops.section('Fiber', SEC_BM, '-GJ', G_S * J_IPE)
    d2 = BM_D / 2.0
    ops.patch('rect', MAT_STEEL, 2, 8, d2 - BM_TF, -BM_BF / 2, d2, BM_BF / 2)    # top flange
    ops.patch('rect', MAT_STEEL, 2, 8, -d2, -BM_BF / 2, -d2 + BM_TF, BM_BF / 2)  # bottom flange
    ops.patch('rect', MAT_STEEL, 8, 1, -d2 + BM_TF, -BM_TW / 2, d2 - BM_TF, BM_TW / 2)  # web

    # --- geometric transformations: geomTransf(type, tag, *vecxz)
    # local y = vecxz x local x. Columns: P-Delta (gravity on the sway).
    ops.geomTransf('PDelta', TR_COL, 1.0, 0.0, 0.0)
    # X beams: vecxz = -Y  -> local y = +Z (section depth vertical, Wy<0 = downward)
    ops.geomTransf('Linear', TR_BX, 0.0, -1.0, 0.0)
    # Y beams: vecxz = +X  -> local y = +Z
    ops.geomTransf('Linear', TR_BY, 1.0, 0.0, 0.0)

    # --- integration rules: beamIntegration('Lobatto', tag, secTag, nIP)
    ops.beamIntegration('Lobatto', INT_COL, SEC_COL, 3)
    ops.beamIntegration('Lobatto', INT_BM, SEC_BM, 3)

    # --- columns (forceBeamColumn = force-based fiber element)
    cols = []
    for st in range(1, NF + 1):
        for iy in range(1, NY + 1):
            for ix in range(1, NX + 1):
                if removed is not None and (st, ix, iy) == tuple(removed):
                    continue                          # this column is "removed"
                # element('forceBeamColumn', tag, iNode, jNode, transfTag, integTag)
                ops.element('forceBeamColumn', col_tag(st, ix, iy),
                            node_tag(st - 1, ix, iy), node_tag(st, ix, iy),
                            TR_COL, INT_COL)
                cols.append(col_tag(st, ix, iy))

    # --- beams with their uniform gravity load (two-way slab split 1/4 per edge)
    beams = []   # (tag, w [kN/m], iNode, jNode)
    for fl in range(1, NF + 1):
        for iy in range(1, NY + 1):                   # X-direction beams
            n_pan = (iy > 1) + (iy < NY)              # panels on each side
            for ix in range(1, NX):
                tag = beam_tag(fl, ix, iy, 'X')
                ni, nj = node_tag(fl, ix, iy), node_tag(fl, ix + 1, iy)
                ops.element('forceBeamColumn', tag, ni, nj, TR_BX, INT_BM)
                beams.append((tag, n_pan * Q * BAY / 4.0, ni, nj))
        for ix in range(1, NX + 1):                   # Y-direction beams
            n_pan = (ix > 1) + (ix < NX)
            for iy in range(1, NY):
                tag = beam_tag(fl, ix, iy, 'Y')
                ni, nj = node_tag(fl, ix, iy), node_tag(fl, ix, iy + 1)
                ops.element('forceBeamColumn', tag, ni, nj, TR_BY, INT_BM)
                beams.append((tag, n_pan * Q * BAY / 4.0, ni, nj))
    return cols, beams


# %% Robust solution helpers ---------------------------------------------------
ALGOS = [('Newton',), ('ModifiedNewton', '-initial'), ('KrylovNewton',)]


def _try_algorithms(run_one):
    """Try Newton, then ModifiedNewton -initial, then KrylovNewton."""
    for alg in ALGOS:
        ops.algorithm(*alg)                           # algorithm(type, *options)
        if run_one() == 0:                            # 0 = converged
            ops.algorithm('Newton')
            return True
    ops.algorithm('Newton')
    return False


def static_step(dlam, depth=0):
    def one():
        ops.integrator('LoadControl', dlam)           # load-factor increment dlam
        return ops.analyze(1)                         # one static step
    if _try_algorithms(one):
        return True
    if depth >= MAX_HALVING:
        return False
    return static_step(dlam / 2, depth + 1) and static_step(dlam / 2, depth + 1)


def dyn_step(dt, depth=0):
    if _try_algorithms(lambda: ops.analyze(1, dt)):   # analyze(nSteps, dt)
        return True
    if depth >= MAX_HALVING:
        return False
    return dyn_step(dt / 2, depth + 1) and dyn_step(dt / 2, depth + 1)


# %% Ground motion sequence ----------------------------------------------------
ec = np.loadtxt(os.path.join(REC_DIR, 'elCentro_NS_dt0.02_g.txt'))          # dt 0.02, g
nr = np.loadtxt(os.path.join(REC_DIR, 'Northridge_LOS270_dt0.01_g.txt'))    # dt 0.01, g
main = SF_MAIN * ec[:int(round(MAIN_DUR / DT))]
t_nr = np.arange(len(nr)) * 0.01
t_after = np.arange(0.0, AFTER_DUR, DT)
after = SF_AFTER * np.interp(t_after, t_nr, nr)          # resample 0.01 -> 0.02 s
ag = np.concatenate([main, np.zeros(int(round(GAP / DT))), after,
                     np.zeros(int(round(TAIL / DT)))])
N_STEPS = len(ag)
t_ag = np.arange(N_STEPS) * DT
T_MAIN_END = len(main) * DT                  # end of mainshock
T_GAP_END = T_MAIN_END + GAP                 # end of the free-vibration gap
print(f'Sequence: {N_STEPS} steps, {t_ag[-1] + DT:.1f} s, PGA main = '
      f'{abs(main).max():.2f} g, PGA after = {abs(after).max():.2f} g')

LINES = [(ix, iy) for iy in range(1, NY + 1) for ix in range(1, NX + 1)]


# %% 2.-4. Loop over the removal cases ------------------------------------------
results = {}
for label, removed in CASES:
    t0 = time.time()
    cols, beams = build_model(removed)

    # %% 2. Gravity (D + 0.25L) - load control in 10 steps
    ops.timeSeries('Linear', 1)                       # load factor = pseudo-time
    ops.pattern('Plain', 1, 1)                        # pattern(Plain, tag, tsTag)
    for tag, w, ni, nj in beams:
        # eleLoad: uniform load per length in local (y, z); local y = up -> -w
        ops.eleLoad('-ele', tag, '-type', '-beamUniform', -w, 0.0)
    ops.constraints('Transformation')                 # needed for rigidDiaphragm
    ops.numberer('RCM')                               # bandwidth-reducing numbering
    ops.system('BandGeneral')                         # banded solver
    ops.test('NormDispIncr', 1.0e-8, 50)              # convergence test (tol, maxIter)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.1)
    ops.analysis('Static')
    ok = all(static_step(0.1) for _ in range(10))
    if not ok:
        print(f'[{label}] gravity analysis failed - case skipped')
        results[label] = None
        continue
    ops.reactions()                                   # update nodal reactions
    W_base = sum(ops.nodeReaction(node_tag(0, ix, iy), 3) for ix, iy in LINES)
    dz_rem = 0.0
    if removed is not None:                           # sag of the joint above the gap
        st, ix, iy = removed
        dz_rem = -ops.nodeDisp(node_tag(st, ix, iy), 3)   # sag, positive downward

    # eigen(nModes) returns eigenvalues omega^2 (uses the Transformation handler)
    lam = ops.eigen(3)
    T = [2 * np.pi / np.sqrt(v) for v in lam]
    ops.loadConst('-time', 0.0)                       # freeze gravity, reset time to 0
    ops.wipeAnalysis()                                # clear the static analysis objects

    # %% 3. Damping and ground motion
    w1 = 2 * np.pi / T[0]
    w2 = 3.0 * w1                                     # ~ 2nd translational mode
    a0 = ZETA * 2 * w1 * w2 / (w1 + w2)
    a1 = ZETA * 2 / (w1 + w2)
    # rayleigh(alphaM, betaK, betaKinit, betaKcomm): mass + committed stiffness
    ops.rayleigh(a0, 0.0, 0.0, a1)
    # timeSeries('Path', tag, '-dt', dt, '-values', *a, '-factor', f): ag(t) in m/s^2
    ops.timeSeries('Path', 2, '-dt', DT, '-values', *ag.tolist(), '-factor', G)
    # pattern('UniformExcitation', tag, dir, '-accel', tsTag): base motion along X
    ops.pattern('UniformExcitation', 2, 1, '-accel', 2)

    # %% 4. Multi-shock time history
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-6, 50)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)              # average acceleration
    ops.analysis('Transient')

    tags_all = [[node_tag(fl, ix, iy) for ix, iy in LINES] for fl in range(NF + 1)]
    ux = np.zeros((N_STEPS + 1, NF + 1, len(LINES)))  # X displacement history
    tt = np.zeros(N_STEPS + 1)
    n_done, status = 0, 'complete'
    for k in range(N_STEPS):
        if not dyn_step(DT):
            status = f'no convergence at t={ops.getTime():.2f}s'
            break
        n_done = k + 1
        tt[k + 1] = ops.getTime()                     # current analysis time
        for fl in range(1, NF + 1):
            # nodeDisp(node, dof): committed displacement, dof 1 = UX
            ux[k + 1, fl] = [ops.nodeDisp(n, 1) for n in tags_all[fl]]
        if time.time() - t0 > CASE_BUDGET:            # never let a case run away
            status = f'time budget hit at t={ops.getTime():.2f}s'
            break

    tt, ux = tt[:n_done + 1], ux[:n_done + 1]
    drift = (ux[:, 1:, :] - ux[:, :-1, :]) / H        # [time, story, line]
    peak = np.abs(drift).max(axis=(0, 2))             # per story, max over plan
    win_end = tt >= tt[-1] - 1.0                      # last second of the run
    win_gap = (tt >= T_GAP_END - 1.0) & (tt <= T_GAP_END)
    res_end = np.abs(drift[win_end].mean(axis=0)).max()
    res_gap = np.abs(drift[win_gap].mean(axis=0)).max() if win_gap.any() else np.nan
    i_c = LINES.index((2, 2))
    results[label] = dict(T=T, W=W_base, dz=dz_rem, t=tt, roof=ux[:, NF, i_c],
                          peak_story=peak, peak=peak.max(), res=res_end,
                          res_gap=res_gap, status=status)
    print(f'[{label}] T1={T[0]:.3f}s  sum Rz={W_base:.0f} kN  peak drift='
          f'{100 * peak.max():.2f}%  residual={100 * res_end:.3f}%  ({status}, '
          f'{time.time() - t0:.1f} s)')


# %% 5. Figures ---------------------------------------------------------------
# Figure 1 - model sketch from node coordinates, removed columns in red
fig = plt.figure(figsize=(7, 4.2))
ax = fig.add_subplot(111, projection='3d')
for fl in range(1, NF + 1):
    for iy in range(1, NY + 1):
        for ix in range(1, NX):
            p, q = xyz(fl, ix, iy), xyz(fl, ix + 1, iy)
            ax.plot(*zip(p, q), color='tab:blue', lw=1.2)
    for ix in range(1, NX + 1):
        for iy in range(1, NY):
            p, q = xyz(fl, ix, iy), xyz(fl, ix, iy + 1)
            ax.plot(*zip(p, q), color='tab:blue', lw=1.2)
removed_set = {r for _, r in CASES if r is not None}
for st in range(1, NF + 1):
    for ix, iy in LINES:
        p, q = xyz(st - 1, ix, iy), xyz(st, ix, iy)
        if (st, ix, iy) in removed_set:
            ax.plot(*zip(p, q), color='red', lw=2.5, ls='--')
            ax.text(p[0] - 1.2, p[1] - 1.2, 0.5 * (p[2] + q[2]), f'S{st}', color='red')
        else:
            ax.plot(*zip(p, q), color='k', lw=1.6)
for fl, ix, iy in [(0, 1, 1), (NF, 1, 1), (NF, 2, 2), (NF, 3, 3)]:
    ax.text(*xyz(fl, ix, iy), f' {node_tag(fl, ix, iy)}', fontsize=7)
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_zlabel('Z [m]')
ax.set_title('3D frame: red dashed = column removed in case S1/S2/S3')
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'ch09_model.png'), dpi=110)
plt.close(fig)

# Figure 2 - ground-motion sequence and roof displacement
fig, (a1_, a2_) = plt.subplots(2, 1, figsize=(7, 4.2), sharex=True)
a1_.plot(t_ag, ag, 'k', lw=0.6)
a1_.axvspan(T_MAIN_END, T_GAP_END, color='0.9')
a1_.set_ylabel('a_g [g]')
a1_.set_title('Mainshock - gap - aftershock (one Path series)')
for label, _ in CASES:
    r = results.get(label)
    if r is not None:
        a2_.plot(r['t'], 1000 * r['roof'], lw=0.8, label=label)
a2_.axvspan(T_MAIN_END, T_GAP_END, color='0.9')
a2_.set_xlabel('Time [s]')
a2_.set_ylabel('Roof u_X [mm]')
a2_.legend(fontsize=7, ncol=4)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'ch09_sequence.png'), dpi=110)
plt.close(fig)

# Figure 3 - peak and residual drift per case
labels = [c for c, _ in CASES]
pk = [100 * results[c]['peak'] if results.get(c) else np.nan for c in labels]
rs_end = [100 * results[c]['res'] if results.get(c) else np.nan for c in labels]
rs_gap = [100 * results[c]['res_gap'] if results.get(c) else np.nan for c in labels]
x = np.arange(len(labels))
fig, (b1, b2) = plt.subplots(1, 2, figsize=(7, 4.2))
b1.bar(x, pk, color='tab:blue')
b1.set_xticks(x)
b1.set_xticklabels(labels)
b1.set_ylabel('Peak story drift [%]')
b2.bar(x - 0.2, rs_gap, 0.4, label='after mainshock')
b2.bar(x + 0.2, rs_end, 0.4, label='after aftershock')
b2.set_xticks(x)
b2.set_xticklabels(labels)
b2.set_ylabel('Residual story drift [%]')
b2.legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'ch09_drifts.png'), dpi=110)
plt.close(fig)

# Figure 4 - peak drift profile over the height
fig, ax = plt.subplots(figsize=(7, 4.2))
for c in labels:
    r = results.get(c)
    if r is not None:
        ax.plot(100 * r['peak_story'], np.arange(1, NF + 1), 'o-', label=c)
ax.set_yticks(range(1, NF + 1))
ax.set_xlabel('Peak story drift (max over plan) [%]')
ax.set_ylabel('Story')
ax.grid(alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'ch09_profile.png'), dpi=110)
plt.close(fig)


# %% 6. Key results -------------------------------------------------------------
def _get(label, key, scale=1.0, fmt='{:.3f}'):
    r = results.get(label)
    if r is None or r.get(key) is None or not np.isfinite(r[key]):
        return 'nan'
    return fmt.format(scale * r[key])


for c in ('Intact', 'S1'):
    r = results.get(c)
    print(f'RESULT T1_{c}_s', '{:.3f}'.format(r['T'][0]) if r else 'nan')
print('RESULT dz_grav_S1_mm', _get('S1', 'dz', 1000.0, '{:.1f}'))
for c in labels:
    print(f'RESULT peak_drift_{c}_pct', _get(c, 'peak', 100.0, '{:.2f}'))
for c in labels:
    print(f'RESULT residual_drift_{c}_pct', _get(c, 'res', 100.0, '{:.3f}'))
for c in labels:
    if results.get(c) and results[c]['status'] != 'complete':
        print(f'WARNING {c}: {results[c]["status"]}')
print(f'Total run time {time.time() - T_START:.1f} s')
