# -*- coding: utf-8 -*-
"""
Chapter 10 - Progressive collapse: alternate path (GSA 2016 / UFC 4-023-03 idea)
================================================================================
Units: kN, m, s  (stress kN/m^2, mass t = kN*s^2/m).  g = 9.81 m/s^2

What the script does
  * 2D steel moment frame, 3 stories x 4 bays, fiber sections (Steel02),
    beams with Corotational transformation (large deflections, catenary action).
  * Load combination for alternate-path analysis: 1.2D + 0.5L on all beams.
  * A) Nonlinear DYNAMIC column removal (GSA nonlinear dynamic procedure):
       gravity on the intact frame -> read the forces of the first-story centre
       column -> remove the column and replace it by those forces -> ramp the
       forces to zero in t_r = T_vert/10 -> follow the vertical displacement.
  * B) Nonlinear STATIC "push-down" of the damaged frame: loads on the two bays
       above the removed column are increased from 0 to 2.0 x (1.2D+0.5L)
       (dynamic amplification factor 2), the rest of the frame stays at 1.0.
  * Compares A and B and checks the beam chord rotation against a limit.

Tag convention (chapter 03, 2D version)
  node   = 100*floor + axis               (floor 0 = base, e.g. 3 = base of axis 3)
  column = 1000 + 100*story + axis        (1xxx, 1103 = story 1, axis 3)
  beam   = 2000 + 100*floor + axis        (2xxx, axis = left end of the span)
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

# %% 0. Parameters ------------------------------------------------------------
G = 9.81
N_ST = 3                      # stories (= floors above the base)
N_AX = 5                      # column axes -> 4 bays
L_BAY = 6.0                   # bay length [m]
H = 3.5                       # story height [m]
TRIB = 6.0                    # frame spacing = tributary width of a beam [m]
D_AREA, L_AREA = 5.0, 2.5     # dead (incl. self weight) and live load [kN/m^2]
W_D, W_L = D_AREA * TRIB, L_AREA * TRIB
W_AP = 1.2 * W_D + 0.5 * W_L  # alternate-path combination 1.2D + 0.5L [kN/m]

E_S, FY, B_H = 2.0e8, 355.0e3, 0.01      # Steel S355 [kN/m^2], hardening ratio
# Columns HEB300 (strong axis in the frame plane), beams IPE550 [m]
COL = dict(d=0.300, bf=0.300, tf=0.019, tw=0.011)
BM = dict(d=0.550, bf=0.210, tf=0.0172, tw=0.0111)

REM_AXIS = 3                  # removed column: story 1, axis 3 (centre column)
COL_REM = 1000 + 100 * 1 + REM_AXIS      # element tag 1103
NODE_REM = 100 * 1 + REM_AXIS            # node above the removed column: 103
AFFECTED_AXES = (REM_AXIS - 1, REM_AXIS)  # beams that start at axes 2 and 3

ZETA = 0.05                   # Rayleigh damping ratio
T_HOLD = 0.1                  # quiet time before the removal starts [s]
DAF = 2.0                     # classic dynamic amplification factor (static case)
MAX_HALVING = 4
TIME_BUDGET = 120.0           # wall-clock limit for the dynamic loop [s]

MAT_STEEL = 1
SEC_COL, SEC_BM = 1, 2
TR_COL, TR_BM = 1, 2
INT_COL, INT_BM = 1, 2


def node_tag(floor, axis):
    return 100 * floor + axis


def col_tag(story, axis):
    return 1000 + 100 * story + axis


def beam_tag(floor, axis):
    return 2000 + 100 * floor + axis


def i_section_props(s):
    """Return A, I, Zpl of a doubly symmetric I section (strong axis)."""
    d, bf, tf, tw = s['d'], s['bf'], s['tf'], s['tw']
    hw = d - 2 * tf
    A = 2 * bf * tf + hw * tw
    I = bf * d**3 / 12 - (bf - tw) * hw**3 / 12
    Z = bf * tf * (d - tf) + tw * hw**2 / 4
    return A, I, Z


A_BM, I_BM, Z_BM = i_section_props(BM)
MP_BM = Z_BM * FY                          # plastic moment of the beam [kNm]


# %% 1. Model ------------------------------------------------------------------
def i_fiber_section(tag, s):
    """2D fiber I-section: depth along local y. patch rect corners (y, z)."""
    d2, bf, tf, tw = s['d'] / 2, s['bf'], s['tf'], s['tw']
    ops.section('Fiber', tag)                                   # 2D fiber section
    # patch('rect', mat, nFibY, nFibZ, yI, zI, yJ, zJ)
    ops.patch('rect', MAT_STEEL, 3, 1, d2 - tf, -bf / 2, d2, bf / 2)        # top flange
    ops.patch('rect', MAT_STEEL, 3, 1, -d2, -bf / 2, -d2 + tf, bf / 2)      # bottom flange
    ops.patch('rect', MAT_STEEL, 12, 1, -d2 + tf, -tw / 2, d2 - tf, tw / 2)  # web


MASS = {}      # node -> mass [t], kept in Python for the modal post-processing


def build_model(with_removed_column=True):
    """Build the 2D frame. with_removed_column=False omits column 1103."""
    ops.wipe()                                          # clear the domain
    ops.model('basic', '-ndm', 2, '-ndf', 3)            # 2D, 3 DOF per node
    MASS.clear()
    for fl in range(N_ST + 1):
        for ax in range(1, N_AX + 1):
            tag = node_tag(fl, ax)
            ops.node(tag, (ax - 1) * L_BAY, fl * H)     # node(tag, X, Y)
            if fl == 0:
                ops.fix(tag, 1, 1, 1)                   # fixed base
            else:
                trib = L_BAY / 2 if ax in (1, N_AX) else L_BAY
                m = W_AP * trib / G                     # mass from the applied load
                ops.mass(tag, m, m, 0.0)                # horizontal AND vertical mass
                MASS[tag] = m

    # Steel02(tag, Fy, E, b, R0, cR1, cR2)
    ops.uniaxialMaterial('Steel02', MAT_STEEL, FY, E_S, B_H, 18.0, 0.925, 0.15)
    i_fiber_section(SEC_COL, COL)
    i_fiber_section(SEC_BM, BM)

    ops.geomTransf('PDelta', TR_COL)                    # columns: P-Delta
    ops.geomTransf('Corotational', TR_BM)               # beams: large rotations/catenary
    ops.beamIntegration('Lobatto', INT_COL, SEC_COL, 4)  # 4 points per column
    ops.beamIntegration('Lobatto', INT_BM, SEC_BM, 5)    # 5 points per beam (hinges)

    for st in range(1, N_ST + 1):
        for ax in range(1, N_AX + 1):
            if not with_removed_column and col_tag(st, ax) == COL_REM:
                continue
            # element('forceBeamColumn', tag, iNode, jNode, transfTag, integTag)
            ops.element('forceBeamColumn', col_tag(st, ax), node_tag(st - 1, ax),
                        node_tag(st, ax), TR_COL, INT_COL)
    for fl in range(1, N_ST + 1):
        for ax in range(1, N_AX):
            ops.element('forceBeamColumn', beam_tag(fl, ax), node_tag(fl, ax),
                        node_tag(fl, ax + 1), TR_BM, INT_BM)


def beams(affected=None):
    """All beam tags; affected=True -> only the bays over the removed column,
    affected=False -> all the others."""
    out = []
    for fl in range(1, N_ST + 1):
        for ax in range(1, N_AX):
            is_aff = ax in AFFECTED_AXES
            if affected is None or affected == is_aff:
                out.append(beam_tag(fl, ax))
    return out


def add_beam_load(pattern_tag, ts_tag, tags, w):
    ops.pattern('Plain', pattern_tag, ts_tag)            # pattern(Plain, tag, tsTag)
    for t in tags:
        # local y of a left-to-right beam is global +Y -> -w acts downward
        ops.eleLoad('-ele', t, '-type', '-beamUniform', -w)


# %% Robust solution helpers ---------------------------------------------------
ALGOS = [('Newton',), ('ModifiedNewton', '-initial'), ('KrylovNewton',)]


def _try_algorithms(run_one):
    for alg in ALGOS:
        ops.algorithm(*alg)
        if run_one() == 0:
            ops.algorithm('Newton')
            return True
    ops.algorithm('Newton')
    return False


def static_step(dlam, depth=0):
    def one():
        ops.integrator('LoadControl', dlam)              # load-factor increment
        return ops.analyze(1)
    if _try_algorithms(one):
        return True
    if depth >= MAX_HALVING:
        return False
    return static_step(dlam / 2, depth + 1) and static_step(dlam / 2, depth + 1)


def dyn_step(dt, depth=0):
    if _try_algorithms(lambda: ops.analyze(1, dt)):
        return True
    if depth >= MAX_HALVING:
        return False
    return dyn_step(dt / 2, depth + 1) and dyn_step(dt / 2, depth + 1)


def static_setup(tol=1.0e-8):
    ops.constraints('Plain')                             # only single-point fixes
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', tol, 50)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.1)
    ops.analysis('Static')


# %% 2. Gravity on the intact frame (dynamic case) -------------------------------
build_model(with_removed_column=True)
ops.timeSeries('Linear', 1)
add_beam_load(1, 1, beams(), W_AP)                       # 1.2D + 0.5L everywhere
static_setup()
if not all(static_step(0.1) for _ in range(10)):
    raise SystemExit('Gravity analysis of the intact frame failed')
delta_grav0 = -ops.nodeDisp(NODE_REM, 2)                 # sag before removal

# forces(ele) = element resisting forces in global axes [Fx_i, Fy_i, M_i, Fx_j, Fy_j, M_j]
f_col = ops.eleResponse(COL_REM, 'forces')
F_REP = [-f_col[3], -f_col[4], -f_col[5]]                # what the column does to node 103
P_REMOVED = F_REP[1]                                     # upward support force [kN]
ops.loadConst('-time', 0.0)                              # gravity stays, time -> 0

# %% 3. Column removal and vertical period ---------------------------------------
ops.remove('element', COL_REM)                           # take column 1103 out
lam = ops.eigen(6)                                       # 6 modes of the damaged frame
best, frac_best = 0, -1.0
for k in range(1, len(lam) + 1):
    num, den = 0.0, 0.0
    for n, m in MASS.items():
        fx = ops.nodeEigenvector(n, k, 1)                # mode k, DOF 1 (UX)
        fy = ops.nodeEigenvector(n, k, 2)                # mode k, DOF 2 (UY)
        num += m * fy
        den += m * (fx**2 + fy**2)
    frac = num**2 / den / sum(MASS.values())             # vertical effective-mass ratio
    if frac > frac_best:
        best, frac_best = k, frac
T_LAT = 2 * np.pi / np.sqrt(lam[0])
T_V = 2 * np.pi / np.sqrt(lam[best - 1])
T_R = T_V / 10.0                                         # GSA: removal in <= T_v/10
print(f'T1 (lateral) = {T_LAT:.3f} s, vertical mode {best}: T_v = {T_V:.3f} s '
      f'(eff. vertical mass {100 * frac_best:.0f}%), t_r = {T_R:.4f} s')
ops.wipeAnalysis()

# Replacement forces: factor 1 until T_HOLD, linear drop to 0 at T_HOLD + t_r
T_END = min(T_HOLD + T_R + 8.0 * T_V, 4.0)
DT = min(T_R / 10.0, T_V / 50.0)
n_steps = int(np.ceil(T_END / DT))
if n_steps > 3000:
    DT = T_END / 3000.0
    n_steps = 3000
# timeSeries('Path', tag, '-time', *t, '-values', *f): piecewise-linear factor
ops.timeSeries('Path', 2, '-time', 0.0, T_HOLD, T_HOLD + T_R, T_END + 1.0,
               '-values', 1.0, 1.0, 0.0, 0.0)
ops.pattern('Plain', 2, 2)
ops.load(NODE_REM, *F_REP)                               # load(node, Fx, Fy, Mz)

w_i, w_j = 2 * np.pi / T_V, 3 * 2 * np.pi / T_V           # damping fixed at T_v, T_v/3
a0 = ZETA * 2 * w_i * w_j / (w_i + w_j)
a1 = ZETA * 2 / (w_i + w_j)
ops.rayleigh(a0, 0.0, 0.0, a1)                           # mass + committed stiffness

# %% 4. Nonlinear dynamic analysis --------------------------------------------------
ops.constraints('Plain')
ops.numberer('RCM')
ops.system('BandGeneral')
ops.test('NormDispIncr', 1.0e-6, 50)
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)                     # average acceleration
ops.analysis('Transient')

ALL_NODES = [node_tag(fl, ax) for fl in range(N_ST + 1) for ax in range(1, N_AX + 1)]
t_hist, v_hist = [0.0], [ops.nodeDisp(NODE_REM, 2)]
disp_at_peak, v_min = None, 0.0
status, t_wall = 'complete', time.time()
for k in range(n_steps):
    if not dyn_step(DT):
        status = f'no convergence at t={ops.getTime():.3f}s'
        break
    v = ops.nodeDisp(NODE_REM, 2)                        # vertical displacement (UY)
    t_hist.append(ops.getTime())
    v_hist.append(v)
    if v < v_min:
        v_min = v
        disp_at_peak = {n: (ops.nodeDisp(n, 1), ops.nodeDisp(n, 2)) for n in ALL_NODES}
    if -v > 0.25 * L_BAY:
        status = 'collapse (sag > L/4)'
        break
    if time.time() - t_wall > TIME_BUDGET:
        status = 'time budget hit'
        break
t_hist, v_hist = np.array(t_hist), np.array(v_hist)
DELTA_DYN = -v_min                                       # peak sag [m]
print(f'Dynamic run: {len(t_hist) - 1} steps, dt={DT:.4f} s, {status}')

# %% 5. Nonlinear static push-down with DAF = 2 -----------------------------------
build_model(with_removed_column=False)
ops.timeSeries('Linear', 1)
add_beam_load(1, 1, beams(affected=False), W_AP)         # unaffected bays at 1.0
static_setup()
ok = all(static_step(0.1) for _ in range(10))
ops.loadConst('-time', 0.0)
ops.timeSeries('Linear', 2)
add_beam_load(2, 2, beams(affected=True), W_AP)          # affected bays: factor lambda
lam_curve, d_curve = [0.0], [-ops.nodeDisp(NODE_REM, 2)]
d_lam = 0.05
for _ in range(int(round(DAF / d_lam))):
    if not ok or not static_step(d_lam):
        ok = False
        break
    lam_curve.append(ops.getTime())                      # Linear series: factor = time
    d_curve.append(-ops.nodeDisp(NODE_REM, 2))
lam_curve, d_curve = np.array(lam_curve), np.array(d_curve)
LAM_MAX = lam_curve[-1]
print(f'Push-down reached lambda = {LAM_MAX:.2f} ({"ok" if ok else "stopped"})')


def d_at(lmb):
    return float(np.interp(lmb, lam_curve, d_curve)) if LAM_MAX >= lmb - 1e-9 else np.nan


DELTA_S1 = d_at(1.0)
DELTA_S2 = d_at(DAF)
DAF_DYN = DELTA_DYN / DELTA_S1 if DELTA_S1 > 0 else np.nan
if np.all(np.diff(d_curve) > 0) and d_curve[0] <= DELTA_DYN <= d_curve[-1]:
    LAM_EQ = float(np.interp(DELTA_DYN, d_curve, lam_curve))
else:
    LAM_EQ = np.nan

# %% 6. Acceptance check (chord rotation) -----------------------------------------
# Chord rotation of the beams next to the removed column: theta = sag / span.
# Limit (illustrative): theta_y + theta_pra, with
#   theta_y   = Zpl*Fy*L / (6 E I)            (beam yield rotation, ASCE 41 form)
#   theta_pra = 0.0284 - 0.0004*d[in]         (WUF connection, primary component,
#                                              UFC 4-023-03 (2009) Table 5-2)
THETA_Y = MP_BM * L_BAY / (6.0 * E_S * I_BM)
THETA_PRA = 0.0284 - 0.0004 * (BM['d'] / 0.0254)
THETA_LIM = THETA_Y + THETA_PRA
THETA_DYN = DELTA_DYN / L_BAY
THETA_S2 = DELTA_S2 / L_BAY
CHECK = 'PASS' if (status == 'complete' and THETA_DYN <= THETA_LIM) else 'FAIL'

# %% 7. Figures ----------------------------------------------------------------
# Figure 1 - model sketch
fig, ax = plt.subplots(figsize=(7, 4.2))
for st in range(1, N_ST + 1):
    for a in range(1, N_AX + 1):
        x = (a - 1) * L_BAY
        rem = col_tag(st, a) == COL_REM
        ax.plot([x, x], [(st - 1) * H, st * H], 'r--' if rem else 'k-', lw=2)
        ax.text(x + 0.1, (st - 0.5) * H, str(col_tag(st, a)), fontsize=6,
                color='r' if rem else '0.3')
for fl in range(1, N_ST + 1):
    for a in range(1, N_AX):
        x0 = (a - 1) * L_BAY
        aff = a in AFFECTED_AXES
        ax.plot([x0, x0 + L_BAY], [fl * H, fl * H], color='tab:orange' if aff else 'tab:blue', lw=2)
        for xa in np.linspace(x0 + 0.5, x0 + L_BAY - 0.5, 5):
            ax.annotate('', xy=(xa, fl * H + 0.05), xytext=(xa, fl * H + 0.6),
                        arrowprops=dict(arrowstyle='->', color='0.5', lw=0.6))
for fl in range(N_ST + 1):
    for a in range(1, N_AX + 1):
        ax.plot((a - 1) * L_BAY, fl * H, 'ko', ms=3)
ax.text((REM_AXIS - 1) * L_BAY + 0.15, H + 0.1, str(NODE_REM), color='r', fontsize=8)
for x0 in range(N_AX):
    ax.plot(x0 * L_BAY, 0, 'ks', ms=6)
ax.set_title(f'Removed column {COL_REM} (red dashed); orange = bays loaded x{DAF:.0f} '
             f'in the static case', fontsize=9)
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_aspect('equal')
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'ch10_model.png'), dpi=110)
plt.close(fig)

# Figure 2 - vertical displacement history
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(t_hist, -1000 * v_hist, 'k', lw=1.2, label='nonlinear dynamic (1.2D+0.5L)')
if np.isfinite(DELTA_S1):
    ax.axhline(1000 * DELTA_S1, color='tab:blue', ls='--', label='static, factor 1.0')
if np.isfinite(DELTA_S2):
    ax.axhline(1000 * DELTA_S2, color='tab:red', ls='--', label=f'static, DAF = {DAF:.0f}')
ax.axvspan(T_HOLD, T_HOLD + T_R, color='0.85', label='removal ramp t_r')
ax.set_xlabel('Time [s]')
ax.set_ylabel(f'Sag at node {NODE_REM} [mm]')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'ch10_history.png'), dpi=110)
plt.close(fig)

# Figure 3 - push-down curve
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(1000 * d_curve, lam_curve, 'b.-', label='static push-down (affected bays)')
ax.axvline(1000 * DELTA_DYN, color='k', ls='--', label='dynamic peak sag')
if np.isfinite(LAM_EQ):
    ax.plot(1000 * DELTA_DYN, LAM_EQ, 'ro', label=f'equivalent static factor {LAM_EQ:.2f}')
ax.set_xlabel(f'Sag at node {NODE_REM} [mm]')
ax.set_ylabel('Load factor on affected bays')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'ch10_pushdown.png'), dpi=110)
plt.close(fig)

# Figure 4 - deformed shape at the dynamic peak
if disp_at_peak is not None:
    sf = 0.1 * L_BAY / max(DELTA_DYN, 1e-6)              # scale: peak sag drawn as L/10
    fig, ax = plt.subplots(figsize=(7, 4.2))

    def pos(n, dsp):
        fl, a = divmod(n, 100)
        x, y = (a - 1) * L_BAY, fl * H
        if dsp:
            x += sf * disp_at_peak[n][0]
            y += sf * disp_at_peak[n][1]
        return x, y
    for dsp, sty in [(False, dict(color='0.8', lw=1)), (True, dict(color='k', lw=1.8))]:
        for st in range(1, N_ST + 1):
            for a in range(1, N_AX + 1):
                if col_tag(st, a) == COL_REM:
                    continue
                p, q = pos(node_tag(st - 1, a), dsp), pos(node_tag(st, a), dsp)
                ax.plot([p[0], q[0]], [p[1], q[1]], **sty)
        for fl in range(1, N_ST + 1):
            for a in range(1, N_AX):
                p, q = pos(node_tag(fl, a), dsp), pos(node_tag(fl, a + 1), dsp)
                ax.plot([p[0], q[0]], [p[1], q[1]], **sty)
    ax.set_title(f'Deformed shape at peak sag (x{sf:.0f}); straight chords between joints',
                 fontsize=9)
    ax.set_aspect('equal')
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'ch10_deformed.png'), dpi=110)
    plt.close(fig)


# %% 8. Key results ---------------------------------------------------------------
def f(v, fmt):
    return fmt.format(v) if np.isfinite(v) else 'nan'


print(f'Hand check: tributary column load 3*w*L = {N_ST * W_AP * L_BAY:.0f} kN; '
      f'elastic ramp DLF = 1 + sin(pi*0.1)/(pi*0.1) = '
      f'{1 + np.sin(np.pi * 0.1) / (np.pi * 0.1):.3f}')
print('RESULT P_removed_kN', f(P_REMOVED, '{:.0f}'))
print('RESULT T_vert_s', f(T_V, '{:.3f}'))
print('RESULT t_removal_s', f(T_R, '{:.4f}'))
print('RESULT delta_dyn_mm', f(1000 * DELTA_DYN, '{:.1f}'))
print('RESULT delta_static1_mm', f(1000 * DELTA_S1, '{:.1f}'))
print('RESULT delta_staticDAF2_mm', f(1000 * DELTA_S2, '{:.1f}'))
print('RESULT DAF_dynamic', f(DAF_DYN, '{:.2f}'))
print('RESULT lambda_equiv', f(LAM_EQ, '{:.2f}'))
print('RESULT theta_dyn_rad', f(THETA_DYN, '{:.4f}'))
print('RESULT theta_limit_rad', f(THETA_LIM, '{:.4f}'))
print('RESULT AP_check', CHECK)
if status != 'complete':
    print('WARNING dynamic run:', status)
print(f'Total run time {time.time() - T_START:.1f} s')
