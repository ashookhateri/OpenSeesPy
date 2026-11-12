# -*- coding: utf-8 -*-
"""
Chapter 07 - RC moment frame with fiber sections
=================================================
A 4-story x 3-bay reinforced-concrete moment-resisting frame (2D).
Columns and beams are forceBeamColumn elements with fiber sections
(Concrete02 cover + confined core, Steel02 rebars).

Analyses:  gravity  ->  modal (eigen)  ->  pushover  ->  nonlinear time history
           (El Centro 1940 NS)  ->  peak story-drift profile.

UNITS:  force = kN,  length = m,  time = s
        stress = kN/m^2 (= kPa),  mass = kN*s^2/m (= tonne),  g = 9.81 m/s^2

TAGS (convention of chapter 03):
  node     = 100*floor + axis          (floor 0 = base, axis 1..4 from left)
  column   = 1000 + 100*story + axis   (= top node + 1000; story s: floor s-1 -> s)
  beam     = 2000 + 100*floor + bay    (= left node + 2000; bay b: axis b -> b+1)
  (3xxx braces and 4xxx dampers are not used here; dampers come in chapter 08)
"""
import os
import time
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')                       # no screen needed, only PNG files
import matplotlib.pyplot as plt
import openseespy.opensees as ops

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, 'figs')
os.makedirs(FIG, exist_ok=True)
REC = os.path.join(HERE, '..', 'records', 'elCentro_NS_dt0.02_g.txt')
T_START = time.time()

# %% 0. Parameters ------------------------------------------------------------
g = 9.81                                    # m/s^2
N_STORY, N_BAY = 4, 3
BAY = 5.0                                   # bay width (m)
H_STORY = [3.5, 3.2, 3.2, 3.2]              # story heights, bottom to top (m)
Y = [0.0] + list(np.cumsum(H_STORY))        # floor elevations (m)
H_TOT = Y[-1]

# gravity load on the beams (seismic-weight load: dead + part of live), kN/m
W_FLOOR, W_ROOF = 40.0, 30.0

# concrete (kPa) - Concrete02
fc = 30000.0                                # f'c = 30 MPa
fcc = 1.3 * fc                              # confined core strength (about +30 %)
ft = 0.10 * fc                              # tensile strength
Ec = 2.0 * fc / 0.002                       # initial modulus used by Concrete02
Ets = 0.10 * Ec                             # tension-softening slope
# steel (kPa) - Steel02
fy, Es, b_h = 400000.0, 2.0e8, 0.01         # fy = 400 MPa, Es = 200 GPa, 1 % hardening

A20 = math.pi * 0.020**2 / 4.0              # area of one 20 mm bar (m^2)
COL = dict(b=0.50, h=0.50, cover=0.05, n_top=4, n_bot=4, n_mid=2)   # 12 bars
BEAM = dict(b=0.35, h=0.50, cover=0.05, n_top=4, n_bot=3, n_mid=0)  # 4 top, 3 bottom

# tag tables (chapter 03): materials 1-9 steel, 10-19 concrete; sections 1xx columns,
# 2xx beams; transformations 1 columns, 2 beams
MAT_STEEL, MAT_COVER, MAT_CORE = 1, 10, 11
SEC_COL, SEC_BEAM = 101, 201
TR_COL, TR_BEAM = 1, 2                      # PDelta for columns, Linear for beams
INT_COL, INT_BEAM = 1, 2                    # beamIntegration tags
N_IP = 5                                    # Lobatto integration points per element

# dynamic analysis
DT_REC, DT_ANA, T_END = 0.02, 0.01, 25.0    # record step, analysis step, duration (s)
SF = 1.0                                    # scale factor on the record
XI = 0.05                                   # Rayleigh damping ratio (modes 1 and 3)
TIME_LIMIT = 160.0                          # wall-clock guard (s)


def ntag(floor, axis):
    """Node tag = 100*floor + axis."""
    return 100 * floor + axis


def col_tag(story, axis):
    return 1000 + 100 * story + axis


def beam_tag(floor, bay):
    return 2000 + 100 * floor + bay


# %% 1. Model -----------------------------------------------------------------
def rc_fiber_section(sec_tag, b, h, cover, n_top, n_bot, n_mid):
    """Rectangular RC fiber section. Local y = depth (bending direction), z = width.
    cover = distance from the face to the bar centre line (= edge of the core)."""
    yc, zc = h / 2.0 - cover, b / 2.0 - cover
    # section('Fiber', tag): starts a fiber section; the patch/layer lines below fill it
    ops.section('Fiber', sec_tag)
    # patch('rect', mat, nDivY, nDivZ, yI, zI, yJ, zJ): rectangle of fibers between
    # corner I (bottom-left) and J (top-right). In 2D only the y division matters.
    ops.patch('rect', MAT_CORE, 10, 1, -yc, -zc, yc, zc)            # confined core
    ops.patch('rect', MAT_COVER, 2, 1, yc, -b / 2, h / 2, b / 2)     # top cover
    ops.patch('rect', MAT_COVER, 2, 1, -h / 2, -b / 2, -yc, b / 2)   # bottom cover
    ops.patch('rect', MAT_COVER, 10, 1, -yc, -b / 2, yc, -zc)        # side cover
    ops.patch('rect', MAT_COVER, 10, 1, -yc, zc, yc, b / 2)          # side cover
    # layer('straight', mat, nBars, areaBar, yStart, zStart, yEnd, zEnd): a row of bars
    ops.layer('straight', MAT_STEEL, n_top, A20, yc, zc, yc, -zc)    # top bars
    ops.layer('straight', MAT_STEEL, n_bot, A20, -yc, zc, -yc, -zc)  # bottom bars
    if n_mid > 0:                                                    # web bars, both sides
        s = 2.0 * yc / (n_mid + 1)
        ops.layer('straight', MAT_STEEL, n_mid, A20, yc - s, zc, -yc + s, zc)
        ops.layer('straight', MAT_STEEL, n_mid, A20, yc - s, -zc, -yc + s, -zc)


def build_model():
    """Builds the whole frame. Returns dicts used for plotting and post-processing."""
    ops.wipe()                                      # clear any previous model
    ops.model('basic', '-ndm', 2, '-ndf', 3)        # 2D, 3 DOF per node (ux, uy, rz)

    coords, cols, beams = {}, {}, {}
    # nodes: node(tag, x, y)
    for f in range(N_STORY + 1):
        for a in range(1, N_BAY + 2):
            coords[ntag(f, a)] = ((a - 1) * BAY, Y[f])
            ops.node(ntag(f, a), *coords[ntag(f, a)])
    # supports: fix(tag, ux, uy, rz) with 1 = fixed
    for a in range(1, N_BAY + 2):
        ops.fix(ntag(0, a), 1, 1, 1)
    # rigid floor in its own plane: equalDOF(master, slave, dof) ties ux of each floor
    for f in range(1, N_STORY + 1):
        for a in range(2, N_BAY + 2):
            ops.equalDOF(ntag(f, 1), ntag(f, a), 1)

    # materials -----------------------------------------------------------
    # Concrete02(tag, fpc, epsc0, fpcu, epscu, lambda, ft, Ets): compression negative
    ops.uniaxialMaterial('Concrete02', MAT_COVER, -fc, -0.002, -0.2 * fc, -0.005, 0.1, ft, Ets)
    ops.uniaxialMaterial('Concrete02', MAT_CORE, -fcc, -0.003, -0.2 * fcc, -0.020, 0.1, ft, Ets)
    # Steel02(tag, Fy, E0, b, R0, cR1, cR2): Giuffre-Menegotto-Pinto steel
    ops.uniaxialMaterial('Steel02', MAT_STEEL, fy, Es, b_h, 18.0, 0.925, 0.15)

    rc_fiber_section(SEC_COL, **COL)
    rc_fiber_section(SEC_BEAM, **BEAM)
    # beamIntegration('Lobatto', tag, secTag, nIP): integration points along the member
    ops.beamIntegration('Lobatto', INT_COL, SEC_COL, N_IP)
    ops.beamIntegration('Lobatto', INT_BEAM, SEC_BEAM, N_IP)
    # geomTransf: PDelta adds the P-Delta effect of gravity loads on the columns
    ops.geomTransf('PDelta', TR_COL)
    ops.geomTransf('Linear', TR_BEAM)

    # elements: forceBeamColumn(tag, iNode, jNode, transfTag, integrationTag)
    for s in range(1, N_STORY + 1):
        for a in range(1, N_BAY + 2):
            t = col_tag(s, a)
            cols[t] = (ntag(s - 1, a), ntag(s, a))
            ops.element('forceBeamColumn', t, *cols[t], TR_COL, INT_COL)
    for f in range(1, N_STORY + 1):
        for b in range(1, N_BAY + 1):
            t = beam_tag(f, b)
            beams[t] = (ntag(f, b), ntag(f, b + 1))
            ops.element('forceBeamColumn', t, *beams[t], TR_BEAM, INT_BEAM)

    # lumped horizontal masses from the beam loads (tributary half-bays)
    floor_mass = {}
    for f in range(1, N_STORY + 1):
        w = W_ROOF if f == N_STORY else W_FLOOR
        mf = 0.0
        for a in range(1, N_BAY + 2):
            n_adj = (1 if a > 1 else 0) + (1 if a < N_BAY + 1 else 0)   # beams meeting here
            m = w * BAY / 2.0 * n_adj / g
            # mass(tag, mx, my, mrz): only horizontal mass is needed for the seismic analysis
            ops.mass(ntag(f, a), m, 0.0, 0.0)
            mf += m
        floor_mass[f] = mf
    return coords, cols, beams, floor_mass


def base_shear():
    """Base shear = -(sum of horizontal support reactions)."""
    ops.reactions()                                  # update the nodal reactions
    return -sum(ops.nodeReaction(ntag(0, a), 1) for a in range(1, N_BAY + 2))


# %% 2. Gravity ---------------------------------------------------------------
def gravity(beam_tags):
    # timeSeries('Linear', tag): load factor grows linearly with pseudo-time
    ops.timeSeries('Linear', 1)
    # pattern('Plain', tag, tsTag): an ordinary load pattern
    ops.pattern('Plain', 1, 1)
    for t in beam_tags:
        f = (t - 2000) // 100
        w = W_ROOF if f == N_STORY else W_FLOOR
        # eleLoad uniform load in local y (local y of a left-to-right beam points up)
        ops.eleLoad('-ele', t, '-type', '-beamUniform', -w)
    ops.constraints('Transformation')               # needed because of equalDOF
    ops.numberer('RCM')                             # reduce the bandwidth
    ops.system('BandGeneral')                       # banded solver
    ops.test('NormDispIncr', 1.0e-8, 50)            # convergence test: |du| < tol
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.1)              # 10 steps of 0.1 = full load
    ops.analysis('Static')
    ok = ops.analyze(10)
    if ok != 0:
        print('WARNING: gravity analysis did not converge')
    ops.loadConst('-time', 0.0)                     # freeze gravity, reset time to 0
    return ok


def modal(n=3):
    """Eigen analysis on the current (post-gravity) tangent stiffness."""
    try:
        lam = ops.eigen(n)                          # default solver: genBandArpack
    except Exception:
        lam = []
    if len(lam) < n or min(lam) <= 0.0:             # fall back to the dense solver
        lam = ops.eigen('-fullGenLapack', n)
    w = [math.sqrt(abs(x)) for x in lam]
    return w, [2 * math.pi / x for x in w]


# robust single step: Newton -> ModifiedNewton -initial -> KrylovNewton
ALGOS = [('Newton',), ('ModifiedNewton', '-initial'), ('KrylovNewton',)]
MAX_HALVING = 4


def try_algos(*analyze_args):
    for alg in ALGOS:
        ops.algorithm(*alg)
        if ops.analyze(*analyze_args) == 0:
            ops.algorithm('Newton')
            return True
    ops.algorithm('Newton')
    return False


def dyn_step(dt, depth=0):
    """One transient step of size dt; on failure split it in two halves (recursively)."""
    if try_algos(1, dt):
        return True
    if depth >= MAX_HALVING:
        return False
    return dyn_step(dt / 2.0, depth + 1) and dyn_step(dt / 2.0, depth + 1)


def push_step(node, du, depth=0):
    """One displacement-controlled step; halves du on failure."""
    # integrator('DisplacementControl', node, dof, du): impose du on (node, dof)
    ops.integrator('DisplacementControl', node, 1, du)
    if try_algos(1):
        return True
    if depth >= MAX_HALVING:
        return False
    return push_step(node, du / 2.0, depth + 1) and push_step(node, du / 2.0, depth + 1)


# =============================================================================
# RUN A: gravity + modal + pushover
# =============================================================================
coords, col_ele, beam_ele, floor_mass = build_model()
W_TOTAL = sum(floor_mass.values()) * g
gravity(sorted(beam_ele))
ops.reactions()
grav_sum = sum(ops.nodeReaction(ntag(0, a), 2) for a in range(1, N_BAY + 2))
N_col = {a: ops.nodeReaction(ntag(0, a), 2) for a in range(1, N_BAY + 2)}  # axial at base

# %% 3. Modal analysis --------------------------------------------------------
omega, T = modal(3)
print('Periods (s):', ['%.3f' % x for x in T])

# %% 4. Pushover --------------------------------------------------------------
ROOF = ntag(N_STORY, 1)
F_shape = np.array([floor_mass[f] * Y[f] for f in range(1, N_STORY + 1)])
F_shape = F_shape / F_shape.sum()                   # inverted-triangle, sums to 1 kN
ops.timeSeries('Linear', 2)
ops.pattern('Plain', 2, 2)
for f in range(1, N_STORY + 1):
    # load(node, Fx, Fy, Mz): reference lateral load, scaled by the load factor
    ops.load(ntag(f, 1), float(F_shape[f - 1]), 0.0, 0.0)

ops.wipeAnalysis()                                  # new analysis objects for pushover
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('BandGeneral')
ops.test('NormDispIncr', 1.0e-7, 100)
ops.algorithm('Newton')
ops.integrator('DisplacementControl', ROOF, 1, 0.002)
ops.analysis('Static')

DU, U_TARGET = 0.002, 0.03 * H_TOT                  # push to 3 % roof drift
push_u, push_V, push_drifts = [0.0], [0.0], []
while ops.nodeDisp(ROOF, 1) < U_TARGET:
    if time.time() - T_START > 0.4 * TIME_LIMIT:
        print('WARNING: pushover stopped by the time guard')
        break
    if not push_step(ROOF, DU):
        print('WARNING: pushover did not converge at roof u = %.3f m' % ops.nodeDisp(ROOF, 1))
        break
    push_u.append(ops.nodeDisp(ROOF, 1))
    push_V.append(base_shear())
    u_fl = [ops.nodeDisp(ntag(f, 1), 1) for f in range(N_STORY + 1)]
    push_drifts.append([(u_fl[s] - u_fl[s - 1]) / H_STORY[s - 1] for s in range(1, N_STORY + 1)])
push_u, push_V = np.array(push_u), np.array(push_V)
i_max = int(np.argmax(push_V))
V_max = push_V[i_max]

# hand check: beam-sway mechanism with simple rectangular-block capacities
def m_beam(As, b, d):
    a = As * fy / (0.85 * fc * b)
    return As * fy * (d - a / 2.0)

d_b = BEAM['h'] - BEAM['cover']
Mb = m_beam(BEAM['n_top'] * A20, BEAM['b'], d_b) + m_beam(BEAM['n_bot'] * A20, BEAM['b'], d_b)
d_c = COL['h'] - COL['cover']
P0 = 0.85 * fc * COL['b'] * COL['h']
Mc = sum(COL['n_top'] * A20 * fy * (d_c - COL['cover'])
         + N * COL['h'] / 2.0 * (1.0 - N / P0) for N in N_col.values())
h_eff = float(np.sum(F_shape * np.array(Y[1:])))   # height of the resultant lateral force
V_mech = (N_STORY * N_BAY * Mb + Mc) / h_eff

# =============================================================================
# RUN B: gravity + nonlinear time history (new, undamaged model)
# =============================================================================
build_model()
gravity(sorted(beam_ele))
omega, T = modal(3)
# %% 5. Damping and ground motion --------------------------------------------
wi, wj = omega[0], omega[2]
a0 = XI * 2.0 * wi * wj / (wi + wj)                 # mass-proportional coefficient
a1 = XI * 2.0 / (wi + wj)                           # stiffness-proportional coefficient
# rayleigh(alphaM, betaK, betaKinit, betaKcomm): here proportional to committed stiffness
ops.rayleigh(a0, 0.0, 0.0, a1)

acc = np.loadtxt(REC)                               # El Centro NS, in g
# timeSeries('Path', tag, '-dt', dt, '-values', ..., '-factor', f): record in m/s^2
ops.timeSeries('Path', 3, '-dt', DT_REC, '-values', *acc.tolist(), '-factor', SF * g)
# pattern('UniformExcitation', tag, dir, '-accel', tsTag): base acceleration in X
ops.pattern('UniformExcitation', 3, 1, '-accel', 3)

ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('BandGeneral')
ops.test('NormDispIncr', 1.0e-6, 50)
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)                # average acceleration, unconditionally stable
ops.analysis('Transient')

# %% 6. Time-history loop ----------------------------------------------------
n_steps = int(round(T_END / DT_ANA))
th_t, th_roof, th_V = [0.0], [0.0], [0.0]
drift_max = np.zeros(N_STORY)
th_ok = True
for k in range(n_steps):
    if time.time() - T_START > TIME_LIMIT:
        print('WARNING: time history stopped by the time guard at t = %.2f s' % ops.getTime())
        th_ok = False
        break
    if not dyn_step(DT_ANA):
        print('WARNING: time history failed at t = %.2f s' % ops.getTime())
        th_ok = False
        break
    u_fl = [ops.nodeDisp(ntag(f, 1), 1) for f in range(N_STORY + 1)]
    dr = np.abs([(u_fl[s] - u_fl[s - 1]) / H_STORY[s - 1] for s in range(1, N_STORY + 1)])
    drift_max = np.maximum(drift_max, dr)
    th_t.append(ops.getTime())
    th_roof.append(u_fl[-1])
    th_V.append(base_shear())
th_t, th_roof, th_V = np.array(th_t), np.array(th_roof), np.array(th_V)

# %% 7. Results ---------------------------------------------------------------
def res(name, value):
    print('RESULT', name, value)


res('W_total_kN', int(round(float(W_TOTAL))))
res('gravity_reaction_sum_kN', int(round(float(grav_sum))))
res('T1_s', round(float(T[0]), 3))
res('T2_s', round(float(T[1]), 3))
res('push_Vmax_kN', int(round(float(V_max))))
res('push_Vmax_over_W', round(float((V_max / W_TOTAL)), 3))
res('push_V_mechanism_hand_kN', int(round(float(V_mech))))
res('push_roof_drift_at_Vmax_pct', round(float((100 * push_u[i_max] / H_TOT)), 2))
res('th_roof_disp_max_m', round(float(np.max(np.abs(th_roof))), 3))
res('th_drift_max_pct', round(float((100 * drift_max.max())), 2))
res('th_drift_max_story', int(np.argmax(drift_max)) + 1)
res('th_base_shear_max_kN', int(round(float(np.max(np.abs(th_V))))))
if not th_ok:
    print('NOTE time history incomplete, results cover t <= %.2f s' % th_t[-1])

# %% 8. Figures ---------------------------------------------------------------
def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, name), dpi=110)
    plt.close(fig)

# (a) model sketch with tags
fig, ax = plt.subplots(figsize=(7, 4.2))
for t, (i, j) in col_ele.items():
    (x1, y1), (x2, y2) = coords[i], coords[j]
    ax.plot([x1, x2], [y1, y2], color='tab:blue', lw=2.5)
    ax.text(x1 + 0.12, (y1 + y2) / 2, str(t), fontsize=6, color='tab:blue', rotation=90, va='center')
for t, (i, j) in beam_ele.items():
    (x1, y1), (x2, y2) = coords[i], coords[j]
    ax.plot([x1, x2], [y1, y2], color='tab:orange', lw=2.5)
    ax.text((x1 + x2) / 2, y1 + 0.12, str(t), fontsize=6, color='tab:orange', ha='center')
for n, (x, y) in coords.items():
    ax.plot(x, y, 'ko', ms=3)
    ax.text(x - 0.1, y - (0.55 if y == 0 else 0.3), str(n), fontsize=6, ha='right')
for a in range(1, N_BAY + 2):
    x = (a - 1) * BAY
    ax.plot([x - 0.3, x + 0.3], [-0.05, -0.05], 'k', lw=3)
ax.set_title('RC frame: %d stories x %d bays, columns %.0fx%.0f cm, beams %.0fx%.0f cm'
             % (N_STORY, N_BAY, COL['b'] * 100, COL['h'] * 100, BEAM['b'] * 100, BEAM['h'] * 100),
             fontsize=9)
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
ax.set_aspect('equal'); ax.set_xlim(-1.2, N_BAY * BAY + 1.0); ax.set_ylim(-0.8, H_TOT + 0.8)
save(fig, 'ch07_model.png')

# (b) pushover curve
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(100 * push_u / H_TOT, push_V, 'b-', lw=2, label='fiber model')
ax.axhline(V_mech, color='gray', ls='--', label='hand mechanism estimate')
ax.plot(100 * push_u[i_max] / H_TOT, V_max, 'ro', label='peak base shear')
ax.set_xlabel('roof drift (%)'); ax.set_ylabel('base shear (kN)')
ax.set_title('Pushover, inverted-triangle load', fontsize=10)
ax.grid(alpha=0.3); ax.legend(fontsize=8)
save(fig, 'ch07_pushover.png')

# (c) ground acceleration and roof displacement
fig, (a1x, a2x) = plt.subplots(2, 1, figsize=(7, 4.2), sharex=True)
t_rec = np.arange(len(acc)) * DT_REC
a1x.plot(t_rec, SF * acc, 'k', lw=0.7)
a1x.set_ylabel('ag (g)'); a1x.grid(alpha=0.3)
a2x.plot(th_t, th_roof, 'b', lw=1)
a2x.set_ylabel('roof u (m)'); a2x.set_xlabel('time (s)'); a2x.grid(alpha=0.3)
a2x.set_xlim(0, T_END)
a1x.set_title('El Centro 1940 NS x %.1f: input and roof response' % SF, fontsize=10)
save(fig, 'ch07_roof_history.png')

# (d) peak story-drift profile
fig, ax = plt.subplots(figsize=(7, 4.2))
stories = np.arange(1, N_STORY + 1)
ax.plot(100 * drift_max, stories, 'b-o', lw=2, label='time history, peak')
if push_drifts:
    pd = 100 * np.abs(np.array(push_drifts[max(i_max - 1, 0)]))
    ax.plot(pd, stories, 'r--s', label='pushover at peak base shear')
ax.set_yticks(stories); ax.set_xlabel('story drift ratio (%)'); ax.set_ylabel('story')
ax.set_title('Story drift profile', fontsize=10)
ax.grid(alpha=0.3); ax.legend(fontsize=8)
save(fig, 'ch07_drift_profile.png')

print('Done in %.1f s' % (time.time() - T_START))
