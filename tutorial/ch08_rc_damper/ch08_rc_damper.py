# -*- coding: utf-8 -*-
"""
Chapter 08 - RC frame with nonlinear viscous dampers
=====================================================
The SAME 4-story x 3-bay RC moment frame as chapter 07 (same geometry, fiber
sections, materials, loads and masses; the code is copied so that this script is
self-contained). One diagonal ViscousDamper (twoNodeLink) is added in the middle
bay of every story.

Analyses:  for the bare frame and for the damped frame:
           gravity -> eigen -> nonlinear time history (El Centro 1940 NS)
           then compare roof displacement, peak story drifts and base shear,
           and plot the damper force-deformation loops.

UNITS:  force = kN,  length = m,  time = s
        stress = kN/m^2 (= kPa),  mass = kN*s^2/m (= tonne),  g = 9.81 m/s^2

TAGS (convention of chapter 03):
  node     = 100*floor + axis          (floor 0 = base, axis 1..4 from left)
  column   = 1000 + 100*story + axis   (= top node + 1000; story s: floor s-1 -> s)
  beam     = 2000 + 100*floor + bay    (= left node + 2000; bay b: axis b -> b+1)
  damper   = 4000 + 100*story + bay    (diagonal in bay b of story s)
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
# 20-29 devices
MAT_STEEL, MAT_COVER, MAT_CORE, MAT_DAMPER = 1, 10, 11, 21
SEC_COL, SEC_BEAM = 101, 201
TR_COL, TR_BEAM = 1, 2                      # PDelta for columns, Linear for beams
INT_COL, INT_BEAM = 1, 2                    # beamIntegration tags
N_IP = 5                                    # Lobatto integration points per element

# dynamic analysis
DT_REC, DT_ANA, T_END = 0.02, 0.01, 25.0    # record step, analysis step, duration (s)
SF = 1.0                                    # scale factor on the record
XI = 0.05                                   # Rayleigh damping ratio (modes 1 and 3)
TIME_LIMIT = 160.0                          # wall-clock guard for the whole script (s)

# viscous damper (one per story, bay 2)
DAMP_BAY = 2
KD = 1.0e5                                  # series (brace + device) stiffness, kN/m
CD = 1000.0                                 # damping coefficient, kN*(s/m)^alpha
ALPHA = 0.5                                 # velocity exponent (1 = linear damper)


def ntag(floor, axis):
    """Node tag = 100*floor + axis."""
    return 100 * floor + axis


def col_tag(story, axis):
    return 1000 + 100 * story + axis


def beam_tag(floor, bay):
    return 2000 + 100 * floor + bay


def damper_tag(story, bay):
    return 4000 + 100 * story + bay


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


def build_model(with_dampers=False):
    """Builds the chapter-07 frame; with_dampers=True adds one ViscousDamper diagonal
    per story in bay DAMP_BAY. Returns dicts used for plotting and post-processing."""
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

    # viscous dampers ----------------------------------------------------
    dampers = {}
    if with_dampers:
        # ViscousDamper(tag, K, Cd, alpha): Maxwell model, a spring K in series with a
        # nonlinear dashpot F = Cd * sign(v) * |v|^alpha  (v = deformation rate of the link)
        ops.uniaxialMaterial('ViscousDamper', MAT_DAMPER, KD, CD, ALPHA)
        for s in range(1, N_STORY + 1):
            # diagonal from bottom-right to top-left of the bay (as in chapter 03)
            t = damper_tag(s, DAMP_BAY)
            dampers[t] = (ntag(s - 1, DAMP_BAY + 1), ntag(s, DAMP_BAY))
            # twoNodeLink(tag, iNode, jNode, '-mat', matTag, '-dir', 1): a zero-mass link
            # whose local x axis runs from iNode to jNode; dir 1 = axial deformation.
            # It adds no Rayleigh damping (no '-doRayleigh' flag), so only the device damps.
            ops.element('twoNodeLink', t, *dampers[t], '-mat', MAT_DAMPER, '-dir', 1)

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
    return coords, cols, beams, dampers, floor_mass


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


def damper_state(tag, i, j):
    """Axial deformation (from node displacements) and axial force of one damper link."""
    (xi, yi), (xj, yj) = coords[i], coords[j]
    L = math.hypot(xj - xi, yj - yi)
    ex, ey = (xj - xi) / L, (yj - yi) / L
    du = [ops.nodeDisp(j, 1) - ops.nodeDisp(i, 1), ops.nodeDisp(j, 2) - ops.nodeDisp(i, 2)]
    defo = du[0] * ex + du[1] * ey                  # elongation positive
    f = None
    if FORCE_KEY[0] is not None:
        r = ops.eleResponse(tag, FORCE_KEY[0])      # basic (axial) force of the link
        if r is not None and len(r) > 0:
            f = r[0]
    if f is None:                                   # fallback: project the global end force
        r = ops.eleResponse(tag, 'globalForce')     # [Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j]
        f = r[3] * ex + r[4] * ey
    return defo, f


def pick_force_key(tag):
    """Finds which eleResponse keyword returns the link's basic force in this version."""
    for key in ('basicForce', 'basicForces'):
        try:
            r = ops.eleResponse(tag, key)
        except Exception:
            r = []
        if r is not None and len(r) > 0:
            return key
    return None


def run_time_history(with_dampers, a0, a1, time_budget):
    """Builds the model, applies gravity, El Centro, and returns the response histories."""
    global coords
    coords, cols, beams, dampers, floor_mass = build_model(with_dampers)
    gravity(sorted(beams))
    omega, T = modal(3)

    # Rayleigh damping: SAME coefficients for both runs (taken from the bare frame), so the
    # only difference between the runs is the dampers.
    # rayleigh(alphaM, betaK, betaKinit, betaKcomm)
    ops.rayleigh(a0, 0.0, 0.0, a1)
    # timeSeries('Path', tag, '-dt', dt, '-values', ..., '-factor', f): record in m/s^2
    ops.timeSeries('Path', 3, '-dt', DT_REC, '-values', *ACC.tolist(), '-factor', SF * g)
    # pattern('UniformExcitation', tag, dir, '-accel', tsTag): base acceleration in X
    ops.pattern('UniformExcitation', 3, 1, '-accel', 3)

    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-6, 50)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)            # average acceleration
    ops.analysis('Transient')

    d_tags = sorted(dampers)
    if d_tags:
        FORCE_KEY[0] = pick_force_key(d_tags[0])
    t0 = time.time()
    out = dict(t=[0.0], roof=[0.0], V=[0.0], drift_max=np.zeros(N_STORY),
               d_def={t: [0.0] for t in d_tags}, d_force={t: [0.0] for t in d_tags},
               T=T, ok=True, dampers=dampers, cols=cols, beams=beams)
    for k in range(int(round(T_END / DT_ANA))):
        if time.time() - t0 > time_budget:
            print('WARNING: time guard stopped the run at t = %.2f s' % ops.getTime())
            out['ok'] = False
            break
        if not dyn_step(DT_ANA):
            print('WARNING: no convergence at t = %.2f s' % ops.getTime())
            out['ok'] = False
            break
        u_fl = [ops.nodeDisp(ntag(f, 1), 1) for f in range(N_STORY + 1)]
        dr = np.abs([(u_fl[s] - u_fl[s - 1]) / H_STORY[s - 1] for s in range(1, N_STORY + 1)])
        out['drift_max'] = np.maximum(out['drift_max'], dr)
        out['t'].append(ops.getTime())
        out['roof'].append(u_fl[-1])
        out['V'].append(base_shear())               # includes the damper force at the base
        for t in d_tags:
            dd, ff = damper_state(t, *dampers[t])
            out['d_def'][t].append(dd)
            out['d_force'][t].append(ff)
    for key in ('t', 'roof', 'V'):
        out[key] = np.array(out[key])
    print('run with_dampers=%s finished in %.1f s' % (with_dampers, time.time() - t0))
    return out


# =============================================================================
# %% 3. Bare frame: modal properties, Rayleigh coefficients, mode shape
# =============================================================================
ACC = np.loadtxt(REC)                               # El Centro NS, in g
FORCE_KEY = [None]
coords, col_ele, beam_ele, _, floor_mass = build_model(False)
gravity(sorted(beam_ele))
omega, T_bare = modal(3)
wi, wj = omega[0], omega[2]
A0 = XI * 2.0 * wi * wj / (wi + wj)                 # mass-proportional coefficient
A1 = XI * 2.0 / (wi + wj)                           # stiffness-proportional coefficient
# first mode shape at the floors (ux of axis 1), normalised to 1 at the roof
phi = np.array([ops.nodeEigenvector(ntag(f, 1), 1, 1) for f in range(N_STORY + 1)])
phi = phi / phi[-1]
m_fl = np.array([floor_mass[f] for f in range(1, N_STORY + 1)])
print('Bare-frame periods (s):', ['%.3f' % x for x in T_bare])

# %% 4. Time history WITHOUT dampers ------------------------------------------
budget = 0.45 * TIME_LIMIT
bare = run_time_history(False, A0, A1, budget)

# %% 5. Time history WITH dampers ---------------------------------------------
budget = max(30.0, TIME_LIMIT - (time.time() - T_START) - 10.0)
damp = run_time_history(True, A0, A1, budget)
damper_ele = damp['dampers']

# %% 6. Hand check: added damping ratio of the dampers (energy method) --------
# beta_add = sum_j lam*C*(phi_r,j*cos(theta_j))^(1+a) / (2*pi*A^(1-a)*w^(2-a)*sum_i m_i*phi_i^2)
# lam = 2^(2+a)*Gamma(1+a/2)^2/Gamma(2+a)  (= pi for a linear damper, a = 1)
lam = 2.0 ** (2 + ALPHA) * math.gamma(1 + ALPHA / 2) ** 2 / math.gamma(2 + ALPHA)
A_roof = np.max(np.abs(damp['roof']))               # amplitude: peak roof displacement
num = 0.0
for s in range(1, N_STORY + 1):
    cos_t = BAY / math.hypot(BAY, H_STORY[s - 1])
    num += lam * CD * (abs(phi[s] - phi[s - 1]) * cos_t) ** (1 + ALPHA)
xi_add = num / (2 * math.pi * A_roof ** (1 - ALPHA) * omega[0] ** (2 - ALPHA)
                * np.sum(m_fl * phi[1:] ** 2))

# %% 7. Results ---------------------------------------------------------------
def res(name, value):
    print('RESULT', name, value)


r_bare, r_damp = np.max(np.abs(bare['roof'])), np.max(np.abs(damp['roof']))
f_max = max(np.max(np.abs(v)) for v in damp['d_force'].values())
res('T1_bare_s', round(float(T_bare[0]), 3))
res('T1_dampers_locked_s', round(float(damp['T'][0]), 3))
res('roof_max_bare_m', round(float(r_bare), 3))
res('roof_max_damped_m', round(float(r_damp), 3))
res('roof_reduction_pct', round(float(100 * (1 - r_damp / r_bare)), 1))
res('drift_max_bare_pct', round(float(100 * bare['drift_max'].max()), 2))
res('drift_max_damped_pct', round(float(100 * damp['drift_max'].max()), 2))
res('Vb_max_bare_kN', int(round(float(np.max(np.abs(bare['V']))))))
res('Vb_max_damped_kN', int(round(float(np.max(np.abs(damp['V']))))))
res('damper_force_max_kN', int(round(float(f_max))))
res('xi_added_est_pct', round(float(100 * xi_add), 1))
if not (bare['ok'] and damp['ok']):
    print('NOTE at least one run is incomplete: bare to t = %.2f s, damped to t = %.2f s'
          % (bare['t'][-1], damp['t'][-1]))

# %% 8. Figures ---------------------------------------------------------------
def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, name), dpi=110)
    plt.close(fig)

# (a) model sketch with tags
fig, ax = plt.subplots(figsize=(7, 4.2))
for t, (i, j) in damp['cols'].items():
    (x1, y1), (x2, y2) = coords[i], coords[j]
    ax.plot([x1, x2], [y1, y2], color='tab:blue', lw=2.5)
    ax.text(x1 + 0.12, (y1 + y2) / 2, str(t), fontsize=6, color='tab:blue', rotation=90, va='center')
for t, (i, j) in damp['beams'].items():
    (x1, y1), (x2, y2) = coords[i], coords[j]
    ax.plot([x1, x2], [y1, y2], color='tab:orange', lw=2.5)
    ax.text((x1 + x2) / 2, y1 + 0.12, str(t), fontsize=6, color='tab:orange', ha='center')
for t, (i, j) in damper_ele.items():
    (x1, y1), (x2, y2) = coords[i], coords[j]
    ax.plot([x1, x2], [y1, y2], color='tab:purple', lw=1.5, ls='--')
    xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
    ax.plot(xm, ym, 's', color='tab:purple', ms=7)
    ax.text(xm + 0.25, ym, str(t), fontsize=6, color='tab:purple', va='center')
for n, (x, y) in coords.items():
    ax.plot(x, y, 'ko', ms=3)
    ax.text(x - 0.1, y - (0.55 if y == 0 else 0.3), str(n), fontsize=6, ha='right')
for a in range(1, N_BAY + 2):
    x = (a - 1) * BAY
    ax.plot([x - 0.3, x + 0.3], [-0.05, -0.05], 'k', lw=3)
ax.set_title('Chapter 07 frame + ViscousDamper diagonals (bay %d): C = %.0f, alpha = %.2f'
             % (DAMP_BAY, CD, ALPHA), fontsize=9)
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
ax.set_aspect('equal'); ax.set_xlim(-1.2, N_BAY * BAY + 1.0); ax.set_ylim(-0.8, H_TOT + 0.8)
save(fig, 'ch08_model.png')

# (b) roof displacement, with and without dampers
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(bare['t'], bare['roof'], color='gray', lw=1, label='without dampers')
ax.plot(damp['t'], damp['roof'], color='tab:purple', lw=1.2, label='with dampers')
ax.set_xlabel('time (s)'); ax.set_ylabel('roof displacement (m)')
ax.set_title('El Centro 1940 NS x %.1f: roof displacement' % SF, fontsize=10)
ax.set_xlim(0, T_END); ax.grid(alpha=0.3); ax.legend(fontsize=8)
save(fig, 'ch08_roof_compare.png')

# (c) peak story drift profiles
fig, ax = plt.subplots(figsize=(7, 4.2))
stories = np.arange(1, N_STORY + 1)
ax.plot(100 * bare['drift_max'], stories, '-o', color='gray', lw=2, label='without dampers')
ax.plot(100 * damp['drift_max'], stories, '-s', color='tab:purple', lw=2, label='with dampers')
ax.set_yticks(stories); ax.set_xlabel('peak story drift ratio (%)'); ax.set_ylabel('story')
ax.set_title('Peak story drifts', fontsize=10)
ax.grid(alpha=0.3); ax.legend(fontsize=8)
save(fig, 'ch08_drift_compare.png')

# (d) damper hysteresis loops, one panel per story
fig, axs = plt.subplots(2, 2, figsize=(7, 4.2), sharex=True, sharey=True)
for k, t in enumerate(sorted(damper_ele)):
    axk = axs.flat[k]
    axk.plot(1000 * np.array(damp['d_def'][t]), damp['d_force'][t], color='tab:purple', lw=0.6)
    axk.set_title('story %d (link %d)' % ((t - 4000) // 100, t), fontsize=8)
    axk.grid(alpha=0.3)
for axk in axs[1, :]:
    axk.set_xlabel('deformation (mm)', fontsize=8)
for axk in axs[:, 0]:
    axk.set_ylabel('force (kN)', fontsize=8)
save(fig, 'ch08_damper_hysteresis.png')

print('Done in %.1f s' % (time.time() - T_START))
