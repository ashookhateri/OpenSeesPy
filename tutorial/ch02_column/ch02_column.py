# -*- coding: utf-8 -*-
"""
Chapter 02 - One column: P-Delta, buckling, fiber section, pushover
=====================================================================
Units used everywhere in this script:  kN, m, s
  -> force kN, length m, stress kN/m^2 (= kPa), moment kN*m, curvature 1/m
  g = 9.81 m/s^2

Part (a)  Elastic steel cantilever column (HEB 200) with the 'PDelta'
          geometric transformation:
            - amplification of the lateral tip drift by the axial load,
              compared with 1/(1 - P/Pcr) and with the exact solution,
            - buckling load found by bisection (lateral stiffness -> 0),
              compared with Euler, Pcr = pi^2 EI / (2L)^2.
Part (b)  Reinforced-concrete column with a fiber section
          (Concrete02 core + cover, Steel02 bars):
            - moment-curvature at constant axial load,
            - lateral pushover of the cantilever column (forceBeamColumn,
              PDelta) under the same axial load.
"""
import os
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')                     # draw straight to PNG files, no screen needed
import matplotlib.pyplot as plt
import openseespy.opensees as ops

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, 'figs')
os.makedirs(FIGDIR, exist_ok=True)

def savefig(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'ch02_' + name + '.png'), dpi=110)
    plt.close(fig)

g = 9.81                                  # m/s^2 (not needed for static work, kept for consistency)

# Robust step-by-step static solver used in part (b) ---------------------------------
ALGOS = (('Newton',), ('ModifiedNewton', '-initial'), ('KrylovNewton',))

def disp_step(node, dof, du, depth=0, max_depth=5):
    """One DisplacementControl step of size du. If every algorithm fails,
    split the step into two halves (recursively, at most max_depth times)."""
    for algo in ALGOS:
        ops.algorithm(*algo)                               # change the solution algorithm
        ops.integrator('DisplacementControl', node, dof, du)  # control dof 'dof' of 'node' by du per step
        if ops.analyze(1) == 0:                            # 0 = converged
            ops.algorithm('Newton')
            return True
    if depth >= max_depth:
        return False
    return (disp_step(node, dof, 0.5 * du, depth + 1, max_depth) and
            disp_step(node, dof, 0.5 * du, depth + 1, max_depth))

def static_objects(tol=1.0e-8, iters=30):
    """Analysis objects for a static analysis (integrator is set later)."""
    ops.wipeAnalysis()                            # clear previous analysis objects
    ops.constraints('Plain')                      # simple fixities only
    ops.numberer('RCM')                           # band-reducing DOF numbering
    ops.system('BandGeneral')                     # banded, non-symmetric solver (works also when K is not positive definite)
    ops.test('NormDispIncr', tol, iters)          # convergence: norm of displacement increment < tol
    ops.algorithm('Newton')                       # Newton-Raphson iterations


# =====================================================================================
# PART (a): elastic column with PDelta transformation
# =====================================================================================
# %% 0. Parameters, part (a)
H_col = 4.0                    # column height [m]
nEle_a = 8                     # elements along the height (PDelta works per element chord)
E_s = 200.0e6                  # steel modulus [kN/m^2]
A_s = 78.1e-4                  # HEB 200 area [m^2]
I_s = 2003.0e-8                # HEB 200 weak-axis I [m^4] (the column buckles about its weak axis)
EI_s = E_s * I_s
Hlat = 1.0                     # small lateral tip load [kN]
P_euler = math.pi**2 * EI_s / (2.0 * H_col)**2      # cantilever: effective length K L = 2 L


# %% 1. Model, part (a)
def build_elastic_column():
    ops.wipe()                                    # remove everything from the previous run
    ops.model('basic', '-ndm', 2, '-ndf', 3)      # 2D, 3 DOFs per node (ux, uy, rz)
    for i in range(nEle_a + 1):
        ops.node(i + 1, 0.0, i * H_col / nEle_a)  # node(tag, x, y): vertical column along global Y
    ops.fix(1, 1, 1, 1)                           # fix(tag, ux, uy, rz): fixed base
    # 'PDelta' transformation: small strains, but adds the geometric stiffness -N/L of
    # each element chord, so the axial force reduces the lateral stiffness.
    ops.geomTransf('PDelta', 1)
    for e in range(nEle_a):
        # elasticBeamColumn(tag, iNode, jNode, A, E, Iz, transfTag)
        ops.element('elasticBeamColumn', e + 1, e + 1, e + 2, A_s, E_s, I_s, 1)
    return nEle_a + 1                             # tag of the top node


def elastic_tip_drift(P):
    """Apply axial P (held constant) and then Hlat; return the lateral tip drift [m].
    Returns None if the solver fails (the stiffness matrix is singular)."""
    top = build_elastic_column()
    static_objects()
    # 'Linear' algorithm: one solve per step with the current tangent K = Ke + Kg(N).
    # The problem is linear here (axial force first, then a small lateral load), and
    # it avoids any convergence-tolerance trouble when P is close to Pcr.
    ops.algorithm('Linear')
    # Stage 1: axial load only
    ops.timeSeries('Linear', 1)                   # load factor = pseudo-time
    ops.pattern('Plain', 1, 1)                    # pattern(type, tag, tsTag)
    ops.load(top, 0.0, -P, 0.0)                   # load(node, Fx, Fy, Mz): compression at the top
    ops.integrator('LoadControl', 1.0)            # full axial load in one step (axial response is linear)
    ops.analysis('Static')
    if ops.analyze(1) != 0:
        return None
    ops.loadConst('-time', 0.0)                   # freeze the axial load and reset the time to 0
    # Stage 2: small lateral load with the axial force already present
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    ops.load(top, Hlat, 0.0, 0.0)
    ops.integrator('LoadControl', 1.0)
    if ops.analyze(1) != 0:
        return None
    return ops.nodeDisp(top, 1)                   # nodeDisp(node, dof): dof 1 = ux


# Model sketch data for part (a) (drawn later together with the section)
build_elastic_column()
xy_a = np.array([ops.nodeCoord(n) for n in ops.getNodeTags()])


# %% 2. P-Delta amplification, part (a)
d0 = elastic_tip_drift(0.0)                       # first-order drift (no axial load)
d0_th = Hlat * H_col**3 / (3.0 * EI_s)            # H L^3 / (3 E I)
alphas = np.arange(0.0, 0.951, 0.05)              # P / P_euler
AF_ops = []
for a in alphas:
    d = elastic_tip_drift(a * P_euler)
    AF_ops.append(d / d0 if d is not None else np.nan)
AF_ops = np.array(AF_ops)

def AF_exact(a):
    """Exact amplification for a cantilever with axial P and tip load H:
    delta = H (tan u - u) / (P k), u = kL = (pi/2) sqrt(P/Pe)."""
    if a <= 0:
        return 1.0
    u = 0.5 * math.pi * math.sqrt(a)
    return 3.0 * (math.tan(u) - u) / u**3

AF_exact_v = np.array([AF_exact(a) for a in alphas])
AF_approx_v = 1.0 / (1.0 - alphas)


# %% 3. Buckling load by bisection, part (a)
# Below Pcr the column still pushes back (drift has the sign of H); above Pcr the
# tangent lateral stiffness is negative and the linearised solve returns a drift of
# the opposite sign. Bisection on that sign change gives the numerical Pcr.
P_lo, P_hi = 0.5 * P_euler, 1.5 * P_euler
for _ in range(40):
    Pm = 0.5 * (P_lo + P_hi)
    d = elastic_tip_drift(Pm)
    if d is not None and d > 0:
        P_lo = Pm
    else:
        P_hi = Pm
    if (P_hi - P_lo) / P_lo < 1.0e-6:
        break
P_cr_ops = 0.5 * (P_lo + P_hi)

fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(alphas, AF_approx_v, '--', color='0.4', label='approx. 1/(1 - P/Pcr)')
ax.plot(alphas, AF_exact_v, '-', color='tab:blue', label='exact (beam-column theory)')
ax.plot(alphas, AF_ops, 'o', color='tab:red', ms=5, label='OpenSees, PDelta, %d elements' % nEle_a)
ax.set_xlabel('P / P_euler')
ax.set_ylabel('drift amplification  delta / delta_0')
ax.set_title('Part (a): P-Delta amplification of the tip drift')
ax.set_ylim(0, 22)
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
savefig(fig, 'pdelta')


# =====================================================================================
# PART (b): RC column with fiber section
# =====================================================================================
# %% 4. Parameters, part (b)
b_sec = 0.50                   # section width (out of plane, local z) [m]
h_sec = 0.50                   # section depth (in the bending plane, local y) [m]
cov = 0.05                     # cover measured to the bar centres [m]
fc = 30.0e3                    # concrete strength [kN/m^2] (30 MPa)
fy = 400.0e3                   # steel yield stress [kN/m^2] (400 MPa)
Es = 200.0e6                   # steel modulus [kN/m^2]
d_bar = 0.020                  # bar diameter [m]
A_bar = math.pi * d_bar**2 / 4.0
L_col = 3.0                    # height of the RC cantilever column [m]
P_rc = 0.15 * fc * b_sec * h_sec   # constant axial load = 0.15 Ag f'c [kN]

# Material tags
COVER, CORE, STEEL = 1, 2, 3
SEC = 1


# %% 5. Fiber section, part (b)
def define_rc_section():
    # Concrete02(tag, fpc, epsc0, fpcu, epsU, lambda, ft, Ets)
    #   fpc/epsc0: peak stress/strain (negative = compression); fpcu/epsU: crushing stress/strain;
    #   lambda: unloading slope ratio; ft: tensile strength; Ets: tension softening slope.
    ops.uniaxialMaterial('Concrete02', COVER, -fc, -0.002, -0.2 * fc, -0.005, 0.1, 0.1 * fc, 0.1 * fc / 0.002)
    fcc = 1.3 * fc                                # confined core: higher strength and much more ductile
    ops.uniaxialMaterial('Concrete02', CORE, -fcc, -0.0045, -0.8 * fcc, -0.02, 0.1, 0.1 * fc, 0.1 * fc / 0.002)
    # Steel02(tag, Fy, E0, b, R0, cR1, cR2): Giuffre-Menegotto-Pinto steel, b = hardening ratio
    ops.uniaxialMaterial('Steel02', STEEL, fy, Es, 0.01, 18.0, 0.925, 0.15)

    yc, zc = h_sec / 2 - cov, b_sec / 2 - cov     # core edge = line through the bar centres
    yt, zt = h_sec / 2, b_sec / 2                 # outer edge
    ops.section('Fiber', SEC)                     # start a fiber section with tag SEC
    # patch('rect', matTag, nFibY, nFibZ, yI, zI, yJ, zJ): rectangle from corner I to corner J.
    # In 2D only the y subdivision matters (bending about z).
    ops.patch('rect', CORE, 20, 1, -yc, -zc, yc, zc)          # confined core
    ops.patch('rect', COVER, 2, 1, yc, -zt, yt, zt)           # top cover strip
    ops.patch('rect', COVER, 2, 1, -yt, -zt, -yc, zt)         # bottom cover strip
    ops.patch('rect', COVER, 20, 1, -yc, zc, yc, zt)          # side cover (+z)
    ops.patch('rect', COVER, 20, 1, -yc, -zt, yc, -zc)        # side cover (-z)
    # layer('straight', matTag, nBars, areaBar, yStart, zStart, yEnd, zEnd): a row of bars
    ops.layer('straight', STEEL, 4, A_bar, yc, zc, yc, -zc)    # 4 bars at the top
    ops.layer('straight', STEEL, 4, A_bar, -yc, zc, -yc, -zc)  # 4 bars at the bottom
    ops.layer('straight', STEEL, 2, A_bar, 0.0, zc, 0.0, -zc)  # 2 bars at mid-depth

n_bars = 10
rho_l = n_bars * A_bar / (b_sec * h_sec)


# %% 6. Moment-curvature at constant axial load, part (b)
# A force-based element of length Lm = 1 m loaded by an end moment carries a uniform
# moment, so every section has the same curvature kappa = theta_tip / Lm and the same
# axial strain eps0 = u_tip / Lm. Controlling the tip rotation gives M-kappa directly,
# including the softening branch.
Lm = 1.0
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)
ops.node(1, 0.0, 0.0)
ops.node(2, Lm, 0.0)
ops.fix(1, 1, 1, 1)
define_rc_section()
ops.geomTransf('Linear', 1)                       # no P-Delta: we want the section response only
ops.beamIntegration('Lobatto', 1, SEC, 3)         # beamIntegration('Lobatto', tag, secTag, nIP)
ops.element('forceBeamColumn', 1, 1, 2, 1, 1)     # forceBeamColumn(tag, iNode, jNode, transfTag, integrationTag)

static_objects()
ops.timeSeries('Linear', 1)                       # load factor = pseudo-time
ops.pattern('Plain', 1, 1)
ops.load(2, -P_rc, 0.0, 0.0)                      # axial compression along local x
ops.integrator('LoadControl', 0.1)                # 10 steps of 0.1 -> full axial load
ops.analysis('Static')
ok_axial = ops.analyze(10)
ops.loadConst('-time', 0.0)                       # keep the axial load constant from now on

ops.timeSeries('Linear', 2)
ops.pattern('Plain', 2, 2)
ops.load(2, 0.0, 0.0, 1.0)                        # reference moment of 1 kN*m (scaled by the load factor)
kappa_max = 0.08                                  # target curvature [1/m]
n_mk = 400
dk = kappa_max / n_mk
eps_y = fy / Es
kap, mom, eps_t = [0.0], [0.0], [0.0]
mk_ok = (ok_axial == 0)
for i in range(n_mk if mk_ok else 0):
    if not disp_step(2, 3, dk * Lm):              # rotation increment = d(kappa) * Lm
        mk_ok = False
        print('WARNING: moment-curvature stopped at kappa = %.4f 1/m' % kap[-1])
        break
    ops.reactions()                               # update reactions
    th = ops.nodeDisp(2, 3)
    u = ops.nodeDisp(2, 1)
    k = th / Lm
    kap.append(k)
    mom.append(-ops.nodeReaction(1, 3))           # base reaction moment = - section moment
    eps_t.append(u / Lm + (h_sec / 2 - cov) * k)  # strain of the tension bars: eps0 - y*kappa at y = -(h/2 - cov)
kap, mom, eps_t = np.array(kap), np.array(mom), np.array(eps_t)

# first yield of the tension steel (linear interpolation)
iy = int(np.argmax(eps_t >= eps_y)) if np.any(eps_t >= eps_y) else -1
if iy > 0:
    t = (eps_y - eps_t[iy - 1]) / (eps_t[iy] - eps_t[iy - 1])
    kap_y = kap[iy - 1] + t * (kap[iy] - kap[iy - 1])
    M_y = mom[iy - 1] + t * (mom[iy] - mom[iy - 1])
else:
    kap_y, M_y = float('nan'), float('nan')
M_max = float(np.max(mom))


# %% 7. Pushover of the RC cantilever column, part (b)
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)
ops.node(1, 0.0, 0.0)
ops.node(2, 0.0, L_col)
ops.fix(1, 1, 1, 1)
define_rc_section()
ops.geomTransf('PDelta', 1)                       # include the P-Delta moment P * drift
ops.beamIntegration('Lobatto', 1, SEC, 5)         # 5 integration points: base and top sections are sampled
ops.element('forceBeamColumn', 1, 1, 2, 1, 1)

static_objects()
# Gravity: axial load in 10 load-controlled steps
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, 0.0, -P_rc, 0.0)
ops.integrator('LoadControl', 0.1)                # 10 steps of 0.1
ops.analysis('Static')
ok_grav = ops.analyze(10)
ops.loadConst('-time', 0.0)

# Lateral: displacement control of the top node, reference load 1 kN
ops.timeSeries('Linear', 2)
ops.pattern('Plain', 2, 2)
ops.load(2, 1.0, 0.0, 0.0)
target = 0.05 * L_col                             # 5 % drift
du = 0.001                                        # 1 mm per step
n_po = int(round(target / du))
drift, Vb = [0.0], [0.0]
po_ok = (ok_grav == 0)
for i in range(n_po if po_ok else 0):
    if not disp_step(2, 1, du):
        po_ok = False
        print('WARNING: pushover stopped at drift = %.2f %%' % (100 * drift[-1]))
        break
    ops.reactions()
    drift.append(ops.nodeDisp(2, 1) / L_col)
    Vb.append(-ops.nodeReaction(1, 1))            # base shear = - horizontal reaction
drift, Vb = np.array(drift), np.array(Vb)
ip = int(np.argmax(Vb))
V_max = float(Vb[ip])
drift_Vmax = float(drift[ip])
# hand check: the base moment at peak is V L + P Delta and is limited by M_max of the section
V_hand = (M_max - P_rc * drift_Vmax * L_col) / L_col


# %% 8. Figures
# Figure 1: model sketch (column of part a) and the fiber section of part b
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 4.2), gridspec_kw={'width_ratios': [1, 1.4]})
ax1.plot(xy_a[:, 0], xy_a[:, 1], '-', color='0.3', lw=3)
ax1.plot(xy_a[:, 0], xy_a[:, 1], 'o', color='tab:blue', ms=5)
for n, (x, y) in enumerate(xy_a, start=1):
    ax1.text(x + 0.08, y, str(n), fontsize=7, va='center')
ax1.plot([-0.4, 0.4], [0, 0], 'k-', lw=2)
for xx in np.linspace(-0.4, 0.3, 8):
    ax1.plot([xx, xx + 0.1], [-0.12, 0], 'k-', lw=0.8)
ax1.annotate('', xy=(0, H_col + 0.05), xytext=(0, H_col + 0.9),
             arrowprops=dict(arrowstyle='->', color='tab:red', lw=2))
ax1.text(0.08, H_col + 0.6, 'P', color='tab:red')
ax1.annotate('', xy=(0.05, H_col), xytext=(-0.8, H_col),
             arrowprops=dict(arrowstyle='->', color='tab:green', lw=2))
ax1.text(-0.8, H_col + 0.12, 'H', color='tab:green')
ax1.set_xlim(-1.2, 1.0)
ax1.set_ylim(-0.4, H_col + 1.1)
ax1.set_aspect('equal')
ax1.set_title('(a) cantilever column', fontsize=9)
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
# fiber section: z horizontal, y vertical (y = bending direction)
yc, zc, yt, zt = h_sec / 2 - cov, b_sec / 2 - cov, h_sec / 2, b_sec / 2
ax2.add_patch(plt.Rectangle((-zt, -yt), b_sec, h_sec, fc='0.85', ec='k', label='cover (Concrete02)'))
ax2.add_patch(plt.Rectangle((-zc, -yc), 2 * zc, 2 * yc, fc='0.6', ec='k', label='core (Concrete02)'))
for yy in np.linspace(-yc, yc, 21):
    ax2.plot([-zc, zc], [yy, yy], color='0.45', lw=0.5)
bars = [(z, yc) for z in np.linspace(zc, -zc, 4)] + [(z, -yc) for z in np.linspace(zc, -zc, 4)] + \
       [(zc, 0.0), (-zc, 0.0)]
for k, (z, y) in enumerate(bars):
    ax2.add_patch(plt.Circle((z, y), d_bar / 2 * 1.5, fc='tab:red', ec='k',
                             label='bars (Steel02)' if k == 0 else None))
ax2.set_xlim(-0.33, 0.33)
ax2.set_ylim(-0.33, 0.33)
ax2.set_aspect('equal')
ax2.set_xlabel('z [m]')
ax2.set_ylabel('y [m]')
ax2.set_title('(b) RC fiber section %.0fx%.0f cm, %d bars d%.0f' %
              (100 * b_sec, 100 * h_sec, n_bars, 1000 * d_bar), fontsize=9)
ax2.legend(fontsize=7, loc='upper right', framealpha=0.9)
savefig(fig, 'model')

# Figure 3: moment-curvature
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(kap, mom, '-', color='tab:blue', lw=1.5, label='OpenSees fiber section, P = %.0f kN' % P_rc)
if iy > 0:
    ax.plot([kap_y], [M_y], 'o', color='tab:red', label='first yield of tension steel')
ax.set_xlabel('curvature kappa [1/m]')
ax.set_ylabel('moment M [kN m]')
ax.set_title('Part (b): moment-curvature at constant axial load')
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
savefig(fig, 'mphi')

# Figure 4: pushover
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(100 * drift, Vb, '-', color='tab:blue', lw=1.5, label='OpenSees pushover (PDelta)')
ax.plot([100 * drift_Vmax], [V_max], 'o', color='tab:red', label='peak base shear')
ax.axhline(M_max / L_col, color='0.4', ls='--', lw=1, label='M_max / L (no P-Delta)')
ax.set_xlabel('drift ratio [%]')
ax.set_ylabel('base shear V [kN]')
ax.set_title('Part (b): lateral pushover of the RC cantilever column')
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
savefig(fig, 'pushover')


# %% 9. Results
i_half = int(np.argmin(np.abs(alphas - 0.5)))
print('RESULT', 'Pcr_euler_kN', round(P_euler, 1))
print('RESULT', 'Pcr_ops_kN', round(P_cr_ops, 1))
print('RESULT', 'AF_half_ops', round(float(AF_ops[i_half]), 3))
print('RESULT', 'AF_half_exact', round(float(AF_exact_v[i_half]), 3))
print('RESULT', 'AF_half_approx', round(float(AF_approx_v[i_half]), 3))
print('RESULT', 'P_axial_rc_kN', round(P_rc, 1))
print('RESULT', 'My_kNm', round(float(M_y), 1))
print('RESULT', 'kappa_y_1perm', round(float(kap_y), 5))
print('RESULT', 'Mmax_kNm', round(M_max, 1))
print('RESULT', 'Vmax_kN', round(V_max, 1))
print('RESULT', 'drift_at_Vmax_pct', round(100 * drift_Vmax, 2))
print('RESULT', 'Vmax_hand_kN', round(V_hand, 1))
ops.wipe()
