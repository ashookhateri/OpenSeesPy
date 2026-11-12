# -*- coding: utf-8 -*-
"""
Chapter 04 - One-bay, one-story 2D portal frame
===============================================
Units: kN, m, s  (mass in t = kN*s^2/m, E in kPa = kN/m^2), g = 9.81 m/s^2

Steps
  1. Model: two fixed-base columns and one beam (elasticBeamColumn).
  2. Gravity: uniform load w on the beam -> beam end moments vs slope-deflection.
  3. Lateral: horizontal force F at the roof -> drift, stiffness, base moment
     vs the closed-form portal-frame solution.
  4. Modal: mass lumped at the roof -> period vs T = 2*pi*sqrt(m/K).
Tags follow chapter 03: node = 100*floor + axis, column = 1000 + top node,
beam = 2000 + left node.
"""
import os
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')                        # write PNG files, no screen
import matplotlib.pyplot as plt
import openseespy.opensees as ops

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, 'figs')
os.makedirs(FIGS, exist_ok=True)
g = 9.81

# %% 0. Parameters
# ---------------------------------------------------------------------------
H = 3.5                    # story height (m)
L = 6.0                    # bay width (m)
E = 25.0e6                 # concrete Young's modulus (kPa)
bc, hc = 0.40, 0.40        # column section b x h (m)
bb, hb = 0.30, 0.50        # beam section b x h (m)
Ac, Ic = bc * hc, bc * hc**3 / 12.0
Ab, Ib = bb * hb, bb * hb**3 / 12.0
w = 30.0                   # gravity load on the beam (kN/m), also used for the mass
F = 50.0                   # lateral force at the roof (kN)

# tags (chapter 03 convention)
N1, N2 = 1, 2              # base nodes, floor 0, axes 1 and 2
N101, N102 = 101, 102      # roof nodes, floor 1
C1, C2 = 1101, 1102        # columns (1000 + top node)
B1 = 2101                  # beam (2000 + left node)
T_COL, T_BEAM = 1, 2       # geometric transformation tags


def run_static(n_steps, d_lambda):
    """Robust static loop: Newton, then ModifiedNewton -initial, then
    KrylovNewton; if all fail the step is halved (up to 4 times)."""
    algos = [('Newton',), ('ModifiedNewton', '-initial'), ('KrylovNewton',)]

    def one_step(dl, depth):
        for alg in algos:
            ops.algorithm(*alg)                   # solution algorithm
            ops.integrator('LoadControl', dl)     # load-factor increment
            if ops.analyze(1) == 0:               # 0 = converged
                ops.algorithm('Newton')
                return True
        if depth >= 4:
            return False
        return one_step(dl / 2.0, depth + 1) and one_step(dl / 2.0, depth + 1)

    for k in range(n_steps):
        if not one_step(d_lambda, 0):
            print('WARNING: static analysis failed at step', k + 1)
            return False
    return True


def static_objects():
    """Analysis objects for a static, load-controlled analysis."""
    ops.constraints('Plain')                  # only fixities: Plain handler
    ops.numberer('RCM')                       # reverse Cuthill-McKee numbering
    ops.system('BandGeneral')                 # banded general solver
    ops.test('NormDispIncr', 1.0e-10, 20)     # tolerance, max iterations
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.1)
    ops.analysis('Static')


# %% 1. Model
# ---------------------------------------------------------------------------
ops.wipe()                                    # clear memory
ops.model('basic', '-ndm', 2, '-ndf', 3)      # 2D, 3 DOF per node (ux, uy, rz)

ops.node(N1, 0.0, 0.0)                        # node(tag, x, y)
ops.node(N2, L, 0.0)
ops.node(N101, 0.0, H)
ops.node(N102, L, H)

ops.fix(N1, 1, 1, 1)                          # fixed base: ux, uy, rz restrained
ops.fix(N2, 1, 1, 1)

m_node = w * L / 2.0 / g                      # half the beam load as mass at each roof node
ops.mass(N101, m_node, m_node, 0.0)           # mass(tag, mx, my, m_rz)
ops.mass(N102, m_node, m_node, 0.0)

ops.geomTransf('Linear', T_COL)               # first-order columns (to match the hand check)
ops.geomTransf('Linear', T_BEAM)

# elasticBeamColumn(tag, iNode, jNode, A, E, Iz, transfTag)
ops.element('elasticBeamColumn', C1, N1, N101, Ac, E, Ic, T_COL)
ops.element('elasticBeamColumn', C2, N2, N102, Ac, E, Ic, T_COL)
ops.element('elasticBeamColumn', B1, N101, N102, Ab, E, Ib, T_BEAM)

# %% 2. Gravity analysis
# ---------------------------------------------------------------------------
ops.timeSeries('Linear', 1)                   # factor = pseudo-time
ops.pattern('Plain', 1, 1)                    # pattern 1 uses time series 1
ops.eleLoad('-ele', B1, '-type', '-beamUniform', -w)   # w downward (local y = +Y)
static_objects()
run_static(10, 0.1)                           # 10 steps up to the full gravity load
ops.loadConst('-time', 0.0)                   # freeze gravity, reset time

# local end forces of the beam: [N_i, V_i, M_i, N_j, V_j, M_j]
fb_g = ops.eleResponse(B1, 'localForce')
M_end_ops = abs(fb_g[2])                      # hogging moment at the beam end
ops.reactions()                               # update support reactions
Rx_g = [ops.nodeReaction(n, 1) for n in (N1, N2)]
Rm_g = [ops.nodeReaction(n, 3) for n in (N1, N2)]
ux_g = ops.nodeDisp(N101, 1)
disp_g = {n: ops.nodeDisp(n) for n in (N1, N2, N101, N102)}   # [ux, uy, rz] after gravity
fc_g = {t: ops.eleResponse(t, 'localForce') for t in (C1, C2)}

# hand check (slope-deflection, no sway, axial deformations neglected)
kc, kb = E * Ic / H, E * Ib / L               # EI/length of column and beam
rho = kb / kc
M_end_hand = w * L**2 / 12.0 * 2.0 / (2.0 + rho)
print('gravity: beam end moment OpenSees %.2f, hand %.2f kNm' % (M_end_ops, M_end_hand))
print('gravity: base thrust %.2f kN (hand %.2f)' % (abs(Rx_g[0]), 1.5 * M_end_hand / H))
print('RESULT M_end_grav_ops %.2f' % M_end_ops)
print('RESULT M_end_grav_hand %.2f' % M_end_hand)

# %% 3. Lateral load analysis (on top of gravity)
# ---------------------------------------------------------------------------
ops.pattern('Plain', 2, 1)                    # new pattern, same Linear series
ops.load(N101, F, 0.0, 0.0)                   # load(node, Fx, Fy, Mz)
run_static(10, 0.1)                           # F applied in 10 equal steps (linear model)

ux_ops = 0.5 * (ops.nodeDisp(N101, 1) + ops.nodeDisp(N102, 1)) - ux_g
K_ops = F / ux_ops
ops.reactions()
Mbase_ops = abs(ops.nodeReaction(N1, 3) - Rm_g[0])   # lateral part of the base moment
Vbase = -(ops.nodeReaction(N1, 1) + ops.nodeReaction(N2, 1) - sum(Rx_g))
fc_l = {t: ops.eleResponse(t, 'localForce') for t in (C1, C2)}
fb_l = ops.eleResponse(B1, 'localForce')
disp_l = {n: ops.nodeDisp(n) for n in (N1, N2, N101, N102)}

# hand check: fixed-base portal, beam stiffness ratio rho = (Ib/L)/(Ic/H)
K_hand = 24.0 * E * Ic / H**3 * (1.0 + 6.0 * rho) / (4.0 + 6.0 * rho)
ux_hand = F / K_hand
psi = ux_hand / H                                      # chord rotation of the columns
theta = 3.0 * kc * psi / (2.0 * kc + 3.0 * kb)         # joint rotation
Mbase_hand = 2.0 * kc * (3.0 * psi - theta)
print('lateral: base shear %.2f kN (applied %.2f)' % (Vbase, F))
print('RESULT ux_roof_ops_mm %.3f' % (1000.0 * ux_ops))
print('RESULT ux_roof_hand_mm %.3f' % (1000.0 * ux_hand))
print('RESULT K_ops %.0f' % K_ops)
print('RESULT K_hand %.0f' % K_hand)
print('RESULT M_base_lat_ops %.2f' % Mbase_ops)
print('RESULT M_base_lat_hand %.2f' % Mbase_hand)

# %% 4. Modal analysis
# ---------------------------------------------------------------------------
ops.wipeAnalysis()                            # remove the static analysis objects
lam = ops.eigen(1)                            # first eigenvalue omega^2 (sway mode)
T1_ops = 2.0 * math.pi / math.sqrt(lam[0])
phi = {n: [ops.nodeEigenvector(n, 1, d) for d in (1, 2, 3)] if n > 100 else [0.0, 0.0, 0.0]
       for n in (N1, N2, N101, N102)}
T1_hand = 2.0 * math.pi * math.sqrt(2.0 * m_node / K_hand)
print('periods (s):', ', '.join('%.4f' % (2 * math.pi / math.sqrt(x)) for x in lam))
print('RESULT T1_ops %.4f' % T1_ops)
print('RESULT T1_hand %.4f' % T1_hand)

# %% 5. Figures
# ---------------------------------------------------------------------------
XY = {N1: (0.0, 0.0), N2: (L, 0.0), N101: (0.0, H), N102: (L, H)}
ELE = {C1: (N1, N101), C2: (N2, N102), B1: (N101, N102)}


def hermite_shape(i, j, ui, uj, scale, n=21):
    """Deformed shape of a 2D beam element from its end displacements
    (ux, uy, rz) using cubic Hermite functions for the transverse part."""
    (xi, yi), (xj, yj) = XY[i], XY[j]
    Le = math.hypot(xj - xi, yj - yi)
    c, s = (xj - xi) / Le, (yj - yi) / Le
    ua_i, ua_j = c * ui[0] + s * ui[1], c * uj[0] + s * uj[1]      # axial
    v_i, v_j = -s * ui[0] + c * ui[1], -s * uj[0] + c * uj[1]      # transverse
    xs = np.linspace(0.0, 1.0, n)
    N1_ = 1 - 3 * xs**2 + 2 * xs**3
    N2_ = Le * (xs - 2 * xs**2 + xs**3)
    N3_ = 3 * xs**2 - 2 * xs**3
    N4_ = Le * (-xs**2 + xs**3)
    v = N1_ * v_i + N2_ * ui[2] + N3_ * v_j + N4_ * uj[2]
    u = (1 - xs) * ua_i + xs * ua_j
    X = xi + xs * Le * c + scale * (u * c - v * s)
    Y = yi + xs * Le * s + scale * (u * s + v * c)
    return X, Y


def draw_frame(ax, color='0.75', lw=1.0):
    for t, (i, j) in ELE.items():
        ax.plot([XY[i][0], XY[j][0]], [XY[i][1], XY[j][1]], color=color, lw=lw)


# (a) model sketch with tags, dimensions and loads
fig, ax = plt.subplots(figsize=(7, 4.2))
draw_frame(ax, 'k', 2.0)
for n, (x, y) in XY.items():
    ax.plot(x, y, 'ko', ms=4)
    ax.text(x + 0.12, y + 0.12, str(n), fontsize=8, color='tab:blue')
for t, (i, j) in ELE.items():
    xm, ym = 0.5 * (XY[i][0] + XY[j][0]), 0.5 * (XY[i][1] + XY[j][1])
    ax.text(xm + (0.15 if t < 2000 else 0.0), ym - (0.0 if t < 2000 else 0.35), str(t),
            fontsize=8, color='tab:red', ha='left' if t < 2000 else 'center')
for x in (0.0, L):
    ax.plot([x - 0.35, x + 0.35], [0, 0], 'k-', lw=4)          # fixed supports
for x in np.linspace(0.3, L - 0.3, 9):                          # distributed load arrows
    ax.annotate('', xy=(x, H + 0.05), xytext=(x, H + 0.6),
                arrowprops=dict(arrowstyle='->', color='tab:green', lw=1))
ax.plot([0, L], [H + 0.6, H + 0.6], color='tab:green', lw=1)
ax.text(L / 2, H + 0.7, 'w = %.0f kN/m' % w, color='tab:green', ha='center', fontsize=8)
ax.annotate('', xy=(-0.05, H), xytext=(-1.4, H),
            arrowprops=dict(arrowstyle='->', color='tab:orange', lw=2))
ax.text(-1.4, H + 0.15, 'F = %.0f kN' % F, color='tab:orange', fontsize=8)
ax.annotate('', xy=(0, -0.6), xytext=(L, -0.6), arrowprops=dict(arrowstyle='<->', lw=0.8))
ax.text(L / 2, -0.55, 'L = %.1f m' % L, ha='center', va='bottom', fontsize=8)
ax.annotate('', xy=(L + 0.8, 0), xytext=(L + 0.8, H), arrowprops=dict(arrowstyle='<->', lw=0.8))
ax.text(L + 0.9, H / 2, 'H = %.1f m' % H, va='center', fontsize=8)
ax.text(L / 2, H / 2, 'columns %.2fx%.2f m\nbeam %.2fx%.2f m\nE = %.0f GPa'
        % (bc, hc, bb, hb, E / 1e6), ha='center', va='center', fontsize=8, color='0.3')
ax.set_xlim(-2.0, L + 2.0); ax.set_ylim(-1.0, H + 1.2)
ax.set_aspect('equal'); ax.axis('off')
ax.set_title('Portal frame: tags, geometry and loads', fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'ch04_model.png'), dpi=110)
plt.close(fig)


# (b) bending moment diagrams, gravity and lateral parts (drawn on the tension side)
def moment_along(q, wy, Le, n=21):
    """Internal moment (sagging +) along a member from its local end forces q."""
    xs = np.linspace(0.0, Le, n)
    return xs, -q[2] + q[1] * xs + wy * xs**2 / 2.0


fig, axs = plt.subplots(1, 2, figsize=(7, 4.2))
cases = [('gravity', {C1: fc_g[C1], C2: fc_g[C2], B1: fb_g}, -w),
         ('lateral F only', {t: [a - b for a, b in zip(fc_l[t], fc_g[t])] for t in (C1, C2)}, 0.0)]
cases[1][1][B1] = [a - b for a, b in zip(fb_l, fb_g)]
for ax, (name, Q, wy) in zip(axs, cases):
    draw_frame(ax, 'k', 1.5)
    Mmax = max(np.max(np.abs(moment_along(Q[t], wy if t == B1 else 0.0,
                                          L if t == B1 else H)[1])) for t in ELE)
    sc = 1.2 / Mmax if Mmax > 0 else 1.0
    for t, (i, j) in ELE.items():
        (xi, yi), (xj, yj) = XY[i], XY[j]
        Le = math.hypot(xj - xi, yj - yi)
        c, s = (xj - xi) / Le, (yj - yi) / Le
        xs, M = moment_along(Q[t], wy if t == B1 else 0.0, Le)
        # sagging (+) puts tension on the local -y side
        X = xi + xs * c + sc * M * s
        Y = yi + xs * s - sc * M * c
        ax.fill(np.concatenate([[xi], X, [xj]]), np.concatenate([[yi], Y, [yj]]),
                color='tab:blue' if name == 'gravity' else 'tab:orange', alpha=0.35)
        picks = [(0, 'left'), (len(xs) - 1, 'right')]
        if t == B1 and wy != 0.0:
            picks.append((len(xs) // 2, 'center'))     # midspan value under gravity
        for k, ha in picks:
            if t != B1:
                ha = 'right' if X[k] < XY[i][0] else 'left'
            ax.text(X[k], Y[k], '%.1f' % abs(M[k]), fontsize=6.5, ha=ha,
                    va='bottom' if Y[k] >= (yi + yj) / 2 or t != B1 else 'top')
    ax.set_title('M (kNm), %s' % name, fontsize=9)
    ax.set_aspect('equal'); ax.axis('off')
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'ch04_moments.png'), dpi=110)
plt.close(fig)

# (c) deformed shape under F and mode 1 shape
fig, axs = plt.subplots(1, 2, figsize=(7, 4.2))
# lateral part only = (gravity + lateral) - gravity
d_lat = {n: [a - b for a, b in zip(disp_l[n], disp_g[n])] for n in disp_l}
for ax, (name, D) in zip(axs, [('deformed under F (x%d)', d_lat), ('mode 1, T = %.3f s' % T1_ops, phi)]):
    umax = max(abs(D[n][0]) for n in (N101, N102)) or 1.0
    sc = 0.6 / umax
    draw_frame(ax)
    for t, (i, j) in ELE.items():
        X, Y = hermite_shape(i, j, D[i], D[j], sc)
        ax.plot(X, Y, color='tab:red', lw=2)
    ttl = name % sc if '%d' in name else name
    ax.set_title(ttl, fontsize=9)
    ax.set_aspect('equal'); ax.axis('off')
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'ch04_deformed_mode.png'), dpi=110)
plt.close(fig)

ops.wipe()
print('done, figures in', FIGS)
