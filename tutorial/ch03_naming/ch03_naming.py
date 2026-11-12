# -*- coding: utf-8 -*-
"""
Chapter 03 - A tag (naming) convention for frame models
=======================================================
Units: kN, m, s  (so mass in t = kN*s^2/m, stress in kPa = kN/m^2), g = 9.81 m/s^2

What this script does
  1. Defines a tag convention as small Python functions (node, element,
     material, section, transformation tags).
  2. Generates an N-story x M-bay 2D frame from parameters with these tags
     (columns, beams, braces in one bay, linear viscous dampers in another bay).
  3. Runs a gravity analysis and an eigen analysis on it, using the tag
     formulas to find the base nodes, the roof nodes and so on.
  4. Generates a small 3D grid frame with the 3D convention and runs eigen.
  5. Plots both models with node and element labels.

Tag convention (2D)
  node     = 100*floor + axis               floor 0..9, axis 1..99   e.g. 302
  element  = 1000*type + 100*story + axis   type 1 col, 2 beam, 3 brace, 4 damper
             column  -> top node tag  + 1000 (column 1302 ends at node 302)
             beam    -> left node tag + 2000 (beam 2302 starts at node 302)
Tag convention (3D)
  node     = 10000*floor + 100*axisY + axisX                         e.g. 20301
  element  = 1000000*type + 100000*dir + node tag (top / start node)
             dir 0 = vertical or along X, dir 1 = along Y
             column 1020301, X-beam 2020301, Y-beam 2120301
"""
import os
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')                      # no screen needed: write PNG files
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D    # noqa: F401  (registers the 3d projection)
import openseespy.opensees as ops

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, 'figs')
os.makedirs(FIGS, exist_ok=True)
g = 9.81

# %% 0. The tag convention as functions
# ---------------------------------------------------------------------------
# Element type digit (first digit of an element tag)
COL, BEAM, BRACE, DAMPER = 1, 2, 3, 4
TYPE_NAME = {COL: 'column', BEAM: 'beam', BRACE: 'brace', DAMPER: 'damper'}

# Material tags: 1-9 steel, 10-19 concrete, 20-29 devices (dampers, links)
MAT = {'steel_brace': 1, 'damper_linear': 21}
# Section tags: 1xx columns, 2xx beams, 3xx braces (same digit as element type)
SEC = {'col': 101, 'beam': 201}
# Geometric transformation tags: 2D -> 1 columns, 2 beams
#                                 3D -> 11 columns, 12 X-beams, 13 Y-beams
TRANSF = {'col2d': 1, 'beam2d': 2, 'col3d': 11, 'beamX3d': 12, 'beamY3d': 13}


def node2d(floor, axis):
    """2D node tag: floor 0 = base, axis 1 = left-most column line."""
    assert 0 <= floor <= 9 and 1 <= axis <= 99, 'outside the 2D tag range'
    return 100 * floor + axis


def ele2d(etype, story, axis):
    """2D element tag: type digit, then story, then axis (or bay) number."""
    assert etype in TYPE_NAME and 1 <= story <= 9 and 1 <= axis <= 99
    return 1000 * etype + 100 * story + axis


def decode_ele2d(tag):
    """Inverse of ele2d: tag -> (type name, story, axis)."""
    return TYPE_NAME[tag // 1000], (tag % 1000) // 100, tag % 100


def node3d(floor, iy, ix):
    """3D node tag: floor, axis number along Y, axis number along X."""
    assert 0 <= floor <= 9 and 1 <= iy <= 99 and 1 <= ix <= 99
    return 10000 * floor + 100 * iy + ix


def ele3d(etype, direction, node_tag):
    """3D element tag: type digit, direction digit, then 5-digit node tag."""
    assert etype in TYPE_NAME and direction in (0, 1) and node_tag < 100000
    return 1000000 * etype + 100000 * direction + node_tag


# quick self-check of the formulas (round trip)
assert decode_ele2d(ele2d(COL, 3, 2)) == ('column', 3, 2)
assert ele2d(COL, 3, 2) == 1000 + node2d(3, 2)      # column tag = 1000 + top node


def run_static(n_steps, d_lambda):
    """Robust static loop: Newton, then ModifiedNewton -initial, then
    KrylovNewton; if all fail, the step is cut in half (up to 4 times)."""
    algos = [('Newton',), ('ModifiedNewton', '-initial'), ('KrylovNewton',)]

    def one_step(dl, depth):
        for alg in algos:
            ops.algorithm(*alg)                   # choose the solution algorithm
            ops.integrator('LoadControl', dl)     # load factor increment dl
            if ops.analyze(1) == 0:               # 0 means converged
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


# %% 1. Parameters of the 2D frame
# ---------------------------------------------------------------------------
N_STORY = 4                       # number of stories (<= 9 with this convention)
BAYS = [5.0, 5.0, 5.0]            # bay widths (m) -> M = 3 bays
HEIGHTS = [3.6, 3.2, 3.2, 3.2]    # story heights (m), bottom to top
E_c = 25.0e6                      # concrete modulus (kPa)
b_c, h_c = 0.40, 0.40             # column section (m)
b_b, h_b = 0.30, 0.50             # beam section (m)
A_c, I_c = b_c * h_c, b_c * h_c**3 / 12.0
A_b, I_b = b_b * h_b, b_b * h_b**3 / 12.0
E_s, A_br = 200.0e6, 2.0e-3       # brace steel modulus (kPa) and area (m2)
eta_d, A_d = 2.0e4, 1.0           # damper: stress per strain rate (kPa*s), area
w_grav = 25.0                     # gravity load on every beam (kN/m)
BRACED_BAY = 1                    # braces in bay 1
DAMPED_BAY = 3                    # dampers in bay 3
M_BAY = len(BAYS)
assert len(HEIGHTS) == N_STORY

x_axis = np.concatenate([[0.0], np.cumsum(BAYS)])       # x of axes 1..M+1
y_floor = np.concatenate([[0.0], np.cumsum(HEIGHTS)])   # y of floors 0..N

# %% 2. Build the 2D model
# ---------------------------------------------------------------------------
ops.wipe()                                    # clear any previous model
ops.model('basic', '-ndm', 2, '-ndf', 3)      # 2 dimensions, 3 DOF per node (ux, uy, rz)

# nodes: node(tag, x, y)
for fl in range(N_STORY + 1):
    for ax in range(1, M_BAY + 2):
        ops.node(node2d(fl, ax), float(x_axis[ax - 1]), float(y_floor[fl]))

# supports: fix(tag, ux, uy, rz) with 1 = fixed -> base nodes are floor 0
for ax in range(1, M_BAY + 2):
    ops.fix(node2d(0, ax), 1, 1, 1)

# lumped masses from the beam gravity load (tributary half-bays), X and Y only
for fl in range(1, N_STORY + 1):
    for ax in range(1, M_BAY + 2):
        trib = 0.0
        if ax > 1:
            trib += 0.5 * BAYS[ax - 2]
        if ax <= M_BAY:
            trib += 0.5 * BAYS[ax - 1]
        m = w_grav * trib / g
        ops.mass(node2d(fl, ax), m, m, 0.0)   # mass(tag, mx, my, m_rot)

# materials: uniaxialMaterial('Elastic', tag, E, eta)
ops.uniaxialMaterial('Elastic', MAT['steel_brace'], E_s)          # brace steel
ops.uniaxialMaterial('Elastic', MAT['damper_linear'], 0.0, eta_d)  # zero stiffness, viscous eta
#   (a linear viscous placeholder; chapter 08 uses the ViscousDamper material)

# sections: section('Elastic', tag, E, A, Iz) and an integration rule with
# the same number: beamIntegration('Lobatto', tag, secTag, nIP)
ops.section('Elastic', SEC['col'], E_c, A_c, I_c)
ops.section('Elastic', SEC['beam'], E_c, A_b, I_b)
ops.beamIntegration('Lobatto', SEC['col'], SEC['col'], 5)
ops.beamIntegration('Lobatto', SEC['beam'], SEC['beam'], 5)

# geometric transformations: PDelta for columns, Linear for beams
ops.geomTransf('PDelta', TRANSF['col2d'])
ops.geomTransf('Linear', TRANSF['beam2d'])

elements2d = {}   # tag -> (iNode, jNode, type)  kept for plotting / bookkeeping

# columns: forceBeamColumn(tag, iNode, jNode, transfTag, integrationTag)
for st in range(1, N_STORY + 1):
    for ax in range(1, M_BAY + 2):
        t, i, j = ele2d(COL, st, ax), node2d(st - 1, ax), node2d(st, ax)
        ops.element('forceBeamColumn', t, i, j, TRANSF['col2d'], SEC['col'])
        elements2d[t] = (i, j, COL)

# beams: numbered by floor and bay (= left axis)
for fl in range(1, N_STORY + 1):
    for bay in range(1, M_BAY + 1):
        t, i, j = ele2d(BEAM, fl, bay), node2d(fl, bay), node2d(fl, bay + 1)
        ops.element('forceBeamColumn', t, i, j, TRANSF['beam2d'], SEC['beam'])
        elements2d[t] = (i, j, BEAM)

# braces: Truss(tag, iNode, jNode, A, matTag), bottom-left to top-right
for st in range(1, N_STORY + 1):
    b = BRACED_BAY
    t, i, j = ele2d(BRACE, st, b), node2d(st - 1, b), node2d(st, b + 1)
    ops.element('Truss', t, i, j, A_br, MAT['steel_brace'])
    elements2d[t] = (i, j, BRACE)

# dampers: Truss with the viscous material, bottom-right to top-left
for st in range(1, N_STORY + 1):
    b = DAMPED_BAY
    t, i, j = ele2d(DAMPER, st, b), node2d(st - 1, b + 1), node2d(st, b)
    ops.element('Truss', t, i, j, A_d, MAT['damper_linear'])
    elements2d[t] = (i, j, DAMPER)

n_nodes_2d = (N_STORY + 1) * (M_BAY + 1)
count = {k: sum(1 for v in elements2d.values() if v[2] == k) for k in TYPE_NAME}
print('2D model: %d nodes, %d elements' % (n_nodes_2d, len(elements2d)))
for k, name in TYPE_NAME.items():
    print('   %-7s %3d   tags %s ... %s' % (name, count[k],
          min(t for t in elements2d if t // 1000 == k),
          max(t for t in elements2d if t // 1000 == k)))

# %% 3. Gravity analysis of the 2D frame
# ---------------------------------------------------------------------------
ops.timeSeries('Linear', 1)                   # load factor grows linearly with pseudo-time
ops.pattern('Plain', 1, 1)                    # pattern(Plain, patternTag, tsTag)
for fl in range(1, N_STORY + 1):
    for bay in range(1, M_BAY + 1):
        # uniform load in local y; beams run left->right so local y = global +Y
        ops.eleLoad('-ele', ele2d(BEAM, fl, bay), '-type', '-beamUniform', -w_grav)

ops.constraints('Plain')                      # only single-point fixities: Plain is enough
ops.numberer('RCM')                           # reverse Cuthill-McKee renumbering
ops.system('BandGeneral')                     # banded solver
ops.test('NormDispIncr', 1.0e-8, 25)          # convergence test (tol, max iterations)
ops.algorithm('Newton')
ops.integrator('LoadControl', 0.1)            # 10 steps of 0.1
ops.analysis('Static')
ok_grav = run_static(10, 0.1)
ops.loadConst('-time', 0.0)                   # keep gravity, reset pseudo-time

ops.reactions()                               # compute support reactions
sum_Ry = sum(ops.nodeReaction(node2d(0, ax), 2) for ax in range(1, M_BAY + 2))
W_total = w_grav * sum(BAYS) * N_STORY
print('RESULT sumRy_2d %.2f' % sum_Ry)
print('RESULT W_total_2d %.2f' % W_total)
print('RESULT n_nodes_2d %d' % n_nodes_2d)
print('RESULT n_elements_2d %d' % len(elements2d))

# the tag formula finds the axial force of any column directly
col_tag = ele2d(COL, 1, 2)                    # first-story column on axis 2
# eleResponse(tag, 'localForce') -> [N_i, V_i, M_i, N_j, V_j, M_j] in local axes;
# N_i > 0 means the support pushes the column up, i.e. compression
f_loc = ops.eleResponse(col_tag, 'localForce')
if len(f_loc) >= 6:
    print('axial force in column %d = %.1f kN (compression)' % (col_tag, f_loc[0]))

# %% 4. Eigen analysis of the 2D frame
# ---------------------------------------------------------------------------
ops.wipeAnalysis()                            # drop the static analysis objects
lam2d = ops.eigen(3)                          # first 3 eigenvalues (omega^2)
T2d = [2.0 * math.pi / math.sqrt(l) for l in lam2d]
print('2D periods (s):', ', '.join('%.3f' % t for t in T2d))
print('RESULT T1_2d %.3f' % T2d[0])

# collect coordinates (from OpenSees itself) and mode-1 displacements
coords2d = {}
phi1 = {}
for fl in range(N_STORY + 1):
    for ax in range(1, M_BAY + 2):
        n = node2d(fl, ax)
        coords2d[n] = ops.nodeCoord(n)            # [x, y]
        phi1[n] = (ops.nodeEigenvector(n, 1, 1), ops.nodeEigenvector(n, 1, 2))

# %% 5. Plot the 2D model with labels
# ---------------------------------------------------------------------------
COLORS = {COL: 'tab:blue', BEAM: 'tab:green', BRACE: 'tab:red', DAMPER: 'tab:purple'}
fig, ax_ = plt.subplots(figsize=(7, 4.2))
for t, (i, j, k) in elements2d.items():
    xi, yi = coords2d[i]
    xj, yj = coords2d[j]
    ls = '--' if k == DAMPER else '-'
    ax_.plot([xi, xj], [yi, yj], ls, color=COLORS[k], lw=1.6 if k <= BEAM else 1.0)
    xm, ym = 0.5 * (xi + xj), 0.5 * (yi + yj)
    if k == BRACE:
        xm, ym = xi + 0.65 * (xj - xi), yi + 0.65 * (yj - yi)   # off the bay centre
    elif k == DAMPER:
        xm, ym = xi + 0.35 * (xj - xi), yi + 0.35 * (yj - yi)
    ax_.text(xm, ym, str(t), fontsize=5.5, color=COLORS[k], ha='center', va='center',
             bbox=dict(fc='white', ec='none', pad=0.4))
for n, (x, y) in coords2d.items():
    ax_.plot(x, y, 'ko', ms=2.5)
    ax_.text(x - 0.12, y + 0.12, str(n), fontsize=5.5, ha='right', va='bottom',
             bbox=dict(fc='lightyellow', ec='0.6', lw=0.4, pad=0.5))
for ax in range(1, M_BAY + 2):
    x0 = coords2d[node2d(0, ax)][0]
    ax_.plot([x0 - 0.25, x0 + 0.25], [0, 0], 'k-', lw=2.5)   # fixed support mark
for k, name in TYPE_NAME.items():
    ax_.plot([], [], '--' if k == DAMPER else '-', color=COLORS[k],
             label='%s  %dxxx' % (name, k))
ax_.legend(fontsize=6.5, loc='upper left', bbox_to_anchor=(1.0, 1.0))
ax_.set_title('2D frame: node tag = 100*floor + axis', fontsize=9)
ax_.set_xlabel('x (m)'); ax_.set_ylabel('y (m)')
ax_.set_xlim(-1.2, x_axis[-1] + 0.8); ax_.set_ylim(-0.6, y_floor[-1] + 0.6)
ax_.set_aspect('equal')
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'ch03_frame2d_tags.png'), dpi=110, bbox_inches='tight')
plt.close(fig)

# mode 1 of the 2D frame (straight lines between displaced nodes)
umax = max(abs(v[0]) for v in phi1.values()) or 1.0
sc = 0.12 * y_floor[-1] / umax
fig, ax_ = plt.subplots(figsize=(7, 4.2))
for t, (i, j, k) in elements2d.items():
    xi, yi = coords2d[i]; xj, yj = coords2d[j]
    ax_.plot([xi, xj], [yi, yj], color='0.8', lw=1.0)
    ax_.plot([xi + sc * phi1[i][0], xj + sc * phi1[j][0]],
             [yi + sc * phi1[i][1], yj + sc * phi1[j][1]], color=COLORS[k], lw=1.4)
ax_.set_title('2D frame, mode 1: T1 = %.3f s' % T2d[0], fontsize=9)
ax_.set_xlabel('x (m)'); ax_.set_ylabel('y (m)')
ax_.set_aspect('equal')
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'ch03_mode1_2d.png'), dpi=110)
plt.close(fig)

# %% 6. Build a small 3D grid with the 3D convention
# ---------------------------------------------------------------------------
N3 = 2                       # stories
BX = [5.0, 5.0]              # bays along X (m)
BY = [4.0]                   # bays along Y (m)
H3 = [3.2, 3.2]              # story heights (m)
q_floor = 8.0                # floor load used for the masses (kN/m2)
G_c = E_c / 2.4              # shear modulus (kPa), nu = 0.2
J_c = 0.141 * b_c**4         # torsion constant of a square section (m4)
J_b = 0.229 * h_b * b_b**3   # torsion constant, rectangle h/b = 1.67 (approx.)
I_b_weak = h_b * b_b**3 / 12.0

xs = np.concatenate([[0.0], np.cumsum(BX)])
ys = np.concatenate([[0.0], np.cumsum(BY)])
zs = np.concatenate([[0.0], np.cumsum(H3)])
NX, NY = len(xs), len(ys)

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)      # 3D, 6 DOF per node

for fl in range(N3 + 1):
    for iy in range(1, NY + 1):
        for ix in range(1, NX + 1):
            ops.node(node3d(fl, iy, ix), float(xs[ix - 1]), float(ys[iy - 1]), float(zs[fl]))
            if fl == 0:
                ops.fix(node3d(fl, iy, ix), 1, 1, 1, 1, 1, 1)   # all 6 DOFs fixed
            else:
                tx = 0.5 * ((BX[ix - 2] if ix > 1 else 0.0) + (BX[ix - 1] if ix <= len(BX) else 0.0))
                ty = 0.5 * ((BY[iy - 2] if iy > 1 else 0.0) + (BY[iy - 1] if iy <= len(BY) else 0.0))
                m = q_floor * tx * ty / g
                ops.mass(node3d(fl, iy, ix), m, m, 0.0, 0.0, 0.0, 0.0)

# geomTransf(type, tag, vecxz...): vecxz is any vector in the local x-z plane
ops.geomTransf('PDelta', TRANSF['col3d'], 1.0, 0.0, 0.0)    # columns (local x = Z)
ops.geomTransf('Linear', TRANSF['beamX3d'], 0.0, 1.0, 0.0)  # X-beams: local y vertical
ops.geomTransf('Linear', TRANSF['beamY3d'], -1.0, 0.0, 0.0) # Y-beams: local y vertical

elements3d = {}
# elasticBeamColumn(tag, i, j, A, E, G, J, Iy, Iz, transfTag)
for st in range(1, N3 + 1):
    for iy in range(1, NY + 1):
        for ix in range(1, NX + 1):
            j = node3d(st, iy, ix)
            t = ele3d(COL, 0, j)
            ops.element('elasticBeamColumn', t, node3d(st - 1, iy, ix), j,
                        A_c, E_c, G_c, J_c, I_c, I_c, TRANSF['col3d'])
            elements3d[t] = (node3d(st - 1, iy, ix), j, COL)
for fl in range(1, N3 + 1):
    for iy in range(1, NY + 1):
        for ix in range(1, NX):                       # beams along X
            i, j = node3d(fl, iy, ix), node3d(fl, iy, ix + 1)
            t = ele3d(BEAM, 0, i)
            ops.element('elasticBeamColumn', t, i, j, A_b, E_c, G_c, J_b,
                        I_b_weak, I_b, TRANSF['beamX3d'])
            elements3d[t] = (i, j, BEAM)
    for iy in range(1, NY):                           # beams along Y
        for ix in range(1, NX + 1):
            i, j = node3d(fl, iy, ix), node3d(fl, iy + 1, ix)
            t = ele3d(BEAM, 1, i)
            ops.element('elasticBeamColumn', t, i, j, A_b, E_c, G_c, J_b,
                        I_b_weak, I_b, TRANSF['beamY3d'])
            elements3d[t] = (i, j, BEAM)

n_nodes_3d = (N3 + 1) * NX * NY
print('3D model: %d nodes, %d elements' % (n_nodes_3d, len(elements3d)))
print('RESULT n_nodes_3d %d' % n_nodes_3d)
print('RESULT n_elements_3d %d' % len(elements3d))

# %% 7. Eigen analysis of the 3D grid
# ---------------------------------------------------------------------------
lam3d = ops.eigen(3)
T3d = [2.0 * math.pi / math.sqrt(l) for l in lam3d]
# which direction does each mode move? compare roof-corner X and Y components
roof = node3d(N3, 1, 1)
for k in range(3):
    ux = ops.nodeEigenvector(roof, k + 1, 1)
    uy = ops.nodeEigenvector(roof, k + 1, 2)
    print('3D mode %d: T = %.3f s, roof corner |ux| = %.3f, |uy| = %.3f'
          % (k + 1, T3d[k], abs(ux), abs(uy)))
print('RESULT T1_3d %.3f' % T3d[0])
print('RESULT T2_3d %.3f' % T3d[1])

coords3d = {}
for fl in range(N3 + 1):
    for iy in range(1, NY + 1):
        for ix in range(1, NX + 1):
            n = node3d(fl, iy, ix)
            coords3d[n] = ops.nodeCoord(n)

# %% 8. Plot the 3D grid with labels
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(7, 4.2))
a3 = fig.add_subplot(111, projection='3d')
for t, (i, j, k) in elements3d.items():
    pi, pj = coords3d[i], coords3d[j]
    a3.plot([pi[0], pj[0]], [pi[1], pj[1]], [pi[2], pj[2]], color=COLORS[k], lw=1.4)
    pm = [0.5 * (pi[d] + pj[d]) for d in range(3)]
    a3.text(pm[0], pm[1], pm[2], str(t), fontsize=4.5, color=COLORS[k])
for n, p in coords3d.items():
    a3.scatter([p[0]], [p[1]], [p[2]], color='k', s=4)
    a3.text(p[0], p[1], p[2] + 0.15, str(n), fontsize=5, color='k')
a3.set_xlabel('X (m)', fontsize=7); a3.set_ylabel('Y (m)', fontsize=7)
a3.set_zlabel('Z (m)', fontsize=7)
a3.tick_params(labelsize=6)
a3.view_init(elev=22, azim=-58)
a3.set_title('3D grid: node = 10000*floor + 100*axisY + axisX', fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'ch03_frame3d_tags.png'), dpi=110)
plt.close(fig)

ops.wipe()
print('done, figures in', FIGS)
