# -*- coding: utf-8 -*-
"""
Chapter 01 - Elastic cantilever beam
=====================================
Units used everywhere in this script:  kN, m, s
  -> force kN, length m, time s, stress kN/m^2 (= kPa), mass t (= kN*s^2/m)
  g = 9.81 m/s^2

What the script does
  1. Builds a 2D elastic cantilever (fixed at x = 0, free at x = L) made of
     elasticBeamColumn elements and draws a sketch of the model.
  2. Case A: point load P at the free end  -> tip deflection vs PL^3/(3EI).
  3. Case B: uniform load w on the whole span -> tip deflection vs wL^4/(8EI).
  4. Case C: modal analysis (eigen) -> first two bending frequencies vs
     Euler-Bernoulli theory.
Every OpenSees command is commented with what it does and what its
arguments mean, so this file doubles as a first "anatomy of a model".
"""
import os
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')                     # no screen needed: draw straight to PNG files
import matplotlib.pyplot as plt
import openseespy.opensees as ops

# Folder for figures: figs/ next to this script
HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, 'figs')
os.makedirs(FIGDIR, exist_ok=True)

def savefig(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'ch01_' + name + '.png'), dpi=110)
    plt.close(fig)

# %% 0. Parameters (kN, m, s)
g = 9.81                     # m/s^2
L = 4.0                      # span of the cantilever [m]
nEle = 10                    # number of elements along the span
E = 200.0e6                  # Young's modulus of steel [kN/m^2] (200 GPa)
A = 53.8e-4                  # cross-section area of an IPE 300 [m^2]
Iz = 8356.0e-8               # second moment of area (strong axis) of an IPE 300 [m^4]
rho = 7.85                   # steel density [t/m^3]
mbar = rho * A               # mass per unit length [t/m]
P = 10.0                     # tip point load, acting downwards [kN]
w = 5.0                      # uniform distributed load, acting downwards [kN/m]
EI = E * Iz                  # flexural stiffness [kN*m^2]

TIP = nEle + 1               # tag of the free-end node (nodes are numbered 1..nEle+1)
BASE = 1                     # tag of the fixed node


# %% 1. Model
def build_model(with_mass=False):
    """Build the cantilever from scratch: model -> nodes -> fix -> transf -> elements."""
    ops.wipe()                                    # delete any previous model, loads and analysis objects
    ops.model('basic', '-ndm', 2, '-ndf', 3)      # 2D model (ndm=2) with 3 DOFs per node: ux, uy, rotation

    # Nodes: node(tag, x, y). Equal spacing from x = 0 to x = L along the global X axis.
    for i in range(nEle + 1):
        ops.node(i + 1, i * L / nEle, 0.0)

    # Boundary condition: fix(nodeTag, ux, uy, rz) with 1 = restrained, 0 = free.
    ops.fix(BASE, 1, 1, 1)                        # fully fixed support at x = 0

    # Geometric transformation: maps element local axes to global axes.
    # 'Linear' = small displacements, no P-Delta (correct for this linear problem).
    ops.geomTransf('Linear', 1)                   # geomTransf(type, transfTag)

    # Elements: elasticBeamColumn(tag, iNode, jNode, A, E, Iz, transfTag, ...)
    # Optional '-mass', mbar gives the mass per unit length and '-cMass' asks for a
    # consistent (not lumped) mass matrix, which is more accurate for few elements.
    for e in range(nEle):
        if with_mass:
            ops.element('elasticBeamColumn', e + 1, e + 1, e + 2, A, E, Iz, 1,
                        '-mass', mbar, '-cMass')
        else:
            ops.element('elasticBeamColumn', e + 1, e + 1, e + 2, A, E, Iz, 1)


def static_analysis_objects():
    """The 'analysis objects': how the equations are numbered, stored, solved and stepped."""
    ops.wipeAnalysis()                            # remove analysis objects of a previous run
    ops.constraints('Plain')                      # how fixities are enforced: 'Plain' is enough for simple fix commands
    ops.numberer('RCM')                           # DOF numbering: Reverse Cuthill-McKee keeps the matrix band narrow
    ops.system('BandGeneral')                     # store/solve K as a banded general matrix
    ops.test('NormDispIncr', 1.0e-10, 20)         # convergence test: |du| < 1e-10 m, max 20 iterations
    ops.algorithm('Newton')                       # solution algorithm: Newton-Raphson (1 iteration for linear models)
    ops.integrator('LoadControl', 1.0)            # each analysis step adds 1.0 to the load factor
    ops.analysis('Static')                        # static analysis (no inertia, no damping)


def analyze_static(nsteps=1):
    """Run nsteps; if Newton fails, retry the step with other algorithms. Returns True on success."""
    for _ in range(nsteps):
        ok = ops.analyze(1)                       # analyze(numSteps): returns 0 on success
        if ok != 0:
            for algo in (('ModifiedNewton', '-initial'), ('KrylovNewton',)):
                ops.algorithm(*algo)
                ok = ops.analyze(1)
                if ok == 0:
                    break
            ops.algorithm('Newton')
        if ok != 0:
            print('WARNING: static analysis failed')
            return False
    return True


# ---- Model sketch (figure 1), drawn from the node coordinates stored in OpenSees
build_model()
xy = np.array([ops.nodeCoord(n) for n in ops.getNodeTags()])   # nodeCoord(tag) -> [x, y]
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(xy[:, 0], xy[:, 1], '-', color='0.3', lw=3, zorder=1)
ax.plot(xy[:, 0], xy[:, 1], 'o', color='tab:blue', ms=6, zorder=2)
for n in ops.getNodeTags():
    x, y = ops.nodeCoord(n)
    ax.annotate(str(n), (x, y), textcoords='offset points', xytext=(0, -15), ha='center', fontsize=8)
for e in ops.getEleTags():
    ni, nj = ops.eleNodes(e)
    xm = 0.5 * (ops.nodeCoord(ni)[0] + ops.nodeCoord(nj)[0])
    ax.annotate('[%d]' % e, (xm, 0), textcoords='offset points', xytext=(0, -28), ha='center',
                fontsize=7, color='tab:red')
# fixed support: wall with hatching
ax.plot([0, 0], [-0.5, 0.5], 'k-', lw=2)
for yy in np.linspace(-0.5, 0.4, 8):
    ax.plot([-0.15, 0], [yy, yy + 0.1], 'k-', lw=0.8)
# tip load (case A) and distributed load (case B)
ax.annotate('', xy=(L, 0.05), xytext=(L, 0.9), arrowprops=dict(arrowstyle='->', color='tab:red', lw=2))
ax.text(L + 0.06, 0.75, 'P = %.0f kN (case A)' % P, color='tab:red', fontsize=9)
for xx in np.linspace(0.2, L - 0.2, 9):
    ax.annotate('', xy=(xx, 0.05), xytext=(xx, 0.45), arrowprops=dict(arrowstyle='->', color='tab:green', lw=1))
ax.plot([0.2, L - 0.2], [0.45, 0.45], color='tab:green', lw=1)
ax.text(0.4, 0.52, 'w = %.0f kN/m (case B)' % w, color='tab:green', fontsize=9)
ax.set_xlim(-0.4, L + 1.3)
ax.set_ylim(-0.9, 1.1)
ax.set_aspect('equal')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title('Cantilever: %d elasticBeamColumn elements, L = %.1f m (IPE 300)\nnode tags in black, [element tags] in red' % (nEle, L), fontsize=10)
savefig(fig, 'model')


# %% 2. Case A - tip point load
build_model()
ops.timeSeries('Linear', 1)                       # load factor grows linearly with pseudo-time: lambda = t
ops.pattern('Plain', 1, 1)                        # pattern(type, patternTag, tsTag): a group of loads scaled by series 1
ops.load(TIP, 0.0, -P, 0.0)                       # load(nodeTag, Fx, Fy, Mz): downward force at the free end
static_analysis_objects()
analyze_static(1)                                 # one step of LoadControl 1.0 -> full load applied

xA = np.array([ops.nodeCoord(n)[0] for n in range(1, nEle + 2)])
vA = np.array([ops.nodeDisp(n, 2) for n in range(1, nEle + 2)])   # nodeDisp(tag, dof): dof 2 = uy
tipA = -ops.nodeDisp(TIP, 2)                      # downward tip deflection [m], made positive
ops.reactions()                                   # ask OpenSees to compute support reactions
MbaseA = ops.nodeReaction(BASE, 3)                # nodeReaction(tag, dof): dof 3 = moment at the support
tipA_th = P * L**3 / (3.0 * EI)                   # Euler-Bernoulli: delta = P L^3 / (3 E I)


# %% 3. Case B - uniformly distributed load
build_model()
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
# eleLoad: a load applied along elements. '-beamUniform', Wy acts in the LOCAL y direction.
# For an element drawn from left to right, local y points up (global +Y), so a downward
# load is Wy = -w.
ops.eleLoad('-ele', *range(1, nEle + 1), '-type', '-beamUniform', -w)
static_analysis_objects()
analyze_static(1)

vB = np.array([ops.nodeDisp(n, 2) for n in range(1, nEle + 2)])
tipB = -ops.nodeDisp(TIP, 2)
tipB_th = w * L**4 / (8.0 * EI)                   # delta = w L^4 / (8 E I)

# Figure 2: deflected shapes against the closed-form elastic curves
xs = np.linspace(0, L, 200)
vA_th = P * xs**2 * (3 * L - xs) / (6 * EI)
vB_th = w * xs**2 * (6 * L**2 - 4 * L * xs + xs**2) / (24 * EI)
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(xs, -vA_th * 1e3, '-', color='tab:red', lw=1.2, label='theory, tip load P')
ax.plot(xA, vA * 1e3, 'o', color='tab:red', label='OpenSees nodes, tip load P')
ax.plot(xs, -vB_th * 1e3, '-', color='tab:green', lw=1.2, label='theory, uniform load w')
ax.plot(xA, vB * 1e3, 's', color='tab:green', label='OpenSees nodes, uniform load w')
ax.axhline(0, color='0.5', lw=0.8)
ax.set_xlabel('x [m]')
ax.set_ylabel('vertical displacement uy [mm]')
ax.set_title('Deflected shape of the cantilever')
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
savefig(fig, 'deflection')


# %% 4. Case C - modal analysis
build_model(with_mass=True)
ops.wipeAnalysis()
nModes = 4                                        # ask for a few modes: an axial mode sits among the bending ones
try:
    lam = ops.eigen(nModes)                       # eigen(numModes): returns eigenvalues lambda = omega^2 [rad^2/s^2]
except Exception:
    lam = []
if len(lam) < nModes or min(lam) <= 0:
    lam = ops.eigen('-fullGenLapack', nModes)     # fallback solver: dense, slower but very robust

modes = []                                        # (frequency, shape) of the bending modes only
for k in range(1, nModes + 1):
    ux = np.array([ops.nodeEigenvector(n, k, 1) for n in range(1, nEle + 2)])  # nodeEigenvector(tag, mode, dof)
    uy = np.array([ops.nodeEigenvector(n, k, 2) for n in range(1, nEle + 2)])
    f = math.sqrt(lam[k - 1]) / (2 * math.pi)     # f = omega / (2 pi) [Hz]
    if np.max(np.abs(uy)) > np.max(np.abs(ux)):   # bending mode: vertical motion dominates
        modes.append((f, uy / uy[-1]))            # normalise to +1 at the tip

betaL = [1.875104, 4.694091]                      # roots of 1 + cos(bL) cosh(bL) = 0
f_th = [(b**2) / (2 * math.pi) * math.sqrt(EI / (mbar * L**4)) for b in betaL]

# Figure 3: first two bending mode shapes against the analytical ones
def phi_th(bl, x):
    b = bl / L
    s = (math.cosh(bl) + math.cos(bl)) / (math.sinh(bl) + math.sin(bl))
    return (np.cosh(b * x) - np.cos(b * x)) - s * (np.sinh(b * x) - np.sin(b * x))

fig, ax = plt.subplots(figsize=(7, 4.2))
cols = ['tab:blue', 'tab:orange']
for k in range(min(2, len(modes))):
    ph = phi_th(betaL[k], xs)
    ax.plot(xs, ph / ph[-1], '-', color=cols[k], lw=1.2,
            label='theory mode %d, f = %.1f Hz' % (k + 1, f_th[k]))
    ax.plot(xA, modes[k][1], 'o', color=cols[k],
            label='OpenSees mode %d, f = %.1f Hz' % (k + 1, modes[k][0]))
ax.axhline(0, color='0.5', lw=0.8)
ax.set_xlabel('x [m]')
ax.set_ylabel('normalised mode shape')
ax.set_title('Bending mode shapes of the cantilever')
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
savefig(fig, 'modes')


# %% 5. Results
print('RESULT', 'tipA_mm', round(float(tipA * 1e3), 3))
print('RESULT', 'tipA_theory_mm', round(float(tipA_th * 1e3), 3))
print('RESULT', 'Mbase_A_kNm', round(float(MbaseA), 2))
print('RESULT', 'tipB_mm', round(float(tipB * 1e3), 3))
print('RESULT', 'tipB_theory_mm', round(float(tipB_th * 1e3), 3))
if len(modes) >= 2:
    print('RESULT', 'f1_Hz', round(float(modes[0][0]), 2))
    print('RESULT', 'f1_theory_Hz', round(float(f_th[0]), 2))
    print('RESULT', 'f2_Hz', round(float(modes[1][0]), 2))
    print('RESULT', 'f2_theory_Hz', round(float(f_th[1]), 2))
else:
    print('WARNING: fewer than two bending modes were found')
ops.wipe()
