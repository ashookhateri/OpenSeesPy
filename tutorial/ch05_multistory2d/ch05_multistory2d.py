# -*- coding: utf-8 -*-
"""
Chapter 05 - Multi-story 2D elastic frame (5 stories x 3 bays)
================================================================
Gravity -> modal analysis -> equivalent lateral force (ELF) -> El Centro
linear time history with Rayleigh damping -> story drifts.

UNITS: kN, m, s  (stress in kN/m^2 = kPa, mass in tonne = kN.s^2/m)
g = 9.81 m/s^2

Tag convention (chapter 03):
  node    = 100*floor + axis            (floor 0 = ground, axes 1..4 left->right)
  element = 1000*type + 100*story + axis   (type 1 column, 2 beam)
  column  = 1000 + top node tag         (1xxx, e.g. 1101 = story 1, axis 1)
  beam    = 2000 + left node tag        (2xxx, e.g. 2503 = roof, bay 3)
  transformation tags: 1 columns, 2 beams
"""
import os
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')                      # draw to files, no screen needed
import matplotlib.pyplot as plt
import openseespy.opensees as ops

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, 'figs')
os.makedirs(FIGS, exist_ok=True)
REC = os.path.join(HERE, '..', 'records', 'elCentro_NS_dt0.02_g.txt')
T_START = time.time()
WALL_LIMIT = 150.0                         # s, safety stop so the script never runs > 3 min

# %% 0. Parameters ------------------------------------------------------------
g = 9.81                 # m/s^2
nStory = 5               # number of stories
nBay = 3                 # number of bays
H = 3.2                  # story height (m)
L = 5.0                  # bay width (m)
E = 25.0e6               # concrete modulus (kN/m^2) = 25 GPa
bc, hc = 0.50, 0.50      # column section b x h (m)
bb, hb = 0.35, 0.50      # beam section b x h (m)
Ac, Ic = bc*hc, bc*hc**3/12.0
Ab, Ib = bb*hb, bb*hb**3/12.0
w = 30.0                 # gravity line load on every beam (kN/m), also the seismic weight
C_elf = 0.12             # assumed seismic coefficient V = C*W (design value, illustrative)
zeta = 0.05              # 5 % Rayleigh damping, fixed at modes 1 and 3
nModes = 5               # number of modes to compute
nAxis = nBay + 1

NODES = {}               # tag -> (x, y), kept in Python for plotting
ELES = {}                # tag -> (iNode, jNode, kind)
MASS = {}                # tag -> horizontal mass (t)


COL, BEAM = 1, 2         # element type digits of the chapter-03 convention


def node2d(floor, axis):
    """Node tag = 100*floor + axis (same helper as chapter 03)."""
    return 100*floor + axis


def ele2d(etype, story, axis):
    """Element tag = 1000*type + 100*story + axis (or bay)."""
    return 1000*etype + 100*story + axis


def build_model():
    """Create the elastic frame. Called again before every independent analysis."""
    ops.wipe()                                   # clear any previous model
    ops.model('basic', '-ndm', 2, '-ndf', 3)     # 2D model, 3 DOF per node (ux, uy, rz)
    NODES.clear(); ELES.clear(); MASS.clear()

    # nodes: every floor level x every axis
    for f in range(nStory + 1):
        for a in range(1, nAxis + 1):
            t = node2d(f, a)
            x, y = (a - 1)*L, f*H
            ops.node(t, x, y)                    # node(tag, x, y)
            NODES[t] = (x, y)
    # supports: fixed base (ux, uy, rz restrained = 1)
    for a in range(1, nAxis + 1):
        ops.fix(node2d(0, a), 1, 1, 1)             # fix(tag, ux, uy, rz)

    # lumped horizontal mass from the beam load w (tributary half-bays)
    for f in range(1, nStory + 1):
        for a in range(1, nAxis + 1):
            trib = L/2.0 if a in (1, nAxis) else L          # tributary beam length
            m = w*trib/g                                     # t
            ops.mass(node2d(f, a), m, 0.0, 0.0)                # mass(tag, mx, my, mrz): horizontal only
            MASS[node2d(f, a)] = m

    # geometric transformations: map local element axes to global ones
    ops.geomTransf('Linear', 1)                  # tag 1: columns (chapter-03 table) (linear, no P-Delta in this elastic chapter)
    ops.geomTransf('Linear', 2)                  # tag 2: beams

    # columns 1xxx
    for s in range(1, nStory + 1):
        for a in range(1, nAxis + 1):
            t = ele2d(COL, s, a)            # e.g. 1302 = column under node 302
            i, j = node2d(s - 1, a), node2d(s, a)
            # elasticBeamColumn(tag, iNode, jNode, A, E, Iz, transfTag)
            ops.element('elasticBeamColumn', t, i, j, Ac, E, Ic, 1)
            ELES[t] = (i, j, 'col')
    # beams 2xxx
    for f in range(1, nStory + 1):
        for b in range(1, nBay + 1):
            t = ele2d(BEAM, f, b)           # e.g. 2302 = beam starting at node 302
            i, j = node2d(f, b), node2d(f, b + 1)
            ops.element('elasticBeamColumn', t, i, j, Ab, E, Ib, 2)
            ELES[t] = (i, j, 'beam')


def gravity():
    """Apply w on all beams in 10 load steps and keep it (loadConst)."""
    ops.timeSeries('Linear', 1)                  # load factor grows linearly with pseudo-time
    ops.pattern('Plain', 1, 1)                   # pattern(type, patternTag, tsTag)
    for t, (i, j, k) in ELES.items():
        if k == 'beam':
            # uniform load in the element's local y (local y points up for a left->right beam)
            ops.eleLoad('-ele', t, '-type', '-beamUniform', -w)
    ops.constraints('Plain')                     # only single-point constraints (fix) are used
    ops.numberer('RCM')                          # reverse Cuthill-McKee DOF numbering (small bandwidth)
    ops.system('BandGeneral')                    # banded solver
    ops.test('NormDispIncr', 1.0e-8, 20)         # convergence test: |du| < 1e-8 within 20 iterations
    ops.algorithm('Newton')                      # Newton-Raphson
    ops.integrator('LoadControl', 0.1)           # 10 steps of 0.1
    ops.analysis('Static')                       # static analysis
    ok = ops.analyze(10)                         # run the 10 steps
    ops.loadConst('-time', 0.0)                  # freeze gravity, reset pseudo-time to 0
    return ok


def floor_ux(f):
    """Horizontal displacement of floor f (node on axis 1)."""
    return ops.nodeDisp(node2d(f, 1), 1)


# %% 1. Model ------------------------------------------------------------------
build_model()
W_total = w*nBay*L*nStory                        # total seismic weight (kN)
print('INFO nodes', len(NODES), 'elements', len(ELES), 'W_total', W_total)

# %% 2. Gravity ----------------------------------------------------------------
ok = gravity()
ops.reactions()                                  # compute support reactions
Ry = sum(ops.nodeReaction(node2d(0, a), 2) for a in range(1, nAxis + 1))
print('RESULT Ry_gravity_kN', round(Ry, 1))      # hand check: = w * 15 m * 5 floors = 2250 kN

# %% 3. Modal analysis --------------------------------------------------------
lam = np.array(ops.eigen(nModes))                # eigenvalues omega^2 (default genBandArpack solver)
omega = np.sqrt(lam)
T = 2*np.pi/omega
print('INFO periods (s):', np.round(T, 4))
mtags = sorted(MASS)
mvec = np.array([MASS[t] for t in mtags])
M_tot = mvec.sum()
Gamma, Meff, phi_roof = [], [], []
for k in range(1, nModes + 1):
    phi = np.array([ops.nodeEigenvector(t, k, 1) for t in mtags])   # x-component of mode k
    G = (mvec @ phi)/(mvec @ phi**2)            # participation factor for ground motion in x
    Gamma.append(G)
    Meff.append((mvec @ phi)**2/(mvec @ phi**2)/M_tot)
    phi_roof.append(ops.nodeEigenvector(node2d(nStory, 1), k, 1))
print('INFO effective mass ratios:', np.round(Meff, 3), 'sum', round(sum(Meff), 3))
# mode shapes of all nodes for plotting (first 3 modes)
shapes = {k: {t: (ops.nodeEigenvector(t, k, 1), ops.nodeEigenvector(t, k, 2)) for t in NODES}
          for k in (1, 2, 3)}
print('RESULT T1_s', round(T[0], 3))
print('RESULT T2_s', round(T[1], 3))
print('RESULT T3_s', round(T[2], 3))
print('RESULT Meff1_pct', round(100*Meff[0], 1))

# %% 4. Equivalent lateral force ---------------------------------------------
V = C_elf*W_total                                # design base shear
k_exp = 1.0 if T[0] <= 0.5 else (2.0 if T[0] >= 2.5 else 1.0 + (T[0] - 0.5)/2.0)
Wf = np.full(nStory, w*nBay*L)                   # floor weights
hf = H*np.arange(1, nStory + 1)                  # floor heights
Fx = V*Wf*hf**k_exp/np.sum(Wf*hf**k_exp)         # vertical distribution of V
u_grav = np.array([floor_ux(f) for f in range(nStory + 1)])
ops.timeSeries('Linear', 2)
ops.pattern('Plain', 2, 2)
for f in range(1, nStory + 1):
    for a in range(1, nAxis + 1):
        ops.load(node2d(f, a), Fx[f - 1]/nAxis, 0.0, 0.0)   # load(node, Fx, Fy, Mz): share among 4 nodes
ops.integrator('LoadControl', 1.0)               # linear model: one full step is exact
ops.analyze(1)
u_elf = np.array([floor_ux(f) for f in range(nStory + 1)]) - u_grav
drift_elf = np.diff(u_elf)/H                     # interstory drift ratios
ops.reactions()
V_check = -sum(ops.nodeReaction(node2d(0, a), 1) for a in range(1, nAxis + 1))
print('INFO ELF k =', round(k_exp, 3), 'forces', np.round(Fx, 1), 'base shear check', round(V_check, 1))
print('RESULT V_ELF_kN', round(V, 1))
print('RESULT roof_ELF_mm', round(1000*u_elf[-1], 2))
print('RESULT drift_ELF_max_pct', round(100*np.max(np.abs(drift_elf)), 3))

# %% 5. Linear time history (El Centro, 5 % Rayleigh) -------------------------
acc = np.loadtxt(REC)                            # ground acceleration in g
dt_rec = 0.02
nPts = len(acc)
build_model()                                    # fresh model
gravity()
ops.wipeAnalysis()                               # remove the static analysis objects
w1, w3 = omega[0], omega[2]
a0 = zeta*2.0*w1*w3/(w1 + w3)                    # mass-proportional coefficient
a1 = zeta*2.0/(w1 + w3)                          # stiffness-proportional coefficient
ops.rayleigh(a0, a1, 0.0, 0.0)                   # rayleigh(alphaM, betaK(current), betaKinit, betaKcomm)
ops.timeSeries('Path', 3, '-dt', dt_rec, '-values', *acc.tolist(), '-factor', g)  # record in m/s^2
ops.pattern('UniformExcitation', 3, 1, '-accel', 3)   # base acceleration in direction 1 (x)
ops.constraints('Plain')
ops.numberer('RCM')
ops.system('BandGeneral')
ops.test('NormDispIncr', 1.0e-8, 20)
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)             # average acceleration: unconditionally stable
ops.analysis('Transient')

ALGOS = [('Newton',), ('ModifiedNewton', '-initial'), ('KrylovNewton',)]


def step(dt, level=0):
    """One step of size dt; try 3 algorithms, then split into two half steps (max 4 levels)."""
    for alg in ALGOS:
        ops.algorithm(*alg)
        if ops.analyze(1, dt) == 0:
            if alg is not ALGOS[0]:
                ops.algorithm(*ALGOS[0])
            return True
    if level >= 4:
        return False
    return step(dt/2.0, level + 1) and step(dt/2.0, level + 1)


dt = dt_rec                                      # 0.02 s: T1/dt is large, fine for a linear model
nSteps = nPts                                    # 1559 steps = 31.2 s
base_cols = [ele2d(COL, 1, a) for a in range(1, nAxis + 1)]
t_hist, roof_hist, vb_hist = [0.0], [0.0], [0.0]
peak_drift = np.zeros(nStory)
u0 = np.array([floor_ux(f) for f in range(nStory + 1)])   # gravity offset (tiny)
status = 'completed'
for n in range(nSteps):
    if not step(dt):
        status = 'failed at t=%.2f s' % ops.getTime(); break
    if time.time() - T_START > WALL_LIMIT:
        status = 'stopped (wall-clock limit) at t=%.2f s' % ops.getTime(); break
    u = np.array([floor_ux(f) for f in range(nStory + 1)]) - u0
    peak_drift = np.maximum(peak_drift, np.abs(np.diff(u))/H)
    # base shear = sum of global x end-forces of the first-story columns at their base node
    vb = sum(ops.eleForce(t)[0] for t in base_cols)
    t_hist.append(ops.getTime()); roof_hist.append(u[-1]); vb_hist.append(vb)
print('INFO time history', status)
roof_hist = np.array(roof_hist); vb_hist = np.array(vb_hist)

# hand check: first-mode SDOF estimate u_roof = Gamma1*phi_roof1*Sd(T1, 5 %)


def sdof_peak(Tn, z, ag, dt):
    """Peak relative displacement of a linear SDOF (Newmark average acceleration)."""
    wn = 2*np.pi/Tn
    k, c, m = wn**2, 2*z*wn, 1.0
    beta, gam = 0.25, 0.5
    u = v = 0.0
    a = -ag[0]
    khat = k + gam/(beta*dt)*c + m/(beta*dt**2)
    umax = 0.0
    for p in -ag[1:]:
        ph = (p + m*(u/(beta*dt**2) + v/(beta*dt) + a*(1/(2*beta) - 1))
              + c*(gam/(beta*dt)*u + (gam/beta - 1)*v + dt*(gam/(2*beta) - 1)*a))
        un = ph/khat
        vn = gam/(beta*dt)*(un - u) + (1 - gam/beta)*v + dt*(1 - gam/(2*beta))*a
        an = (un - u)/(beta*dt**2) - v/(beta*dt) - (1/(2*beta) - 1)*a
        u, v, a = un, vn, an
        umax = max(umax, abs(u))
    return umax


Sd1 = sdof_peak(T[0], zeta, acc*g, dt_rec)
roof_est = abs(Gamma[0]*phi_roof[0])*Sd1
print('INFO Sd(T1) =', round(Sd1*1000, 2), 'mm, Gamma1*phi_roof1 =', round(abs(Gamma[0]*phi_roof[0]), 3))
print('RESULT roof_TH_max_mm', round(1000*np.max(np.abs(roof_hist)), 2))
print('RESULT roof_SDOF_estimate_mm', round(1000*roof_est, 2))
print('RESULT Vbase_TH_max_kN', round(np.max(np.abs(vb_hist)), 1))
print('RESULT drift_TH_max_pct', round(100*np.max(peak_drift), 3))

# %% 6. Figures ---------------------------------------------------------------


def draw_frame(ax, disp=None, scale=0.0, labels=False, color='k', lw=1.5):
    for t, (i, j, kind) in ELES.items():
        (xi, yi), (xj, yj) = NODES[i], NODES[j]
        if disp is not None:
            xi += scale*disp[i][0]; yi += scale*disp[i][1]
            xj += scale*disp[j][0]; yj += scale*disp[j][1]
        ax.plot([xi, xj], [yi, yj], color=color, lw=lw)
        if labels:
            ax.text((xi + xj)/2, (yi + yj)/2, str(t), fontsize=5.5, ha='center', va='center',
                    color='tab:blue' if kind == 'col' else 'tab:red',
                    bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.8))
    if labels:
        for t, (x, y) in NODES.items():
            ax.plot(x, y, 'ko', ms=2.5)
            ax.text(x + 0.12, y + 0.12, str(t), fontsize=5.5, color='k')


fig, ax = plt.subplots(figsize=(7, 4.2))
draw_frame(ax, labels=True)
for a in range(1, nAxis + 1):
    x, y = NODES[node2d(0, a)]
    ax.plot([x - 0.3, x + 0.3], [y, y], 'k-', lw=3)
ax.set_title('Ch05 model: 5 stories x 3 bays (node, column 1xxx, beam 2xxx tags)', fontsize=9)
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)'); ax.set_aspect('equal')
ax.set_xlim(-1, nBay*L + 1.0)
fig.tight_layout(); fig.savefig(os.path.join(FIGS, 'ch05_model.png'), dpi=110); plt.close(fig)

fig, axs = plt.subplots(1, 3, figsize=(7, 4.2), sharey=True)
for k, ax in zip((1, 2, 3), axs):
    sh = shapes[k]
    mx = max(abs(v[0]) for v in sh.values()) or 1.0
    draw_frame(ax, color='0.8', lw=1.0)
    draw_frame(ax, disp=sh, scale=2.5/mx, color='tab:blue', lw=1.6)
    ax.set_title('Mode %d, T = %.3f s' % (k, T[k - 1]), fontsize=9)
    ax.set_aspect('equal'); ax.set_xlabel('x (m)')
axs[0].set_ylabel('y (m)')
fig.tight_layout(); fig.savefig(os.path.join(FIGS, 'ch05_modes.png'), dpi=110); plt.close(fig)

fig, (a1x, a2x) = plt.subplots(2, 1, figsize=(7, 4.2), sharex=True)
a1x.plot(np.arange(nPts)*dt_rec, acc, 'k', lw=0.7)
a1x.set_ylabel('ag (g)'); a1x.set_title('El Centro 1940 NS and roof response', fontsize=9)
a2x.plot(t_hist, 1000*roof_hist, 'tab:blue', lw=0.9)
a2x.axhline(1000*roof_est, color='tab:red', ls='--', lw=0.8, label='SDOF estimate')
a2x.axhline(-1000*roof_est, color='tab:red', ls='--', lw=0.8)
a2x.set_xlabel('time (s)'); a2x.set_ylabel('roof u (mm)'); a2x.legend(fontsize=7)
for a_ in (a1x, a2x):
    a_.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIGS, 'ch05_roof_th.png'), dpi=110); plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4.2))
st = np.arange(1, nStory + 1)
ax.plot(100*np.abs(drift_elf), st, 'o-', label='ELF (C = %.2f)' % C_elf)
ax.plot(100*peak_drift, st, 's-', label='El Centro, peak (elastic)')
ax.set_xlabel('interstory drift ratio (%)'); ax.set_ylabel('story')
ax.set_yticks(st); ax.grid(alpha=0.3); ax.legend(fontsize=8)
ax.set_title('Story drift profile', fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIGS, 'ch05_drifts.png'), dpi=110); plt.close(fig)

print('INFO run time %.1f s' % (time.time() - T_START))
ops.wipe()
