# -*- coding: utf-8 -*-
"""
Chapter 06 - Four-story steel concentrically braced frame (CBF)
================================================================
Fiber W-section columns and beams (Steel02, forceBeamColumn), buckling braces
made of two corotational fiber elements with an initial out-of-straightness,
a leaning (P-Delta) column for the gravity load carried elsewhere in the
building. Gravity -> modal -> pushover -> nonlinear time history (Northridge).

UNITS: kN, m, s  (stress in kN/m^2 = kPa, mass in tonne = kN.s^2/m)
g = 9.81 m/s^2

Tag convention (chapter 03):
  node     = 100*floor + axis               floor 0 = ground, axes 1..4 frame, axis 5 leaning column
  element  = 1000*type + 100*story + axis   type 1 column, 2 beam, 3 brace
             column 1xxx = 1000 + top node tag; beam / link 2xxx = 2000 + left node tag
             brace  3xxx = 3000 + 100*story + bay   (lower half of the brace)
                    3xxx + 50                        (upper half of the brace)
  extra nodes (this chapter): brace mid-point = 9000 + 100*story + bay
  materials: 1 frame steel, 2 brace steel, 23 rigid elastic (links, leaning column)
  sections : 101/102 columns, 201 beams, 301/302 braces
  transformations: 1 columns (PDelta), 2 beams (Linear), 3 braces (Corotational)
"""
import os
import time
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')                      # draw to files, no screen needed
import matplotlib.pyplot as plt
import openseespy.opensees as ops

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, 'figs')
os.makedirs(FIGS, exist_ok=True)
REC = os.path.join(HERE, '..', 'records', 'Northridge_LOS270_dt0.01_g.txt')
T_START = time.time()
WALL_PUSH = 70.0                           # s, wall-clock budget for the pushover
WALL_TOTAL = 165.0                         # s, total budget: the script never runs > 3 min

# %% 0. Parameters ------------------------------------------------------------
g = 9.81
nStory, nBay = 4, 3
H, L = 3.5, 5.0                            # story height, bay width (m)
BRACED_BAYS = (1, 3)                       # single diagonals in the two outer bays (opposite directions)
LEAN = 5                                   # axis number of the leaning column
Es = 2.0e8                                 # steel modulus (kPa)
Fy_frame, Fy_brace = 345.0e3, 317.0e3      # yield stresses (kPa)
b_hard = 0.01                              # strain-hardening ratio
R0, cR1, cR2 = 18.0, 0.925, 0.15           # Steel02 transition parameters (Menegotto-Pinto)
# W sections (d, bf, tf, tw) in m
W360x134 = (0.356, 0.369, 0.0180, 0.0112)  # columns, stories 1-2
W310x97 = (0.308, 0.305, 0.0154, 0.0099)   # columns, stories 3-4
W410x60 = (0.407, 0.178, 0.0128, 0.0077)   # beams, all floors
# square HSS braces (b, t) in m
HSS160x10 = (0.160, 0.010)                 # braces, stories 1-2
HSS140x10 = (0.140, 0.010)                 # braces, stories 3-4
w_beam = 20.0                              # gravity line load on the frame beams (kN/m)
P_lean = 1300.0                            # gravity carried by the leaning column, per floor (kN)
W_floor = w_beam*nBay*L + P_lean           # seismic weight per floor = 1600 kN
W_total = W_floor*nStory                   # 6400 kN
imp = 1.0/500.0                            # brace initial out-of-straightness (fraction of length)
zeta = 0.02                                # 2 % Rayleigh damping (modes 1 and 3)
SF = 1.0                                   # ground-motion scale factor
roof_target = 0.015*nStory*H               # pushover target: 1.5 % roof drift
dU = 0.001                                 # pushover displacement step (m)
nAxis = nBay + 1

NODES, ELES, MASS = {}, {}, {}             # Python copies for plotting / post-processing
BR = {}                                    # brace data: lower tag -> dict


def node2d(floor, axis):
    return 100*floor + axis


def ele2d(etype, story, axis):
    return 1000*etype + 100*story + axis


def add_node(t, x, y):
    ops.node(t, x, y)                      # node(tag, x, y)
    NODES[t] = (x, y)


def W_section(tag, mat, d, bf, tf, tw):
    """Fiber W section, bending about the strong axis (local y = depth)."""
    ops.section('Fiber', tag)              # start a fiber section with this tag
    # patch('rect', matTag, nFibersY, nFibersZ, yI, zI, yJ, zJ): yI,zI lower-left, yJ,zJ upper-right
    ops.patch('rect', mat, 3, 1, d/2 - tf, -bf/2, d/2, bf/2)         # top flange
    ops.patch('rect', mat, 3, 1, -d/2, -bf/2, -d/2 + tf, bf/2)       # bottom flange
    ops.patch('rect', mat, 10, 1, -d/2 + tf, -tw/2, d/2 - tf, tw/2)  # web


def HSS_section(tag, mat, b, t):
    """Fiber square hollow section (sharp corners)."""
    ops.section('Fiber', tag)
    ops.patch('rect', mat, 2, 1, b/2 - t, -b/2, b/2, b/2)            # top wall
    ops.patch('rect', mat, 2, 1, -b/2, -b/2, -b/2 + t, b/2)          # bottom wall
    ops.patch('rect', mat, 10, 1, -b/2 + t, -b/2, b/2 - t, -b/2 + t)  # left wall
    ops.patch('rect', mat, 10, 1, -b/2 + t, b/2 - t, b/2 - t, b/2)    # right wall


def hss_props(b, t):
    A = b**2 - (b - 2*t)**2
    I = (b**4 - (b - 2*t)**4)/12.0
    return A, I


def build_model():
    ops.wipe()                                       # clear memory
    ops.model('basic', '-ndm', 2, '-ndf', 3)         # 2D, 3 DOF per node (ux, uy, rz)
    for d in (NODES, ELES, MASS, BR):
        d.clear()

    # ---- nodes, supports and masses
    for f in range(nStory + 1):
        for a in range(1, nAxis + 1):
            add_node(node2d(f, a), (a - 1)*L, f*H)
        add_node(node2d(f, LEAN), nBay*L + L, f*H)   # leaning column line, one bay to the right
    for a in range(1, nAxis + 1):
        ops.fix(node2d(0, a), 1, 1, 1)               # fix(tag, ux, uy, rz): fixed column bases
    ops.fix(node2d(0, LEAN), 1, 1, 1)                # leaning column base (a truss: acts as a pin)
    for f in range(1, nStory + 1):
        ops.fix(node2d(f, LEAN), 0, 0, 1)            # truss-only nodes: restrain the free rotation
        for a in range(1, nAxis + 1):
            m = w_beam*(L/2 if a in (1, nAxis) else L)/g
            ops.mass(node2d(f, a), m, 0.0, 0.0)      # mass(tag, mx, my, mrz): horizontal mass only
            MASS[node2d(f, a)] = m
        ops.mass(node2d(f, LEAN), P_lean/g, 0.0, 0.0)
        MASS[node2d(f, LEAN)] = P_lean/g

    # ---- materials
    # Steel02(tag, Fy, E0, b, R0, cR1, cR2): Giuffre-Menegotto-Pinto steel with kinematic hardening
    ops.uniaxialMaterial('Steel02', 1, Fy_frame, Es, b_hard, R0, cR1, cR2)
    ops.uniaxialMaterial('Steel02', 2, Fy_brace, Es, b_hard, R0, cR1, cR2)
    ops.uniaxialMaterial('Elastic', 23, Es)          # Elastic(tag, E): links and leaning column

    # ---- sections and integration rules
    W_section(101, 1, *W360x134)
    W_section(102, 1, *W310x97)
    W_section(201, 1, *W410x60)
    HSS_section(301, 2, *HSS160x10)
    HSS_section(302, 2, *HSS140x10)
    for s in (101, 102, 201, 301, 302):
        # beamIntegration('Lobatto', tag, secTag, nIP): 5 Gauss-Lobatto points (ends included)
        ops.beamIntegration('Lobatto', s, s, 5)

    # ---- geometric transformations
    ops.geomTransf('PDelta', 1)                      # columns: P-Delta (large axial force)
    ops.geomTransf('Linear', 2)                      # beams
    ops.geomTransf('Corotational', 3)                # braces: large displacements -> buckling

    # ---- columns 1xxx (forceBeamColumn: exact for the given section behaviour, few elements)
    for s in range(1, nStory + 1):
        sec = 101 if s <= 2 else 102
        for a in range(1, nAxis + 1):
            t, i, j = ele2d(1, s, a), node2d(s - 1, a), node2d(s, a)
            # forceBeamColumn(tag, iNode, jNode, transfTag, integrationTag)
            ops.element('forceBeamColumn', t, i, j, 1, sec)
            ELES[t] = (i, j, 'col')
        # leaning column: corotTruss(tag, iNode, jNode, A, matTag), carries P_lean, adds P-Delta
        t, i, j = ele2d(1, s, LEAN), node2d(s - 1, LEAN), node2d(s, LEAN)
        ops.element('corotTruss', t, i, j, 0.1, 23)
        ELES[t] = (i, j, 'lean')

    # ---- beams 2xxx and rigid links to the leaning column
    for f in range(1, nStory + 1):
        for b in range(1, nBay + 1):
            t, i, j = ele2d(2, f, b), node2d(f, b), node2d(f, b + 1)
            ops.element('forceBeamColumn', t, i, j, 2, 201)
            ELES[t] = (i, j, 'beam')
        t, i, j = ele2d(2, f, nAxis), node2d(f, nAxis), node2d(f, LEAN)
        ops.element('Truss', t, i, j, 0.1, 23)       # Truss(tag, iNode, jNode, A, matTag): axially stiff link
        ELES[t] = (i, j, 'link')

    # ---- braces 3xxx: two corotational fiber elements with a bowed mid node
    for s in range(1, nStory + 1):
        sec = 301 if s <= 2 else 302
        for b in BRACED_BAYS:
            if b == 1:                                  # bay 1: bottom-left -> top-right
                i, j = node2d(s - 1, b), node2d(s, b + 1)
            else:                                       # bay 3: bottom-right -> top-left
                i, j = node2d(s - 1, b + 1), node2d(s, b)
            (xi, yi), (xj, yj) = NODES[i], NODES[j]
            Lb = math.hypot(xj - xi, yj - yi)
            cx, cy = (xj - xi)/Lb, (yj - yi)/Lb
            mid = 9000 + 100*s + b
            e0 = imp*Lb                                 # out-of-straightness at mid-length
            add_node(mid, (xi + xj)/2 - e0*cy, (yi + yj)/2 + e0*cx)   # offset normal to the chord
            t1, t2 = ele2d(3, s, b), ele2d(3, s, b) + 50
            ops.element('forceBeamColumn', t1, i, mid, 3, sec)
            ops.element('forceBeamColumn', t2, mid, j, 3, sec)
            ELES[t1] = (i, mid, 'brace'); ELES[t2] = (mid, j, 'brace')
            BR[t1] = dict(i=i, j=j, mid=mid, L0=Lb)


def gravity():
    ops.timeSeries('Linear', 1)                      # load factor = pseudo-time
    ops.pattern('Plain', 1, 1)                       # pattern(type, tag, tsTag)
    for t, (i, j, k) in ELES.items():
        if k == 'beam':
            ops.eleLoad('-ele', t, '-type', '-beamUniform', -w_beam)   # downward line load
    for f in range(1, nStory + 1):
        ops.load(node2d(f, LEAN), 0.0, -P_lean, 0.0)  # load(node, Fx, Fy, Mz)
    ops.constraints('Plain')                         # only fix() constraints in this model
    ops.numberer('RCM')                              # bandwidth-reducing DOF numbering
    ops.system('BandGeneral')                        # banded general solver
    ops.test('NormDispIncr', 1.0e-8, 50)             # convergence: |du| < 1e-8 m
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.1)               # 10 load steps
    ops.analysis('Static')
    ok = run_steps('load', 10, 0.1)
    ops.loadConst('-time', 0.0)                      # keep gravity, reset time
    return ok


ALGOS = [('Newton',), ('ModifiedNewton', '-initial'), ('KrylovNewton',)]
CTRL = node2d(nStory, 1)                             # roof node on axis 1 (control node)


def one_step(kind, h, level=0):
    """One increment h; try Newton, ModifiedNewton -initial, KrylovNewton,
    then cut the increment in half (up to 1/16)."""
    for alg in ALGOS:
        ops.algorithm(*alg)
        if kind == 'load':
            ops.integrator('LoadControl', h)                 # load-factor increment h
            ok = ops.analyze(1)
        elif kind == 'disp':
            ops.integrator('DisplacementControl', CTRL, 1, h)  # roof ux increment h
            ok = ops.analyze(1)
        else:
            ok = ops.analyze(1, h)                           # transient step of size h
        if ok == 0:
            if alg is not ALGOS[0]:
                ops.algorithm(*ALGOS[0])
            return True
    if level >= 4:
        return False
    return one_step(kind, h/2, level + 1) and one_step(kind, h/2, level + 1)


def run_steps(kind, n, h):
    for _ in range(n):
        if not one_step(kind, h):
            return False
    return True


def ux(t):
    return ops.nodeDisp(t, 1)


def brace_state(t1):
    """Axial force N (tension +) and axial deformation of the brace whose lower element is t1."""
    d = BR[t1]
    xi, yi = NODES[d['i']]; xm, ym = NODES[d['mid']]; xj, yj = NODES[d['j']]
    ui, um, uj = ops.nodeDisp(d['i']), ops.nodeDisp(d['mid']), ops.nodeDisp(d['j'])
    # current direction of the lower element (end i -> mid node)
    vx, vy = (xm + um[0]) - (xi + ui[0]), (ym + um[1]) - (yi + ui[1])
    n = math.hypot(vx, vy)
    f = ops.eleForce(t1)                              # global resisting forces [Fx_i, Fy_i, M_i, Fx_j, ...]
    N = -(f[0]*vx + f[1]*vy)/n                        # end-i force of a member in tension points to j
    Lcur = math.hypot((xj + uj[0]) - (xi + ui[0]), (yj + uj[1]) - (yi + ui[1]))
    return N, Lcur - d['L0']


def base_shear():
    ops.reactions()                                   # update support reactions
    return -sum(ops.nodeReaction(node2d(0, a), 1) for a in list(range(1, nAxis + 1)) + [LEAN])


def floor_ux():
    return np.array([ux(node2d(f, 1)) for f in range(nStory + 1)])


# %% 1. Model ------------------------------------------------------------------
build_model()
print('INFO nodes', len(NODES), 'elements', len(ELES), 'W_total', W_total)

# %% 2. Gravity ----------------------------------------------------------------
ok_g = gravity()
print('INFO gravity', 'ok' if ok_g else 'FAILED')

# %% 3. Modal analysis --------------------------------------------------------
lam = np.array(ops.eigen(3))                          # first 3 eigenvalues (tangent stiffness after gravity)
omega = np.sqrt(lam); T = 2*np.pi/omega
print('INFO periods', np.round(T, 4))
print('RESULT T1_s', round(T[0], 3))
print('RESULT T2_s', round(T[1], 3))

# hand values for the story-1 brace (HSS160x10, length Lb)
A1, I1 = hss_props(*HSS160x10)
Lb1 = math.hypot(L, H)
Py1 = A1*Fy_brace
Pe1 = math.pi**2*Es*I1/Lb1**2                         # Euler load, pinned ends (K = 1)
print('RESULT brace_Py_kN', round(Py1, 0))
print('RESULT brace_Pe_pinned_kN', round(Pe1, 0))

# %% 4. Pushover ---------------------------------------------------------------
ops.timeSeries('Linear', 2)
ops.pattern('Plain', 2, 2)
hf = H*np.arange(1, nStory + 1)
Fx = W_floor*hf/np.sum(W_floor*hf)                    # inverted-triangle pattern, sum = 1
for f in range(1, nStory + 1):
    for a in range(1, nAxis + 1):
        ops.load(node2d(f, a), Fx[f - 1]/nAxis, 0.0, 0.0)   # spread over the 4 frame nodes
ops.test('NormDispIncr', 1.0e-6, 50)
nPush = int(round(roof_target/dU))
c_br, t_br = ele2d(3, 1, 3), ele2d(3, 1, 1)           # story-1 compression (bay 3) / tension (bay 1) braces
push = {'roof': [0.0], 'V': [0.0], 'Nc': [0.0], 'Nt': [0.0]}
u_g = ux(CTRL)
status_p = 'completed'
for k in range(nPush):
    if not one_step('disp', dU):
        status_p = 'stopped by non-convergence at roof drift %.3f %%' % (100*(ux(CTRL) - u_g)/(nStory*H)); break
    if time.time() - T_START > WALL_PUSH:
        status_p = 'stopped by wall-clock limit'; break
    push['roof'].append(ux(CTRL) - u_g); push['V'].append(base_shear())
    push['Nc'].append(brace_state(c_br)[0]); push['Nt'].append(brace_state(t_br)[0])
print('INFO pushover', status_p, 'steps', len(push['roof']) - 1)
for kk in push:
    push[kk] = np.array(push[kk])
drift_p = 100*push['roof']/(nStory*H)
iV = int(np.argmax(push['V'])); iB = int(np.argmin(push['Nc']))
print('INFO first brace buckling at roof drift %.3f %%, V = %.0f kN' % (drift_p[iB], push['V'][iB]))
print('INFO story-1 strength estimate (Py + Pc)*cos(theta) = %.0f kN'
      % ((Py1 - push['Nc'][iB])*L/Lb1))
print('RESULT V_push_max_kN', round(push['V'][iV], 0))
print('RESULT V_push_max_over_W', round(push['V'][iV]/W_total, 3))
print('RESULT drift_at_Vmax_pct', round(drift_p[iV], 3))
print('RESULT brace_Pc_model_kN', round(-push['Nc'][iB], 0))

# %% 5. Nonlinear time history (Northridge LOS270) ----------------------------
acc = np.loadtxt(REC)*SF                               # in g
dt_rec = 0.01
t_free = 3.0                                           # free vibration after the record (residual drift)
build_model(); gravity()
ops.wipeAnalysis()
w1, w3 = omega[0], omega[2]
a0 = zeta*2*w1*w3/(w1 + w3); a1 = zeta*2/(w1 + w3)
ops.rayleigh(a0, 0.0, 0.0, a1)                         # rayleigh(alphaM, betaK, betaKinit, betaKcomm): committed K
ops.timeSeries('Path', 3, '-dt', dt_rec, '-values', *acc.tolist(), '-factor', g)   # m/s^2, zero after the end
ops.pattern('UniformExcitation', 3, 1, '-accel', 3)    # horizontal base acceleration
ops.constraints('Plain'); ops.numberer('RCM'); ops.system('BandGeneral')
ops.test('NormDispIncr', 1.0e-6, 50)
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)                   # average acceleration
ops.analysis('Transient')
nTH = int(round((len(acc)*dt_rec + t_free)/dt_rec))    # 2300 steps
u0 = floor_ux()
th = {'t': [0.0], 'roof': [0.0], 'V': [0.0], 'N': [0.0], 'd': [0.0]}
peak_sd = np.zeros(nStory)
status_t = 'completed'
for k in range(nTH):
    if not one_step('trans', dt_rec):
        status_t = 'non-convergence'; break
    u = floor_ux() - u0
    sd = np.abs(np.diff(u))/H
    peak_sd = np.maximum(peak_sd, sd)
    N, dl = brace_state(t_br)
    th['t'].append(ops.getTime()); th['roof'].append(u[-1]); th['V'].append(base_shear())
    th['N'].append(N); th['d'].append(dl)
    if sd.max() > 0.10:
        status_t = 'collapse (story drift > 10 %)'; break
    if time.time() - T_START > WALL_TOTAL:
        status_t = 'wall-clock limit'; break
for kk in th:
    th[kk] = np.array(th[kk])
t_end = th['t'][-1]
print('INFO time history', status_t, 'at t = %.2f s' % t_end)
print('RESULT TH_time_reached_s', round(t_end, 2))
print('RESULT roof_drift_TH_max_pct', round(100*np.max(np.abs(th['roof']))/(nStory*H), 3))
print('RESULT story_drift_TH_max_pct', round(100*np.max(peak_sd), 3))
nlast = min(len(th['roof']), 100)                      # last 1 s of free vibration
print('RESULT residual_roof_drift_pct', round(100*abs(np.mean(th['roof'][-nlast:]))/(nStory*H), 3))
print('INFO peak story drifts (%):', np.round(100*peak_sd, 3))

# %% 6. Figures ---------------------------------------------------------------
STY = {'col': ('tab:blue', '-', 2.0), 'beam': ('tab:green', '-', 2.0), 'brace': ('tab:red', '-', 1.4),
       'lean': ('0.5', '--', 1.2), 'link': ('0.5', ':', 1.2)}
fig, ax = plt.subplots(figsize=(7, 4.2))
for t, (i, j, kind) in ELES.items():
    (xi, yi), (xj, yj) = NODES[i], NODES[j]
    c, ls, lw = STY[kind]
    ax.plot([xi, xj], [yi, yj], color=c, ls=ls, lw=lw)
    if kind in ('col', 'beam', 'lean') or (kind == 'brace' and t % 100 < 50):
        ax.text((xi + xj)/2, (yi + yj)/2, str(t), fontsize=5.5, ha='center', va='center', color=c,
                bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.8))
for t, (x, y) in NODES.items():
    ax.plot(x, y, 'ko', ms=2.0)
    ax.text(x + 0.1, y + 0.12, str(t), fontsize=5, color='k')
for a in list(range(1, nAxis + 1)) + [LEAN]:
    x, y = NODES[node2d(0, a)]
    ax.plot([x - 0.3, x + 0.3], [y, y], 'k-', lw=3)
for kind, lab in (('col', 'columns 1xxx'), ('beam', 'beams 2xxx'), ('brace', 'braces 3xxx'),
                  ('lean', 'leaning column'), ('link', 'rigid link')):
    c, ls, lw = STY[kind]
    ax.plot([], [], color=c, ls=ls, lw=lw, label=lab)
ax.legend(fontsize=6.5, loc='upper left', bbox_to_anchor=(1.0, 1.0))
ax.set_aspect('equal'); ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
ax.set_title('Ch06 model: 4-story CBF + leaning column (brace mid nodes 9xxx)', fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIGS, 'ch06_model.png'), dpi=110); plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(100*th['roof']/(nStory*H), th['V'], color='0.7', lw=0.6, label='time history (Northridge)')
ax.plot(drift_p, push['V'], 'tab:blue', lw=2, label='pushover')
ax.plot(drift_p[iB], push['V'][iB], 'rv', ms=8, label='story-1 brace buckles')
ax.plot(drift_p[iV], push['V'][iV], 'ko', ms=6, label='peak base shear')
ax.set_xlabel('roof drift (%)'); ax.set_ylabel('base shear (kN)')
ax.grid(alpha=0.3); ax.legend(fontsize=7)
ax2 = ax.twinx(); ax2.set_ylim(np.array(ax.get_ylim())/W_total); ax2.set_ylabel('V / W')
ax.set_title('Base shear vs roof drift', fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIGS, 'ch06_pushover.png'), dpi=110); plt.close(fig)

fig, (a1x, a2x) = plt.subplots(2, 1, figsize=(7, 4.2), sharex=True)
a1x.plot(np.arange(len(acc))*dt_rec, acc, 'k', lw=0.6)
a1x.set_ylabel('ag (g)'); a1x.set_title('Northridge 1994 LOS270 (x%.2f) and roof drift' % SF, fontsize=9)
a2x.plot(th['t'], 100*th['roof']/(nStory*H), 'tab:blue', lw=0.9)
a2x.set_xlabel('time (s)'); a2x.set_ylabel('roof drift (%)')
for a_ in (a1x, a2x):
    a_.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIGS, 'ch06_roof_th.png'), dpi=110); plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(1000*th['d'], th['N'], 'tab:red', lw=0.8)
ax.axhline(Py1, color='k', ls='--', lw=0.8, label='Py = A Fy')
ax.axhline(-Pe1, color='tab:purple', ls=':', lw=1.0, label='-Pe (pinned, K = 1)')
ax.axhline(push['Nc'][iB], color='tab:orange', ls='-.', lw=0.8, label='peak compression in pushover')
ax.set_xlabel('brace axial deformation (mm)'); ax.set_ylabel('brace axial force (kN, tension +)')
ax.set_title('Story-1 brace %d hysteresis during the time history' % t_br, fontsize=9)
ax.grid(alpha=0.3); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIGS, 'ch06_brace_hysteresis.png'), dpi=110); plt.close(fig)

print('INFO run time %.1f s' % (time.time() - T_START))
ops.wipe()
