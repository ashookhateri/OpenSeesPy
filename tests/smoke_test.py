# Run on the phone (Pydroid 3 or Termux) after: pip install openseespyandroid
import openseespyandroid.opensees as ops

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)
ops.node(1, 0.0, 0.0); ops.node(2, 0.0, 3.0)
ops.fix(1, 1, 1, 1)
ops.geomTransf('Linear', 1)
ops.element('elasticBeamColumn', 1, 1, 2, 0.09, 2.0e8, 6.75e-4, 1)
ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)
ops.load(2, 10.0, 0.0, 0.0)
ops.system('BandGeneral'); ops.numberer('RCM'); ops.constraints('Plain')
ops.integrator('LoadControl', 1.0); ops.algorithm('Linear'); ops.analysis('Static')
assert ops.analyze(1) == 0
d = ops.nodeDisp(2, 1)
exact = 10.0 * 3.0**3 / (3 * 2.0e8 * 6.75e-4)
print('tip disp =', d, 'exact =', exact)
assert abs(d - exact) / exact < 1e-6
print('OK - openseespyandroid works')
