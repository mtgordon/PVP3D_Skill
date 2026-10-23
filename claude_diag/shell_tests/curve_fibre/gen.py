"""Single-hex uniaxial tests for the curve material (2026-10-23; the user: the wall should not be stiff in compression, it folds).
A 1 mm cube, symmetry BCs on x=0 / y=0 / z=0, the x=1 face pulled / pushed in x (laterally free = uniaxial). Reaction Rx on
x=1 = nominal stress (1 mm2). Materials: 'yeoh' = Curve_soft (the model's curves now); 'mix' = uncoupled solid mixture of a soft
Mooney-Rivlin matrix + tension-only fibre-pow-linear-uncoupled fibres along x (the canal) and y (around the U).
usage: py -3.10 gen.py NAME DISP MAT [c1m E beta lam0 [k]]  -> NAME.feb, run with febio4, Rx in NAME.log"""
import sys, subprocess, os, re
name, disp, mat = sys.argv[1], float(sys.argv[2]), sys.argv[3]
if mat == 'yeohw':
    M = ('<material id="1" name="m" type="Yeoh"><density>1.06e-09</density><c1>0.002087444</c1><c2>0.006592414</c2>'
         '<k>0.35</k></material>')
elif mat == 'cfd':
    c1m, E, beta, lam0 = map(float, sys.argv[4:8]); kk = float(sys.argv[8]) if len(sys.argv) > 8 else 0.035
    M = (f'<material id="1" name="m" type="uncoupled solid mixture"><density>1.06e-09</density><k>{kk}</k>'
         f'<solid type="Mooney-Rivlin"><c1>{c1m}</c1><c2>0</c2></solid>'
         f'<solid type="continuous fiber distribution uncoupled"><fibers type="fiber-pow-linear-uncoupled"><E>{E}</E><beta>{beta}</beta><lam0>{lam0}</lam0></fibers>'
         f'<distribution type="spherical"/><scheme type="fibers-3d-gkt"><nph>7</nph><nth>31</nth></scheme></solid></material>')
elif mat == 'mixa':
    # (10-23 L442c review) in-plane tension-only fibres at the given angles (deg, from the pull direction x, in the x-y plane)
    c1m, E, beta, lam0, kk = map(float, sys.argv[4:9]); angs = [float(a) for a in sys.argv[9].split(',')]
    import math
    fib = lambda a: (f'<solid type="fiber-pow-linear-uncoupled"><fiber type="vector">{math.cos(math.radians(a)):.6f},'
                     f'{math.sin(math.radians(a)):.6f},0</fiber><E>{E}</E><beta>{beta}</beta><lam0>{lam0}</lam0></solid>')
    M = (f'<material id="1" name="m" type="uncoupled solid mixture"><density>1.06e-09</density><k>{kk}</k>'
         f'<solid type="Mooney-Rivlin"><c1>{c1m}</c1><c2>0</c2></solid>{"".join(fib(a) for a in angs)}</material>')
elif mat == 'yeoh':
    M = ('<material id="1" name="m" type="Yeoh"><density>1.06e-09</density><c1>0.001252466</c1><c2>0.003955448</c2>'
         '<k>0.21</k></material>')
else:
    c1m, E, beta, lam0 = map(float, sys.argv[4:8]); kk = float(sys.argv[8]) if len(sys.argv) > 8 else 0.21
    fib = lambda v: (f'<solid type="fiber-pow-linear-uncoupled"><fiber type="vector">{v}</fiber><E>{E}</E><beta>{beta}</beta>'
                     f'<lam0>{lam0}</lam0></solid>')
    M = (f'<material id="1" name="m" type="uncoupled solid mixture"><density>1.06e-09</density><k>{kk}</k>'
         f'<solid type="Mooney-Rivlin"><c1>{c1m}</c1><c2>0</c2></solid>{fib("1,0,0")}{fib("0,1,0")}</material>')
P = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
nodes = ''.join(f'<node id="{i+1}">{x},{y},{z}</node>' for i, (x, y, z) in enumerate(P))
t = f'''<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="4.0"><Module type="solid"/>
<Control><analysis>STATIC</analysis><time_steps>20</time_steps><step_size>0.05</step_size><time_stepper type="default"><max_retries>10</max_retries><opt_iter>15</opt_iter><dtmin>1e-5</dtmin><dtmax>0.05</dtmax></time_stepper><solver type="solid"><max_refs>25</max_refs></solver></Control>
<Material>{M}</Material>
<Mesh><Nodes name="all">{nodes}</Nodes><Elements type="hex8" name="S"><elem id="1">1,2,3,4,5,6,7,8</elem></Elements>
<NodeSet name="x0">1,4,5,8</NodeSet><NodeSet name="y0">1,2,5,6</NodeSet><NodeSet name="z0">1,2,3,4</NodeSet><NodeSet name="x1">2,3,6,7</NodeSet></Mesh>
<MeshDomains><SolidDomain name="S" mat="m"/></MeshDomains>
<Boundary><bc name="a" node_set="x0" type="zero displacement"><x_dof>1</x_dof><y_dof>0</y_dof><z_dof>0</z_dof></bc>
<bc name="b" node_set="y0" type="zero displacement"><x_dof>0</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc>
<bc name="c" node_set="z0" type="zero displacement"><x_dof>0</x_dof><y_dof>0</y_dof><z_dof>1</z_dof></bc>
<bc name="p" node_set="x1" type="prescribed displacement"><dof>x</dof><value lc="1">{disp}</value><relative>0</relative></bc></Boundary>
<LoadData><load_controller id="1" type="loadcurve"><interpolate>LINEAR</interpolate><points><pt>0,0</pt><pt>1,1</pt></points></load_controller></LoadData>
<Output><logfile><node_data data="Rx" node_set="x1" delim=","/></logfile></Output></febio_spec>'''
open(name + '.feb', 'w').write(t)
env = dict(os.environ, OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
subprocess.run([os.environ.get('FEBIO4', r'C:\Program Files\FEBioStudio\bin\febio4.exe'),'-i', name + '.feb', '-silent'], env=env, capture_output=True)
log = open(name + '.log', errors='ignore').read()
out = []
for m in re.finditer(r'Time = ([\d.e+-]+)\s*\n.*?Data = Rx\s*\n((?:\d+,[^\n]*\n)+)', log, re.S):
    T = float(m.group(1)); R = sum(float(l.split(',')[1]) for l in m.group(2).strip().split('\n'))
    out.append((T * disp, R))
print(name, 'ok' if 'N O R M A L' in log else 'FAILED', ' '.join(f'{u:+.2f}:{r:.5f}' for u, r in out[::4] + out[-1:]))
