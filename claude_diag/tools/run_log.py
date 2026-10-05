"""Run log spreadsheet: one row per FEBio run, rebuilt from the run folders (the user's request, 2026-09-26; laid out
per the user's feedback the same evening, then extended to every run folder).

Every folder in claude_diag/runs (except _*) gets a row, in the order the runs were built (the .feb's time). The text
of each row comes from claude_diag/run_log_notes.csv (run, why, result, label, focus, source); the rest from the run's
own files. Columns, left to right:
  Run               the file name (a link to the .feb), its line (springs / lofts, from the name or its ancestors), the
                    run or file it was built from;
  Why, Result       why it was run and what came of it, in words, and whether that text was written at the time or
                    drafted later (the notes file); status, t reached, failed attempts, wall time, start and threads
                    (the .log, runs/_autostop.txt); then the measurements, filled only for the runs whose question they
                    answer (the notes file's focus: rotation = the perineal body's rotation about -x as
                    tools/peb_table.py; pvw / la = the PVW / LA displacement at the end, as tools/wall_disp.py);
  What was changed  from the run's .feb.changes.txt, highlighted, with a label (NOT IN SOURCE, ...);
  then the model, read from each .feb: the loads; one column per property of every material in use (the tissues; each
  connective-tissue family, first whether it is a loft or connectors, then each loft's material and each connector
  set's spring material; the arcus chain; the stabilisation springs, collapsed into one group; the display body);
  contacts and constraints; run control; BCs and mesh.
Yellow = differs from what the run was built from (changed in this run; a blank yellow cell = removed in this run);
light orange = as in the base but different from the line's reference (springs: L19_pen5, "the Abaqus replica";
lofts: L21_lofts4_pen_r2), i.e. changed in an earlier run since that reference (only for the reference's descendants);
grey = not in this run's model.
usage: py -3.10 run_log.py [--out XLSX]   (default: <project>/PVP3D_run_log.xlsx; if that file is open in Excel, a
       timestamped copy is written next to it instead)
"""
import csv
import datetime
import math
import os
import re
import subprocess
import sys
import urllib.parse
import xml.etree.ElementTree as ET
from collections import OrderedDict

import numpy as np
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import DIAG_DIR, JOBS_DIR, RUNS_DIR  # noqa: E402
from xplt import Xplt  # noqa: E402
from batch18_table import log_info  # noqa: E402
from peb_table import rot_series  # noqa: E402
from peb_rotation import region_nodes  # noqa: E402

NOTES = os.path.join(DIAG_DIR, 'run_log_notes.csv')
OUT = os.path.join(JOBS_DIR, 'PVP3D_run_log.xlsx')
REFERENCE = {'springs': 'L19_pen5', 'lofts': 'L21_lofts4_pen_r2'}
# The base models (the user's request, 2026-09-28: mark when a new base is adopted and what changed from the previous one).
# (line, run, adopted, the decision, previous base: a run or a .feb path relative to the project folder). The differences
# on the "Base models" sheet are computed from the two .feb files. Add a row here whenever a new base is adopted.
TRUNK = '(before the two lines)'
BASES = (
    (TRUNK, 'LP2DM_m10', '2026-09-23',
     'Base for further work after the 2026-09-23 sessions (README_2026-09-23b): the user\'s TRUSS_fixed2 file with '
     'Load-LA, the Yeoh LA (the source\'s Marlow law refit), a dynamic analysis and the arcus chain mass x10.',
     'PVP3DModel_v21_v3_HYBRID_withLA_TRUSS_fixed2.feb'),
    (TRUNK, 'L6_ls', '2026-09-24',
     'The user: L5X_all + the line search became the base (copied to Box).', 'LP2DM_m10'),
    (TRUNK, 'L9_fitall_edge', '2026-09-24',
     'The user: the fitted lofts with the PM-PeB anchor-edge BC; "the reference", on the fast base (arcus chain mass x10).',
     'L6_ls'),
    ('springs', 'L19_pen5', '2026-09-25',
     'The user: the springs line\'s reference, "the Abaqus replica": every Abaqus connector as a spring, no lofts.',
     'L9_fitall_edge'),
    ('lofts', 'L21_lofts4_pen_r2', '2026-09-25',
     'The user: the lofts line\'s reference: the user\'s four fitted lofts (AVW-Para, CL, USL, PM), the soft lofts as '
     'their connectors.', 'L9_fitall_edge'),
    ('springs', 'L26_springs_la3_pm', '2026-09-25 (evening)',
     'The current springs line and the base of the 2026-09-26 work: Load-LA x 1/3 (the user\'s choice until the LA\'s '
     'large displacement is explained; NOT IN SOURCE), the time stepper\'s opt_iter 25 -> 10 (the fastest setting), the '
     'source\'s display body PM_Plane.', 'L19_pen5'),
    ('lofts', 'L26_lofts_la3_pm', '2026-09-25 (evening)',
     'The current lofts line and the base of the 2026-09-26 work: the same three changes as on the springs line.',
     'L21_lofts4_pen_r2'),
    ('springs', 'L81_springs_la3_nopinch_pvwsegup2_pen05', '2026-09-27 ~20:00',
     'The user\'s decision (1): the PVW-LA contact = the 30 pinch facets out (NOT IN SOURCE) + penalty 0.5 + seg_up 2.',
     'L26_springs_la3_pm'),
    ('lofts', 'L82_lofts_la3_nopinch_pvwpen05_segup2', '2026-09-27 ~20:00',
     'The user\'s decision (1) on the lofts line: the same PVW-LA contact.', 'L26_lofts_la3_pm'),
    ('springs', 'L83_springs_la3_la100_pebyeoh', '2026-09-27 ~20:00',
     'The user\'s decisions (2)-(3): the perineal-body refit (the source\'s Marlow data fit to Yeoh) and the healthy LA '
     '(the source\'s PCM-LA_Yamada100% table; NOT IN SOURCE: the source uses the 50 % table) with Load-LA x 1/3. The new '
     'springs line, the comparison line.', 'L81_springs_la3_nopinch_pvwsegup2_pen05'),
    ('lofts', 'L83_lofts_la3_la100_pebyeoh', '2026-09-27 ~20:00',
     'The user\'s decisions (2)-(3) on the lofts line: the new lofts line (for faster screens).',
     'L82_lofts_la3_nopinch_pvwpen05_segup2'),
    ('springs', 'L87_springs_newline_rhoi0', '2026-09-28',
     'The user: "go ahead and change rhoi": solver rhoi 0.5 -> 0 (the generalized-alpha integrator\'s full damping of '
     'high-frequency, step-to-step jitter; a solver setting, not a model change) on both lines. The same solution '
     '(Ba / Bp / C within 0.1 mm, LA x1.09) in 22 % fewer iterations, 28 % less time. The walls stay Ogden (the user: '
     '"I\'m fine with having ogden as the material for now when testing other things").',
     'L83_springs_la3_la100_pebyeoh'),
    ('lofts', 'L87_lofts_newline_rhoi0', '2026-09-28',
     'The user: rhoi 0.5 -> 0 on the lofts line too (this run is the check that its solution does not move).',
     'L83_lofts_la3_la100_pebyeoh'),
    ('springs', 'L91_newline_vwyeoh_rhoi0_seamspr001', '2026-09-29 ~11:45',
     'The user: "adjust the main file to use the 0.01 springs for now and then go back to the non-ogden material for '
     'now": the source-faithful vaginal walls (the Marlow data refit as Yeoh: c1 0.005964127, c2 0.01883547, k 1; AVW, '
     'PVW, cervix) and the canal seam springs (NOT IN SOURCE: the 116 AVW / cervix lateral-edge nodes, the whole length '
     'of both sides, each joined to the PVW point it faces by a zero-length spring of 0.01 N/mm in every direction; '
     '3-30x softer per node than the wall tissue itself would be). With the free seam the faithful walls stall at '
     't 0.78; with the springs t = 1 (3 failed, 12073 iterations). Ba / Bp / C -9.0 / -16.6 / -39.8 mm (the Ogden '
     'line -3.3 / -19.2 / -38.3). The lofts line is not changed (it stalls with any softer wall).',
     'L87_springs_newline_rhoi0'),
    ('springs', 'L127_tube_r30s', '2026-09-30',
     'The user adopted the tube (faithful walls kept): the canal seam springs removed; the AVW / cervix and PVW joined '
     'along both sides by a U of wall tissue (Vagina_AVW) with a 3 mm inner radius round the lumen corner, the lumen-edge '
     'nodes left separate on the base\'s contact (NOT IN SOURCE: the source joins the walls there by contact only), so '
     'the walls can peel apart (build_batch122.py / 123); SlidingElastic1 seg_up 0 -> 2 (solver-side: the answer is '
     'unchanged, and it stops the lumen-edge chatter that made the paper case P1 take 4:22 instead of 0:49). t = 1, '
     '0 failed, 40 min; Ba / Bp / C -9.3 / -15.6 / -40.2 mm (L91 -9.0 / -16.5 / -39.8); seam opening p90 3.4 mm '
     '(L91 2.35). README_2026-09-30.md.',
     'L91_newline_vwyeoh_rhoi0_seamspr001'),
    ('springs', 'L128_tube_r30s_pm', '2026-09-30',
     'The user: "make a new base model that includes that as well" (the PM as a structure): L127_tube_r30s + PM_Plane '
     'as a deformable 2 mm shell (PM_Yeoh, fitted to OPAL325_PM_mid\'s PARAVAG_H_highdensity data), its outer arc '
     'clamped, the rigid constraint removed, and the PM_PeB, PM_avw_bottom, PM and PM-LA anchors (47) tied to it '
     '(penalty 100, maxaug 0) with their fixed BCs removed (NOT IN SOURCE; build_batch128.py, as batch 113 step 2). '
     't = 1, 0 failed, 26 min; Ba / Bp / C -9.3 / -15.5 / -40.2 mm (L127_tube_r30s -9.3 / -15.6 / -40.2); the PM moves '
     'up to 2.8 mm.',
     'L127_tube_r30s'),
)
BASE_OF = {b[1]: b for b in BASES}
LOADS =('Load-AVW', 'Load-PVW', 'Load-PeB-top', 'Load-top', 'Load-LA')
LA_DOMAINS = ('LA_PCMPRM', 'LA_PCM', 'LA_ICM', 'LA_ICM_tri')
TISSUES = (('AVW', ('_PickedSet347',)), ('PVW', ('_PickedSet64',)), ('Cervix', ('_PickedSet346',)),
           ('PeB', ('_PickedSet66',)), ('LA', LA_DOMAINS))
FAMILIES = (   # (family, name prefixes of its loft domains (<name>_fan) and connector sets)
    ('AVW-Para', ('AVW-Para',)),
    ('P-arcus', ('P-arcus', 'Parcus')),
    ('CL', ('CL-L', 'CL-R')),
    ('USL', ('USL-L', 'USL-R')),
    ('PM', ('PM_fan', 'PM_conn')),
    ('PM-PeB', ('PM_PeB',)),
    ('PM-AVW bottom', ('PM_avw_bottom',)),
    ('PeB-constrain', ('PeB-constrin',)),
    ('PM-LA', ('PM-LA',)),
    ('Sphincter side (PeB to LA)', ('LA_sphincter_side', 'Sphincter-L', 'Sphincter-R')),
    ('Sphincter posterior (PeB to LA)', ('LA_sphincter_post', 'Sphincter-P')),
)
CHAIN_PREFIXES = ('ATLA', 'Posterior_Arcus')
STAB = 'Stabilisation springs (not in the source)'
CURVE = 'force curve (elongation mm, force N)'
PARAM_ORDER = ('c1', 'm1', 'c2', 'm2', 'c3', 'm3', 'c4', 'm4', 'c5', 'm5', 'c6', 'm6', 'E', 'v', 'G', 'k', 'density',
               'shell thickness', 'shell formulation', 'shell_normal_nodal', 'sets', 'springs', 'scale', 'measure', CURVE)
OLD_JOBS = ()  # the project's old absolute run paths; not needed elsewhere


# ---- reading a .feb (each file once; only its cells are kept) --------------------------------------------------------
def txt(el, tag, default=None):
    e = el.find(tag) if el is not None else None
    return e.text.strip() if e is not None and e.text is not None else default


def num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def val(s):
    """A number if the text is one, else the text."""
    v = num(s)
    return v if v is not None else s


class Feb:
    def __init__(self, path):
        self.root = r = ET.parse(path).getroot()
        self.partial = r.find('Material') is None or r.find('Mesh') is None   # e.g. a restart input: steps only
        if self.partial:
            return
        self.materials = {m.get('name'): m for m in r.find('Material')}
        md = r.find('MeshDomains')
        self.domains = OrderedDict((d.get('name'), d) for d in md) if md is not None else OrderedDict()
        self.dsets, self.nodesets, self.nnodes, self.nelems = OrderedDict(), {}, 0, 0
        self.surfaces, self.pairs = {}, {}
        for ch in r.find('Mesh'):
            if ch.tag == 'Surface':
                self.surfaces[ch.get('name')] = len(ch)
            elif ch.tag == 'SurfacePair':
                self.pairs[ch.get('name')] = (txt(ch, 'primary'), txt(ch, 'secondary'))
            elif ch.tag == 'DiscreteSet':
                self.dsets[ch.get('name')] = len(ch.findall('delem'))
            elif ch.tag == 'NodeSet':
                ids = [v for v in re.split(r'[,\s]+', ch.text or '') if v] + [n.get('id') for n in ch.findall('n')]
                self.nodesets[ch.get('name')] = len(ids)
            elif ch.tag == 'Nodes':
                self.nnodes += len(ch.findall('node'))
            elif ch.tag == 'Elements':
                self.nelems += len(ch.findall('elem'))
        disc = r.find('Discrete')
        dmats = disc.findall('discrete_material') if disc is not None else []
        self.set_mat = {}   # FEBio reads <discrete dmat="k"> as the k-th discrete_material in the list, not by id
        for d in (disc.findall('discrete') if disc is not None else []):
            k = int(d.get('dmat'))
            self.set_mat[d.get('discrete_set')] = dmats[k - 1] if 0 < k <= len(dmats) else None
        ld = r.find('LoadData')
        self.curves = {lc.get('id'): lc for lc in ld} if ld is not None else {}


_CELLS = {}


def cells_of(path):
    """The model cells of a .feb ({column key: value}), cached by path; None if there is no such file."""
    if not path:
        return None
    if path not in _CELLS:
        try:
            _CELLS[path] = cells(Feb(path)) if os.path.exists(path) else None
        except Exception as e:   # an unreadable file: keep the row, compare nothing
            print(f'  cannot read {path}: {type(e).__name__}: {e}', flush=True)
            _CELLS[path] = {'partial': True}
    return _CELLS[path]


def run_feb(run):
    return os.path.join(RUNS_DIR, run, run + '.feb')


def curve_desc(f, lc):
    c = f.curves.get(lc)
    if c is None:
        return f'curve {lc}'
    pts = ' '.join('(' + p.text.strip() + ')' for p in c.iter('pt'))
    return f"{txt(c, 'interpolate', '').lower()}, extend {txt(c, 'extend', '').lower()}: {pts}"


def spring_curve(dm):
    pts = [p.text.strip() for p in dm.iter('pt')]
    if not pts:
        return None
    force = dm.find('force')
    how = [txt(force, k, '') for k in ('interpolate', 'extend')] if force is not None else []
    note = '' if [h.lower() for h in how] == ['linear', 'constant'] else f"[{', '.join(how)}] "
    return note + '; '.join(pts)


def family_of(name):
    for fam, prefixes in FAMILIES:
        if any(name.startswith(p) for p in prefixes):
            return fam
    return None


def side_of(name):
    if re.search(r'-L(\b|_)|_Left|_left', name):
        return 'L'
    if re.search(r'-R(\b|_)|_Right|_right', name):
        return 'R'
    return ''


def mat_params(m):
    p = OrderedDict()
    if m is None:
        return p
    for ch in m:
        if ch.tag in ('pressure_model',) or len(ch):
            continue
        p[ch.tag] = val(ch.text.strip()) if ch.text else None
    return p


def spring_params(f, s, dm):
    p = OrderedDict([('springs', f.dsets[s])])
    p.update(mat_params(dm))
    c = spring_curve(dm)
    if c is not None:
        p[CURVE] = c
    return p


def natural(names):
    return sorted(names, key=lambda s: [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', s)])


def materials(f):
    """{(section, item): (material name, law, {property: value})} for every material a domain or spring set uses."""
    out = OrderedDict()
    used, used_sets = set(), set()

    def shell_settings(p, d):
        if d.tag == 'ShellDomain':
            p['shell formulation'] = d.get('type') or 'default'
            p['shell_normal_nodal'] = val(txt(d, 'shell_normal_nodal', '1'))   # FEBio's default is 1 (averaged normals)

    def dom_item(dname, sec, item):
        d = f.domains[dname]
        m = f.materials.get(d.get('mat'))
        p = mat_params(m)
        if txt(d, 'shell_thickness') is not None:
            p['shell thickness'] = val(txt(d, 'shell_thickness'))
        shell_settings(p, d)
        out[(sec, item)] = (d.get('mat'), m.get('type') if m is not None else '?', p)
        used.add(dname)
    for region, doms in TISSUES:
        present = [d for d in doms if d in f.domains]
        if not present:
            continue
        d = f.domains[present[0]]
        m = f.materials.get(d.get('mat'))
        p = mat_params(m)
        thick = sorted({txt(f.domains[x], 'shell_thickness') for x in present} - {None})
        if thick:
            p['shell thickness'] = val(thick[0]) if len(thick) == 1 else ', '.join(thick)
        shell_settings(p, d)
        out[('Tissue materials', region)] = (d.get('mat'), m.get('type') if m is not None else '?', p)
        used.update(present)
    for fam, _ in FAMILIES:
        for dname in f.domains:
            if dname.endswith('_fan') and family_of(dname) == fam:
                dom_item(dname, 'Connective tissue: ' + fam, dname + ' (loft)')
        for s in natural(f.dsets):
            dm = f.set_mat.get(s)
            if dm is not None and family_of(s) == fam:
                out[('Connective tissue: ' + fam, s + ' (connectors)')] = (dm.get('name'), dm.get('type'), spring_params(f, s, dm))
                used_sets.add(s)
    for dname in f.domains:
        if dname.startswith(CHAIN_PREFIXES) and dname.endswith('_tube'):
            dom_item(dname, 'Arcus chain', dname + ' (tube)')
        elif dname == 'arcus_beam':
            dom_item(dname, 'Arcus chain', 'arcus_beam (beam along the chain)')
        elif dname == 'chain_mass':
            dom_item(dname, 'Arcus chain', 'chain_mass (mass-only trusses)')
    for s in natural(f.dsets):
        dm = f.set_mat.get(s)
        if dm is not None and s.startswith(CHAIN_PREFIXES):
            out[('Arcus chain', s + ' (springs)')] = (dm.get('name'), dm.get('type'), spring_params(f, s, dm))
            used_sets.add(s)
    stab = [s for s in f.dsets if s.startswith('stab_') and f.set_mat.get(s) is not None]
    if stab:    # hundreds of one-spring sets, each with its own material: one group
        p = OrderedDict([('sets', len(stab)), ('springs', sum(f.dsets[s] for s in stab))])
        keys = OrderedDict()
        for s in stab:
            for k, v in mat_params(f.set_mat[s]).items():
                keys.setdefault(k, set()).add(v)
        for k, vs in keys.items():
            p[k] = next(iter(vs)) if len(vs) == 1 else ', '.join(sorted(str(v) for v in vs))
        curves = {spring_curve(f.set_mat[s]) for s in stab} - {None}
        if curves:
            p[CURVE] = next(iter(curves)) if len(curves) == 1 else f'{len(curves)} different curves'
        types = {f.set_mat[s].get('type') for s in stab}
        out[(STAB, 'stab_* ground springs')] = ('stab_<node>_mat', ', '.join(sorted(types)), p)
        used_sets.update(stab)
    if 'PM_Plane' in f.domains:
        dom_item('PM_Plane', 'Display body', 'PM_Plane')
    for dname in f.domains:
        if dname not in used:
            dom_item(dname, 'Other materials', dname)
    for s in natural(f.dsets):
        dm = f.set_mat.get(s)
        if dm is not None and s not in used_sets:
            out[('Other materials', s + ' (springs)')] = (dm.get('name'), dm.get('type'), spring_params(f, s, dm))
    return out


def family_forms(f):
    forms = {}
    for fam, _ in FAMILIES:
        sides = OrderedDict()
        for dname in f.domains:
            if dname.endswith('_fan') and family_of(dname) == fam:
                sides.setdefault(side_of(dname), set()).add('loft')
        for s in f.dsets:
            if f.set_mat.get(s) is not None and family_of(s) == fam:
                sides.setdefault(side_of(s), set()).add('connectors')
        per = {k: ' and '.join(sorted(v)) for k, v in sides.items()}
        if not per:
            forms[fam] = 'none'
        elif len(set(per.values())) == 1:
            forms[fam] = next(iter(per.values()))
        else:
            forms[fam] = ', '.join(f'{k or "centre"} {v}' for k, v in sorted(per.items()))
    return forms


def model_props(f):
    """The non-material columns (keys as in the MODEL columns)."""
    P = {}
    r = f.root
    controls = r.findall('.//Control')   # one, or one per step in a multi-step file (the first step's settings shown)
    c = controls[0]
    tsr, sol = c.find('time_stepper'), c.find('solver')
    P['analysis'] = ' then '.join(dict.fromkeys(txt(x, 'analysis', '') for x in controls))
    P['steps'] = len(controls)
    P['end_t'] = round(sum(int(txt(x, 'time_steps')) * float(txt(x, 'step_size')) for x in controls), 9)
    P['step'] = val(txt(c, 'step_size'))
    P['plot_level'] = txt(c, 'plot_level', '')
    for k in ('dtmax', 'opt_iter', 'max_retries', 'aggressiveness', 'cutback', 'dtmin'):
        P[k] = val(txt(tsr, k))
    for k in ('max_refs', 'max_ups', 'lstol', 'dtol', 'etol', 'rtol', 'symmetric_stiffness', 'rhoi'):
        P[k] = val(txt(sol, k))
    loads, other = r.find('Loads'), []
    psym = set()
    for L in (loads if loads is not None else []):
        nm = L.get('name')
        if L.tag == 'surface_load' and L.get('type') == 'pressure' and nm in LOADS:
            p = L.find('pressure')
            P['load:' + nm] = val(p.text.strip())
            psym.add(txt(L, 'symmetric_stiffness', 'default'))
            if p.get('lc'):
                P['lc:' + nm] = curve_desc(f, p.get('lc'))
        elif L.tag == 'body_load' and L.get('type') == 'mass damping':
            C = L.find('C')
            P['damp_C'] = val(C.text.strip())
            P['damp_curve'] = curve_desc(f, C.get('lc')) if C.get('lc') else 'always on'
        else:
            other.append(f"{nm} ({L.get('type')})")
    curves = {v for k, v in P.items() if k.startswith('lc:')}
    P['load_curve'] = next(iter(curves)) if len(curves) == 1 else '; '.join(
        f"{k[3:]}: {v}" for k, v in P.items() if k.startswith('lc:'))
    for k in [k for k in P if k.startswith('lc:')]:
        del P[k]
    P['load_other'] = '; '.join(other) or 'none'
    P['p_sym'] = ', '.join(sorted(psym)) if psym else None
    cont = r.find('Contact')
    cont = list(cont) if cont is not None else []
    for key, name in (('pvwla', 'PVW_LA'), ('canal', 'SlidingElastic1')):
        cc = next((x for x in cont if x.get('name') == name), None)
        P[key + '_type'] = cc.get('type') if cc is not None else 'none'
        P[key + '_pen'] = val(txt(cc, 'penalty')) if cc is not None else None
        P[key + '_auto'] = val(txt(cc, 'auto_penalty')) if cc is not None else None
        prim, sec = f.pairs.get(cc.get('surface_pair'), (None, None)) if cc is not None else (None, None)
        P[key + '_f1'] = f.surfaces.get(prim)
        P[key + '_f2'] = f.surfaces.get(sec)
    tied = [x for x in cont if x.get('type', '').startswith('tied')]
    ties = [x.get('name') for x in tied]
    P['ties_n'] = len(ties)
    P['ties'] = ', '.join(ties) or 'none'
    pens = sorted({txt(x, 'penalty', '') for x in tied} - {''}, key=lambda s: num(s) or 0)
    P['ties_pen'] = (val(pens[0]) if len(pens) == 1 else ', '.join(pens)) if pens else None
    known = {'PVW_LA', 'SlidingElastic1'} | set(ties)
    P['c_other'] = ', '.join(f"{x.get('name')} ({x.get('type')})" for x in cont if x.get('name') not in known) or 'none'
    cons = r.find('Constraints')
    P['constraints'] = ', '.join(f"{x.get('name')} ({x.get('type')})" for x in cons) if cons is not None and len(cons) else 'none'
    bnd = r.find('Boundary')
    bcs = list(bnd) if bnd is not None else []
    P['bc_n'] = len(bcs)
    P['bc_nodes'] = sum(f.nodesets.get(b.get('node_set'), 0) for b in bcs)
    P['nodes'], P['elements'], P['domains'] = f.nnodes, f.nelems, len(f.domains)
    return P


def contact_settings(f):
    """Every setting of every contact and constraint set, {('c', name, setting): value}, for the Base models sheet only
    (the Runs sheet shows the main ones): contact parameters, the facet counts of each side, and for a constraint set its
    parameters and how many linear constraints it holds."""
    S = {}
    cont = f.root.find('Contact')
    for x in (list(cont) if cont is not None else []):
        name = f"{x.get('name')} ({x.get('type')})"
        for p in x:
            if len(p) == 0:
                S[('c', name, p.tag)] = val(p.text)
        prim, sec = f.pairs.get(x.get('surface_pair'), (None, None))
        S[('c', name, 'primary facets')] = f.surfaces.get(prim)
        S[('c', name, 'secondary facets')] = f.surfaces.get(sec)
    cons = f.root.find('Constraints')
    for x in (list(cons) if cons is not None else []):
        name = f"{x.get('name')} ({x.get('type')})"
        for p in x:
            if len(p) == 0:
                S[('c', name, p.tag)] = val(p.text)
        n = len(x.findall('linear_constraint'))
        if n:
            S[('c', name, 'linear constraints')] = n
    return S


def cells(f):
    """Every model cell of a run: {column key: value} ({'partial': True} for a file that is not a whole model)."""
    if f.partial:
        return {'partial': True}
    C = {('m', k): v for k, v in model_props(f).items()}
    C.update(contact_settings(f))
    for fam, form in family_forms(f).items():
        C[('form', fam)] = form
    for (sec, item), (mname, law, params) in materials(f).items():
        C[('mat', sec, item, 'material')] = f'{mname} ({law})'
        for p, v in params.items():
            C[('mat', sec, item, p)] = v
    return C


# ---- the run's history ---------------------------------------------------------------------------------------------
_BASE = {}


def run_folders():
    return [d for d in os.listdir(RUNS_DIR) if os.path.isdir(os.path.join(RUNS_DIR, d)) and not d.startswith('_')]


def base_info(run):
    """(shown as, base run or None, base .feb path or None, note, change lines) from <run>.feb.changes.txt."""
    if run in _BASE:
        return _BASE[run]
    p = os.path.join(RUNS_DIR, run, run + '.feb.changes.txt')
    if not os.path.exists(p):
        _BASE[run] = ('', None, None, '', [])
        return _BASE[run]
    lines = [l.rstrip() for l in open(p, encoding='latin-1') if l.strip()]
    b = next((l[5:].strip() for l in lines if l.startswith('base:')), '')
    body = [l for l in lines if not l.startswith('base:')]
    shown, brun, bpath, note = b, None, None, ''
    m = re.search(r'([A-Za-z]:[\\/][^()]*?\.feb|[^\s()]+\.feb)', b)
    if m:
        path = m.group(1).strip()
        name = re.sub(r'\.feb$', '', re.split(r'[\\/]', path)[-1])
        note = b[m.end():].strip()
        if os.path.isdir(os.path.join(RUNS_DIR, name)):
            brun, bpath, shown = name, run_feb(name), name
        else:
            low = path.lower()
            for old in OLD_JOBS:
                if low.startswith(old):
                    path = os.path.join(JOBS_DIR, path[len(old):])
            if not os.path.isabs(path):
                path = os.path.join(JOBS_DIR, path)
            bpath, shown = (path if os.path.exists(path) else None), os.path.basename(path)
            if b[:m.start()].strip():
                shown = b.split(' (')[0]
    else:
        first = b.split(' (')[0].strip()
        tok = first.split()[0] if first else ''
        if tok and os.path.isdir(os.path.join(RUNS_DIR, tok)):
            brun, bpath = tok, run_feb(tok)
        shown, note = first, b[len(first):].strip()
    _BASE[run] = (shown, brun, bpath, note, body)
    return _BASE[run]


def compare_path(run):
    """The .feb this run is compared with: its base (the base's base for an identical copy)."""
    shown, brun, bpath, note, _ = base_info(run)
    if brun and 'identical copy' in note:
        b2 = base_info(brun)
        if b2[2]:
            return b2[2]
    return bpath


def ancestors(run):
    out, r = [], run
    for _ in range(80):
        brun = base_info(r)[1]
        if not brun or brun in out or brun == run:
            break
        out.append(brun)
        r = brun
    return out


def line_of(run):
    for line, ref in REFERENCE.items():
        if run == ref:
            return line
    for r in [run] + ancestors(run):
        if '_lofts' in r or r in (REFERENCE['lofts'],):
            return 'lofts'
        if '_springs' in r or 'allconn' in r or r == REFERENCE['springs']:
            return 'springs'
    return ''


def change_text(note, body, limit=280):
    out = [note] if note else []
    for l in body:
        out.append(l if len(l) <= limit else l[:limit].rstrip() + ' [the list continues in the .feb.changes.txt]')
    return '\n'.join(out)


def labels(note, body, extra):
    t = ' '.join(body) + ' ' + note
    labs = [extra] if extra else []
    for key, lab in (('NOT IN SOURCE', 'NOT IN SOURCE'), ('NOT FAITHFUL', 'NOT FAITHFUL'), ('DIAGNOSTIC', 'DIAGNOSTIC'),
                     ('identical copy', 'repeat'), ('display body', 'display only'),
                     ('back to the source', 'back to the source')):
        if key in t and lab not in labs:
            labs.append(lab)
    if re.search(r'run control|run length|solver', t):
        labs.append('run control / solver')
    return ', '.join(dict.fromkeys(labs))


def run_pattern(run):
    return re.compile(r'(?<![\w-])' + re.escape(run) + r'(?![\w-])')


def strip_brackets(s):
    while True:
        s2 = re.sub(r'\([^()]*\)', '', s)
        if s2 == s:
            return s
        s = s2


def started(run, lines):
    pat = run_pattern(run)
    for i, l in enumerate(lines):
        if pat.search(l) and 'start' in l:
            try:
                when = datetime.datetime.strptime(l[:19], '%Y-%m-%d %H:%M:%S')
            except ValueError:
                continue
            thr = None
            for j in range(i, max(-1, i - 4), -1):
                m = re.search(r'\((\d+) threads', lines[j])
                if m:
                    thr = int(m.group(1))
                    break
            return when, thr
    return None, None


def status(run, st, procs, lines):
    if st == 'NORMAL':
        return 'finished normally'
    if st == 'ERROR':
        return 'error termination'
    if any(run + '.feb' in p for p in procs):
        return 'running'
    pat = run_pattern(run)
    for l in lines:   # a run named only inside brackets ("answered by ...") was not the one stopped
        if pat.search(strip_brackets(l)) and 'stopped' in l:
            low = l.lower()
            if 'no new converged step' in low or 'stalled 30' in low:
                return 'stopped (stalled 30 min)'
            if 'answered' in low:
                return 'stopped (question answered)'
            if 'superseded' in low:
                return 'stopped (superseded)'
            if 'user' in low:
                return 'stopped (by the user)'
            return 'stopped by hand'
    return 'ended without a normal finish'


def results(run, focus, procs, lines):
    R = {}
    try:
        t, _, fails, _, st, wall = log_info(run)
    except FileNotFoundError:
        R['status'] = 'not run'
        return R
    R.update(t=t, failed=fails, wall=wall, status=status(run, st, procs, lines))
    R['started'], R['threads'] = started(run, lines)
    xp = os.path.join(RUNS_DIR, run, run + '.xplt')
    if not focus or not os.path.exists(xp):
        return R
    try:
        if 'rotation' in focus:
            s = rot_series(run, every=1)
            if s is not None:
                tt, rr, _, _ = s
                if tt[-1] >= 1 - 1e-3:
                    R['rot1'] = rr[int(np.argmin(np.abs(tt - 1.0)))]
                ip = int(np.argmax(rr))
                R.update(peak=rr[ip], peak_t=tt[ip], rot_end=rr[-1])
        if 'pvw' in focus or 'la' in focus:
            x = Xplt(xp)
            conv = [i for i, s in enumerate(x.states) if s[1] == 0]
            if conv:
                u = x.var(conv[-1], 'displacement')
                if 'pvw' in focus:
                    un = u[region_nodes(run, x, ('_PickedSet64',))]
                    R.update(pvw_med=float(np.median(np.linalg.norm(un, axis=1))), pvw_dy=float(un[:, 1].mean()),
                             pvw_dz=float(un[:, 2].mean()))
                if 'la' in focus:
                    a = np.linalg.norm(u[region_nodes(run, x, LA_DOMAINS)], axis=1)
                    if len(a):
                        R.update(la_med=float(np.median(a)), la_p90=float(np.percentile(a, 90)), la_max=float(a.max()))
    except Exception as e:   # a damaged or partial plot file: keep the row, say so
        R['measure_note'] = f'measurement failed: {type(e).__name__}'
    return R


# ---- columns -------------------------------------------------------------------------------------------------------
# (section, header, sub-header or None, key, width, format); format: text | wrap | gen | num1 | sign1 | t3 | int | date |
# link | curve
STATIC = [
    ('Run', 'File name (click to open the .feb)', None, 'run', 46, 'link'),
    ('Run', 'Line', None, 'line', 9, 'text'),
    ('Run', 'Base model (see the Base models sheet)', None, 'basemark', 16, 'wrap'),
    ('Run', 'Built from (base run or file)', None, 'base', 30, 'wrap'),
    ('Why', 'Why this run was done', None, 'why', 60, 'wrap'),
    ('Result', 'Result in words', None, 'result', 55, 'wrap'),
    ('Result', 'Why / result text', None, 'source', 22, 'wrap'),
    ('Result', 'Status', None, 'status', 16, 'wrap'),
    ('Result', 't reached', None, 't', 8, 't3'),
    ('Result', 'Failed attempts', None, 'failed', 9, 'int'),
    ('Result', 'Wall time', None, 'wall', 10, 'text'),
    ('Result', 'Started', None, 'started', 16, 'date'),
    ('Result', 'Threads', None, 'threads', 8, 'int'),
    ('Result', 'PeB rotation (deg; + = cranial end forward and down)', 'at t = 1', 'rot1', 9, 'sign1'),
    ('Result', 'PeB rotation (deg; + = cranial end forward and down)', 'peak', 'peak', 9, 'sign1'),
    ('Result', 'PeB rotation (deg; + = cranial end forward and down)', 'peak at t', 'peak_t', 9, 't3'),
    ('Result', 'PeB rotation (deg; + = cranial end forward and down)', 'at the end', 'rot_end', 9, 'sign1'),
    ('Result', 'PVW displacement at the end (mm)', '|u| median', 'pvw_med', 9, 'num1'),
    ('Result', 'PVW displacement at the end (mm)', 'mean u_y (+ forward)', 'pvw_dy', 9, 'sign1'),
    ('Result', 'PVW displacement at the end (mm)', 'mean u_z (+ up)', 'pvw_dz', 9, 'sign1'),
    ('Result', 'LA displacement at the end (mm)', '|u| median', 'la_med', 9, 'num1'),
    ('Result', 'LA displacement at the end (mm)', '|u| p90', 'la_p90', 9, 'num1'),
    ('Result', 'LA displacement at the end (mm)', '|u| max', 'la_max', 9, 'num1'),
    ('What was changed', 'What was changed from the base', None, 'changed', 60, 'wrap'),
    ('What was changed', 'Label', None, 'label', 18, 'wrap'),
]
LOAD_COLS = [('Loads', 'Pressure (MPa)', L, ('m', 'load:' + L), 10, 'gen') for L in LOADS] + [
    ('Loads', 'Pressure tangent (symmetric_stiffness)', None, ('m', 'p_sym'), 11, 'gen'),
    ('Loads', 'Load curve (time, factor)', None, ('m', 'load_curve'), 20, 'wrap'),
    ('Loads', 'Settle damping (mass damping)', 'C (1/s)', ('m', 'damp_C'), 9, 'gen'),
    ('Loads', 'Settle damping (mass damping)', 'switch-on curve', ('m', 'damp_curve'), 20, 'wrap'),
    ('Loads', 'Other loads', None, ('m', 'load_other'), 12, 'wrap'),
]
TAIL_COLS = [
    ('Contacts and constraints', 'PVW-LA contact', 'type', ('m', 'pvwla_type'), 13, 'wrap'),
    ('Contacts and constraints', 'PVW-LA contact', 'penalty', ('m', 'pvwla_pen'), 9, 'gen'),
    ('Contacts and constraints', 'PVW-LA contact', 'auto_penalty', ('m', 'pvwla_auto'), 9, 'gen'),
    ('Contacts and constraints', 'PVW-LA contact', 'primary facets', ('m', 'pvwla_f1'), 8, 'gen'),
    ('Contacts and constraints', 'PVW-LA contact', 'secondary facets', ('m', 'pvwla_f2'), 9, 'gen'),
    ('Contacts and constraints', 'Canal contact (SlidingElastic1)', 'type', ('m', 'canal_type'), 13, 'wrap'),
    ('Contacts and constraints', 'Canal contact (SlidingElastic1)', 'penalty', ('m', 'canal_pen'), 9, 'gen'),
    ('Contacts and constraints', 'Canal contact (SlidingElastic1)', 'auto_penalty', ('m', 'canal_auto'), 9, 'gen'),
    ('Contacts and constraints', 'Canal contact (SlidingElastic1)', 'primary facets', ('m', 'canal_f1'), 8, 'gen'),
    ('Contacts and constraints', 'Canal contact (SlidingElastic1)', 'secondary facets', ('m', 'canal_f2'), 9, 'gen'),
    ('Contacts and constraints', 'Tied contacts (lofts to the walls)', 'count', ('m', 'ties_n'), 8, 'gen'),
    ('Contacts and constraints', 'Tied contacts (lofts to the walls)', 'penalty', ('m', 'ties_pen'), 9, 'gen'),
    ('Contacts and constraints', 'Tied contacts (lofts to the walls)', 'names', ('m', 'ties'), 26, 'wrap'),
    ('Contacts and constraints', 'Other contacts', None, ('m', 'c_other'), 16, 'wrap'),
    ('Contacts and constraints', 'Constraints', None, ('m', 'constraints'), 20, 'wrap'),
    ('Run control', 'Analysis', None, ('m', 'analysis'), 10, 'text'),
    ('Run control', 'Steps', None, ('m', 'steps'), 6, 'gen'),
    ('Run control', 'End time (s)', None, ('m', 'end_t'), 8, 'gen'),
    ('Run control', 'Step size (s)', None, ('m', 'step'), 8, 'gen'),
] + [('Run control', 'Time stepper', k, ('m', k), 9, 'gen')
     for k in ('dtmax', 'opt_iter', 'max_retries', 'aggressiveness', 'cutback', 'dtmin')] + [
    ('Run control', 'Solver', k, ('m', k), 8 if k != 'symmetric_stiffness' else 11, 'gen')
    for k in ('max_refs', 'max_ups', 'lstol', 'dtol', 'etol', 'rtol', 'symmetric_stiffness', 'rhoi')] + [
    ('Run control', 'Plot level', None, ('m', 'plot_level'), 14, 'wrap')] + [
    ('Boundary conditions and mesh', 'Boundary conditions', 'count', ('m', 'bc_n'), 8, 'gen'),
    ('Boundary conditions and mesh', 'Boundary conditions', 'nodes fixed', ('m', 'bc_nodes'), 8, 'gen'),
    ('Boundary conditions and mesh', 'Mesh', 'nodes', ('m', 'nodes'), 8, 'gen'),
    ('Boundary conditions and mesh', 'Mesh', 'elements', ('m', 'elements'), 9, 'gen'),
    ('Boundary conditions and mesh', 'Mesh', 'domains', ('m', 'domains'), 8, 'gen'),
]
SECTION_ORDER = (['Tissue materials'] + ['Connective tissue: ' + f for f, _ in FAMILIES] +
                 ['Arcus chain', STAB, 'Display body', 'Other materials'])
TISSUE_ORDER = [t for t, _ in TISSUES]
PARAM_WIDTH = {'material': 22, 'measure': 10, CURVE: 34, 'springs': 8, 'sets': 7}
GROUP_COLOR = {'Run': '44546A', 'Why': '2F5597', 'Result': '548235', 'What was changed': 'BF8F00', 'Loads': '7F6000',
               'Tissue materials': '843C0C', 'Arcus chain': '7030A0', STAB: '5B5B5B', 'Display body': '404040',
               'Other materials': '404040', 'Contacts and constraints': '1F4E79', 'Run control': '595959',
               'Boundary conditions and mesh': '404040'}
CT_COLORS = ('C55A11', '9C5700')   # connective-tissue families alternate
CHANGED = PatternFill('solid', fgColor='FFE699')     # changed in this run (vs its base); blank = removed in this run
INHERITED = PatternFill('solid', fgColor='FBE2D5')   # differs from the line's reference, changed in an earlier run
ABSENT = PatternFill('solid', fgColor='EDEDED')      # not in this run's model
BASEFILL = PatternFill('solid', fgColor='C6EFCE')    # a base model (the Base models sheet)
NUMFMT = {'num1': '0.0', 'sign1': '+0.0;-0.0;0.0', 't3': '0.000', 'int': '0', 'date': 'yyyy-mm-dd hh:mm',
          'gen': 'General'}
FONT = Font(name='Arial', size=10)
THIN = Side(style='thin', color='BFBFBF')
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)


def material_columns(all_cells):
    """The material columns present in any run, in section order; a property that is zero or absent in every run is
    left out (with its Ogden exponent: a term whose c is 0 in every run does nothing)."""
    items = OrderedDict()
    for C in all_cells:
        for k in C:
            if k[0] == 'mat':
                items.setdefault((k[1], k[2]), set()).add(k[3])

    def nonzero(sec, item, p):
        return any(C.get(('mat', sec, item, p)) not in (None, 0, 0.0, '') for C in all_cells)

    def item_order(key):
        sec, item = key
        si = SECTION_ORDER.index(sec) if sec in SECTION_ORDER else 99
        if sec == 'Tissue materials':
            return (si, TISSUE_ORDER.index(item), [])
        rank = 0 if item.endswith('(loft)') else (1 if item.endswith(('(tube)', 'trusses)', 'chain)')) else 2)
        return (si, rank, [(0, int(t)) if t.isdigit() else (1, t) for t in re.split(r'(\d+)', item)])
    cols, fams_done = [], set()
    for sec, item in sorted(items, key=item_order):
        if sec.startswith('Connective tissue: ') and sec not in fams_done:
            fams_done.add(sec)
            cols.append((sec, 'Loft or connectors', None, ('form', sec[len('Connective tissue: '):]), 12, 'wrap'))
        names = {C.get(('mat', sec, item, 'material')) for C in all_cells} - {None}
        head = item + (': ' + next(iter(names)) if len(names) == 1 else '')
        params = [p for p in PARAM_ORDER if p in items[(sec, item)] and nonzero(sec, item, p)
                  and not (p[0] == 'm' and p[1:].isdigit() and not nonzero(sec, item, 'c' + p[1:]))]
        params += sorted(p for p in items[(sec, item)] - set(PARAM_ORDER) - {'material'} if nonzero(sec, item, p))
        if len(names) > 1:
            params = ['material'] + params
        for p in params:
            fmt = 'curve' if p == CURVE else ('text' if p in ('measure', 'material') else 'gen')
            cols.append((sec, head, p, ('mat', sec, item, p), PARAM_WIDTH.get(p, 9), fmt))
    return cols


def columns(all_cells):
    return STATIC + LOAD_COLS + material_columns(all_cells) + TAIL_COLS


# ---- the workbook --------------------------------------------------------------------------------------------------
def lighten(hex6, f=0.78):
    rgb = [int(hex6[i:i + 2], 16) for i in (0, 2, 4)]
    return ''.join(f'{int(v + (255 - v) * f):02X}' for v in rgb)


def section_color(sec):
    fams = [f for f, _ in FAMILIES]
    if sec.startswith('Connective tissue: '):
        fam = sec[len('Connective tissue: '):]
        return CT_COLORS[fams.index(fam) % 2] if fam in fams else CT_COLORS[0]
    return GROUP_COLOR.get(sec, '404040')


def file_url(path):
    return 'file:///' + urllib.parse.quote(path.replace('\\', '/'), safe=':/')


def build_rows(notes):
    procs = subprocess.run(['powershell', '-NoProfile', '-Command',
                            "Get-CimInstance Win32_Process -Filter \"Name='febio4.exe'\" | ForEach-Object { $_.CommandLine }"],
                           capture_output=True, text=True).stdout.splitlines()
    ap = os.path.join(RUNS_DIR, '_autostop.txt')
    lines = open(ap, encoding='latin-1').read().splitlines() if os.path.exists(ap) else []
    by_run = {n['run']: n for n in notes}
    runs = run_folders()
    runs.sort(key=lambda r: (os.path.getmtime(run_feb(r)) if os.path.exists(run_feb(r)) else 0, r))
    rows = []
    for i, run in enumerate(runs, 1):
        n = by_run.get(run, {})
        shown, brun, bpath, note, body = base_info(run)
        focus = {w.strip() for w in (n.get('focus') or '').split(',') if w.strip()}
        line = line_of(run)
        has_log = os.path.exists(os.path.join(RUNS_DIR, run, run + '.feb.changes.txt'))
        bm = BASE_OF.get(run)
        row = {'run': run, 'line': line, 'base': shown or ('' if has_log else '(no change log)'),
               'basemark': f"BASE: {bm[0].strip('()')}, adopted {bm[2]}" if bm else '',
               'why': n.get('why', ''), 'result': n.get('result', ''),
               'source': n.get('source', ''), 'changed': change_text(note, body),
               'label': labels(note, body, n.get('label', ''))}
        row.update(results(run, focus, procs, lines))
        if row.get('measure_note'):
            row['result'] = (row['result'] + ' ' if row['result'] else '') + f"[{row['measure_note']}]"
        row['cells'] = cells_of(run_feb(run)) or {}
        row['base_cells'] = cells_of(compare_path(run))
        ref = REFERENCE.get(line)
        row['ref_cells'] = cells_of(run_feb(ref)) if ref and ref != run and ref in ancestors(run) else None
        if row['cells'].get('partial'):   # not a whole model (a restart input): nothing to compare
            row['base_cells'] = row['ref_cells'] = None
            row['label'] = ', '.join(x for x in (row['label'], "restart input: the model is its base's") if x)
        if row['base_cells'] and row['base_cells'].get('partial'):
            row['base_cells'] = None
        rows.append(row)
        if i % 25 == 0:
            print(f'  {i} / {len(runs)} runs read', flush=True)
    return rows


def widen(cols):
    """Widen a header group whose columns are too narrow for its row-2 text in 3 lines."""
    out, i = [list(c) for c in cols], 0
    while i < len(out):
        j = i
        while j + 1 < len(out) and out[j + 1][:2] == out[i][:2] and out[j + 1][2] is not None and out[i][2] is not None:
            j += 1
        if out[i][2] is not None:
            need = len(out[i][1]) / 3 + 2
            have = sum(c[4] for c in out[i:j + 1])
            if have < need:
                for c in out[i:j + 1]:
                    c[4] = math.ceil(c[4] * need / have)
        i = j + 1
    return [tuple(c) for c in out]


def write(rows, out):
    cols = widen(columns([r['cells'] for r in rows]))
    wb = Workbook()
    ws = wb.active
    ws.title = 'Runs'
    hdr_font = Font(name='Arial', size=10, bold=True)
    # row 1: sections; row 2: the item a group of sub-columns belongs to; row 3: the column's own header (the filter row)
    j0 = 1
    for j in range(1, len(cols) + 1):
        sec = cols[j - 1][0]
        if j == len(cols) or cols[j][0] != sec:
            if j > j0:
                ws.merge_cells(start_row=1, start_column=j0, end_row=1, end_column=j)
            c = ws.cell(row=1, column=j0, value=sec)
            c.font = Font(name='Arial', size=11, bold=True, color='FFFFFF')
            for jj in range(j0, j + 1):
                ws.cell(row=1, column=jj).fill = PatternFill('solid', fgColor=section_color(sec))
            c.alignment = Alignment(horizontal='left', vertical='center')
            j0 = j + 1
    j0 = 1
    for j in range(1, len(cols) + 1):
        sec, head, sub = cols[j - 1][:3]
        light = PatternFill('solid', fgColor=lighten(section_color(sec)))
        if sub is None:
            ws.cell(row=2, column=j).fill = light
            c = ws.cell(row=3, column=j, value=head)
            c.font, c.fill, c.alignment = hdr_font, light, Alignment(wrap_text=True, vertical='bottom')
            j0 = j + 1
            continue
        c3 = ws.cell(row=3, column=j, value=sub)
        c3.font, c3.fill, c3.alignment = hdr_font, light, Alignment(wrap_text=True, vertical='bottom')
        if j == len(cols) or cols[j][:2] != (sec, head) or cols[j][2] is None:
            if j > j0:
                ws.merge_cells(start_row=2, start_column=j0, end_row=2, end_column=j)
            c = ws.cell(row=2, column=j0, value=head)
            c.font = hdr_font
            c.alignment = Alignment(wrap_text=True, vertical='bottom')
            for jj in range(j0, j + 1):
                ws.cell(row=2, column=jj).fill = light
            j0 = j + 1
    for j, col in enumerate(cols, 1):
        ws.column_dimensions[get_column_letter(j)].width = col[4]
        for rr in (2, 3):
            ws.cell(row=rr, column=j).border = BORDER
    ws.row_dimensions[1].height = 18
    ws.row_dimensions[2].height = 42
    ws.row_dimensions[3].height = 42
    for i, row in enumerate(rows, 4):
        nlines = 1
        for j, (sec, head, sub, key, width, fmt) in enumerate(cols, 1):
            model = isinstance(key, tuple)
            v = row['cells'].get(key) if model else row.get(key)
            if isinstance(v, (np.floating, np.integer)):
                v = v.item()
            c = ws.cell(row=i, column=j, value=v if v != '' else None)
            c.font, c.border = FONT, BORDER
            c.alignment = Alignment(wrap_text=fmt in ('wrap', 'link'), vertical='top')
            if fmt in NUMFMT:
                c.number_format = NUMFMT[fmt]
            if fmt == 'link':
                c.hyperlink = file_url(run_feb(v))
                c.font = Font(name='Arial', size=10, color='0563C1', underline='single')
            if key == 'changed':
                c.fill = CHANGED
            elif key == 'basemark' and v:
                c.fill = BASEFILL
                c.font = Font(name='Arial', size=10, bold=True)
            elif model:
                b = row['base_cells'].get(key) if row['base_cells'] is not None else v
                rf = row['ref_cells'].get(key) if row['ref_cells'] is not None else v
                if b != v:
                    c.fill = CHANGED
                elif rf != v:
                    c.fill = INHERITED
                elif v is None:
                    c.fill = ABSENT
            if isinstance(v, str) and fmt in ('wrap', 'link'):
                per = max(1, int(width * 1.15) - 1)
                nlines = max(nlines, sum(max(1, math.ceil(len(p) / per)) for p in v.split('\n')))
        ws.row_dimensions[i].height = min(409, 13.5 * nlines + 3)
    ws.freeze_panes = 'B4'
    ws.auto_filter.ref = f'A3:{get_column_letter(len(cols))}{3 + len(rows)}'
    base_sheet(wb, rows, cols)
    legend(wb)
    try:
        wb.save(out)
        return out, len(cols)
    except PermissionError:
        alt = out.replace('.xlsx', f'_{datetime.datetime.now():%Y%m%d_%H%M%S}.xlsx')
        wb.save(alt)
        return alt, len(cols)


def base_chain(run, prev):
    """How a base was built from the previous one: the runs in between, or, if it was not built from it, from the
    nearest ancestor the two share."""
    anc = ancestors(run)
    if prev in anc:
        return ' -> '.join(reversed(anc[:anc.index(prev)])) or '(built directly from it)'
    pa = [prev] + ancestors(prev)
    common = next((a for a in anc if a in pa), None)
    if common is None:
        return 'not built from it'
    between = list(reversed(anc[:anc.index(common)]))
    return f'not built from it; from {common} (its ancestor too) through ' + (' -> '.join(between) or '(directly)')


def base_sheet(wb, rows, cols):
    """One block per base model: its line, date and decision, the runs it was built through from the previous base, then
    one row per model property that differs from the previous base (both read from the .feb files, as the Runs sheet)."""
    ws = wb.create_sheet('Base models', 1)
    heads = [('Line', 12), ('Base model (click to open the .feb)', 40), ('Adopted', 13), ('Decision', 58),
             ('Previous base', 34), ('Built from the previous base through', 36), ('Section', 24), ('Item', 34),
             ('Property', 16), ("Previous base's value", 28), ("This base's value", 28)]
    hdr = Font(name='Arial', size=10, bold=True, color='FFFFFF')
    for j, (h, w) in enumerate(heads, 1):
        c = ws.cell(row=1, column=j, value=h)
        c.font, c.fill, c.border = hdr, PatternFill('solid', fgColor='375623'), BORDER
        c.alignment = Alignment(wrap_text=True, vertical='bottom')
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.row_dimensions[1].height = 30
    by_run = {r['run']: r for r in rows}
    label = {key: (sec, head, sub) for sec, head, sub, key, _, _ in cols if isinstance(key, tuple)}
    order = [key for *_, key, _, _ in cols if isinstance(key, tuple)]

    def model(run_or_file):
        if os.path.isdir(os.path.join(RUNS_DIR, run_or_file)):
            C = by_run.get(run_or_file, {}).get('cells') or cells_of(run_feb(run_or_file))
        else:
            C = cells_of(os.path.join(JOBS_DIR, run_or_file))
        return {} if not C or C.get('partial') else C

    def plain(v):
        return v.item() if isinstance(v, (np.floating, np.integer)) else v

    r = 2
    for line, run, adopted, decision, prev in BASES:
        cur, old = model(run), model(prev)
        diffs = [k for k in order if cur.get(k) != old.get(k)]
        # every contact / constraint setting (not columns on the Runs sheet); the ones the Runs sheet shows are left out
        def on_runs_sheet(k):
            return (k[1].split(' (')[0] in ('PVW_LA', 'SlidingElastic1')
                    and k[2] in ('penalty', 'auto_penalty', 'primary facets', 'secondary facets'))
        extra = sorted((k for k in set(cur) | set(old)
                        if k[0] == 'c' and cur.get(k) != old.get(k) and not on_runs_sheet(k)),
                       key=lambda k: (k[1], k[2]))
        for k in extra:
            label[k] = ('Contacts and constraints (every setting)', k[1], k[2])
        diffs += extra
        top = [line, run, adopted, decision, prev,
               base_chain(run, prev) if os.path.isdir(os.path.join(RUNS_DIR, prev)) else '(the user\'s file)',
               (f'{len(diffs)} {"property differs" if len(diffs) == 1 else "properties differ"} from the previous base:'
                if old else '(previous base not readable)')]
        for j, v in enumerate(top, 1):
            c = ws.cell(row=r, column=j, value=v)
            c.font, c.fill, c.border = Font(name='Arial', size=10, bold=True), BASEFILL, BORDER
            c.alignment = Alignment(wrap_text=True, vertical='top')
        for j in range(len(top) + 1, len(heads) + 1):
            ws.cell(row=r, column=j).fill = BASEFILL
            ws.cell(row=r, column=j).border = BORDER
        link = ws.cell(row=r, column=2)
        link.hyperlink = file_url(run_feb(run))
        link.font = Font(name='Arial', size=10, bold=True, color='0563C1', underline='single')
        ws.row_dimensions[r].height = min(409, 13.5 * max(math.ceil(len(decision) / 64),
                                                          math.ceil(len(top[5]) / 40)) + 4)
        r += 1
        for k in diffs:
            sec, head, sub = label[k]
            vals = [sec, head, sub or '', plain(old.get(k)), plain(cur.get(k))]
            for j, v in enumerate(vals, 7):
                c = ws.cell(row=r, column=j, value='(none)' if v is None else v)
                c.font, c.border = FONT, BORDER
                c.alignment = Alignment(wrap_text=True, vertical='top')
            longest = max(len(str(v)) for v in vals)
            ws.row_dimensions[r].height = min(409, 13.5 * max(1, math.ceil(longest / 30)) + 3)
            r += 1
        r += 1
    ws.freeze_panes = 'C2'


LEGEND = [
    ('How the sheet is built', 'One row per FEBio run folder, in the order the runs were built, rebuilt by '
     'claude_diag/tools/run_log.py. The "why", "result in words" and focus of each run come from '
     'claude_diag/run_log_notes.csv; everything else is read from the run\'s own files: <run>.feb (the model), '
     '<run>.feb.changes.txt (what was changed and the base), <run>.log (status, t reached, failed attempts, wall time), '
     '<run>.xplt (the measurements) and runs/_autostop.txt (start time, threads, stops). Runs live in '
     '<project>/claude_diag/runs/<run>/.'),
    ('Why / result text', '"written ..." = written when the run was discussed; "drafted 2026-09-26 from ..." = drafted '
     'afterwards from the notes named (the READMEs, SESSION_SUMMARY.md, the build script, the change log). Drafted '
     'text may miss what was said in conversation; the numbers in the other columns are read from the files.'),
    ('Colours', 'Yellow = differs from what the run was built from, i.e. changed in this run (a blank yellow cell = '
     'removed in this run). Light orange = the same as the base but different from the line\'s reference (springs '
     'line: L19_pen5, "the Abaqus replica"; lofts line: L21_lofts4_pen_r2), i.e. changed in an earlier run since that '
     'reference; only for runs built from the reference. Grey = not in this run\'s model. The "What was changed" column '
     'is always yellow. Runs built from a file outside the runs folder (the 2026-09-22/23 runs, from the original '
     'PVP3DModel_*.feb files) are compared with that file.'),
    ('Header rows', 'Row 1: section. Row 2: the item a group of columns belongs to (for materials: where it is used and '
     'the material\'s name). Row 3: the column itself (the property); the filter buttons are on this row.'),
    ('Line', 'springs = the Abaqus replica line (every connector a spring), lofts = the lofts line; from the run\'s name '
     'or its ancestors. Blank for the early runs made before the two lines.'),
    ('Base models', 'The "Base model" column (green) marks the runs adopted as a base: the model the next runs were built '
     'from and compared with, with the date it was adopted. The "Base models" sheet lists them in order per line: the '
     'decision behind each, the runs it was built through from the previous base, and every model property that differs '
     'from the previous base (read from the two .feb files, as the columns of this sheet). The list itself is kept in '
     'claude_diag/tools/run_log.py (BASES).'),
    ('Measurements', 'Filled only for the runs whose question they answer (the focus column of the notes file). '
     'PeB rotation: the perineal body\'s (_PickedSet66) best-fit rigid rotation about -x, + = cranial end forward (+y) '
     'and down (-z) (tools/peb_table.py, every converged state). PVW and LA displacement: over the PVW nodes '
     '(_PickedSet64) and the LA shells at the last converged state (tools/wall_disp.py). The end = the last converged '
     'state (t = 1.5 for a finished hold run).'),
    ('Status', 'finished normally; running; stopped (question answered / superseded / stalled 30 min / by the user / by '
     'hand, from runs/_autostop.txt); ended without a normal finish (the process died or was stopped without a note in '
     'the log); error termination (FEBio gave up: max retries); not run (built only).'),
    ('Shell settings', 'shell formulation = the ShellDomain type (default or three-field-shell); shell_normal_nodal = 1 '
     '(FEBio\'s default: averaged nodal normals) or 0 (element normals, the fix for shells pinched at rest).'),
    ('Label', 'NOT IN SOURCE = not in the Abaqus source model; NOT FAITHFUL = deliberately unlike the source (a test); '
     'DIAGNOSTIC = a test, not a proposed fix; run control / solver = no model change; display only; repeat; back to '
     'the source; control; line reference; session base; reference / saved as.'),
    ('Loads', 'The five pressures (MPa; the source value is 0.014 MPa = 140 cm H2O) and the curve they follow; the '
     'settle damping (mass damping, C in 1/s, times its switch-on curve).'),
    ('Materials', 'One column per property of every material the model uses, grouped by where it is used. Units: '
     'mm, N, s, tonne, so stresses and moduli (c1, c2, k, E) are in MPa and density is in tonne/mm3 (1.06e-9 = '
     '1060 kg/m3). Ogden: c_i, m_i (the i-th term; terms whose c is 0 in every run are left out), k (bulk modulus). '
     'Yeoh: c1, c2, k. Isotropic elastic and linear truss: E, v, density. Shells also list their thickness (mm). When '
     'an item\'s material changes between runs, its first column gives the material\'s name.'),
    ('Connective tissue', 'Each family (AVW-Para, P-arcus, CL, USL, PM, ...) first says whether it is a loft, '
     'connectors or none in this run, then lists the lofts\' materials and the connector sets\' spring materials. '
     'A loft is a shell lofted over the source\'s connector lines, with its own material fitted to those connectors '
     '(the early runs used a shared beam material, Beam-CL-USL, E = 21). Connectors are the source\'s connectors as '
     'nonlinear springs: springs = how many, scale = the multiplier on the force curve (1 = the source), measure = '
     'what the curve is a function of, force curve = its points (elongation mm, force N). PM-PeB, PM-AVW bottom and '
     'PeB-constrain are already impaired to almost nothing in the source.'),
    ('Arcus chain', 'The arcus tendons: tubes in the earliest runs, then chains of springs, with mass-only trusses '
     '(chain_mass) and in a few runs a stabilising beam. The chain mass density is 0.00011 in the source and 0.0011 '
     '(x10) in the fast base used for screening.'),
    ('Stabilisation springs', 'The converted file\'s zero-length ground springs (stab_*; not in the Abaqus source), one '
     'set and one material per node, shown as one group: sets, springs, and the scale and curve they share. The 328 '
     'midline ones were removed on 2026-09-24 in the morning (batch 19), the other 145 that evening (L11_stab).'),
    ('Other settings', 'Settings that are the same in every run and are not material properties (e.g. Ogden '
     'pressure_model = default, the springs\' linear interpolation and constant extension) are not listed; a spring '
     'curve that differs from linear/constant says so in its cell.'),
]


def legend(wb):
    lg = wb.create_sheet('How to read')
    lg.column_dimensions['A'].width = 22
    lg.column_dimensions['B'].width = 120
    c = lg.cell(row=1, column=1, value='PVP3D run log: how to read it')
    c.font = Font(name='Arial', size=14, bold=True)
    c = lg.cell(row=2, column=1, value=f'Built {datetime.datetime.now():%Y-%m-%d %H:%M}.')
    c.font = FONT
    r = 4
    for key, text in LEGEND:
        a = lg.cell(row=r, column=1, value=key)
        a.font = Font(name='Arial', size=10, bold=True)
        a.alignment = Alignment(vertical='top')
        b = lg.cell(row=r, column=2, value=text)
        b.font = FONT
        b.alignment = Alignment(wrap_text=True, vertical='top')
        lg.row_dimensions[r].height = 13.5 * math.ceil(len(text) / 135) + 4
        r += 1
    r += 1
    for fill, text in ((CHANGED, 'changed in this run'), (INHERITED, 'changed in an earlier run (differs from the '
                                                                     'line\'s reference)'), (ABSENT, 'not in this model')):
        a = lg.cell(row=r, column=1, value='')
        a.fill = fill
        b = lg.cell(row=r, column=2, value=text)
        b.font = FONT
        r += 1


def main():
    args = sys.argv[1:]
    out = args[args.index('--out') + 1] if '--out' in args else OUT
    notes = list(csv.DictReader(open(NOTES, encoding='utf-8-sig')))
    rows = build_rows(notes)
    path, ncols = write(rows, out)
    print(f'wrote {path}: {len(rows)} runs, {ncols} columns')


if __name__ == '__main__':
    main()
