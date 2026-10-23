"""Normalized structural summary of a .feb, for diffing model versions section by section.

Mesh blocks are reduced to name / type / count / id-range plus a short content hash,
so a coordinate-level change shows up as one changed hash instead of thousands of
diff lines. Every other section (Control, Material, MeshDomains, Boundary, Loads,
Contact, Constraints, Discrete, LoadData, Output) is flattened to one sorted
"tag attrs = value" line per parameter.

usage: py -3 summarize_feb.py model.feb [--collapse PREFIX ...] > model.summary.txt
       diff old.summary.txt new.summary.txt
  --collapse PREFIX   summarize DiscreteSets / discrete materials / bindings whose name
                      starts with PREFIX as a single count+hash line (useful for hundreds
                      of generated stabilization springs, e.g. --collapse stab_)
Renumbered FEBioStudio re-exports change every hash; compare counts and the non-Mesh
sections first in that case.
Needs Python 3 (standard library only).
"""
import argparse
import hashlib
import xml.etree.ElementTree as ET


def h(s):
    return hashlib.md5(s.encode()).hexdigest()[:8]


def canon(el):
    """Whitespace-insensitive serialization for hashing. ET.tostring() includes the element's tail and the
    indentation between children, so re-indenting, or appending a sibling after the last block, changed the
    hash of an unchanged block."""
    attrs = ','.join(f'{k}={v}' for k, v in sorted(el.attrib.items()))
    return f'<{el.tag} {attrs}>{(el.text or "").strip()}' + ''.join(canon(c) for c in el) + f'</{el.tag}>'


def flat(el, indent=''):
    attrs = ' '.join(f'{k}={v}' for k, v in sorted(el.attrib.items()))
    if len(el):
        out = [f'{indent}<{el.tag} {attrs}>']
        for c in el:
            out.extend(flat(c, indent + '  '))
        return out
    return [f'{indent}{el.tag} {attrs} = {(el.text or "").strip()}']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('feb')
    ap.add_argument('--collapse', action='append', default=[])
    a = ap.parse_args()

    def collapsed(name):
        return any((name or '').startswith(p) for p in a.collapse)

    root = ET.parse(a.feb).getroot()
    for sec in root:
        if sec.tag == 'Mesh':
            coll = []
            for child in sec:
                name = child.get('name')
                items = list(child)
                raw = canon(child)
                if child.tag == 'Nodes':
                    ids = [int(n.get('id')) for n in items]
                    print(f'Mesh/Nodes {name}: n={len(ids)} ids={min(ids)}-{max(ids)} hash={h(raw)}')
                elif child.tag == 'Elements':
                    ids = [int(n.get('id')) for n in items]
                    print(f'Mesh/Elements {name} type={child.get("type")}: n={len(ids)} '
                          f'ids={min(ids)}-{max(ids)} hash={h(raw)}')
                elif child.tag == 'NodeSet':
                    txt = (child.text or '').replace('\n', ',')
                    n = len([v for v in txt.split(',') if v.strip()]) + len(items)
                    print(f'Mesh/NodeSet {name}: n={n} hash={h(raw)}')
                elif child.tag == 'Surface':
                    print(f'Mesh/Surface {name}: n={len(items)} types={sorted({i.tag for i in items})} hash={h(raw)}')
                elif child.tag == 'SurfacePair':
                    print(f'Mesh/SurfacePair {name}: {child.find("primary").text} | {child.find("secondary").text}')
                elif child.tag == 'DiscreteSet' and collapsed(name):
                    coll.append(raw)
                elif child.tag == 'DiscreteSet':
                    print(f'Mesh/DiscreteSet {name}: n={len(items)} hash={h(raw)}')
                else:
                    print(f'Mesh/{child.tag} {name}: hash={h(raw)}')
            if coll:
                print(f'Mesh/DiscreteSet (collapsed {a.collapse}): count={len(coll)} hash={h("".join(coll))}')
        elif sec.tag == 'Discrete':
            coll = [c for c in sec if collapsed(c.get('name') or c.get('discrete_set'))]
            if coll:
                print(f'Discrete (collapsed {a.collapse}): entries={len(coll)} '
                      f'hash={h("".join(canon(c) for c in coll))}')
            for c in sec:
                if c not in coll:
                    for line in flat(c):
                        print('Discrete/' + line)
        elif sec.tag == 'Constraints':
            for c in sec:
                lcs = c.findall('linear_constraint')
                params = [p for p in c if p.tag != 'linear_constraint']
                print(f'Constraints/<constraint name={c.get("name")} type={c.get("type")}> '
                      f'linear_constraint count={len(lcs)} '
                      f'hash={h("".join(canon(x) for x in lcs))}')
                for p in params:
                    for line in flat(p, '  '):
                        print('Constraints/' + line)
        else:
            for line in flat(sec):
                print(line)


if __name__ == '__main__':
    main()
