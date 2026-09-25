"""Print a normalized structural summary of a .feb (for diffing model versions).

Mesh blocks are summarized by name/count/id-range plus a hash of their content,
so coordinate-level changes show up as a hash change without drowning the diff.
"""
import hashlib
import sys
import xml.etree.ElementTree as ET


def h(s):
    return hashlib.md5(s.encode()).hexdigest()[:8]


def flat(el, indent=''):
    """Flatten an element's params (non-Mesh sections) to sorted lines."""
    out = []
    attrs = ' '.join(f'{k}={v}' for k, v in sorted(el.attrib.items()))
    txt = (el.text or '').strip()
    if len(el):
        out.append(f'{indent}<{el.tag} {attrs}>')
        for c in el:
            out.extend(flat(c, indent + '  '))
    else:
        out.append(f'{indent}{el.tag} {attrs} = {txt}')
    return out


def main(path):
    root = ET.parse(path).getroot()
    for sec in root:
        if sec.tag == 'Mesh':
            for child in sec:
                name = child.get('name')
                items = list(child)
                raw = ET.tostring(child, encoding='unicode')
                if child.tag == 'Nodes':
                    ids = [int(n.get('id')) for n in items]
                    print(f'Mesh/Nodes {name}: n={len(ids)} ids={min(ids)}-{max(ids)} hash={h(raw)}')
                elif child.tag == 'Elements':
                    ids = [int(n.get('id')) for n in items]
                    print(f'Mesh/Elements {name} type={child.get("type")}: n={len(ids)} ids={min(ids)}-{max(ids)} hash={h(raw)}')
                elif child.tag == 'NodeSet':
                    txt = (child.text or '').strip()
                    n = len([v for v in txt.replace('\n', ',').split(',') if v.strip()]) + len(items)
                    print(f'Mesh/NodeSet {name}: n={n} hash={h(raw)}')
                elif child.tag == 'Surface':
                    print(f'Mesh/Surface {name}: n={len(items)} types={sorted(set(i.tag for i in items))} hash={h(raw)}')
                elif child.tag == 'SurfacePair':
                    print(f'Mesh/SurfacePair {name}: {child.find("primary").text} | {child.find("secondary").text}')
                elif child.tag == 'DiscreteSet':
                    if name.startswith('stab_'):
                        continue
                    print(f'Mesh/DiscreteSet {name}: n={len(items)} hash={h(raw)}')
                else:
                    print(f'Mesh/{child.tag} {name}: hash={h(raw)}')
            nst = sum(1 for c in sec if c.tag == 'DiscreteSet' and c.get('name', '').startswith('stab_'))
            print(f'Mesh/DiscreteSet stab_*: count={nst}')
        elif sec.tag == 'Discrete':
            stab = [c for c in sec if 'stab_' in (c.get('name') or c.get('discrete_set') or '')]
            print(f'Discrete: stab_* entries={len(stab)} hash={h("".join(ET.tostring(c, encoding="unicode") for c in stab))}')
            for c in sec:
                if c in stab:
                    continue
                for line in flat(c):
                    print('Discrete/' + line)
        else:
            for line in flat(sec):
                print(line)


if __name__ == '__main__':
    main(sys.argv[1])
