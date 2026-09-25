"""Minimal FEBio .xplt (FEBio 3/4 plot file) reader for post-mortem diagnostics.

Reads the dictionary, the mesh (node coordinates + domains) and an index of every
state, then loads nodal / element variables for any single state on demand, so a
multi-GB per-iteration plot file never has to fit in memory. Safe to run on a
.xplt that FEBio is still writing (a partially written last state is ignored).

Things that are easy to get wrong (all verified against FEBio 4.13 output):
  * Element connectivity in the plot file is 0-based *node indices*; map to
    FEBio node IDs with `x.node_ids[index]`.
  * Element data (`relative volume`, `stress`, ...) comes back as
    {region_id: array}, where region_id is the 1-based domain index in
    `x.domains` order (domain i -> key i+1). Domains appear in <MeshDomains>
    order, followed by one region per <discrete> binding. Shell domain names are
    often missing from the file, so take names from the .feb (feb_postmortem.py
    does this).
  * State status flags: FEBioStudio "debug" runs mark converged states 0 and
    iteration states 2. A plain `febio4 -i` run with PLOT_MINOR_ITRS marks
    *every* state 2, so use `converged_state_indices(x, log_path)`, which
    matches the log's "converged at time" values instead.

Usage:  py -3 xplt_reader.py model.xplt      (prints dictionary, domains, state count)
Needs:  Python 3 + numpy.
"""
import os
import re
import struct
import sys
import zlib

import numpy as np

FEBIO_TAG = 0x00464542
PLT_ROOT = 0x01000000
PLT_HEADER = 0x01010000
PLT_HDR_VERSION = 0x01010001
PLT_HDR_COMPRESSION = 0x01010004
PLT_DICTIONARY = 0x01020000
PLT_DIC_ITEM = 0x01020001
PLT_DIC_ITEM_TYPE = 0x01020002
PLT_DIC_ITEM_FMT = 0x01020003
PLT_DIC_ITEM_NAME = 0x01020004
PLT_DIC_GLOBAL = 0x01021000
PLT_DIC_NODAL = 0x01023000
PLT_DIC_DOMAIN = 0x01024000
PLT_DIC_SURFACE = 0x01025000
PLT_DIC_EDGE = 0x01026000
PLT_MESH = 0x01040000
PLT_NODE_SECTION = 0x01041000
PLT_NODE_HEADER = 0x01041100
PLT_NODE_SIZE = 0x01041101
PLT_NODE_DIM = 0x01041102
PLT_NODE_COORDS = 0x01041200
PLT_DOMAIN_SECTION = 0x01042000
PLT_DOMAIN = 0x01042100
PLT_DOMAIN_HDR = 0x01042101
PLT_DOM_ELEM_TYPE = 0x01042102
PLT_DOM_PART_ID = 0x01042103
PLT_DOM_ELEMS = 0x01032104
PLT_DOM_NAME = 0x01032105
PLT_DOM_ELEM_LIST = 0x01042200
PLT_ELEMENT = 0x01042201
PLT_STATE = 0x02000000
PLT_STATE_HEADER = 0x02010000
PLT_STATE_HDR_TIME = 0x02010002
PLT_STATE_STATUS = 0x02010003
PLT_STATE_DATA = 0x02020000
PLT_STATE_VARIABLE = 0x02020001
PLT_STATE_VAR_ID = 0x02020002
PLT_STATE_VAR_DATA = 0x02020003
PLT_GLOBAL_DATA = 0x02020100
PLT_NODE_DATA = 0x02020300
PLT_ELEMENT_DATA = 0x02020400
PLT_FACE_DATA = 0x02020500

# dictionary item type -> floats per item (FLOAT, VEC3F, MAT3FS, MAT3FD, TENS4FS, MAT3F)
TYPE_SIZE = {0: 1, 1: 3, 2: 6, 3: 3, 4: 21, 5: 9}


def chunks(buf, off=0, end=None):
    end = len(buf) if end is None else end
    while off + 8 <= end:
        cid, size = struct.unpack_from('<II', buf, off)
        yield cid, off + 8, size
        off += 8 + size


def _name(raw):
    """FEBio 4 writes length-prefixed strings; FEBio 3 writes fixed 64-char fields."""
    if len(raw) >= 4:
        n = struct.unpack_from('<I', raw, 0)[0]
        if 0 < n <= len(raw) - 4:
            return raw[4:4 + n].split(b'\0')[0].decode(errors='replace')
    return raw.split(b'\0')[0].decode(errors='replace')


class Xplt:
    def __init__(self, path):
        self.path = path
        self.f = open(path, 'rb')
        tag = struct.unpack('<I', self.f.read(4))[0]
        if tag != FEBIO_TAG:
            raise ValueError(f'{path}: not a FEBio plot file (tag {tag:#x})')
        self.compressed = 0
        self.version = 0
        self.dict = {'global': [], 'nodal': [], 'domain': [], 'surface': []}
        self.domains = []     # dicts: etype, ne, eids (array), conn (0-based node indices), name?
        self.states = []      # (time, status, data_offset, data_size)
        self._read_root()
        self._index_states()

    def _hdr(self):
        b = self.f.read(8)
        return struct.unpack('<II', b) if len(b) == 8 else None

    def _read_root(self):
        cid, size = self._hdr()
        assert cid == PLT_ROOT, hex(cid)
        buf = self.f.read(size)
        for cid, o, s in chunks(buf):
            if cid == PLT_HEADER:
                for c2, o2, s2 in chunks(buf, o, o + s):
                    if c2 == PLT_HDR_VERSION:
                        self.version = struct.unpack_from('<I', buf, o2)[0]
                    elif c2 == PLT_HDR_COMPRESSION:
                        self.compressed = struct.unpack_from('<I', buf, o2)[0]
            elif cid == PLT_DICTIONARY:
                self._read_dict(buf, o, o + s)
        cid, size = self._hdr()          # FEBio 3+: mesh chunk follows the root chunk
        assert cid == PLT_MESH, hex(cid)
        self._read_mesh(self.f.read(size))
        self.first_state_off = self.f.tell()

    def _read_dict(self, buf, o, e):
        keys = {PLT_DIC_GLOBAL: 'global', PLT_DIC_NODAL: 'nodal', PLT_DIC_DOMAIN: 'domain',
                PLT_DIC_SURFACE: 'surface', PLT_DIC_EDGE: 'edge'}
        for cid, o2, s2 in chunks(buf, o, e):
            key = keys.get(cid)
            if key is None:
                continue
            self.dict.setdefault(key, [])
            for c3, o3, s3 in chunks(buf, o2, o2 + s2):
                if c3 != PLT_DIC_ITEM:
                    continue
                item = {}
                for c4, o4, s4 in chunks(buf, o3, o3 + s3):
                    if c4 == PLT_DIC_ITEM_TYPE:
                        item['type'] = struct.unpack_from('<I', buf, o4)[0]
                    elif c4 == PLT_DIC_ITEM_FMT:
                        item['fmt'] = struct.unpack_from('<I', buf, o4)[0]
                    elif c4 == PLT_DIC_ITEM_NAME:
                        item['name'] = _name(buf[o4:o4 + s4])
                self.dict[key].append(item)

    def _read_mesh(self, buf):
        for cid, o, s in chunks(buf):
            if cid == PLT_NODE_SECTION:
                nn = dim = 0
                for c2, o2, s2 in chunks(buf, o, o + s):
                    if c2 == PLT_NODE_HEADER:
                        for c3, o3, s3 in chunks(buf, o2, o2 + s2):
                            if c3 == PLT_NODE_SIZE:
                                nn = struct.unpack_from('<I', buf, o3)[0]
                            elif c3 == PLT_NODE_DIM:
                                dim = struct.unpack_from('<I', buf, o3)[0]
                    elif c2 == PLT_NODE_COORDS:
                        rec = np.frombuffer(buf, dtype=np.float32, count=nn * (dim + 1),
                                            offset=o2).reshape(nn, dim + 1)
                        self.node_ids = rec[:, 0].view(np.int32).copy()
                        self.X = rec[:, 1:].astype(float)
                self.nn = nn
            elif cid == PLT_DOMAIN_SECTION:
                for c2, o2, s2 in chunks(buf, o, o + s):
                    if c2 != PLT_DOMAIN:
                        continue
                    dom = {}
                    for c3, o3, s3 in chunks(buf, o2, o2 + s2):
                        if c3 == PLT_DOMAIN_HDR:
                            for c4, o4, s4 in chunks(buf, o3, o3 + s3):
                                if c4 == PLT_DOM_ELEM_TYPE:
                                    dom['etype'] = struct.unpack_from('<I', buf, o4)[0]
                                elif c4 == PLT_DOM_ELEMS:
                                    dom['ne'] = struct.unpack_from('<I', buf, o4)[0]
                                elif c4 == PLT_DOM_NAME:
                                    dom['name'] = _name(buf[o4:o4 + s4])
                        elif c3 == PLT_DOM_ELEM_LIST:
                            ids, conns = [], []
                            for c4, o4, s4 in chunks(buf, o3, o3 + s3):
                                if c4 == PLT_ELEMENT:
                                    arr = np.frombuffer(buf, dtype=np.int32, count=s4 // 4, offset=o4)
                                    ids.append(int(arr[0]))
                                    conns.append(arr[1:].copy())
                            dom['eids'] = np.array(ids)
                            dom['conn'] = np.array(conns)
                    self.domains.append(dom)

    def _index_states(self):
        f = self.f
        fsize = os.fstat(f.fileno()).st_size
        f.seek(self.first_state_off)
        while True:
            pos = f.tell()
            h = self._hdr()
            if h is None:
                break
            cid, size = h
            if pos + 8 + size > fsize:            # state still being written
                break
            if cid != PLT_STATE:
                f.seek(pos + 8 + size)
                continue
            t = status = off = dsize = None
            end = pos + 8 + size
            while f.tell() + 8 <= end:
                c2, s2 = self._hdr()
                here = f.tell()
                if c2 == PLT_STATE_HEADER:
                    hb = f.read(s2)
                    for c3, o3, s3 in chunks(hb):
                        if c3 == PLT_STATE_HDR_TIME:
                            t = struct.unpack_from('<f', hb, o3)[0]
                        elif c3 == PLT_STATE_STATUS:
                            status = struct.unpack_from('<I', hb, o3)[0]
                elif c2 == PLT_STATE_DATA:
                    off, dsize = here, s2
                f.seek(here + s2)
            self.states.append((t, status, off, dsize))
            f.seek(end)

    def var(self, i, name):
        """Variable `name` at state i: nodal -> (nn, k) array; domain/surface -> {region_id: array}."""
        t, st, off, size = self.states[i]
        self.f.seek(off)
        buf = self.f.read(size)
        for kind, sect in (('nodal', PLT_NODE_DATA), ('domain', PLT_ELEMENT_DATA),
                           ('surface', PLT_FACE_DATA), ('global', PLT_GLOBAL_DATA)):
            items = self.dict.get(kind, [])
            names = [it.get('name') for it in items]
            if name not in names:
                continue
            want = names.index(name) + 1
            k = TYPE_SIZE.get(items[want - 1].get('type'), 1)
            for cid, o, s in chunks(buf):
                if cid != sect:
                    continue
                for c2, o2, s2 in chunks(buf, o, o + s):
                    if c2 != PLT_STATE_VARIABLE:
                        continue
                    vid = None
                    for c3, o3, s3 in chunks(buf, o2, o2 + s2):
                        if c3 == PLT_STATE_VAR_ID:
                            vid = struct.unpack_from('<I', buf, o3)[0]
                        elif c3 == PLT_STATE_VAR_DATA and vid == want:
                            raw = buf[o3:o3 + s3]
                            if self.compressed:
                                raw = zlib.decompress(raw)
                            out, p = {}, 0
                            while p + 8 <= len(raw):
                                rid, rs = struct.unpack_from('<II', raw, p)
                                arr = np.frombuffer(raw, dtype=np.float32, count=rs // 4, offset=p + 8)
                                out[rid] = arr.reshape(-1, k) if k > 1 else arr
                                p += 8 + rs
                            if kind == 'nodal':
                                return out[0] if 0 in out else next(iter(out.values()))
                            return out
            return None
        raise KeyError(f'{name!r} not in plot dictionary {self.dict}')


def converged_state_indices(x, log_path=None):
    """Indices of converged states. Uses status flags when the file has them (FEBioStudio
    debug runs), otherwise the last state at each 'converged at time' value in the log."""
    flagged = [i for i, s in enumerate(x.states) if s[1] == 0]
    if len(flagged) > 1 or log_path is None:
        return flagged
    log = open(log_path, encoding='latin-1').read()
    ts = [float(v) for v in re.findall(r'^------- converged at time : ([0-9.eE+-]+)', log, re.M)]
    times = np.array([s[0] for s in x.states])
    out = []
    for t in ts:
        hit = np.where(np.abs(times - t) <= 1e-6 * max(abs(t), 1e-3))[0]
        if len(hit):
            out.append(int(hit.max()))
    return sorted(set(out)) or flagged


if __name__ == '__main__':
    x = Xplt(sys.argv[1])
    print(f'version {x.version:#x}, compressed={x.compressed}, nodes={x.nn}, domains={len(x.domains)}')
    for k, v in x.dict.items():
        print(f'  {k}: {[it.get("name") for it in v]}')
    flags = {}
    for s in x.states:
        flags[s[1]] = flags.get(s[1], 0) + 1
    print(f'states: {len(x.states)} (status flag counts {flags}); '
          f'time {x.states[0][0]:.6g} .. {x.states[-1][0]:.6g}' if x.states else 'no states')
