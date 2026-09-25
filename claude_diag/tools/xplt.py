"""Minimal FEBio .xplt reader (FEBio 3/4 plot format) for post-mortem diagnostics.

Reads the dictionary, mesh (nodes + domains), and an index of states; can then
load nodal / element variables for any state on demand without reading the whole
(possibly multi-GB) file into memory.
"""
import os
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
PLT_PARTS_SECTION = 0x01046000
PLT_STATE = 0x02000000
PLT_STATE_HEADER = 0x02010000
PLT_STATE_HDR_ID = 0x02010001
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

ELEM_NODES = {0: 8, 1: 6, 2: 4, 3: 4, 4: 3, 5: 2, 6: 20, 7: 10, 8: 6, 9: 8, 10: 15, 11: 27, 12: 3, 13: 6,
              14: 9, 15: 5, 16: 10, 17: 15, 18: 13}
TYPE_SIZE = {0: 1, 1: 3, 2: 6, 3: 9, 4: 9, 5: 6, 6: 3, 7: 9}   # FLOAT VEC3F MAT3FS MAT3FD TENS4FS MAT3F ...


def chunks(buf, off=0, end=None):
    end = len(buf) if end is None else end
    while off + 8 <= end:
        cid, size = struct.unpack_from('<II', buf, off)
        yield cid, off + 8, size
        off += 8 + size


class Xplt:
    def __init__(self, path):
        self.path = path
        self.f = open(path, 'rb')
        tag = struct.unpack('<I', self.f.read(4))[0]
        assert tag == FEBIO_TAG, hex(tag)
        self.compressed = 0
        self.dict = {'global': [], 'nodal': [], 'domain': [], 'surface': []}
        self.domains = []
        self.states = []   # (time, status, offset_of_state_data, size)
        self._read_root()
        self._index_states()

    def _read_chunk_header(self):
        b = self.f.read(8)
        if len(b) < 8:
            return None
        return struct.unpack('<II', b)

    def _read_root(self):
        cid, size = self._read_chunk_header()
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
        # mesh follows the root chunk in FEBio 3+
        cid, size = self._read_chunk_header()
        assert cid == PLT_MESH, hex(cid)
        mbuf = self.f.read(size)
        self._read_mesh(mbuf)
        self.first_state_off = self.f.tell()

    def _read_dict(self, buf, o, e):
        keymap = {PLT_DIC_GLOBAL: 'global', PLT_DIC_NODAL: 'nodal', PLT_DIC_DOMAIN: 'domain',
                  PLT_DIC_SURFACE: 'surface', PLT_DIC_EDGE: 'edge'}
        for cid, o2, s2 in chunks(buf, o, e):
            key = keymap.get(cid)
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
                        raw = buf[o4:o4 + s4]
                        # FEBio 4 stores length-prefixed strings; FEBio 3 stores 64-char
                        if s4 >= 4:
                            n = struct.unpack_from('<I', raw, 0)[0]
                            if 0 < n <= s4 - 4:
                                item['name'] = raw[4:4 + n].split(b'\0')[0].decode(errors='replace')
                            else:
                                item['name'] = raw.split(b'\0')[0].decode(errors='replace')
                        else:
                            item['name'] = raw.split(b'\0')[0].decode(errors='replace')
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
                        rec = np.frombuffer(buf, dtype=np.float32, count=nn * (dim + 1), offset=o2).reshape(nn, dim + 1)
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
                                elif c4 == PLT_DOM_PART_ID:
                                    dom['part'] = struct.unpack_from('<I', buf, o4)[0]
                                elif c4 == PLT_DOM_ELEMS:
                                    dom['ne'] = struct.unpack_from('<I', buf, o4)[0]
                                elif c4 == PLT_DOM_NAME:
                                    raw = buf[o4:o4 + s4]
                                    n = struct.unpack_from('<I', raw, 0)[0] if s4 >= 4 else 0
                                    if 0 < n <= s4 - 4:
                                        dom['name'] = raw[4:4 + n].split(b'\0')[0].decode(errors='replace')
                                    else:
                                        dom['name'] = raw.split(b'\0')[0].decode(errors='replace')
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
            h = self._read_chunk_header()
            if h is None:
                break
            cid, size = h
            if pos + 8 + size > fsize:     # state still being written
                break
            if cid != PLT_STATE:
                f.seek(pos + 8 + size)
                continue
            # read state header only
            t = None; status = None; data_off = None; data_size = None
            end = pos + 8 + size
            while f.tell() + 8 <= end:
                c2, s2 = self._read_chunk_header()
                here = f.tell()
                if c2 == PLT_STATE_HEADER:
                    hb = f.read(s2)
                    for c3, o3, s3 in chunks(hb):
                        if c3 == PLT_STATE_HDR_TIME:
                            t = struct.unpack_from('<f', hb, o3)[0]
                        elif c3 == PLT_STATE_STATUS:
                            status = struct.unpack_from('<I', hb, o3)[0]
                elif c2 == PLT_STATE_DATA:
                    data_off, data_size = here, s2
                f.seek(here + s2)
            self.states.append((t, status, data_off, data_size))
            f.seek(end)

    def _state_data(self, i):
        t, st, off, size = self.states[i]
        self.f.seek(off)
        buf = self.f.read(size)
        return buf

    def var(self, i, name):
        """Return data for a named variable at state i.

        nodal -> array (nn, k); domain -> dict domain_index -> array (ne, k)
        """
        buf = self._state_data(i)
        for kind, sect in (('nodal', PLT_NODE_DATA), ('domain', PLT_ELEMENT_DATA), ('surface', PLT_FACE_DATA),
                           ('global', PLT_GLOBAL_DATA)):
            items = self.dict.get(kind, [])
            names = [it.get('name') for it in items]
            if name not in names:
                continue
            want = names.index(name) + 1
            item = items[want - 1]
            k = TYPE_SIZE.get(item['type'], 1)
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
                            out = {}
                            p = 0
                            while p + 8 <= len(raw):
                                rid, rs = struct.unpack_from('<II', raw, p)
                                arr = np.frombuffer(raw, dtype=np.float32, count=rs // 4, offset=p + 8)
                                out[rid] = arr.reshape(-1, k) if k > 1 else arr
                                p += 8 + rs
                            if kind == 'nodal':
                                return out[0] if 0 in out else next(iter(out.values()))
                            return out
            return None
        raise KeyError(name)


if __name__ == '__main__':
    x = Xplt(sys.argv[1])
    print('version', hex(getattr(x, 'version', 0)), 'compressed', x.compressed)
    for k, v in x.dict.items():
        print(k, [(it.get('name'), it.get('type'), it.get('fmt')) for it in v])
    print('nodes', x.nn, 'domains', len(x.domains))
    for i, d in enumerate(x.domains):
        print(' ', i, d.get('name'), d.get('etype'), d.get('ne'), d['conn'].shape if 'conn' in d else None)
    print('states', len(x.states))
    for s in x.states[:5] + x.states[-5:]:
        print('  t=%.6g status=%s' % (s[0], s[1]))
