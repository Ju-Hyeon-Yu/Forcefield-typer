# AmberTools antechamber 의 GAFF2 원자 타입 판정(bondtype -j part → atomtype -p gaff2)을 옮긴 타이퍼
"""GAFF2 타이퍼 — `antechamber -at gaff2 -j 5` 와 같은 결과를 내도록 처리 순서를 그대로 옮겼다.

  1. 입력 결합 차수(SMILES 의 케쿨레 구조)를 유지한 채 방향족(7, 8)·비편재(9) 결합만 다시 지정
     (bondtype.c 의 finalize, -j part)
  2. 고리 찾기와 방향성 분류 AR1~AR5 (ring.c 의 ringdetect, aromatic)
  3. ATOMTYPE_GFF2.DEF 규칙을 위에서부터 맞춰 첫 번째로 맞는 타입 (atomtype.c 의 jat, apcheck,
     jatspecial, cematch)
  4. 공액 타입 쌍 조정 cc/cd, ce/cf, cg/ch, nc/nd, ne/nf, pc/pd, pe/pf, cp/cq (atadjust, cpadjust)

원본: AmberTools / AmberClassic antechamber (Junmei Wang 외), GPL-2.0-or-later.
규칙 파일 data/ATOMTYPE_GFF2.DEF 는 AmberTools 26 의 것을 그대로 담았다.
웹 이식(web/src/gaff2.js)과 결과를 맞추기 위해 반복 순서까지 원본과 같게 둔다.
"""
import itertools
import json
import os
import re

from rdkit import Chem

from . import TypingError, moltemplate_ff_path

DEF_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "ATOMTYPE_GFF2.DEF")
_PT = Chem.GetPeriodicTable()


# --------------------------------------------------------------------------- 규칙 파일 해석
def _parse_ces(keyword):
    """화학 환경 문자열 → 사슬 목록 (jatspecial 의 문자 단위 해석을 그대로 옮김).

    사슬 = [{name, connum, ap, cesname}, ...] (중심 원자에서 바깥으로)
    """
    size = 22
    atname, atconnum = [""] * size, [0] * size
    ap, cesname = [""] * size, [""] * size
    chains = []
    layer, tmpapindex, tmpap = 0, 0, ""
    cesname_index, tmpcesname = 0, ""
    index0, cea_id = 0, 1

    def alpha(c):
        return ("a" <= c <= "z") or ("A" <= c <= "Z")

    def emit(last, length):
        nodes = [{"name": atname[j + 1], "connum": atconnum[j + 1], "ap": ap[j + 1],
                  "cesname": cesname[j + 1]} for j in range(last + 1)]
        chains.append(nodes[:length])

    n = len(keyword)
    for i in range(n):
        c = keyword[i]
        nxt = keyword[i + 1] if i + 1 < n else ""
        prv = keyword[i - 1] if i > 0 else ""
        if tmpapindex == 0 and cesname_index == 0 and alpha(c) and alpha(nxt):
            continue
        if c == "(":
            layer += 1
        if c == ")":
            layer -= 1
        if tmpapindex == 0 and c == "[":
            tmpapindex, tmpap = 1, "["
            continue
        if tmpapindex == 1 and c == "]":
            ap[layer] = tmpap + "]"
            tmpapindex = 0
            continue
        if tmpapindex == 1:
            tmpap += c
            continue
        if cesname_index == 0 and c == "<":
            cesname_index, tmpcesname = 1, ""
            continue
        if cesname_index == 1 and c == ">":
            cesname[layer] = tmpcesname
            cesname_index = 0
            continue
        if cesname_index == 1:
            tmpcesname += c
            continue
        if c == "," and prv != ")":
            emit(layer, layer)
        if c == ")" and prv != ")":
            emit(layer + 1, layer + 1)
        if alpha(c) and alpha(nxt):
            continue
        if alpha(c):
            index0 = 1
            atname[layer] = (prv + c) if alpha(prv) else c
            ap[layer] = ""
            cesname[layer] = f"cea_{cea_id}"
            cea_id += 1
        if "0" <= c <= "9":
            atconnum[layer] = int(c)
        elif index0 == 1:
            atconnum[layer] = 0
            index0 = 0
    return chains


def _parse_ap(text):
    """원자 성질 문자열 '[a.b,c]' → [[(개수, 이름), ...], ...] (',' 는 그리고, '.' 는 또는)."""
    out = []
    for cons in text.strip("[]").split(","):
        alts = []
        for unit in cons.split("."):
            m = re.match(r"(\d*)(.*)", unit)
            alts.append((int(m.group(1)) if m.group(1) else -1, m.group(2)))
        out.append(alts)
    return out


def parse_def(text):
    """ATOMTYPE_GFF2.DEF → {"wild": {이름: [원자번호...]}, "rules": [...]} (웹 이식에서도 그대로 씀)."""
    wild, rules = {}, []
    for line in text.splitlines():
        if line.startswith("WILDATOM"):
            t = line.split()
            wild[t[1]] = [_PT.GetAtomicNumber(re.sub(r"[^A-Za-z]", "", s)) for s in t[2:]]
        elif line.startswith("ATD"):
            t = line.split()[1:]
            fields = []
            for f in t[1:8]:
                fields.append(f)
                if f == "&":
                    break
            rule = {"type": t[0], "res": fields[0]}
            keys = ("z", "connum", "nh", "new", "ap", "ces")
            for k, f in zip(keys, fields[1:]):
                if f == "&":
                    break
                if k in ("z", "connum", "nh", "new"):
                    rule[k] = None if f == "*" else int(f)
                elif k == "ap":
                    rule[k] = None if f == "*" else _parse_ap(f)
                else:
                    rule[k] = _parse_ces(f)
            rules.append(rule)
    return {"wild": wild, "rules": rules}


_DEF = None


def default_def():
    global _DEF
    if _DEF is None:
        with open(DEF_PATH) as f:
            _DEF = parse_def(f.read())
    return _DEF


# --------------------------------------------------------------------------- 분자 그래프
class _Mol:
    """antechamber 의 ATOM/BOND 배열 대응 — con 은 결합 순서대로, bonds 는 [i, j, 형식]."""

    def __init__(self, z, bonds):
        self.n = len(z)
        self.z = list(z)
        self.bonds = [list(b) for b in bonds]
        self.con = [[] for _ in z]
        for i, j, _ in self.bonds:
            self.con[i].append(j)
            self.con[j].append(i)
        self.connum = [len(c) for c in self.con]

    def bond_type(self, a, b):
        for i, j, t in self.bonds:
            if (i == a and j == b) or (i == b and j == a):
                return t
        return None


def _ring_ok(m, a):
    z, c = m.z[a], m.connum[a]
    if z not in (6, 7, 8, 15, 16):
        return False
    if z == 6 and c <= 2:
        return False
    if z in (8, 16) and c == 1:
        return False
    return True


def _find_rings(m):
    """ring_detect_cycle + purify: 이웃 앞 4개만 따라가는 3~10원 고리 (현이 있는 고리는 뺌)."""
    found, path = [], []

    def cycle(start):
        path.append(start)
        sel = len(path)
        for nb in m.con[start][:4]:
            if not _ring_ok(m, nb) or nb in path:
                continue
            if sel > 10:
                break
            if 2 <= sel <= 9 and nb in m.con[path[0]][:4]:
                found.append(path + [nb])
            cycle(nb)
        path.pop()

    for a in range(m.n):
        if _ring_ok(m, a):
            cycle(a)
    rings, seen = [], set()
    for r in found:
        key = tuple(sorted(r))
        if key in seen:
            continue
        seen.add(key)
        members = set(r)
        if any(sum(1 for x in m.con[a][:6] if x in members) == 3 for a in key):
            continue
        rings.append(list(key))
    return rings


def _ringdetect(m, rings, amber_bcc):
    """ringproperty + aromatic (+ amber_bcc 이면 인돌형 5원 고리 보정) → 원자별 고리·방향성 성질."""
    n = m.n
    rg = [[0] * 11 for _ in range(n)]
    for r in rings:
        for a in r:
            rg[a][0] += 1
            rg[a][len(r)] += 1
    init, ewd = [0] * n, [0] * n
    for a in range(n):
        z, c = m.z[a], m.connum[a]
        if z == 6:
            init[a] = 2 if c == 3 else (-2 if c == 4 else 0)
        elif z == 7:
            ewd[a] = 1
            init[a] = 2 if c <= 3 else 0
        elif z == 8:
            ewd[a] = 1
            init[a] = 1 if c == 2 else 0
        elif z == 15:
            init[a] = 2 if c == 2 else (1 if c == 3 else (-1 if c >= 4 else 0))
        elif z == 16:
            ewd[a] = 1
            init[a] = 1 if c == 2 else (-1 if c >= 3 else 0)
        elif z in (9, 17, 35, 53):
            ewd[a] = 1
    ar = [[0] * 6 for _ in range(n)]          # ar[a][1..5]
    for r in rings:
        size, tot = len(r), sum(init[a] for a in r)
        if tot == -2 * size:
            kind = 5
        elif any(init[a] < 0 for a in r):
            kind = 4
        else:
            kind = None
            if size <= tot <= 2 * size:
                for i, j, t in m.bonds:
                    idx = sum(1 for a in r if i == a and rg[j][0] == 0) + \
                        sum(1 for a in r if j == a and rg[i][0] == 0)
                    if idx == 1 and t in (2, 8):
                        kind = 3
                        break
            if kind is None and tot == 12 and size == 6:
                ok = True
                for a in r:
                    if m.z[a] in (7, 15) and not any(t in (8, 2, 10) and (i == a or j == a)
                                                     for i, j, t in m.bonds):
                        ok = False
                if ok:
                    kind = 1
            if kind is None:
                kind = 2 if tot >= size + 3 else 4
        for a in r:
            ar[a][kind] += 1
    if amber_bcc:
        for a in range(n):
            if rg[a][5] >= 1 and rg[a][6] >= 1 and ar[a][1] + ar[a][2] >= 2:
                for r in rings:
                    if len(r) == 5 and a in r:
                        for b in r:
                            ar[b][2] -= 1
    nr = [0 if any(ar[a][k] > 0 for k in range(1, 6)) else 1 for a in range(n)]
    return {"rg": rg, "ar": ar, "nr": nr, "ewd": ewd}


def _term(m, a, z):
    """원자 a 가 원자번호 z 의 말단(이웃 1개) 원자인지."""
    return m.z[a] == z and m.connum[a] == 1


def _bondtype_part(m, rings):
    """bondtype -j part 의 finalize: 방향족 고리 결합 7/8, 비편재 결합 9, N–O 6."""
    prop = _ringdetect(m, rings, amber_bcc=True)
    ar = prop["ar"]
    conj = [0] * m.n
    for a in range(m.n):
        z, c, con = m.z[a], m.connum[a], m.con[a]
        tO = lambda k: _term(m, con[k], 8) or _term(m, con[k], 16)  # noqa: E731
        if z == 6 and c == 3 and sum(tO(k) for k in range(3)) >= 2:
            conj[a] = 1
        if z == 15 and c == 4 and sum(tO(k) for k in range(4)) >= 2:
            conj[a] = 1
        if z == 16 and c == 4 and sum(tO(k) for k in range(4)) >= 3:
            conj[a] = 1
        if z == 7 and c == 3 and sum(tO(k) for k in range(3)) >= 2:
            conj[a] = 1
    termOS = lambda a: _term(m, a, 8) or _term(m, a, 16)  # noqa: E731
    for b in m.bonds:
        bi, bj = b[0], b[1]
        # part1: 5·6원 방향족 고리 안의 결합
        if (ar[bi][1] > 0 and ar[bj][1] > 0) or (ar[bi][2] > 0 and ar[bj][2] > 0):
            for r in rings:
                if len(r) >= 7 or len(r) <= 4:
                    continue
                if any(ar[a][1] <= 0 and ar[a][2] <= 0 for a in r):
                    continue
                if bi in r and bj in r:
                    if b[2] == 1:
                        b[2] = 7
                    elif b[2] == 2:
                        b[2] = 8
                    break
        # part2: 카복실레이트·니트로·술포네이트 등의 =O → 비편재
        if b[2] == 2:
            if conj[bi] == 1 and termOS(bj):
                b[2] = 9
                continue
            if conj[bj] == 1 and termOS(bi):
                b[2] = 9
                continue
        # part3: N–O(말단)
        if m.z[bi] == 7 and m.connum[bi] in (2, 3) and termOS(bj):
            for x in m.con[bi][:6]:
                if x == bj:
                    continue
                if termOS(x):
                    b[2] = 6
                break
            if b[2] == 6:
                continue
        if m.z[bj] == 7 and m.connum[bj] in (2, 3) and termOS(bi):
            b[2] = 6
            continue
        # part4: 말단 O/S 와의 단일결합 → 비편재
        if b[2] == 1 and (termOS(bi) or termOS(bj)):
            b[2] = 9
            continue
    # part5: 3배위 S 의 말단 O/S 두 개, 4배위 P 의 말단 O/S 두 개 이상
    for a in range(m.n):
        if m.z[a] == 16 and m.connum[a] == 3:
            sel = [x for x in m.con[a][:3] if termOS(x)]
            if len(sel) == 2:
                for b in m.bonds:
                    if (b[0] == a and b[1] in sel) or (b[1] == a and b[0] in sel):
                        b[2] = 9
        if m.z[a] == 15 and m.connum[a] == 4:
            sel = [x for x in m.con[a][:4] if termOS(x)]
            if len(sel) >= 2:
                for b in m.bonds:
                    if (b[0] == a and b[1] in sel) or (b[1] == a and b[0] in sel):
                        b[2] = 9


# --------------------------------------------------------------------------- 원자 타입 판정
_JBOND = {1: (1, 7, 9, 10), -1: (1, 9), 2: (2, 8, 9, 10), -2: (2, 9), 3: (3,), -3: (3,),
          4: (7, 8, 10), 9: (9,)}
_BOND_PROP = {"SB": ("SB", -1), "sb": ("sb", 1), "DB": ("DB", -2), "db": ("db", 2),
              "TB": ("TB", -3), "tb": ("tb", 3), "DL": ("DL", 9), "AB": ("AB", 4)}


class _Typer:
    def __init__(self, m, rules):
        self.m = m
        self.rules = rules
        self.wild = rules["wild"]
        rings = _find_rings(m)
        _bondtype_part(m, rings)
        prop = _ringdetect(m, rings, amber_bcc=False)
        self.rg, self.ar, self.nr, self.ewd = prop["rg"], prop["ar"], prop["nr"], prop["ewd"]
        cnt = {k: [0] * m.n for k in ("sb", "SB", "db", "DB", "tb", "TB", "AB", "DL")}
        add = {1: ("sb", "SB"), 2: ("db", "DB"), 3: ("tb", "TB"), 7: ("AB", "sb"), 8: ("AB", "db"),
               9: ("sb", "SB", "DL"), 10: ("AB",)}
        for i, j, t in m.bonds:
            for k in add.get(t, ()):
                cnt[k][i] += 1
                cnt[k][j] += 1
        self.cnt = cnt

    def jbond(self, a, b, kind):
        t = self.m.bond_type(a, b)
        return t is not None and t in _JBOND[kind]

    def apcheck(self, a, pre, ap):
        for alts in ap:
            suc = False
            for count, name in alts:
                if self._prop(a, pre, count, name):
                    suc = True
                    break
            if not suc:
                return False
        return True

    def _prop(self, a, pre, count, name):
        def cmp(v):
            return count == v or (count == -1 and v > 0)
        if name == "RG":
            return cmp(self.rg[a][0])
        m = re.fullmatch(r"RG(\d+)", name)
        if m:
            return int(m.group(1)) <= 10 and cmp(self.rg[a][int(m.group(1))])
        if name == "NR":
            return cmp(self.nr[a])
        m = re.fullmatch(r"AR([1-5])", name)
        if m:
            return cmp(self.ar[a][int(m.group(1))])
        key = name[:2]
        if key in _BOND_PROP:
            q = 1 if name[2:3] == "'" else 0
            if q == 1 and name[3:4] == "'":
                q = 2
            counter, kind = _BOND_PROP[key]
            if not cmp(self.cnt[counter][a]):
                return False
            if q == 0:
                return True
            if pre == -1:
                return False
            return self.jbond(a, pre, kind) == (q == 1)
        return False

    def name_match(self, a, name):
        sym = _PT.GetElementSymbol(self.m.z[a]) if self.m.z[a] > 0 else ""
        if len(sym) == len(name) and sym == name:
            return True
        return name in self.wild and self.m.z[a] in self.wild[name]

    def ces(self, center, chains):
        maxchain = max(len(c) for c in chains)
        found = [[] for _ in chains]          # 사슬별로 맞는 경로 (나온 순서)
        order = []                            # schain 전체 순서 (id, 원자들)
        path = [center]

        def rec():
            length = len(path) - 1
            for k, ch in enumerate(chains):
                if length != len(ch):
                    continue
                ok = True
                for j in range(1, length + 1):
                    a, node = path[j], ch[j - 1]
                    if node["connum"] and self.m.connum[a] != node["connum"]:
                        ok = False
                        break
                    if node["name"].startswith("EW") and self.ewd[a] != 1:
                        ok = False
                        break
                    if not self.name_match(a, node["name"]):
                        ok = False
                        break
                    if node["ap"] and not self.apcheck(a, path[j - 1], _parse_ap(node["ap"])):
                        ok = False
                        break
                if ok:
                    order.append((k, list(path[1:])))
                    found[k].append(len(order) - 1)
            for nb in self.m.con[path[-1]][:6]:
                if nb in path:
                    continue
                if len(path) > maxchain:
                    return
                path.append(nb)
                rec()
                path.pop()

        rec()
        if any(not f for f in found):
            return False
        sel = [-1] * len(chains)

        def check():
            for i in range(len(chains)):
                for j in range(i + 1, len(chains)):
                    t1, t2 = sel[i], sel[j]
                    if t1 == t2:
                        return False
                    a1, a2 = order[t1][1], order[t2][1]
                    mn = min(len(a1), len(a2))
                    if all(a1[k] == a2[k] for k in range(mn)):
                        return False
                    for k in range(mn):
                        same = chains[i][k]["cesname"] == chains[j][k]["cesname"]
                        if (a1[k] == a2[k]) != same:
                            return False
            return True

        def dc(level):
            for idx in range(len(order)):
                if order[idx][0] != level:
                    continue
                sel[level] = idx
                if level == len(chains) - 1:
                    if check():
                        return True
                elif dc(level + 1):
                    return True
            return False

        return dc(0)

    def match(self, a, rule):
        m = self.m
        if rule["res"] == "&":
            return True
        if rule["res"] != "*":
            return False
        for key in ("z", "connum", "nh", "new", "ap", "ces"):
            if key not in rule:
                return True
            v = rule[key]
            if key == "z" and v != m.z[a]:
                return False
            if key == "connum" and v is not None and v != m.connum[a]:
                return False
            if key == "nh" and v is not None and v != sum(1 for x in m.con[a][:6] if m.z[x] == 1):
                return False
            if key == "new" and v is not None:
                if not m.con[a] or v != sum(1 for x in m.con[m.con[a][0]][:6] if self.ewd[x] == 1):
                    return False
            if key == "ap" and v is not None and not self.apcheck(a, -1, v):
                return False
            if key == "ces" and not self.ces(a, v):
                return False
        return True

    def run(self):
        types = []
        for a in range(self.m.n):
            t = None
            for rule in self.rules["rules"]:
                if self.match(a, rule):
                    t = rule["type"]
                    break
            types.append(t)
        self._atadjust(types)
        self._cpadjust(types)
        return types

    def _atadjust(self, types):
        pair = {"cc": "cd", "ce": "cf", "cg": "ch", "pc": "pd", "pe": "pf", "nc": "nd", "ne": "nf"}
        n = self.m.n
        idx1, idx2 = [0] * n, [0] * n
        num, first = 0, True
        for a in range(n):
            if types[a] in pair:
                idx2[a] = 1
                if first:
                    idx1[a] = 1
                    first = False
                num += 1
        if num == 0:
            return
        for _ in range(num - 1):
            flag = 0
            for i, j, t in self.m.bonds:
                if idx2[i] + idx2[j] != 2:
                    continue
                if flag == 0 and idx1[i] == 0 and idx1[j] == 0:
                    idx1[i] = 1
                if idx1[i] == 0 and idx1[j] != 0:
                    flag = 1
                    if t in (1, 7):
                        idx1[i] = idx1[j]
                    if t in (2, 8, 3):
                        idx1[i] = -idx1[j]
                if idx1[j] == 0 and idx1[i] != 0:
                    flag = 1
                    if t in (1, 7):
                        idx1[j] = idx1[i]
                    if t in (2, 8, 3):
                        idx1[j] = -idx1[i]
        for a in range(n):
            if idx1[a] == -1 and types[a] in pair:
                types[a] = pair[types[a]]

    def _cpadjust(self, types):
        n = self.m.n
        idx1, idx2 = [0] * n, [0] * n
        num, first = 0, True
        for a in range(n):
            if types[a] == "cp":
                idx2[a] = 1
                if first:
                    idx1[a] = 1
                    first = False
                num += 1
        if num == 0:
            return
        for _ in range(num - 1):
            for i, j, t in self.m.bonds:
                if idx2[i] + idx2[j] != 2:
                    continue
                if idx1[i] == 0 and idx1[j] != 0:
                    idx1[i] = idx1[j] if t == 1 else -idx1[j]
                if idx1[j] == 0 and idx1[i] != 0:
                    idx1[j] = idx1[i] if t == 1 else -idx1[i]
        for a in range(n):
            if idx1[a] == -1 and types[a] == "cp":
                types[a] = "cq"

    def consistency_warnings(self, types):
        """atomtype.c 의 error(): 공액 타입 조합이 결합 차수와 맞지 않으면 경고."""
        sign = {"cc": 1, "ce": 1, "nc": 1, "ne": 1, "pc": 1, "pe": 1, "cg": 1,
                "cd": -1, "cf": -1, "nd": -1, "nf": -1, "pd": -1, "pf": -1, "ch": -1}
        out = []
        for i, j, t in self.m.bonds:
            si, sj = sign.get(types[i], 0), sign.get(types[j], 0)
            bad = (si == sj != 0 and t not in (1, 7)) or (si * sj == -1 and t in (1, 7))
            pi, pj = types[i], types[j]
            if (pi == pj == "cp" or pi == pj == "cq") and t != 1:
                bad = True
            if {pi, pj} == {"cp", "cq"} and t == 1:
                bad = True
            if bad:
                out.append(f"원자 {i + 1}({pi})–{j + 1}({pj}): antechamber 도 이 공액 타입 조합이 틀릴 수 있다고 "
                           "경고하는 경우입니다")
        return out


def kekule_graph(mol):
    """RDKit 분자(명시적 H) → (원자번호 목록, [(i, j, 케쿨레 차수)]) — 결합 순서 그대로.

    결합에 "kek" 속성(molio.stamp_kekule)이 있으면 그 차수를 쓴다.
    """
    z = [a.GetAtomicNum() for a in mol.GetAtoms()]
    if all(b.HasProp("kek") for b in mol.GetBonds()):
        return z, [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetIntProp("kek")) for b in mol.GetBonds()]
    kek = Chem.Mol(mol)
    Chem.Kekulize(kek, clearAromaticFlags=True)
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), int(b.GetBondTypeAsDouble())) for b in kek.GetBonds()]
    return z, bonds


def assign_types(z, bonds, rules=None):
    """그래프 → (GAFF2 타입 목록, antechamber 식 경고 목록)."""
    t = _Typer(_Mol(z, bonds), rules or default_def())
    types = t.run()
    return types, t.consistency_warnings(types)


# --------------------------------------------------------------------------- moltemplate gaff2.lt
# --------------------------------------------------------------------------- 대응 타입 파라미터 (parmchk2)
# PARMCHK.DAT(AmberTools 26)의 EQUA 줄 — GAFF2 소문자 타입만, 파일에 적힌 순서 그대로.
# parmchk2 는 파라미터가 없는 항을 이 "같은 타입"으로 바꿔 찾은 값을 벌점 0 으로 그대로 쓴다.
EQUA = {
    "c5": ["c3"], "c6": ["c3"], "cs": ["c"], "hw": ["ho"], "n+": ["n4"],
    "n5": ["n7", "n3"], "n6": ["n7", "n3"], "n7": ["n3"], "n8": ["n3"], "n9": ["n3"],
    "ni": ["n"], "nj": ["n"], "nk": ["nx", "n4", "ny"], "nl": ["nx", "n4", "ny"],
    "nm": ["nh"], "nn": ["nh"], "np": ["n3"], "nq": ["n3"], "ns": ["n"], "nt": ["n"],
    "nu": ["nh"], "nv": ["nh"], "nx": ["n4"], "ny": ["n4"], "nz": ["n4"],
    "o2": ["o"], "op": ["os"], "oq": ["os"], "sp": ["ss"], "sq": ["ss"],
}
# parmchk2 의 중첩 반복 순서 (바깥 → 안쪽 자리). improper 는 중심 원자를 셋째 자리에 둔다.
_EQUA_LOOP = {"bond": (0, 1), "angle": (1, 0, 2), "dihedral": (1, 2, 0, 3), "improper": (0, 1, 2, 3)}


def _equa_combos(kind, types, places=None):
    """types 의 일부를 EQUA 타입으로 바꾼 조합 (parmchk2 반복 순서, 원래 조합 제외).

    places 를 주면 그 자리만 바꾼다 (일반 이면각은 가운데 두 원자만).
    """
    order = _EQUA_LOOP[kind]
    lists = [[types[p]] + (EQUA.get(types[p], []) if places is None or p in places else [])
             for p in order]
    for combo in itertools.product(*lists):
        t = [None] * len(types)
        for p, x in zip(order, combo):
            t[p] = x
        if tuple(t) != tuple(types):
            yield tuple(t)


def _equa_decide(rules, kind, types):
    """한 항의 처리: ("ok", None) 원래 규칙 있음, ("fill", 규칙 번호) 대응 타입 규칙 복사, ("miss", None).

    parmchk2 의 chk_bond·chk_angle·chk_torsion 단계 순서를 따른다. 이면각은 원래 타입의 특정 규칙
    → 대응 타입의 특정 규칙 → 원래 타입의 일반(X) 규칙 → 가운데 원자만 바꾼 일반 규칙 순이다.
    improper 는 moltemplate(gaff_imp.py) 방식대로 바깥 세 원자 순서와 상관없이 같은 순서로 찾는다.
    """
    if kind in ("bond", "angle"):
        if rules.match(kind, types) is not None:
            return "ok", None
        for t in _equa_combos(kind, types):
            i = rules.match(kind, t)
            if i is not None:
                return "fill", i
        return "miss", None
    if rules.match(kind, types, wild=False) is not None:
        return "ok", None
    for t in _equa_combos(kind, types):
        i = rules.match(kind, t, wild=False)
        if i is not None:
            return "fill", i
    if rules.match(kind, types, exact=False) is not None:
        return "ok", None
    for t in _equa_combos(kind, types, places=(1, 2) if kind == "dihedral" else None):
        i = rules.match(kind, t, exact=False)
        if i is not None:
            return "fill", i
    return "miss", None


def _canonical(kind, t, center=2):
    """parmchk2 frcmod 의 이름 순서 (결합·각은 양 끝 알파벳 순, 이면각은 가운데 쌍 기준, improper 는
    중심 원자를 빼고 알파벳 순)."""
    t = tuple(t)
    if kind == "improper":
        outer = sorted(x for i, x in enumerate(t) if i != center)
        outer.insert(center, t[center])
        return tuple(outer)
    if kind == "dihedral":
        return t if (t[1], t[2], t[0], t[3]) <= (t[2], t[1], t[3], t[0]) else t[::-1]
    return min(t, t[::-1])


def equa_fill(mol, labels, rules):
    """gaff2.lt 에 규칙이 없는 항 중 대응 타입(EQUA)으로 값을 그대로 가져올 수 있는 것.

    반환: [(종류, 타입 튜플, 가져온 계수 타입 이름)] — 처음 나온 순서, 같은 타입 튜플은 한 번.
    """
    from .ltrules import topology
    terms = topology(mol, rules.improper_center)
    out, seen = [], set()
    for kind in ("bond", "angle", "dihedral", "improper"):
        if kind == "improper" and rules.improper_center is None:
            continue
        for t in terms[kind]:
            types = tuple(labels[i] for i in t)
            key = _canonical(kind, types, rules.improper_center)
            if (kind, key) in seen:
                continue
            seen.add((kind, key))
            what, i = _equa_decide(rules, kind, types)
            if what == "fill":
                src = rules.rules[kind][i][0]
                if src in rules.coeffs[kind]:
                    out.append((kind, key, src))
    return out


def equa_block(rules, adds, obj="GAFF2"):
    """equa_fill 결과를 역장 객체에 덧붙이는 moltemplate 문장 (gaff2.lt 파일은 그대로 둠)."""
    if not adds:
        return ""
    names = {}
    lines = ["# GAFF2 2.11 에 파라미터가 없는 항 — parmchk2 처럼 대응 타입(EQUA)의 gaff2.lt 값을 그대로 복사",
             f"{obj} {{", '  write_once("In Settings") {']
    for kind, key, src in adds:
        name = "-".join(key)
        if name in rules.coeffs[kind]:            # 역장 파일의 이름과 겹치지 않게
            name += "_equa"
        names[(kind, key)] = name
        lines.append(f"    {kind}_coeff @{kind}:{name} {rules.coeffs[kind][src]}   # = {src}")
    lines.append("  }")
    for kind in ("bond", "angle", "dihedral", "improper"):
        rows = [(key, names[(k, key)]) for k, key, _ in adds if k == kind]
        if rows:
            lines.append(f'  write_once("{rules.sections[kind]}") {{')
            lines += [f"    @{kind}:{name} " + " ".join(f"@atom:{x}" for x in key) for key, name in rows]
            lines.append("  }")
    lines.append("}")
    return "\n".join(lines) + "\n"


class Gaff2Typer:
    """moltemplate gaff2.lt 용 타이퍼.

    moltemplate 의 gaff2.lt 는 GAFF2 2.11 판(gaff2.dat, 2016)으로 만든 것이라 GAFF2 2.2 에서 생긴
    c5·c6(5·6원 고리 sp3 C)이 없다. 2.11 판 규칙에서는 이 원자들이 c3 이므로, 역장 파일에 c5·c6 이
    없으면 c3 로 둔다(규칙 파일에서 c5·c6 두 줄을 뺀 것과 같은 결과).
    2.11 판에는 ns, nt, nx 처럼 질량·LJ 만 있고 결합 파라미터가 없는 타입이 있어서, AmberTools 에서
    parmchk2 가 하는 것처럼 대응 타입(EQUA)의 값을 그대로 복사해 .lt 에 함께 쓴다(ff_extra).
    """
    NAME = "gaff2"
    LABEL = "GAFF2"
    LT_FILE = "gaff2.lt"
    LT_OBJECT = "GAFF2"
    TYPE_CHARGES = False

    def __init__(self, ff_file=None):
        from .ltrules import LtRules
        self.rules = LtRules(ff_file or moltemplate_ff_path("gaff2.lt"))
        if not self.rules.masses:
            raise TypingError(f"{self.rules.path} 에서 Data Masses 타입을 찾지 못했습니다")

    def assign(self, mol):
        from .ltrules import KIND_KO, missing_terms
        z, bonds = kekule_graph(mol)
        types, warnings = assign_types(z, bonds)
        old = [i for i, t in enumerate(types) if t in ("c5", "c6") and t not in self.rules.masses]
        for i in old:
            types[i] = "c3"
        if old:
            warnings.append(f"원자 {', '.join(str(i + 1) for i in old)}: 이 gaff2.lt 는 GAFF2 2.11 판이라 "
                            "5·6원 고리 sp3 C 를 c3 로 둡니다 (AmberTools GAFF2 2.2 의 c5·c6)")
        missing = sorted({t for t in types if t not in self.rules.masses})
        if missing:
            where = [f"{i + 1}({t})" for i, t in enumerate(types) if t in missing]
            raise TypingError(f"gaff2.lt 에 없는 GAFF2 타입이 나왔습니다: {', '.join(missing)} — "
                              f"원자 {', '.join(where)}")
        adds = equa_fill(mol, types, self.rules)
        if adds:
            warnings.append(f"gaff2.lt 에 없는 항 {len(adds)}개는 대응 타입 값을 복사해 .lt 에 함께 씁니다 "
                            f"(parmchk2 EQUA): " + ", ".join(f"{KIND_KO[k]} {'-'.join(key)}←{src}" for k, key, src in adds[:6])
                            + (" …" if len(adds) > 6 else ""))
        covered = {}
        for kind, key, _ in adds:
            covered.setdefault(kind, set()).add(key)
        warnings += missing_terms(mol, types, self.rules, covered)
        return types, warnings

    def ff_extra(self, mol, labels):
        """분자 정의 앞에 쓸 역장 보충 문장 (없으면 빈 문자열)."""
        return equa_block(self.rules, equa_fill(mol, labels, self.rules), self.LT_OBJECT)


def rules_json():
    """웹 이식용 규칙 JSON (화학 환경 사슬과 원자 성질까지 미리 해석)."""
    return json.dumps(default_def(), separators=(",", ":"), ensure_ascii=False)
