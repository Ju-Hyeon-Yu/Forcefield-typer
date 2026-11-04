# foyer(MIT)의 SMARTS 원자 타입 판정 — 문법·고리 찾기·규칙 반복·덮어쓰기를 그대로 옮긴 엔진
"""foyer 1.2.0 (smarts.py, smarts_graph.py, atomtyper.py)과 같은 결과를 내도록 옮겼다.

- 문법: foyer SMARTS (원소 기호, *, #번호, %타입, X연결수, r고리크기, R고리수, ! & , ;, 가지, 고리 닫기).
  결합 기호(- = # : ~)와 $( ) 는 OPLS-AA 규칙에 없어 지원하지 않는다.
- 패턴의 첫 원자가 타입을 받는 원자다. 짝 맞추기는 networkx 의 subgraph_isomorphisms_iter 와 같은
  "유도 부분그래프" 기준이다 (패턴에서 결합이 없는 두 원자는 분자에서도 결합이 없어야 함).
- 규칙은 파일 순서대로, 원자는 번호 순으로 맞춰 보며, 맞으면 바로 그 원자의 후보(whitelist)에 넣고
  덮어쓰는 타입들을 제외 목록(blacklist)에 넣는다. 새로 맞는 것이 없을 때까지 최대 10번 되풀이한다.
- 고리(r, R)는 foyer 의 _find_chordless_cycles(최대 8원)를 그대로 따른다.
웹 이식(web/src/foyer.js)은 여기서 해석한 패턴(JSON)을 그대로 쓴다.
"""
import itertools

from rdkit import Chem

_PT = Chem.GetPeriodicTable()
# foyer smarts.py 의 SYMBOL 정규식과 같은 원소 기호 목록 (긴 것부터 맞춤)
_SYMBOL_RE = (r"C[laroudsemf]?|Os?|N[eaibdpos]?|S[icernbmg]?|P[drmtboau]?|H[eofgas]?|A[lrsgutcm]|B[eraik]?|"
              r"Dy|E[urs]|F[erm]?|G[aed]|I[nr]?|Kr?|L[iaur]|M[gnodt]|R[buhenaf]|T[icebmalh]|U|V|W|Xe|Yb?|Z[nr]")


class SmartsError(ValueError):
    pass


# --------------------------------------------------------------------------- 문법 해석
class _Parser:
    """foyer SMARTS → {"atoms": [식], "edges": [[i, j]], "element": 뿌리 원소 또는 None, "ring": bool}.

    식: ["and", a, b] (& 와 ;), ["or", a, b], ["not", x], ["sym", "C" | "*"], ["num", 6], ["X", n],
        ["r", n], ["R", n], ["label", "opls_145"]
    원자 번호는 문자열에 나온 순서(foyer 의 iter_subtrees_topdown 순서)와 같다.
    """

    def __init__(self, text):
        import re
        self.re_symbol = re.compile(_SYMBOL_RE)
        self.re_num = re.compile(r"\d+")
        self.re_label = re.compile(r"%[A-Za-z_0-9]+")
        self.s = text
        self.i = 0
        self.atoms, self.edges, self.labels = [], [], {}

    def peek(self):
        return self.s[self.i] if self.i < len(self.s) else ""

    def fail(self, msg):
        raise SmartsError(f"SMARTS '{self.s}' {self.i}번째 글자: {msg}")

    def parse(self):
        self.string(None)
        if self.i != len(self.s):
            self.fail("해석하지 못한 글자가 남았습니다")
        for d, ends in self.labels.items():
            if len(ends) != 2:
                self.fail(f"고리 번호 {d} 가 짝이 맞지 않습니다")
            self.edges.append([ends[0], ends[1]])
        ring = any(self._has(a, ("r", "R")) for a in self.atoms)
        return {"atoms": self.atoms, "edges": self.edges, "element": self._root_element(), "ring": ring}

    # _string: _chain _nonlastbranch* _lastbranch?   /  _chain: atom (bond? atom)*   /  branch: bond? _string
    def string(self, trunk):
        self._no_bond()
        cur = self.atom()
        if trunk is not None:
            self.edges.append([trunk, cur])
        while self.peek() and self.peek() not in "()":
            self._no_bond()
            nxt = self.atom()
            self.edges.append([cur, nxt])
            cur = nxt
        while self.peek() == "(":                 # 괄호 가지는 모두 사슬의 마지막 원자에 붙음
            self.i += 1
            self.string(cur)
            if self.peek() != ")":
                self.fail("')' 가 필요합니다")
            self.i += 1
        if self.peek() and self.peek() != ")":    # 괄호 없이 이어지는 마지막 가지
            self.string(cur)

    def _no_bond(self):
        if self.peek() and self.peek() in "-=#:~":
            self.fail("결합 기호는 지원하지 않습니다")

    def atom(self):
        idx = len(self.atoms)
        if self.peek() == "[":
            self.i += 1
            expr = self.weak_and()
            if self.peek() != "]":
                self.fail("']' 가 필요합니다")
            self.i += 1
        else:
            expr = self.symbol()
            if expr is None:
                self.fail("원자가 필요합니다")
        self.atoms.append(expr)
        m = self.re_num.match(self.s, self.i)
        if m:
            for d in m.group(0):
                self.labels.setdefault(d, []).append(idx)
            self.i = m.end()
        return idx

    def symbol(self):
        if self.peek() == "*":
            self.i += 1
            return ["sym", "*"]
        m = self.re_symbol.match(self.s, self.i)
        if not m:
            return None
        self.i = m.end()
        return ["sym", m.group(0)]

    def weak_and(self):
        node = self.or_()
        while self.peek() == ";":
            self.i += 1
            node = ["and", node, self.or_()]
        return node

    def or_(self):
        node = self.and_()
        while self.peek() == ",":
            self.i += 1
            node = ["or", node, self.and_()]
        return node

    def and_(self):
        node = self.unit()
        while self.peek() == "&":
            self.i += 1
            node = ["and", node, self.unit()]
        return node

    def unit(self):
        if self.peek() == "!":
            self.i += 1
            return ["not", self.atom_id()]
        return self.atom_id()

    def atom_id(self):
        c = self.peek()
        if c == "#":
            self.i += 1
            return ["num", self.number()]
        if c == "$":
            self.fail("$( ) 는 지원하지 않습니다")
        if c == "%":
            m = self.re_label.match(self.s, self.i)
            self.i = m.end()
            return ["label", m.group(0)[1:]]
        sym = self.symbol()
        if sym is not None:
            return sym
        if c in ("X", "r", "R"):
            self.i += 1
            return [c, self.number()]
        self.fail("알 수 없는 원자 조건")

    def number(self):
        m = self.re_num.match(self.s, self.i)
        if not m:
            self.fail("숫자가 필요합니다")
        self.i = m.end()
        return int(m.group(0))

    @staticmethod
    def _has(expr, kinds):
        if expr[0] in kinds:
            return True
        return any(isinstance(x, list) and _Parser._has(x, kinds) for x in expr[1:])

    def _root_element(self):
        """atomtyper.find_atomtypes 의 규칙 거르기: 뿌리 원자 기호가 하나뿐이고 부정이 없으면 그 기호."""
        syms = []

        def walk(e):
            if e[0] == "sym":
                syms.append(e[1])
            for x in e[1:]:
                if isinstance(x, list):
                    walk(x)
        walk(self.atoms[0])
        if len(syms) == 1 and not self._has(self.atoms[0], ("not",)):
            return syms[0]
        return None


def parse(text):
    return _Parser(text).parse()


# --------------------------------------------------------------------------- 분자 그래프
class Graph:
    """원자번호 목록 + 결합 목록. 이웃은 결합이 나온 순서(networkx 인접 목록 순서)."""

    def __init__(self, z, bonds):
        self.n = len(z)
        self.z = list(z)
        self.nbrs = [[] for _ in z]
        self.adj = [set() for _ in z]
        for i, j in bonds:
            if j in self.adj[i]:
                continue
            self.nbrs[i].append(j)
            self.nbrs[j].append(i)
            self.adj[i].add(j)
            self.adj[j].add(i)
        self.cycles = None

    def ring_info(self):
        """foyer _find_chordless_cycles(max_cycle_size=8): 원자별 고리 튜플 집합."""
        if self.cycles is None:
            self.cycles = [set(tuple(c) for c in cs) for cs in _chordless_cycles(self, 8)]
        return self.cycles


def _chordless_cycles(g, max_cycle_size):
    cycles = [[] for _ in range(g.n)]
    for i in range(g.n):
        for a, b in itertools.combinations(g.nbrs[i], 2):
            connected = False
            last = a
            possible = [[a, i, b]]
            if b in g.adj[last]:
                cycles[i].append([a, i, b])
                connected = True
            while not connected:
                new = []
                for ring in possible:
                    for nb in g.nbrs[ring[-1]]:
                        if nb != ring[-2]:
                            new.append(ring + [nb])
                possible = new
                for ring in possible:
                    if last in g.adj[ring[-1]]:
                        if not any(x in g.adj[ring[-1]] for x in ring[1:-2]):
                            cycles[i].append(ring)
                            connected = True
                if not possible or len(possible[0]) == max_cycle_size:
                    break
    return cycles


# --------------------------------------------------------------------------- 짝 맞추기
def _atom_ok(expr, g, a, white):
    k = expr[0]
    if k == "and":
        return _atom_ok(expr[1], g, a, white) and _atom_ok(expr[2], g, a, white)
    if k == "or":
        return _atom_ok(expr[1], g, a, white) or _atom_ok(expr[2], g, a, white)
    if k == "not":
        return not _atom_ok(expr[1], g, a, white)
    if k == "sym":
        return expr[1] == "*" or g.z[a] == _PT.GetAtomicNumber(expr[1])
    if k == "num":
        return g.z[a] == expr[1]
    if k == "label":
        return expr[1] in white[a]
    if k == "X":
        return len(g.nbrs[a]) == expr[1]
    if k == "r":
        return any(len(c) == expr[1] for c in g.ring_info()[a])
    if k == "R":
        return len(g.ring_info()[a]) == expr[1]
    raise SmartsError(f"알 수 없는 식 {k}")


def _plan(pat):
    """패턴 원자를 뿌리부터 이어지는 순서로: [(원자, 이미 놓인 이웃 하나)]."""
    n = len(pat["atoms"])
    padj = [set() for _ in range(n)]
    for i, j in pat["edges"]:
        padj[i].add(j)
        padj[j].add(i)
    order, seen = [(0, None)], {0}
    k = 0
    while k < len(order):
        p = order[k][0]
        for q in sorted(padj[p]):
            if q not in seen:
                seen.add(q)
                order.append((q, p))
        k += 1
    if len(order) != n:
        raise SmartsError("이어지지 않은 패턴은 지원하지 않습니다")
    return order, padj


def matches_at(pat, g, root, white):
    """root 원자를 패턴 첫 원자로 하는 유도 부분그래프 짝이 있는지."""
    if "_plan" not in pat:
        pat["_plan"] = _plan(pat)
    order, padj = pat["_plan"]
    atoms = pat["atoms"]
    if not _atom_ok(atoms[0], g, root, white):
        return False
    n = len(order)
    image = {0: root}
    used = {root}

    def rec(k):
        if k == n:
            return True
        p, anchor = order[k]
        for h in g.nbrs[image[anchor]]:
            if h in used or not _atom_ok(atoms[p], g, h, white):
                continue
            if any((q in padj[p]) != (hq in g.adj[h]) for q, hq in image.items()):
                continue
            image[p] = h
            used.add(h)
            if rec(k + 1):
                return True
            del image[p]
            used.discard(h)
        return False

    return rec(1)


# --------------------------------------------------------------------------- 규칙 반복과 결정
def load_rules(rules):
    """[{name, def, overrides}] → 해석한 규칙 목록 (파일 순서)."""
    out = []
    for r in rules:
        pat = parse(r["def"])
        out.append({"name": r["name"], "overrides": set(r.get("overrides", [])), "pat": pat})
    return out


def find_atomtypes(g, rules, max_iter=10):
    """foyer find_atomtypes 와 같은 순서로 판정 → (원자별 후보 집합 whitelist - blacklist 목록)."""
    present = {_PT.GetElementSymbol(z) for z in g.z if z > 0}
    rules = [r for r in rules if r["pat"]["element"] is None or r["pat"]["element"] in present]
    white = [set() for _ in range(g.n)]
    black = [set() for _ in range(g.n)]
    for _ in range(max_iter):
        found = False
        for r in rules:
            for a in range(g.n):
                if matches_at(r["pat"], g, a, white) and r["name"] not in white[a]:
                    white[a].add(r["name"])
                    black[a] |= r["overrides"]
                    found = True
        if not found:
            break
    return [sorted(w - b) for w, b in zip(white, black)]
