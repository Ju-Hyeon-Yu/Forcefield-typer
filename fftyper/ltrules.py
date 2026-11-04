# moltemplate 역장 .lt 의 타입 목록·"By Type" 규칙을 읽어 결합·각·이면각 파라미터가 있는지 확인
import fnmatch
import re


class LtRules:
    """Data Masses 타입과 Data Bonds/Angles/Dihedrals By Type 규칙.

    moltemplate 는 맞는 규칙이 없는 각·이면각을 경고 없이 빼므로, 라벨을 붙인 뒤 미리 확인한다.
    규칙의 @atom 패턴은 와일드카드(*, ?)만 지원한다 (gaff2.lt, oplsaa2024.lt 형식).
    """

    SECTIONS = {"Data Bonds By Type": 2, "Data Angles By Type": 3, "Data Dihedrals By Type": 4}

    def __init__(self, path, atom_alias=None):
        self.path = path
        self.masses = {}
        self.rules = {2: [], 3: [], 4: []}
        section = None
        with open(path) as f:
            for line in f:
                s = line.split("#")[0].strip()
                m = re.match(r'write_once\("([^"]+)"\)', s)
                if m:
                    section = m.group(1)
                    continue
                if s.startswith("}"):
                    section = None
                    continue
                if not s:
                    continue
                t = s.split()
                if section == "Data Masses" and t[0].startswith("@atom:"):
                    self.masses[t[0][6:]] = float(t[1])
                elif section in self.SECTIONS and t[0].startswith("@"):
                    n = self.SECTIONS[section]
                    pats = [x.split(":", 1)[1] for x in t[1:1 + n]]
                    self.rules[n].append(pats)
        self._exact = {n: set() for n in self.rules}
        self._wild = {n: [] for n in self.rules}
        for n, rules in self.rules.items():
            for pats in rules:
                if any(ch in p for p in pats for ch in "*?["):
                    self._wild[n].append(pats)
                else:
                    self._exact[n].add(tuple(pats))
        self._cache = {}

    def has(self, names):
        """타입 이름 튜플(정방향)에 맞는 규칙이 하나라도 있는지 (역방향도 확인)."""
        names = tuple(names)
        if names in self._cache:
            return self._cache[names]
        n = len(names)
        rev = names[::-1]
        ok = names in self._exact[n] or rev in self._exact[n]
        if not ok:
            for pats in self._wild[n]:
                if all(fnmatch.fnmatchcase(a, p) for a, p in zip(names, pats)) or \
                        all(fnmatch.fnmatchcase(a, p) for a, p in zip(rev, pats)):
                    ok = True
                    break
        self._cache[names] = ok
        return ok


def topology(mol):
    """RDKit 분자 → (결합, 결합각, 이면각) 원자 번호 튜플 목록."""
    nb = [[n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    angles = [(i, j, k) for j in range(len(nb)) for x, i in enumerate(nb[j]) for k in nb[j][x + 1:]]
    dihedrals = []
    for j, k in bonds:
        for i in nb[j]:
            if i == k:
                continue
            for l in nb[k]:
                if l == j or l == i:
                    continue
                dihedrals.append((i, j, k, l))
    return bonds, angles, dihedrals


def missing_terms(mol, labels, rules, kind_names=("결합", "결합각", "이면각"), limit=6):
    """파라미터 규칙이 없는 결합·각·이면각을 [경고 문장] 으로 돌려준다."""
    out = []
    for terms, name in zip(topology(mol), kind_names):
        miss = {}
        for t in terms:
            key = tuple(labels[i] for i in t)
            if not rules.has(key):
                k = min(key, key[::-1])
                miss.setdefault(k, []).append(t)
        for key, where in sorted(miss.items()):
            atoms = ", ".join("-".join(str(i + 1) for i in t) for t in where[:3])
            more = f" 외 {len(where) - 3}개" if len(where) > 3 else ""
            out.append(f"{name} {'-'.join(key)} (원자 {atoms}{more}): 역장 파일에 파라미터가 없어 moltemplate 가 "
                       "이 항을 빼고 만듭니다")
            if len(out) >= limit * 3:
                return out
    return out
