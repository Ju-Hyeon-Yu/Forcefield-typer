# moltemplate 역장 .lt 의 타입 목록·"By Type" 규칙·계수 줄을 읽어 결합·각·이면각 파라미터가 있는지 확인
import fnmatch
import itertools
import re

KINDS = ("bond", "angle", "dihedral", "improper")
NATOMS = {"bond": 2, "angle": 3, "dihedral": 4, "improper": 4}
KIND_KO = {"bond": "결합", "angle": "결합각", "dihedral": "이면각", "improper": "improper"}
# improper 대칭 모듈 → 중심 원자 자리 (moltemplate nbody_alt_symmetry). 나머지 셋은 순서 무관.
IMPROPER_CENTER = {"gaff_imp.py": 2}


class LtRules:
    """Data Masses 타입, Data Bonds/Angles/Dihedrals/Impropers By Type 규칙, In Settings 계수 줄.

    moltemplate 는 맞는 규칙이 없는 각·이면각을 경고 없이 빼므로, 라벨을 붙인 뒤 미리 확인한다.
    규칙이 여럿 맞으면 moltemplate 는 파일에서 나중에 나온 규칙을 쓴다(remove_duplicates_nbody.py).
    규칙의 @atom 패턴은 와일드카드(*, ?)만 지원한다 (gaff2.lt, oplsaa2024.lt 형식).
    """

    def __init__(self, path):
        self.path = path
        self.masses = {}
        self.rules = {k: [] for k in KINDS}       # [(계수 타입 이름, [원자 패턴…])] 파일 순서
        self.coeffs = {k: {} for k in KINDS}      # 계수 타입 이름 → "harmonic 232.5 1.5380" 같은 나머지
        self.sections = {}                        # 종류 → By Type 절 이름 (대칭 모듈 이름 포함)
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
                if not s or section is None:
                    continue
                t = s.split()
                m = re.match(r"Data (Bond|Angle|Dihedral|Improper)s By Type", section)
                if section == "Data Masses" and t[0].startswith("@atom:"):
                    self.masses[t[0][6:]] = float(t[1])
                elif m and t[0].startswith("@"):
                    kind = m.group(1).lower()
                    self.sections[kind] = section
                    pats = [x.split(":", 1)[1] for x in t[1:1 + NATOMS[kind]]]
                    self.rules[kind].append((t[0].split(":", 1)[1], pats))
                elif section == "In Settings" and re.match(r"(bond|angle|dihedral|improper)_coeff$", t[0]):
                    kind = t[0][:-6]
                    self.coeffs[kind][t[1].split(":", 1)[1]] = " ".join(t[2:])
        self._exact = {k: {} for k in KINDS}      # 원자 타입 튜플 → 마지막 규칙 번호
        self._wild = {k: [] for k in KINDS}       # [(규칙 번호, 패턴)] 번호 내림차순
        for kind, rules in self.rules.items():
            for i, (_, pats) in enumerate(rules):
                if any(ch in p for p in pats for ch in "*?["):
                    self._wild[kind].append((i, pats))
                else:
                    self._exact[kind][tuple(pats)] = i
            self._wild[kind].reverse()
        sym = re.search(r"\((\S+)\)", self.sections.get("improper", ""))
        self.improper_center = IMPROPER_CENTER.get(sym.group(1)) if sym else None

    def orders(self, kind, names):
        """moltemplate 가 같은 항으로 보는 원자 순서들 (결합·각·이면각은 정·역방향)."""
        names = tuple(names)
        if kind != "improper":
            return [names, names[::-1]]
        c = self.improper_center
        outer = [x for i, x in enumerate(names) if i != c]
        out = []
        for p in itertools.permutations(outer):
            p = list(p)
            p.insert(c, names[c])
            out.append(tuple(p))
        return out

    def match(self, kind, names, exact=True, wild=True):
        """이 타입 조합에 moltemplate 가 줄 규칙 번호(맞는 것 중 마지막), 없으면 None.

        exact/wild 로 와일드카드 없는 규칙·있는 규칙만 볼 수 있다 (parmchk2 의 단계 구분용).
        """
        best = -1
        for t in self.orders(kind, names):
            if exact:
                best = max(best, self._exact[kind].get(t, -1))
            if wild:
                for i, pats in self._wild[kind]:
                    if i <= best:
                        break
                    if all(fnmatch.fnmatchcase(a, p) for a, p in zip(t, pats)):
                        best = i
                        break
        return best if best >= 0 else None


def topology(mol, improper_center=None):
    """RDKit 분자 → {종류: 원자 번호 튜플 목록}. improper 는 중심 원자를 improper_center 자리에 둔다."""
    nb = [sorted(n.GetIdx() for n in a.GetNeighbors()) for a in mol.GetAtoms()]
    bonds = sorted(tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds())
    angles = [(i, j, k) for j in range(len(nb)) for x, i in enumerate(nb[j]) for k in nb[j][x + 1:]]
    dihedrals = [(i, j, k, l) for j, k in bonds for i in nb[j] if i != k
                 for l in nb[k] if l != j and l != i]
    impropers = []
    if improper_center is not None:
        for c in range(len(nb)):
            for trio in itertools.combinations(nb[c], 3):
                t = list(trio)
                t.insert(improper_center, c)
                impropers.append(tuple(t))
    return {"bond": bonds, "angle": angles, "dihedral": dihedrals, "improper": impropers}


def missing_terms(mol, labels, rules, covered=None, limit=6):
    """파라미터 규칙이 없는 결합·각·이면각을 [경고 문장] 으로 돌려준다.

    covered: {종류: {타입 튜플}} — .lt 에 따로 채워 넣는 항 (정·역방향 어느 쪽이든).
    improper 는 원래 일부 원자에만 있으므로 확인하지 않는다.
    """
    covered = covered or {}
    out = []
    terms = topology(mol)
    for kind in ("bond", "angle", "dihedral"):
        miss = {}
        extra = covered.get(kind, set())
        for t in terms[kind]:
            key = tuple(labels[i] for i in t)
            if key in extra or key[::-1] in extra or rules.match(kind, key) is not None:
                continue
            miss.setdefault(min(key, key[::-1]), []).append(t)
        for key, where in sorted(miss.items()):
            atoms = ", ".join("-".join(str(i + 1) for i in t) for t in where[:3])
            more = f" 외 {len(where) - 3}개" if len(where) > 3 else ""
            what = ("moltemplate 가 오류로 멈춥니다 (No bond types defined)" if kind == "bond"
                    else "moltemplate 가 이 항을 빼고 만듭니다")
            out.append(f"{KIND_KO[kind]} {'-'.join(key)} (원자 {atoms}{more}): 역장 파일에 파라미터가 없어 {what}")
            if len(out) >= limit * 3:
                return out
    return out
