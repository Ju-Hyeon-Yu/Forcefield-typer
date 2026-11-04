# foyer OPLS-AA SMARTS 규칙으로 판정해 moltemplate oplsaa2024.lt 타입 번호·타입 전하를 붙이는 타이퍼
"""OPLS-AA 타이퍼.

1. foyer(MIT)의 OPLS-AA 규칙 228개(data/foyer_oplsaa.json)로 foyer 와 똑같이 판정한다(foyer_smarts).
2. foyer 타입(GROMACS 번호, opls_135 …)을 moltemplate oplsaa2024.lt 의 번호(BOSS 번호)로 바꾼다.
   두 번호 체계는 대부분 같지만 뜻이 다른 번호가 있어(NO_2024), 그런 타입이나 2024 파일에 없는 foyer
   추가 타입(opls_151a 처럼 글자가 붙은 것)이 나오면 오류로 멈춘다.
3. oplsaa2024.lt 가 OPLS/2020 으로 전하를 바꾸면서 짝 타입을 새로 둔 곳(OPLS2020_PAIRS)은 그 짝으로 바꾼다.
   foyer 규칙은 고전 OPLS-AA 기준이라 이 짝을 모르고, 그대로 두면 전하 합이 0 이 되지 않는다.
4. 전하는 oplsaa2024.lt 의 타입 전하(In Charges)를 쓴다. moltemplate 도 실행할 때 이 값으로 덮어쓴다.
"""
import json
import os

from rdkit import Chem

from . import TypingError, moltemplate_ff_path
from .foyer_smarts import Graph, find_atomtypes, load_rules

RULES_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "foyer_oplsaa.json")

# foyer 번호와 oplsaa2024.lt 번호가 같아도 가리키는 원자가 다른 것 (대응 없음)
NO_2024 = {
    "opls_732": "2024 의 732 는 아이오도벤젠 I 인데 foyer 는 2-아이오도프로판 I 에 씀",
    "opls_927": "2024 의 927 은 알카인에 붙은 CH3 탄소인데 foyer 는 알카인 sp 탄소에 씀",
    "opls_930": "2024 의 930 은 R3C–C≡ 탄소인데 foyer 는 알카인 옆 탄소의 H 에 씀",
}

# (바꿀 타입, 이웃에 이 타입이 있으면, 바꿀 2024 타입) — 둘 다 oplsaa2024.lt 설명에 OPLS/2020 짝으로 적혀 있음
OPLS2020_PAIRS = [
    ("154", ("159",), "160"),          # 3차 알코올: C 159 "C, tert. alcohols OPLS/2020" 의 O → 160 "mono tert. alcohols"
    ("180", ("189",), "186"),          # OCH2O 아세탈: C 189(0.06) 의 O → 186 "acetal/hemiacetal ether O OPLS/2020"
]

_RULES = None


def default_rules():
    global _RULES
    if _RULES is None:
        with open(RULES_PATH) as f:
            _RULES = json.load(f)
    return _RULES


def to_2024(name, types2024):
    """foyer 타입 이름 → oplsaa2024.lt 타입 번호(문자열) 또는 None."""
    num = name[5:] if name.startswith("opls_") else ""
    if name in NO_2024 or not num.isdigit() or num not in types2024:
        return None
    return num


def graph_of(mol):
    """RDKit 분자(명시적 H) → foyer 판정용 그래프 (원자·결합 순서 그대로)."""
    z = [a.GetAtomicNum() for a in mol.GetAtoms()]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    return Graph(z, bonds)


class OplsTyper:
    NAME = "oplsaa"
    LABEL = "OPLS-AA"
    LT_FILE = "oplsaa2024.lt"
    LT_OBJECT = "OPLSAA"
    TYPE_CHARGES = True

    def __init__(self, ff_file=None):
        from .ltrules import LtRules
        self.rules = LtRules(ff_file or moltemplate_ff_path("oplsaa2024.lt"))
        if not self.rules.charges:
            raise TypingError(f"{self.rules.path} 에서 In Charges 타입 전하를 찾지 못했습니다")
        data = default_rules()
        self.foyer = load_rules(data["rules"])
        self.desc = {r["name"]: r.get("desc", "") for r in data["rules"]}

    def foyer_types(self, mol):
        """foyer 판정 결과 (원자별 타입 이름). 타입이 없거나 여럿이면 TypingError."""
        cands = find_atomtypes(graph_of(mol), self.foyer)
        none = [i for i, c in enumerate(cands) if not c]
        many = [i for i, c in enumerate(cands) if len(c) > 1]
        if none or many:
            sym = [a.GetSymbol() for a in mol.GetAtoms()]
            msg = []
            if none:
                msg.append("OPLS-AA 규칙(foyer)에 맞는 타입이 없는 원자: " +
                           ", ".join(f"{i + 1}({sym[i]})" for i in none))
            if many:
                msg.append("타입이 여럿 맞는 원자: " +
                           ", ".join(f"{i + 1}({sym[i]}: {', '.join(cands[i])})" for i in many))
            raise TypingError("; ".join(msg))
        return [c[0] for c in cands]

    def assign(self, mol):
        from .ltrules import missing_terms
        names = self.foyer_types(mol)
        labels, bad = [], {}
        for i, nm in enumerate(names):
            num = to_2024(nm, self.rules.charges)
            if num is None or num not in self.rules.masses:
                bad.setdefault(nm, []).append(i + 1)
            labels.append(num)
        if not bad:
            nb = [[n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()]
            old = list(labels)
            for i, t in enumerate(old):
                for a, near, b in OPLS2020_PAIRS:
                    if t == a and any(old[j] in near for j in nb[i]) and b in self.rules.masses:
                        labels[i] = b
        if bad:
            raise TypingError("oplsaa2024.lt 에 대응하는 타입이 없습니다 — " + "; ".join(
                f"원자 {', '.join(map(str, at))}: foyer {nm} ({self.desc.get(nm, '')}) — "
                f"{NO_2024.get(nm) or 'oplsaa2024.lt 에 같은 번호가 없음'}" for nm, at in bad.items()))
        return labels, missing_terms(mol, labels, self.rules)

    def charges(self, labels):
        return [self.rules.charges[t] for t in labels]
