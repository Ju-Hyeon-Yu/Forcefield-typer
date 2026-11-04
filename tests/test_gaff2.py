# GAFF2 타이퍼를 antechamber(AmberTools 26) 결과와 원자 단위로 대조하는 테스트
import json
import os
import sys

import pytest
from rdkit import Chem

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from fftyper import TypingError  # noqa: E402
from fftyper import gaff2  # noqa: E402

REF = json.load(open(os.path.join(ROOT, "tests", "data", "gaff2_reference.json")))


def mol(smiles):
    return Chem.AddHs(Chem.MolFromSmiles(smiles))


def test_same_as_antechamber_on_reference_set():
    """tests/antechamber_ref.py 로 만든 정답지(antechamber -at gaff2 -j 5)와 전부 같아야 함."""
    diffs, n = [], 0
    for smi, ref in REF.items():
        z, bonds = gaff2.kekule_graph(mol(smi))
        types, _ = gaff2.assign_types(z, bonds)
        n += len(types)
        if types != ref["types"]:
            diffs.append((smi, [(i + 1, a, b) for i, (a, b) in enumerate(zip(types, ref["types"])) if a != b]))
    assert len(REF) > 350 and n > 5000
    assert not diffs, diffs[:5]


@pytest.mark.parametrize("smiles, expect", [
    ("C=CC=CC=C", ["c2", "ce", "ce", "cf", "cf", "c2"]),              # 공액 사슬: 이중결합 양쪽은 다른 짝
    ("c1ccc(cc1)-c1ccccc1", ["ca", "ca", "ca", "cp", "ca", "ca", "cp"]),
    ("CC(=O)[O-]", ["c3", "c", "o", "o"]),
    ("O=[N+]([O-])c1ccccc1", ["o", "no", "o", "ca"]),
    ("Cn1cc[n+](C)c1", ["c3", "na", "cc", "cd", "na", "c3", "cc"]),
    ("CC[NH+](CC)CC", ["c3", "c3", "nx", "c3", "c3", "c3", "c3"]),
    ("[O-]S(=O)(=O)C(F)(F)F", ["o", "s6", "o", "o", "c3", "f", "f", "f"]),
    ("Ic1ccccc1", ["i", "ca"]),
])
def test_known_types(smiles, expect):
    z, bonds = gaff2.kekule_graph(mol(smiles))
    types, _ = gaff2.assign_types(z, bonds)
    assert types[:len(expect)] == expect


def test_ring_sp3_carbon_is_c3_for_moltemplate_gaff211():
    typer = gaff2.Gaff2Typer()
    types, warns = typer.assign(mol("C1CCCCC1"))
    assert types[:6] == ["c3"] * 6
    assert any("2.11" in w for w in warns)
    z, bonds = gaff2.kekule_graph(mol("C1CCCCC1"))
    assert gaff2.assign_types(z, bonds)[0][:6] == ["c6"] * 6      # AmberTools 2.2 규칙 그대로면 c6


@pytest.mark.parametrize("smiles", ["F[B-](F)(F)F", "[Na+]", "C[Si](C)(C)C"])
def test_type_missing_from_gaff2_lt_raises(smiles):
    with pytest.raises(TypingError, match="gaff2.lt 에 없는"):
        gaff2.Gaff2Typer().assign(mol(smiles))


def test_missing_bonded_parameters_are_reported():
    _, warns = gaff2.Gaff2Typer().assign(mol("CC(N)=O"))      # nt 는 gaff2.lt(2.11)에 결합 파라미터가 없음
    assert any(w.startswith("결합 c-nt") for w in warns)


def test_def_parser_reads_all_rules():
    d = gaff2.default_def()
    assert set(d["wild"]) == {"XX", "XA", "XB", "XC", "XD"}
    assert len(d["rules"]) > 300 and d["rules"][-1]["type"] == "DU"
    cp = next(r for r in d["rules"] if r["type"] == "cp")
    assert len(cp["ces"]) == 3 and all(ch[0]["name"] == "XX" for ch in cp["ces"])
