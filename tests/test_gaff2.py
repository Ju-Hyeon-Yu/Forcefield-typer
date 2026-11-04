# GAFF2 타이퍼를 antechamber(AmberTools 26) 결과와 원자 단위로 대조하는 테스트
import json
import os
import sys

import pytest
from rdkit import Chem

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from fftyper import TypingError  # noqa: E402
from fftyper import cli, gaff2  # noqa: E402

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
    _, warns = gaff2.Gaff2Typer().assign(mol("N#C[N-]C#N"))   # 다이사이안아마이드: cg-n2 결합은 대응 타입도 없음
    assert any(w.startswith("결합 cg-n2") and "오류로 멈춥니다" in w for w in warns)
    assert any(w.startswith("결합각 n1-cg-n2") for w in warns)


EQUA_REF = json.load(open(os.path.join(ROOT, "tests", "data", "gaff2_parmchk2_equa.json")))


def _coeff_values(kind, text):
    """gaff2.lt 계수 줄 → parmchk2 frcmod 와 같은 꼴의 값."""
    t = text.split()
    r = lambda x: round(float(x), 4)                        # noqa: E731  (정답지는 소수 4자리)
    if kind in ("bond", "angle"):
        return [r(t[1]), r(t[2])]
    if kind == "dihedral":                                  # fourier n K1 n1 d1 ...
        return sorted([r(t[2 + 3 * k]), int(t[3 + 3 * k]), r(t[4 + 3 * k])] for k in range(int(t[1])))
    return [r(t[1]), 180.0 if int(t[2]) == -1 else 0.0, int(t[3])]           # cvff K d n


def test_equa_fill_same_as_parmchk2():
    """tests/parmchk2_ref.py 정답지(parmchk2 -p gaff211.dat 의 대응 타입 항)와 항 목록·값이 같아야 함.

    예외 하나: gaff2.dat 2.11 에는 방향만 반대인 같은 이면각이 값이 다르게 두 번 있는 경우가 있다
    (c-n-c3-c3 와 c3-c3-n-c). parmchk2 는 파일 앞쪽 것을, moltemplate 는 뒤쪽 규칙을 쓰므로 우리는
    "대응 타입이었다면 moltemplate 가 쓸 값"(뒤쪽)을 따른다.
    """
    typer = gaff2.Gaff2Typer()
    rules = typer.rules
    exact = {tuple(p): n for n, p in rules.rules["dihedral"] if "*" not in "".join(p)}
    n, diffs = 0, []
    for smi, ref in EQUA_REF.items():
        m = mol(smi)
        types, _ = typer.assign(m)
        ours = {(k, key): src for k, key, src in gaff2.equa_fill(m, types, rules)}
        want = {(k, tuple(names)): vals for k, names, _, vals in ref}
        if set(ours) != set(want):
            diffs.append((smi, sorted(set(ours) ^ set(want))))
            continue
        for (kind, key), src in ours.items():
            n += 1
            exp = sorted(want[kind, key]) if kind == "dihedral" else want[kind, key]
            got = _coeff_values(kind, rules.coeffs[kind][src])
            if got == exp:
                continue
            pats = tuple(dict(rules.rules[kind])[src])
            rev = exact.get(pats[::-1]) if kind == "dihedral" else None
            if rev and rev != src and _coeff_values(kind, rules.coeffs[kind][rev]) == exp:
                continue                                    # 방향만 반대인 중복 이면각 (위 설명)
            diffs.append((smi, kind, key, src, got, exp))
    assert n > 700
    assert not diffs, diffs[:5]


def test_equa_block_is_written_before_the_molecule(tmp_path):
    out = tmp_path / "ace.lt"
    cli.main(["--ff", "gaff2", "--smiles", "CC(N)=O", "-o", str(out)])
    text = out.read_text()
    assert text.index('import "gaff2.lt"') < text.index("\nGAFF2 {\n") < text.index("ace inherits GAFF2")
    assert "bond_coeff @bond:c-nt harmonic 356.2 1.3790   # = c-n" in text
    assert "@bond:c-nt @atom:c @atom:nt" in text
    assert 'write_once("Data Impropers By Type (gaff_imp.py)")' in text
    assert "@improper:c-hn-nt-hn @atom:c @atom:hn @atom:nt @atom:hn" in text


def test_no_equa_block_when_nothing_is_missing(tmp_path):
    out = tmp_path / "bz.lt"
    cli.main(["--ff", "gaff2", "--smiles", "c1ccccc1", "-o", str(out)])
    assert "\nGAFF2 {\n" not in out.read_text()


def test_def_parser_reads_all_rules():
    d = gaff2.default_def()
    assert set(d["wild"]) == {"XX", "XA", "XB", "XC", "XD"}
    assert len(d["rules"]) > 300 and d["rules"][-1]["type"] == "DU"
    cp = next(r for r in d["rules"] if r["type"] == "cp")
    assert len(cp["ces"]) == 3 and all(ch[0]["name"] == "XX" for ch in cp["ces"])
