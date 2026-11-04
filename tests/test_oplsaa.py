# OPLS-AA 타이퍼 시험 — foyer 판정 재현, oplsaa2024.lt 번호·전하 대응, moltemplate·LAMMPS 통합
import json
import math
import os
import re
import shutil
import subprocess
import sys

import pytest
from rdkit import Chem

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from fftyper import TypingError, cli  # noqa: E402
from fftyper import foyer_smarts as F  # noqa: E402
from fftyper.ltrules import topology  # noqa: E402
from fftyper.oplsaa import OplsTyper, default_rules  # noqa: E402

REF = json.load(open(os.path.join(ROOT, "tests", "data", "foyer_oplsaa_reference.json")))
RULES = F.load_rules(default_rules()["rules"])


def engine(z, bonds):
    """foyer_ref.py 와 같은 꼴: 원자별 타입 또는 {"error", "atom"} (처음 문제 원자)."""
    cands = F.find_atomtypes(F.Graph(z, bonds), RULES)
    for i, c in enumerate(cands):
        if len(c) != 1:
            return {"error": "multiple" if c else "none", "atom": i}
    return [c[0] for c in cands]


def mol(smiles):
    return Chem.AddHs(Chem.MolFromSmiles(smiles))


def test_same_as_foyer_on_validation_molecules():
    """foyer/opls_validation 182 분자: 실제 foyer 1.2.0 결과(실패 포함)와 같아야 함."""
    diffs = [name for name, v in REF["validation"].items() if engine(v["z"], v["bonds"]) != v["foyer"]]
    assert len(REF["validation"]) > 180
    assert not diffs, diffs[:5]


def test_same_as_foyer_on_smiles_set():
    """RDKit 원자 순서의 391 분자: 실제 foyer 결과(실패 포함)와 같아야 함."""
    diffs = []
    for smi, want in REF["smiles"].items():
        m = mol(smi)
        got = engine([a.GetAtomicNum() for a in m.GetAtoms()],
                     [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()])
        if got != want:
            diffs.append(smi)
    assert len(REF["smiles"]) > 380
    assert not diffs, diffs[:5]


def test_matches_reference_topologies_foyer_implements():
    """foyer 시험이 '맞게 구현됨'으로 보는 146 분자는 GROMACS 기준 타입과도 같음 (글자 접미사 빼고)."""
    bad = []
    for name in REF["implemented"]:
        v = REF["validation"][name]
        got = engine(v["z"], v["bonds"])
        if not isinstance(got, list) or [g.rstrip("abcde") for g in got] != v["reference"]:
            bad.append(name)
    assert len(REF["implemented"]) == 146
    assert not bad, bad[:5]


@pytest.mark.parametrize("smiles, expect", [
    ("CCO", ["135", "157", "154", "140", "140", "140", "140", "140", "155"]),
    ("c1ccccc1", ["145"] * 6 + ["146"] * 6),
    ("CC(C)(C)O", ["135", "159", "135", "135", "160"]),      # OPLS/2020 3차 알코올 짝: O 154 → 160
    ("COCOC", ["181", "186", "189", "186", "181"]),          # OPLS/2020 아세탈 짝: O 180 → 186
    ("CNC(C)=O", ["242", "238", "235", "135", "236"]),
])
def test_known_2024_types(smiles, expect):
    labels, _ = OplsTyper().assign(mol(smiles))
    assert labels[:len(expect)] == expect


@pytest.mark.parametrize("smiles", ["CCO", "CC(C)(C)O", "COCOC", "C1COCO1", "CNC(C)=O", "CC(=O)O", "c1ccncc1",
                                    "Cc1ccccc1", "CCO.[Li+].[Cl-]"])
def test_type_charges_are_neutral(smiles):
    m = mol(smiles)
    typer = OplsTyper()
    labels, _ = typer.assign(m)
    assert abs(sum(typer.charges(labels)) - Chem.GetFormalCharge(m)) < 1e-6


def test_untypable_and_unmapped_atoms_stop_with_error():
    with pytest.raises(TypingError, match="타입이 없는 원자: 1\\(C\\), 2\\(N\\)"):
        OplsTyper().assign(mol("Cn1cc[n+](C)c1"))               # 이미다졸륨: foyer 규칙 없음
    with pytest.raises(TypingError, match="opls_927"):
        OplsTyper().assign(mol("CC#C"))                         # 번호 뜻이 다른 타입


def test_cli_writes_type_charges_and_warns_non_neutral(tmp_path, capsys):
    out = tmp_path / "ipb.lt"
    cli.main(["--ff", "oplsaa", "--smiles", "CC(C)c1ccccc1", "-o", str(out)])
    text = out.read_text()
    assert 'import "oplsaa2024.lt"' in text and "ipb inherits OPLSAA {" in text
    assert re.search(r"@atom:145\s+-0\.115000", text)
    assert "타입 전하 합 -0.1150 이 SMILES 형식 전하 0 와 다릅니다" in capsys.readouterr().err


# --------------------------------------------------------------------------- moltemplate·LAMMPS
SYSTEM = {"eth": "CCO", "tbu": "CC(C)(C)O", "dmm": "COCOC", "bz": "c1ccccc1", "tol": "Cc1ccccc1",
          "nma": "CNC(C)=O", "aca": "CC(=O)O", "pyr": "c1ccncc1"}


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    if shutil.which("moltemplate.sh") is None:
        pytest.skip("moltemplate 없음")
    d = tmp_path_factory.mktemp("mt_opls")
    for name, smi in SYSTEM.items():
        cli.main(["--ff", "oplsaa", "--smiles", smi, "-o", str(d / f"{name}.lt")])
    with open(d / "system.lt", "w") as f:
        for name in SYSTEM:
            f.write(f'import "{name}.lt"\n')
        for k, name in enumerate(SYSTEM):
            f.write(f"m_{name} = new {name}.move({30 * k},0,0)\n")
        f.write(f'write_once("Data Boundary") {{\n -50 {30 * len(SYSTEM) + 50} xlo xhi\n'
                ' -50 50 ylo yhi\n -50 50 zlo zhi\n}\n')
    r = subprocess.run(["moltemplate.sh", "system.lt"], cwd=d, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    return d


def test_moltemplate_types_every_term(built):
    counts = {}
    for line in open(built / "system.data"):
        m = re.match(r"\s*(\d+)\s+(bonds|angles|dihedrals)$", line)
        if m:
            counts[m.group(2)] = int(m.group(1))
    want = {"bonds": 0, "angles": 0, "dihedrals": 0}
    for smi in SYSTEM.values():
        t = topology(mol(smi))
        want["bonds"] += len(t["bond"])
        want["angles"] += len(t["angle"])
        want["dihedrals"] += len(t["dihedral"])
    assert counts == want


@pytest.mark.skipif(shutil.which("lmp") is None, reason="LAMMPS 없음")
def test_lammps_runs_with_type_charges(built):
    with open(built / "run0.in", "w") as f:
        f.write("include system.in.init\nread_data system.data\ninclude system.in.settings\n"
                "include system.in.charges\nvariable q equal charge(all)\n"
                "thermo_style custom step pe ebond eangle edihed eimp v_q\nrun 0\n")
    r = subprocess.run(["lmp", "-in", "run0.in", "-log", "none"], cwd=built, capture_output=True, text=True)
    assert r.returncode == 0, (r.stdout + r.stderr)[-2000:]
    lines = r.stdout.splitlines()
    k = next(i for i, line in enumerate(lines) if line.split()[:2] == ["Step", "PotEng"])
    vals = [float(x) for x in lines[k + 1].split()[1:]]
    assert all(math.isfinite(v) for v in vals)
    assert abs(vals[-1]) < 1e-6                     # 모든 분자가 중성 → 계 전체 전하 0
