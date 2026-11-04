# GAFF2 .lt(대응 타입 보충 포함)를 moltemplate·LAMMPS 로 실제로 돌려 모든 항이 들어가는지 확인하는 통합 테스트
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
from fftyper import cli  # noqa: E402
from fftyper.ltrules import topology  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("moltemplate.sh") is None, reason="moltemplate 없음")

# 아세트아마이드·요소는 같은 보충 타입(c-nt, hn-nt …)을 함께 쓴다 → 두 파일이 같은 정의를 덧붙여도 되는지 확인
SMILES = {
    "acetamide": "CC(N)=O",
    "urea": "NC(N)=O",
    "nma": "CNC(C)=O",
    "pyrrolidone": "O=C1CCCN1",
    "nh4cl": "[NH4+].[Cl-]",
    "tea_otf": "CC[NH+](CC)CC.[O-]S(=O)(=O)C(F)(F)F",
    "benzene": "c1ccccc1",
}
# 중심 원자가 평면인 improper: C=O 탄소 X-X-c-o, 아마이드 N(nt/ns, 보충), 벤젠 X-X-ca-ha
N_IMPROPERS = 2 + 3 + 2 + 2 + 0 + 0 + 6


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    d = tmp_path_factory.mktemp("mt_gaff2")
    for name, smi in SMILES.items():
        cli.main(["--ff", "gaff2", "--smiles", smi, "-o", str(d / f"{name}.lt")])
    with open(d / "system.lt", "w") as f:
        for name in SMILES:
            f.write(f'import "{name}.lt"\n')
        for k, name in enumerate(SMILES):
            f.write(f"m_{name} = new {name}.move({30 * k},0,0)\n")
        f.write(f'write_once("Data Boundary") {{\n -50 {30 * len(SMILES) + 50} xlo xhi\n'
                ' -50 50 ylo yhi\n -50 50 zlo zhi\n}\n')
    r = subprocess.run(["moltemplate.sh", "system.lt"], cwd=d, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    return d


def read_sections(path):
    secs, cur = {}, None
    for line in open(path):
        s = line.split("#")[0].strip()
        if not s:
            continue
        if re.match(r"^[A-Z][A-Za-z ]+$", s):
            cur = s
            secs[cur] = []
        elif cur:
            secs[cur].append(s.split())
    return secs


def test_every_bond_angle_dihedral_is_typed(built):
    data = read_sections(built / "system.data")
    want = {"bond": 0, "angle": 0, "dihedral": 0}
    for smi in SMILES.values():
        t = topology(Chem.AddHs(Chem.MolFromSmiles(smi)))
        for k in want:
            want[k] += len(t[k])
    assert len(data["Bonds"]) == want["bond"]
    assert len(data["Angles"]) == want["angle"]
    assert len(data["Dihedrals"]) == want["dihedral"]
    assert len(data["Impropers"]) == N_IMPROPERS


def test_copied_coefficients_reach_lammps_input(built):
    data = read_sections(built / "system.data")
    masses = {r[0]: r for r in data["Masses"]}
    atype = {r[0]: r[2] for r in data["Atoms"]}
    names = {}                                  # LAMMPS 원자 타입 번호 → GAFF2 이름 (Masses 주석)
    for line in open(built / "system.data"):
        m = re.match(r"\s*(\d+)\s+\S+\s+#\s*(\S+)", line)
        if m and m.group(1) in masses:
            names.setdefault(m.group(1), m.group(2))
    settings = {}
    for line in open(built / "system.in.settings"):
        t = line.split("#")[0].split()
        if t and t[0] == "bond_coeff":
            settings[t[1]] = t[2:]
    c_nt = [r for r in data["Bonds"] if sorted(names[atype[x]] for x in r[2:4]) == ["c", "nt"]]
    assert len(c_nt) == 3                                   # 아세트아마이드 1 + 요소 2
    assert all(settings[r[1]] == ["harmonic", "356.2", "1.3790"] for r in c_nt)


@pytest.mark.skipif(shutil.which("lmp") is None, reason="LAMMPS 없음")
def test_lammps_reads_and_evaluates(built):
    with open(built / "run0.in", "w") as f:
        f.write("include system.in.init\nread_data system.data\ninclude system.in.settings\n"
                "kspace_modify gewald 0.3\n"            # 전하 없이 만든 계라 Ewald 계수를 직접 줌
                "thermo_style custom step pe ebond eangle edihed eimp\nrun 0\n")
    r = subprocess.run(["lmp", "-in", "run0.in", "-log", "none"], cwd=built, capture_output=True, text=True)
    assert r.returncode == 0, (r.stdout + r.stderr)[-2000:]
    lines = r.stdout.splitlines()
    k = next(i for i, line in enumerate(lines) if line.split()[:2] == ["Step", "PotEng"])
    vals = [float(x) for x in lines[k + 1].split()[1:]]
    assert all(math.isfinite(v) for v in vals) and vals[-1] > 0     # improper 에너지도 계산됨
