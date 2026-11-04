# antechamber(AmberTools) 로 GAFF2 타입 정답지를 만드는 도구 — SMILES 목록 → JSON
"""사용: python3 tests/antechamber_ref.py SMILES_FILE OUT_JSON [--amberhome DIR]

SMILES(+H)의 케쿨레 결합 차수를 mol2 로 쓰고 `antechamber -at gaff2 -j 5`(입력 결합 차수 유지,
방향족·비편재 결합만 재지정)로 타입을 받는다. 원자·결합 순서는 RDKit 순서 그대로이므로
fftyper.gaff2 결과와 원자 단위로 비교할 수 있다.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem


def write_mol2(mol, path):
    kek = Chem.Mol(mol)
    Chem.Kekulize(kek, clearAromaticFlags=True)
    AllChem.Compute2DCoords(kek)
    conf = kek.GetConformer()
    with open(path, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\nMOL\n{kek.GetNumAtoms()} {kek.GetNumBonds()} 1 0 0\nSMALL\nUSER_CHARGES\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for a in kek.GetAtoms():
            p = conf.GetAtomPosition(a.GetIdx())
            el = a.GetSymbol()
            f.write(f"{a.GetIdx() + 1} {el}{a.GetIdx() + 1} {p.x:.4f} {p.y:.4f} 0.0000 {el} 1 MOL "
                    f"{float(a.GetFormalCharge()):.4f}\n")
        f.write("@<TRIPOS>BOND\n")
        for b in kek.GetBonds():
            f.write(f"{b.GetIdx() + 1} {b.GetBeginAtomIdx() + 1} {b.GetEndAtomIdx() + 1} "
                    f"{int(b.GetBondTypeAsDouble())}\n")


def run_antechamber(mol, amberhome):
    env = dict(os.environ, AMBERHOME=amberhome, PATH=os.path.join(amberhome, "bin") + ":" + os.environ["PATH"])
    with tempfile.TemporaryDirectory() as d:
        write_mol2(mol, os.path.join(d, "in.mol2"))
        r = subprocess.run(["antechamber", "-i", "in.mol2", "-fi", "mol2", "-o", "out.mol2", "-fo", "mol2",
                            "-at", "gaff2", "-j", "5", "-dr", "no", "-pf", "y"],
                           cwd=d, env=env, capture_output=True, text=True)
        out = os.path.join(d, "out.mol2")
        if r.returncode != 0 or not os.path.exists(out):
            return {"error": (r.stdout + r.stderr).strip().splitlines()[-1:]}
        types, sec = [], None
        for line in open(out):
            s = line.strip()
            if s.startswith("@<TRIPOS>"):
                sec = s
                continue
            if sec == "@<TRIPOS>ATOM" and s:
                types.append(s.split()[5])
        return {"types": types}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("smiles_file")
    ap.add_argument("out_json")
    ap.add_argument("--amberhome", default=os.environ.get("AMBERHOME"))
    args = ap.parse_args()
    if not args.amberhome:
        sys.exit("AMBERHOME 이 필요합니다")
    out = {}
    for smi in open(args.smiles_file).read().split():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        out[smi] = run_antechamber(Chem.AddHs(m), args.amberhome)
    with open(args.out_json, "w") as f:
        json.dump(out, f, indent=0, ensure_ascii=False)
    print(f"{len(out)} 분자, 오류 {sum('error' in v for v in out.values())}개")


if __name__ == "__main__":
    main()
