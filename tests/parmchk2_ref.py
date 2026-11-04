# parmchk2(AmberTools)로 GAFF2 2.11 대응 타입(EQUA) 보충 항 정답지를 만드는 도구 — SMILES 목록 → JSON
"""사용: python3 tests/parmchk2_ref.py SMILES_FILE TYPES_JSON OUT_JSON [--amberhome DIR]

TYPES_JSON 은 tests/antechamber_ref.py 로 만든 antechamber 타입이다. moltemplate gaff2.lt 가
GAFF2 2.11 판이므로 c5·c6 은 c3 로 바꾸고, gaff2.lt 에 없는 타입이 있는 분자는 뺀다.
`parmchk2 -p gaff211.dat`(moltemplate gaff2.lt 와 같은 원본) 결과 frcmod 에서 대응 타입(EQUA)이나
일반 원자(X)로만 바꿔 찾은 항(벌점 0)을 모은다. 원자 순서는 RDKit 순서 그대로다.
"""
import argparse
import itertools
import json
import os
import re
import subprocess
import sys
import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from fftyper import moltemplate_ff_path  # noqa: E402
from fftyper.gaff2 import EQUA  # noqa: E402
from fftyper.ltrules import LtRules  # noqa: E402

KIND = {"BOND": "bond", "ANGLE": "angle", "DIHE": "dihedral", "IMPROPER": "improper"}
NAT = {"bond": 2, "angle": 3, "dihedral": 4, "improper": 4}


def write_mol2(mol, types, path):
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()
    with open(path, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\nMOL\n{mol.GetNumAtoms()} {mol.GetNumBonds()} 1 0 0\nSMALL\nNO_CHARGES\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for a in mol.GetAtoms():
            p = conf.GetAtomPosition(a.GetIdx())
            f.write(f"{a.GetIdx() + 1} {a.GetSymbol()}{a.GetIdx() + 1} {p.x:.4f} {p.y:.4f} 0.0000 "
                    f"{types[a.GetIdx()]} 1 MOL 0.0\n")
        f.write("@<TRIPOS>BOND\n")
        for b in mol.GetBonds():
            f.write(f"{b.GetIdx() + 1} {b.GetBeginAtomIdx() + 1} {b.GetEndAtomIdx() + 1} 1\n")


def equa_derived(kind, names, src):
    """src 가 names 를 대응 타입으로 한 자리 이상 바꾸고 나머지는 같거나 X 인지 (정·역방향, improper 는
    바깥 세 원자 순서 무관)."""
    cands = [src, src[::-1]]
    if kind == "improper":
        cands = [(p[0], p[1], src[2], p[2]) for p in itertools.permutations((src[0], src[1], src[3]))]
    for c in cands:
        if all(s == o or s == "X" or s in EQUA.get(o, []) for o, s in zip(names, c)) and \
                any(s in EQUA.get(o, []) for o, s in zip(names, c)):
            return True
    return False


def parse_frcmod(path):
    """frcmod → [(종류, 원자 타입 튜플, 출처 튜플, 값)] (대응 타입 유래 항만)."""
    out, kind, rows = [], None, []

    def flush():
        if not rows:
            return
        names, texts = rows[0][0], [r[1] for r in rows]
        m = re.search(r"[Ss]ame as\s+([\w\s+-]+?)\s*(,|$)", texts[-1])
        if m:
            src = tuple(x.strip() for x in m.group(1).split("-"))
            zero = kind == "improper" or "penalty score=  0.0" in texts[-1]
            if zero and equa_derived(kind, names, src):
                nums = [[float(x) for x in re.findall(r"-?\d+\.?\d*", t.split("ame as")[0])] for t in texts]
                if kind == "dihedral":
                    vals = [[round(v[1] / v[0], 4), int(abs(v[3])), round(v[2], 1)] for v in nums]
                elif kind == "improper":
                    vals = [round(nums[0][0], 4), round(nums[0][1], 1), int(nums[0][2])]
                else:
                    vals = [round(x, 4) for x in nums[0][:2]]
                out.append((kind, list(names), list(src), vals))
        rows.clear()

    for line in open(path):
        head = line.strip()
        if head in KIND or head in ("MASS", "NONBON"):
            flush()
            kind = KIND.get(head)
            continue
        if not head or kind is None:
            flush()
            continue
        width = 3 * NAT[kind] - 1
        names = tuple(x.strip() for x in line[:width].split("-"))
        if rows and rows[0][0] != names:
            flush()
        rows.append((names, line[width:].rstrip("\n")))
    flush()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("smiles_file")
    ap.add_argument("types_json")
    ap.add_argument("out_json")
    ap.add_argument("--amberhome", default=os.environ.get("AMBERHOME"))
    args = ap.parse_args()
    if not args.amberhome:
        sys.exit("AMBERHOME 이 필요합니다")
    env = dict(os.environ, AMBERHOME=args.amberhome)
    parm = os.path.join(args.amberhome, "dat", "leap", "parm", "gaff211.dat")
    masses = LtRules(moltemplate_ff_path("gaff2.lt")).masses
    ref = json.load(open(args.types_json))
    out = {}
    for smi in open(args.smiles_file).read().split():
        if "types" not in ref.get(smi, {}):
            continue
        types = ["c3" if t in ("c5", "c6") else t for t in ref[smi]["types"]]
        if any(t not in masses for t in types):
            continue
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        with tempfile.TemporaryDirectory() as d:
            write_mol2(mol, types, os.path.join(d, "in.mol2"))
            subprocess.run([os.path.join(args.amberhome, "bin", "parmchk2"), "-i", "in.mol2", "-f", "mol2",
                            "-o", "out.frcmod", "-p", parm], cwd=d, env=env, check=True, capture_output=True)
            out[smi] = parse_frcmod(os.path.join(d, "out.frcmod"))
    with open(args.out_json, "w") as f:
        json.dump(out, f, indent=0, ensure_ascii=False)
    print(f"{len(out)} 분자, 대응 타입 항 {sum(len(v) for v in out.values())}개")


if __name__ == "__main__":
    main()
