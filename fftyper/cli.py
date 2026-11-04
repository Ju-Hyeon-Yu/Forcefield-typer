# 역장을 골라 SMILES(+mol2)에서 moltemplate 단일 분자 .lt 를 만드는 CLI (fftyper 명령)
"""사용 예
  fftyper --ff gaff2 --smiles "Cn1cc[n+](C)c1.[O-]S(=O)(=O)C(F)(F)F" --mol2 ip.mol2 -o ip.lt
  fftyper --ff oplsaa --smiles "CCO" -o ethanol.lt
  fftyper --ff dreiding --smiles "CCO" -o ethanol.lt

mol2 를 주면 원자 순서·좌표는 mol2 그대로 쓰고, 결합 차수와 화학 정보는 SMILES 에서
가져온다(두 분자를 원소+연결 그래프로 대응). 전하는 DREIDING·GAFF2 는 mol2 의 값,
OPLS-AA 는 원자 타입에 정해진 값을 쓴다.
"""
import argparse
import os
import re
import sys

from rdkit import Chem

from . import FORCE_FIELDS, TypingError, get_typer
from .molio import embed, match_mol2, read_mol2, write_lt


def main(argv=None):
    ap = argparse.ArgumentParser(prog="fftyper",
                                 description="SMILES(+RESP mol2) → moltemplate .lt (DREIDING / GAFF2 / OPLS-AA)")
    ap.add_argument("--ff", required=True, choices=FORCE_FIELDS, help="역장")
    ap.add_argument("--smiles", required=True, help="이온쌍은 '.' 으로 연결 (예: A.B)")
    ap.add_argument("--mol2", help="좌표·전하를 가져올 mol2 (원자 순서 유지)")
    ap.add_argument("-o", "--out", required=True, help="출력 .lt 경로")
    ap.add_argument("--name", help="moltemplate 분자 이름 (기본: 출력 파일 이름)")
    ap.add_argument("--ff-file", default=None,
                    help="라벨·규칙을 읽을 역장 .lt (기본: 설치된 moltemplate 의 파일)")
    ap.add_argument("--seed", type=int, default=1, help="mol2 없을 때 3D 구조 난수 시드")
    args = ap.parse_args(argv)

    name = args.name or re.sub(r"\W", "_", os.path.splitext(os.path.basename(args.out))[0])
    mol0 = Chem.MolFromSmiles(args.smiles)
    if mol0 is None:
        sys.exit("SMILES 를 읽지 못했습니다")
    try:
        typer = get_typer(args.ff, args.ff_file)
    except TypingError as e:
        sys.exit(f"역장 파일을 읽지 못했습니다: {e}")
    template = Chem.AddHs(mol0)
    n = template.GetNumAtoms()
    warnings = []

    if args.mol2:
        atoms, mbonds = read_mol2(args.mol2)
        try:
            t2m = match_mol2(template, atoms, mbonds)
        except ValueError as e:
            sys.exit(f"mol2 대응 실패: {e}")
        m2t = {m: t for t, m in enumerate(t2m)}
        # 원자 번호를 mol2 순서로 바꿔서 이후 출력·경고 번호가 mol2 와 같게 함
        template = Chem.RenumberAtoms(template, [m2t[k] for k in range(n)])
        xyz = [a["xyz"] for a in atoms]
        q = [a["q"] for a in atoms]
        if not typer.TYPE_CHARGES:
            if not all(a["has_q"] for a in atoms):
                warnings.append("mol2 일부 원자에 전하 열이 없어 0 으로 둡니다")
            charge_note = f"mol2 ({os.path.basename(args.mol2)})"
            qsum, formal = sum(q), Chem.GetFormalCharge(template)
            if abs(qsum - formal) > 1e-3:
                warnings.append(f"mol2 전하 합 {qsum:.4f} 이 SMILES 형식 전하 {formal} 와 다릅니다")
    else:
        xyz = embed(template, args.seed)
        q = [0.0] * n
        if not typer.TYPE_CHARGES:
            charge_note = "없음 — 모두 0 (RESP 는 --mol2 로 넣으세요)"
            warnings.append("mol2 가 없어 전하를 모두 0 으로 썼습니다")

    try:
        labels, w = typer.assign(template)
    except TypingError as e:
        sys.exit(f"타이핑 실패: {e}")
    warnings += w
    if typer.TYPE_CHARGES:
        q = typer.charges(labels)
        charge_note = f"{typer.LABEL} 원자 타입 전하 ({typer.LT_FILE})"
        qsum, formal = sum(q), Chem.GetFormalCharge(template)
        if abs(qsum - formal) > 1e-3:
            warnings.append(f"타입 전하 합 {qsum:.4f} 이 SMILES 형식 전하 {formal} 와 다릅니다")

    frags = Chem.GetMolFrags(template)
    frag_of = {i: k for k, fr in enumerate(frags, 1) for i in fr}
    rows = [{"name": f"{a.GetSymbol()}{a.GetIdx() + 1}", "frag": frag_of[a.GetIdx()],
             "label": labels[a.GetIdx()], "q": q[a.GetIdx()], "xyz": xyz[a.GetIdx()]}
            for a in template.GetAtoms()]
    bonds = sorted(tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in template.GetBonds())
    write_lt(args.out, typer, name, args.smiles, rows, bonds, charge_note)

    width = max(len(r["label"]) for r in rows)
    print(f"{args.out}: {typer.LABEL}, 원자 {n}개, 결합 {len(bonds)}개, 조각 {len(frags)}개")
    for k, r in enumerate(rows, 1):
        print(f"  {k:4d} {r['name']:<6} {r['label']:<{width}} {r['q']:>10.6f}")
    for w in dict.fromkeys(warnings):
        print(f"경고: {w}", file=sys.stderr)


if __name__ == "__main__":
    main()
