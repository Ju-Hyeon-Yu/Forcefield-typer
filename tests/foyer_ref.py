# 실제 foyer 로 OPLS-AA 원자 타입 정답지를 만드는 도구 — foyer 검증 분자(.top/.mol2) + SMILES 목록 → JSON
"""사용: python3 tests/foyer_ref.py FOYER_SRC SMILES_FILE OUT_JSON

FOYER_SRC 는 foyer 1.2.0 소스 폴더(github.com/mosdef-hub/foyer)다. parmed, lark, ele, openmm 이 필요하다.
gmso 는 conda 에만 있어서 없으면 빈 대역 모듈을 끼운다(원자 타이핑에는 쓰이지 않음).
foyer 의 find_atomtypes 를 그대로 불러,
- validation: foyer/opls_validation 의 분자를 foyer 시험과 같이 parmed 로 읽은 원자·결합 순서 그대로
- smiles: RDKit AddHs 순서의 원자·결합
각각 원자별 타입 이름 또는 {"error": "multiple"|"none", "atom": 번호} 를 쓴다. 검증 분자는 원소·결합과
.top 에 적힌 기준 타입(reference)도 함께 담아 fftyper 시험이 foyer 없이 돌 수 있게 한다.
"""
import glob
import json
import os
import re
import sys
import tempfile
import types


def load_foyer(src):
    try:
        import gmso  # noqa: F401
    except ImportError:
        stub = types.ModuleType("gmso")
        stub.Topology = type("Topology", (), {})
        sys.modules["gmso"] = stub
    sys.path.insert(0, src)
    import foyer
    from foyer.atomtyper import find_atomtypes
    from foyer.exceptions import FoyerError
    from foyer.topology_graph import TopologyGraph
    return foyer, find_atomtypes, FoyerError, TopologyGraph


def run(find_atomtypes, FoyerError, graph, ff):
    try:
        tm = find_atomtypes(graph, ff)
    except FoyerError as e:
        m = re.search(r"atom (?:numbered )?(\d+)", str(e))
        return {"error": "multiple" if "multiple" in str(e) else "none", "atom": int(m.group(1)) if m else -1}
    return [tm[i]["atomtype"] for i in sorted(tm)]


def main():
    src, smiles_file, out_path = sys.argv[1:4]
    foyer, find_atomtypes, FoyerError, TopologyGraph = load_foyer(src)
    import parmed as pmd
    from rdkit import Chem
    ff = foyer.Forcefield(name="oplsaa")
    out = {"foyer": foyer.__version__, "validation": {}, "smiles": {}}
    vdir = os.path.join(src, "foyer", "opls_validation")
    with open(os.path.join(src, "foyer", "tests", "implemented_opls_tests.txt")) as f:
        out["implemented"] = [x.strip() for x in f if x.strip()]   # foyer 시험이 기준 타입과 맞다고 보는 분자
    for d in sorted(glob.glob(os.path.join(vdir, "*", ""))):
        name = os.path.basename(os.path.dirname(d))
        top = os.path.join(d, f"{name}.top")
        mol2 = glob.glob(os.path.join(d, "*.mol2"))
        if os.path.exists(top):
            cwd = os.getcwd()
            with tempfile.TemporaryDirectory() as tmp:
                os.chdir(tmp)            # foyer 시험처럼 임시 폴더에서 읽음
                try:
                    st = pmd.load_file(top, xyz=os.path.join(d, f"{name}.gro"), parametrize=False)
                finally:
                    os.chdir(cwd)
        elif mol2:
            st = pmd.load_file(mol2[0], structure=True)
        else:
            continue
        g = TopologyGraph.from_parmed(st)
        out["validation"][name] = {
            "z": [a.atomic_number for a in st.atoms],
            "bonds": [[b.atom1.idx, b.atom2.idx] for b in st.bonds],
            "reference": [a.type for a in st.atoms],
            "foyer": run(find_atomtypes, FoyerError, g, ff),
        }
    for smi in open(smiles_file).read().split():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        m = Chem.AddHs(m)
        g = TopologyGraph()
        for a in m.GetAtoms():
            g.add_atom(name=a.GetSymbol(), index=a.GetIdx(), atomic_number=a.GetAtomicNum(), symbol=a.GetSymbol())
        for b in m.GetBonds():
            g.add_bond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        out["smiles"][smi] = run(find_atomtypes, FoyerError, g, ff)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=0, ensure_ascii=False)
    ok = sum(isinstance(v["foyer"], list) for v in out["validation"].values())
    print(f"검증 분자 {len(out['validation'])}개 (foyer 타이핑 성공 {ok}), SMILES {len(out['smiles'])}개 "
          f"(성공 {sum(isinstance(v, list) for v in out['smiles'].values())})")


if __name__ == "__main__":
    main()
