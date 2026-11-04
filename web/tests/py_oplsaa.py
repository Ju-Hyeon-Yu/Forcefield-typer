# SMILES 목록을 Python OPLS-AA 타이퍼로 타이핑해 타입·경고·전하·오류를 JSON 으로 출력 (JS 대조용)
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from rdkit import Chem  # noqa: E402
from fftyper import TypingError  # noqa: E402
from fftyper.oplsaa import OplsTyper  # noqa: E402

typer = OplsTyper()
out = {}
for smi in open(sys.argv[1]).read().split():
    m = Chem.MolFromSmiles(smi)
    if m is None:
        out[smi] = {"parse_error": True}
        continue
    m = Chem.AddHs(m)
    try:
        labels, warnings = typer.assign(m)
        out[smi] = {"types": labels, "warnings": warnings, "charges": typer.charges(labels)}
    except TypingError as e:
        out[smi] = {"error": str(e)}
json.dump(out, sys.stdout, ensure_ascii=False, indent=0)
