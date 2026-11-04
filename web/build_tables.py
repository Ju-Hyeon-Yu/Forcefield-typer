# 역장 .lt 에서 브라우저용 규칙표를 JSON 으로 추출 (DREIDING 라벨표, GAFF2 판정 규칙·gaff2.lt 규칙표)
"""사용: python3 web/build_tables.py [--ff dreiding|gaff2|oplsaa] [역장 .lt 경로] > web/build/<역장>_tables.json

DREIDING: 라벨 쌍별 결합 차수와 이면각 규칙 존재 여부는 fftyper.dreiding.DreidingFF 로 계산하므로
Python 타이퍼와 판정이 같다. r0 와 평형각은 3D 좌표 생성(최소화)에만 쓴다.
GAFF2: ATOMTYPE_GFF2.DEF 해석 결과, gaff2.lt 의 질량·By Type 규칙·계수 줄, 대응 타입(EQUA) 표,
타입 설명을 fftyper 의 것 그대로 담는다.
OPLS-AA: foyer 규칙을 Python 으로 해석한 패턴, 2024 번호 대응 예외·OPLS/2020 짝, oplsaa2024.lt 의
질량·By Type 규칙·긴 이름·타입 전하·설명을 담는다.
"""
import fnmatch
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fftyper.dreiding import DreidingFF  # noqa: E402


def read_sections(path):
    """write_once 블록별 줄 목록."""
    secs, cur = {}, None
    with open(path) as f:
        for line in f:
            s = line.strip()
            m = re.match(r'write_once\("([^"]+)"\)', s)
            if m:
                cur = m.group(1)
                secs.setdefault(cur, [])
                continue
            if s.startswith("}"):
                cur = None
                continue
            if cur and s and not s.startswith("#"):
                secs[cur].append(s.split())
    return secs


def last_match(rules, names):
    """moltemplate By Type 규칙: 나중 규칙이 이긴다."""
    hit = None
    for rule_name, pats in rules:
        if all(fnmatch.fnmatchcase(n, p) for n, p in zip(names, pats)) or \
                all(fnmatch.fnmatchcase(n, p) for n, p in zip(names[::-1], pats)):
            hit = rule_name
    return hit


def gaff2_tables(path):
    from fftyper import gaff2, moltemplate_ff_path
    from fftyper.ltrules import LtRules
    path = path or moltemplate_ff_path("gaff2.lt")
    rules = LtRules(path)
    desc = {}
    for line in open(path):
        m = re.match(r"\s*#\s+@atom:(\S+)\s+(.*\S)", line)
        if m:
            desc.setdefault(m.group(1), m.group(2))
    return {
        "source": os.path.basename(path),
        "moltemplate": _mt_version(),
        "def": gaff2.default_def(),
        "equa": gaff2.EQUA,
        "lt": {"masses": rules.masses, "rules": rules.rules, "coeffs": rules.coeffs,
               "sections": rules.sections, "improperCenter": rules.improper_center},
        "desc": {t: desc.get(t, "") for t in rules.masses},
    }


def oplsaa_tables(path):
    from fftyper import moltemplate_ff_path, oplsaa
    from fftyper.foyer_smarts import parse
    from fftyper.ltrules import LtRules
    path = path or moltemplate_ff_path("oplsaa2024.lt")
    rules = LtRules(path)
    desc = {}
    for line in open(path):
        m = re.match(r"\s*set type @atom:(\S+)\s+charge\s+\S+\s*#\s*(.*\S)", line)
        if m:
            desc[m.group(1)] = re.sub(r"\s+", " ", m.group(2))
    data = oplsaa.default_rules()
    return {
        "source": os.path.basename(path),
        "moltemplate": _mt_version(),
        "foyer": data["meta"],
        "rules": [{"name": r["name"], "overrides": r["overrides"], "pat": parse(r["def"]), "def": r["def"],
                   "desc": r.get("desc", "")} for r in data["rules"]],
        "no2024": oplsaa.NO_2024,
        "pairs": oplsaa.OPLS2020_PAIRS,
        "lt": {"masses": rules.masses, "rules": rules.rules, "coeffs": {k: {} for k in rules.coeffs},
               "sections": rules.sections, "improperCenter": rules.improper_center,
               "full": rules.full, "charges": rules.charges},
        "desc": desc,
    }


def _mt_version():
    try:
        from importlib.metadata import version
        return version("moltemplate")
    except Exception:
        return "unknown"


def main():
    args = sys.argv[1:]
    ffname = "dreiding"
    if args[:1] == ["--ff"]:
        ffname, args = args[1], args[2:]
    path = args[0] if args else None
    if ffname in ("gaff2", "oplsaa"):
        out = gaff2_tables(path) if ffname == "gaff2" else oplsaa_tables(path)
        json.dump(out, sys.stdout, separators=(",", ":"), ensure_ascii=False)
        return
    ff = DreidingFF(path)
    secs = read_sections(ff.path)
    settings = [t for k, v in secs.items() if k == "In Settings" for t in v]
    bond_r0 = {t[1].split(":", 1)[1]: float(t[3]) for t in settings if t[0] == "bond_coeff"}
    angle_t0 = {t[1].split(":", 1)[1]: float(t[3]) for t in settings if t[0] == "angle_coeff"}
    bond_rules = [(t[0].split(":", 1)[1], [x.split(":", 1)[1] for x in t[1:3]])
                  for t in secs["Data Bonds By Type"]]
    angle_rules = [(t[0].split(":", 1)[1], [x.split(":", 1)[1] for x in t[1:4]])
                   for t in secs["Data Angles By Type"]]

    labels = sorted(ff.labels)
    pairs = {}
    for i, a in enumerate(labels):
        for b in labels[i:]:
            order = ff.bond_order(a, b)
            if order is None:
                continue
            r0 = bond_r0[last_match(bond_rules, [a, b])]
            pairs[f"{a}|{b}"] = [round(order, 3), r0, 1 if ff.has_dihedral(a, b) else 0]
    theta0 = {}
    for c in labels:
        rule = last_match(angle_rules, ["X", c, "X"])
        if rule:
            theta0[c] = angle_t0[rule]
    out = {
        "source": os.path.basename(ff.path),
        "moltemplate": _mt_version(),
        "labels": {k: ff.labels[k] for k in labels},
        "pairs": pairs,          # "A|B" (A<=B) → [결합 차수, r0, 이면각 규칙 유무]
        "theta0": theta0,        # 중심 라벨 → 평형각(도)
    }
    json.dump(out, sys.stdout, separators=(",", ":"), ensure_ascii=False)


if __name__ == "__main__":
    main()
