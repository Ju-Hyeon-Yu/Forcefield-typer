# foyer 역장 XML 에서 SMARTS 규칙(이름·정의·덮어쓰기·설명)을 파일 순서대로 뽑아 JSON 으로 저장하는 도구
"""사용: python3 tests/foyer_rules.py FOYER_XML OUT_JSON [foyer 판 번호]

foyer(MIT, © 2015 Vanderbilt University)의 forcefields/xml/oplsaa.xml 에서 def 가 있는 타입만
순서 그대로 담는다. fftyper/data/foyer_oplsaa.json 이 이렇게 만든 파일이다.
"""
import json
import sys
import xml.etree.ElementTree as ET


def main():
    src, out = sys.argv[1], sys.argv[2]
    root = ET.parse(src).getroot()
    rules = []
    for t in root.iter("Type"):
        if not t.get("def"):
            continue
        ov = [x.strip() for x in t.get("overrides", "").split(",") if x.strip()]
        rules.append({"name": t.get("name"), "def": t.get("def"), "overrides": ov,
                      "element": t.get("element"), "class": t.get("class"), "desc": t.get("desc", "")})
    pkg = sys.argv[3] if len(sys.argv) > 3 else ""
    meta = {"source": f"foyer {pkg} forcefields/xml/oplsaa.xml (역장 판 {root.get('version', '')})".replace("  ", " "),
            "name": root.get("name"), "license": "MIT, Copyright (c) 2015 Vanderbilt University"}
    with open(out, "w") as f:
        json.dump({"meta": meta, "rules": rules}, f, indent=0, ensure_ascii=False)
    print(f"{len(rules)} 규칙")


if __name__ == "__main__":
    main()
