# 계획 — DREIDING·GAFF2·OPLS-AA 원자 타입 → moltemplate 단일 분자 `.lt`

## 목적
SMILES(필수)와 mol2(선택)로부터 moltemplate 배포본 역장 파일의 원자 타입을 자동으로 붙이고,
`import "<역장>.lt"` + `<이름> inherits <역장>` 형태의 단일 분자 `.lt` 를 만든다.
역장 파일은 읽기만 하고 바꾸지 않는다. 같은 기능을 오프라인 단일 HTML 로도 제공한다.
`Ju-Hyeon-Yu/Dreiding_typer` 의 DREIDING 타이퍼를 옮겨 와 출발점으로 삼는다.

## 역장별 기준
| 역장 | moltemplate 파일 | 타입 규칙의 기준 | 전하 |
|---|---|---|---|
| DREIDING | `dreiding.lt` | DREIDING 라벨 매뉴얼 본문 규칙 (이전 저장소에서 확정) | mol2 (RESP) |
| GAFF2 | `gaff2.lt` | AmberTools antechamber `-at gaff2` 와 같은 결과 (ATOMTYPE_GFF2.DEF + bondtype/atomtype 처리 이식) | mol2 (RESP) |
| OPLS-AA | `oplsaa2024.lt` | foyer 의 OPLS-AA SMARTS 규칙을 2024 번호에 맞게 옮김 | 원자 타입 전하 |

- 타입을 정할 수 없는 원자는 어느 원자인지 알려 주고 오류로 멈춘다.
- OPLS-AA 타입 전하의 합이 SMILES 형식 전하와 다르면 경고한다.

## 구성
- `fftyper/` — `cli.py`(fftyper 명령), `molio.py`(mol2·원자 대응·3D·.lt), `dreiding.py`, `gaff2.py`, `oplsaa.py`
- `web/` — 같은 타이퍼의 JS 이식과 역장 선택 메뉴가 있는 단일 HTML 빌드
- `tests/` — 역장별 단위 테스트, 정답지 대조, moltemplate 통합 테스트

## 진행 순서 (2026-11-04 모두 완료 — checklist.md, context-notes.md 참고)
1. 저장소 골격 + DREIDING 이전
2. GAFF2 (Python) — antechamber 대조, parmchk2 대응 타입 보충
3. GAFF2 (JS) + HTML 역장 선택
4. OPLS-AA (Python·JS) — foyer·참조 분자 대조
5. 문서·테스트·배포

## 검증 기준
1. GAFF2: 넓은 분자 집합에서 antechamber(AmberTools 26) 결과와 원자 단위로 같다.
2. OPLS-AA: foyer 결과·참조 토폴로지와 같고, 쓰는 타입은 모두 2024 파일에 같은 뜻으로 있다.
3. 생성한 `.lt` 를 moltemplate 로 돌렸을 때 결합·각·이면각 누락이 없고 LAMMPS 가 읽는다.
4. JS 이식 결과가 Python 과 같다 (라벨·경고·오류).
