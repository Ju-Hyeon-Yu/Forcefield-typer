# 작업 메모 — 결정과 근거

## 출발점 (2026-11-04)
- `Ju-Hyeon-Yu/Dreiding_typer`(커밋 6d0934f)의 DREIDING 타이퍼·CLI·HTML 도구를 옮겨 옴.
  DREIDING 라벨 규칙의 결정 근거(매뉴얼 본문 기준, 아마이드 짝 규칙 등)는 그 저장소의
  `context-notes.md` 에 있음. 기존 저장소는 그대로 둠.

## 사용자 결정 (2026-11-04)
- 새 저장소 `Ju-Hyeon-Yu/Forcefield-typer` 에 DREIDING·GAFF2·OPLS-AA 를 모두 넣는다.
- 명령은 하나(`fftyper --ff dreiding|gaff2|oplsaa`), HTML 도 하나에 역장 선택 메뉴.
- 전하: DREIDING·GAFF2 는 mol2 의 RESP, OPLS-AA 는 원자 타입 전하. 타입 전하 합이 형식 전하와
  다르면 경고.
- OPLS-AA: moltemplate `oplsaa2024.lt` 번호 + 타입 전하. 번호를 고르는 규칙은 foyer 규칙표를
  2024 번호에 맞게 옮겨 씀. LigParGen 은 결과가 분자별 새 타입(opls_800…)·CM1A 전하이고 웹에서만
  돌아가 이 목적과 맞지 않음(사용자에게 설명 후 추천안 승인).
- 타입이 없는 원자는 오류로 멈춤.
- GAFF2 먼저, 그다음 OPLS-AA.
- 라이선스 GPL-3.0: GAFF2 를 antechamber 와 같게 만들려면 AmberTools(GPL-2.0-or-later) 규칙 파일과
  처리 코드를 이식해야 하므로 Apache 2.0 에서 바꿈(사용자 승인).

## 조사
- moltemplate 2.22.5 의 `oplsaa.lt` 는 폐기 예정, `oplsaa2024.lt`(Jorgensen 2024, 타입 990개, BOSS 번호) 사용 권고.
  `oplsaa2008.lt` 는 TINKER 번호(135 = 아세탈 C)라 foyer·GROMACS 번호와 다름.
- 2024 파일에는 암모늄(R4N⁺, RNH3⁺, R3NH⁺), 아이오도벤젠·아이오도알케인, I⁻, PF6⁻, 술폰·술폰아마이드는 있으나
  술포네이트 음이온(–SO3⁻), 트리플레이트, TFSI, 이미다졸륨 양이온 타입은 없음.
- foyer `oplsaa.xml`(MIT, 825 타입)은 GROMACS 번호. 2024 와 대부분 같으나 다른 것도 있음
  (732: 2024 아이오도벤젠 I, foyer 아이오도프로판 I) → 설명 대조로 대응표를 만들어야 함.
- GAFF2: moltemplate `gaff2.lt` 타입 98개. AmberTools 26 (conda-forge) 설치해 antechamber 를 정답지로 씀.
  antechamber `-j 5` = bondtype `-j part`(입력 결합 차수 유지, 방향족·비편재 결합만 재지정) → atomtype `-p gaff2`.
  우리 입력은 SMILES 결합 차수이므로 이 경로가 맞음.
