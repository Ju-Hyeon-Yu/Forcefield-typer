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

## GAFF2 결합 파라미터 보충 (2026-11-04)
- moltemplate `gaff2.lt` 는 AmberTools 26 의 `gaff211.dat` 를 amber2lt 로 바꾼 것과 주석 빼고 똑같다
  (직접 다시 만들어 diff 0). 이 판에는 ns, nt, n7–n9, nx–nz, n+, nu, nv, n5, n6, cs 의 결합 파라미터가
  없다. AmberTools 에서는 parmchk2 가 PARMCHK.DAT 의 EQUA(같은 타입, 벌점 0)로 값을 복사해 채운다.
- 사용자 선택(1안): 값을 그대로 복사할 수 있는 EQUA 항만 생성 .lt 에 함께 쓴다. gaff2.lt 는 그대로 두고,
  `import "gaff2.lt"` 뒤에 `GAFF2 { write_once("In Settings") {…} write_once("Data … By Type") {…} }` 로
  역장 객체에 덧붙인다. 같은 정의를 여러 분자 파일이 덧붙여도 moltemplate 는 문제없다(통합 테스트).
  CORR(비슷한 타입, 벌점 > 0)·경험식 값은 넣지 않고 남은 누락은 경고로 알린다.
- 단계는 parmchk2 를 그대로 따른다: 결합·각은 원래 규칙 → EQUA. 이면각은 원래 특정 규칙 → EQUA 특정
  → 원래 일반(X) → 가운데 두 원자만 EQUA 한 일반 규칙. 반복 순서도 parmchk2 와 같다.
- improper 는 moltemplate(gaff_imp.py) 방식으로 판단한다(바깥 세 원자 순서 무관). parmchk2 는 바깥 원자를
  결합 순서대로 비교하는 버릇이 있어 결과 이름이 다를 수 있으나 값(1.1/180/2)은 같다.
- gaff2.dat 2.11 에는 방향만 반대인 같은 이면각이 값이 다르게 두 번 있다(c-n-c3-c3/c3-c3-n-c,
  c3-c3-c-o/o-c-c3-c3, c2-ce-ca-ca/ca-ca-ce-c2, c3-c3-os-c3/c3-os-c3-c3). tleap 은 주기별로 앞의 것을,
  parmchk2 는 앞의 것을, moltemplate 는 뒤의 규칙을 쓴다. EQUA 로 채울 때는 "대응 타입이었다면 moltemplate 가
  쓸 값"(뒤의 것)을 쓴다 — 같은 계 안에서 ns 아마이드와 n 아마이드가 같은 값을 받게 하려는 것.
- 결합 타입이 없으면 moltemplate 는 "No bond types defined" 로 멈추고, 각·이면각은 조용히 빠진다.
  경고 문구를 이에 맞춤.

## GAFF2 검증 (2026-11-04)
- parmchk2 대조: tests/parmchk2_ref.py (`parmchk2 -p gaff211.dat`) 정답지 378 분자 737 항과 항 목록이 같고,
  값은 위 중복 이면각만 다름 → tests/test_gaff2.py.
- AMBER 와 항별 대조(임시 스크립트, 같은 좌표, tleap+gaff211.dat+parmchk2 frcmod vs moltemplate):
  결합 4950 개 전부 같음, 각 8064 같음·37 개는 moltemplate 쪽 누락(모두 경고), 이면각 8933 같음·192 개 누락
  (모두 경고)·35 개는 위 중복 이면각, improper 814 같음·476 개는 parmchk2 가 CORR(벌점)로 더한 것(예 cc/cd/c2
  중심 X-X-ca-ha), 20 개는 카복실산 C 에서 moltemplate X-X-c-o(10.5) vs parmchk2 CORR(1.1). moltemplate 에만
  있는 improper 는 없음. 결합 타입이 없어 moltemplate 가 멈추는 분자 4개(다이사이안아마이드, 트라이사이아노메타나이드,
  FSI, P–H)는 경고로 알림.
