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

## 웹 GAFF2 (2026-11-04)
- web/src/gaff2.js 는 gaff2.py 를 반복 순서까지 그대로 옮겼다. 규칙은 build_tables.py --ff gaff2 가 Python 의
  해석 결과(DEF 규칙, gaff2.lt 규칙·계수 줄, EQUA 표, 타입 설명)를 그대로 JSON 으로 넣는다(약 0.5 MB).
- RDKit 케쿨레 구조는 원자 순서에 따라 달라질 수 있다(391 분자 중 82 개가 번호를 섞으면 달라지고 3 개는 ce/cf·ne/nf
  짝이 바뀜). 웹은 SMILES 순서의 케쿨레 구조를 원자 번호만 바꿔 쓰므로, Python CLI 도 번호를 바꾸기 전에
  stamp_kekule 로 결합 속성 "kek" 에 적어 두고 GAFF2 가 그것을 쓰게 했다. ce/cf 같은 짝은 antechamber 에서도
  원자 순서에 따라 정해지며 둘 다 맞는 배정이다.
- SMILES 만 줄 때 3D 시작 구조는 역장과 상관없이 DREIDING 결합 길이·각 표로 만든다. DREIDING 라벨을 붙일 수
  없는 분자(CO2, 알렌, PF6⁻, ClO4⁻, 아이오도늄 등 18/391)는 원소·혼성으로 가까운 라벨을 골라 쓴다. 5·6배위 중심은
  각 대신 리간드끼리 반발(√2 r)을 넣어 PF6⁻·SF6 가 정팔면체가 되게 했다.

## OPLS-AA (2026-11-04, 사용자 "권장대로")
- 판정: foyer 1.2.0(MIT) OPLS-AA 규칙 228개를 그대로 쓰고, 판정 방식을 fftyper/foyer_smarts.py 로 옮김.
  foyer 의 networkx VF2 는 유도 부분그래프 짝맞춤이고, 규칙은 파일 순서·원자는 번호 순으로 맞추며 맞자마자
  whitelist/blacklist 를 고친다(%타입 참조·부정 때문에 순서가 결과에 영향). 고리는 _find_chordless_cycles(8원,
  같은 고리도 출발 이웃 쌍마다 다른 튜플로 셈 → R 수에 그대로 반영)를 그대로 옮김.
- 검증: tests/foyer_ref.py 로 실제 foyer 를 돌린 정답지(검증 분자 182 개, SMILES 391 개 — 실패와 처음 실패 원자까지)
  와 모두 같음. foyer 는 gmso 가 conda 전용이라 빈 대역 모듈로 import. 우리 엔진이 약 16배 빠름.
- 번호: foyer(GROMACS 번호) → oplsaa2024.lt(BOSS 번호). 228 개 중 146 개는 번호·전하·LJ·결합 타입이 같고,
  나머지도 대부분 같은 화학 환경(2024 의 OPLS/2020 값 갱신, 결합 타입 이름 표기 차이)이라 같은 번호로 연결.
  뜻이 다른 번호 3 개(732 아이오도벤젠 vs 2-아이오도프로판, 927·930 알카인 주변)와 2024 에 없는 foyer 추가 타입
  12 개(151a–153c, 722a)는 대응 없음 → 그 원자가 나오면 오류. 484(설포레인 C)·959/960 처럼 foyer 가 다른 설명으로
  같은 번호의 매개변수를 쓰는 것은 foyer 의 선택을 따라 같은 번호로 둠.
- OPLS/2020 짝: 2024 파일은 159(3차 알코올 C)·189(OCH2O 아세탈 C) 전하를 바꾸면서 짝 타입 160(3차 알코올 O)·
  186(아세탈 O)을 두었다. foyer 는 고전 규칙이라 O 를 154·180 으로 주므로, 이웃이 159·189 이면 160·186 으로 바꿈
  (t-부탄올·다이메톡시메테인·1,3-다이옥솔레인이 중성이 됨).
- 전하: 2024 타입 전하. foyer 규칙 자체가 고전 OPLS-AA 의 분자별 전하 조정(바이페닐·큐멘 ipso C, 아세테이트 CH3 등)을
  모르므로 이런 분자는 합이 형식 전하와 다르고, CLI·웹이 경고한다(391 분자 중 타이핑되는 177 개에서 59 개).
- 결합 파라미터: oplsaa2024.lt 는 짧은 이름을 replace{} 로 결합 타입이 붙은 긴 이름으로 바꾸고 By Type 규칙은 긴 이름의
  와일드카드라, LtRules 가 replace 를 읽어 긴 이름으로 맞춤. improper 대칭은 cenIsortJKL.py(중심이 첫 원자).
