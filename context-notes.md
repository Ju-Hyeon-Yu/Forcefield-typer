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
- 웹 OPLS-AA: 패턴은 build_tables.py 가 Python 해석기(foyer_smarts.parse)로 미리 해석한 JSON 을 넣어 JS 는
  짝맞춤·반복만 한다(해석기를 두 번 짜지 않으려는 것). oplsaa2024.lt 의 By Type 규칙·긴 이름·타입 전하·설명(약 0.4 MB).
  좌표 파일을 넣어도 전하는 타입 전하를 쓰고 mol2 전하 경고는 내지 않는다(Python CLI 와 같음).

## 타입 설명 칸 (2026-11-05, 사용자 "뜻+판정 근거+실제 환경", "원문 + 한국어 풀이", "HTML만")
- 설명은 web/src/explain.js 가 실행 결과(pipeline.run 반환값)만으로 만든다. 판정을 다시 흉내 내지 않도록 타이퍼가
  판정 정보를 함께 돌려준다: GAFF2 는 원자별 처음 맞은 DEF 규칙 번호(ruleIdx)·조정 전 타입(preTypes)·c5/c6 원자·
  antechamber 고리(rg, ar)·결합 형식(gbonds), OPLS-AA 는 foyer 의 whitelist·blacklist·고리 튜플·2024 번호.
  DREIDING 은 baseLabel·hFlag 갈림길을 같은 순서로 문장으로 옮겼고, 테스트가 문장 속 라벨·플래그를 실제 갈림길과
  맞대 본다(6820 원자 일치).
- DEF 규칙 원문 줄은 parse_def 가 "line" 으로 보관하고(ATD … & 앞까지), foyer 규칙 원문 SMARTS 는 build_tables 가
  "def" 로 넣는다. foyer 패턴의 나무 모양은 해석기가 결합을 사슬·가지 순서로 먼저 넣고 고리 닫기를 끝에 넣는 점을
  이용해 앞 n-1 개 결합을 부모-자식으로 쓴다.
- 원문 설명은 역장 파일 그대로(영어) 두고, GAFF2 는 98 개 타입의 한국어 풀이 표를 따로 둔다(원문의 오타·잘못된 말
  — "Phosphate with … connected atoms" 는 P 원자 — 은 풀이에서 바로잡음). OPLS-AA 는 타입이 990 개라 풀이 표 대신
  결합 파라미터 이름·전하·질량과 패턴 풀이로 설명한다.
- DREIDING _b1/_b2 의 뜻은 라벨 쌍 차수 표에서 바로 계산해 보여 준다(같은 혼성의 C 라벨 셋과의 차수). sp2 끼리는
  없음–없음·없음–b1·b2–b2 가 이중, b1–b1·b1–b2·없음–b2 가 단일이다.
- 실제 환경의 "가장 작은 고리"는 원자의 이웃 쌍 사이 최단 경로(그 원자 제외)+2 로 구하며, RDKit MinAtomRingSize 와
  401 분자에서 같다(일회성 대조). GAFF2·OPLS-AA 는 각자의 고리 판정(antechamber AR, foyer 8원 이하 chordless)을 따로 보인다.
- 화면 배치: 넓은 화면(1100px 이상)에서는 설명 칸을 3D 보기·타입 목록 오른쪽 칸에 두고 contain: size + sticky 로
  스크롤을 따라오게 했다(타입 카드를 눌러도, 3D 원자를 눌러도 보임). 좁은 화면에서는 타입 목록 아래에 놓고 3D 에서
  누를 때만 그쪽으로 굴린다. 원자 목록 표가 휴대폰 폭에서 페이지 전체를 넓히던 문제는 .main 열을 minmax(0, 1fr) 로
  막아 고쳤다(표 안에서만 가로 스크롤).
- 3D 보기는 4px 이하로 움직인 왼쪽 누름만 원자 고르기로 본다(회전·이동과 구분). 고른 원자는 점선 고리로 표시.
