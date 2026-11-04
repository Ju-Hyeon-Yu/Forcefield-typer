# 체크리스트

## 1단계 — 저장소 골격 (2026-11-04)
- [x] 라이선스 GPL-3.0 으로 변경 (antechamber 이식 때문, 사용자 승인)
- [x] DREIDING 타이퍼를 `fftyper` 패키지로 옮김, `fftyper --ff` 명령
- [x] 기존 DREIDING 테스트(Python 55개, 웹 7개) 통과

## 2단계 — GAFF2 (Python)
- [x] AmberTools 26 설치(정답지), antechamber 대조 스크립트
- [x] bondtype(-j part) 처리, 고리·방향족(AR1~AR5) 판정 이식
- [x] ATOMTYPE_GFF2.DEF 해석기(원자 성질·화학 환경), cc/cd 등 공액 타입 조정
- [x] 넓은 분자 집합에서 antechamber 와 원자 단위 일치 (391 분자, 5279 원자)
- [x] gaff2.lt 에 없는 타입 → 오류, c5·c6 → c3 (2.11 판)
- [x] 결합 파라미터 없는 타입(ns, nt, nx …): parmchk2 EQUA 값 복사 블록을 .lt 에 함께 씀 (사용자 선택 1안)
- [x] parmchk2 정답지 대조 (378 분자, 737 항), moltemplate·LAMMPS 통합 테스트
- [x] AMBER(tleap+sander)와 항별 대조 — 차이는 경고된 누락·gaff2.lt 중복 이면각·CORR improper 뿐

## 3단계 — GAFF2 (JS) + HTML
- [x] JS 이식(web/src/gaff2.js, ltrules.js), Python 대조 — 391 분자 타입·경고·보충 블록 모두 같음
- [x] 케쿨레 구조를 SMILES 순서로 고정(stamp_kekule) — mol2 입력에서도 Python·웹 결과가 같게
- [x] 역장 선택 메뉴(브라우저에 기억), 3D 시작 구조는 DREIDING 표 + 라벨 못 붙는 분자용 대체, 5·6배위 리간드 반발
- [x] 재빌드, Node 10개·브라우저 9개 테스트

## 4단계 — OPLS-AA
- [x] foyer SMARTS 판정 엔진 (Python) — 실제 foyer 와 검증 분자 182·SMILES 391 모두 같음
- [x] JS 이식(web/src/foyer.js) — 검증 분자 182 개 foyer 와 같음, 391 분자 타입·경고·전하가 Python 과 같음
- [x] 웹: 역장 선택에 OPLS-AA, 타입 전하, 타입 카드에 설명 표시, Node 13개·브라우저 10개 테스트
- [x] foyer 규칙 → oplsaa2024 번호 대응 (설명·매개변수 대조, 뜻이 다른 3 개·없는 12 개는 오류), OPLS/2020 짝 2 개
- [x] 타입 전하, 타입 없는 원자 오류, 전하 합 경고
- [x] moltemplate·LAMMPS 통합 테스트

## 5단계 — 문서·배포
- [x] README·README_kor (역장별 기준·검증 결과·라이선스)
- [x] 깨끗한 가상환경에 pip 설치 후 세 역장 모두 fftyper 실행 확인
- [x] 커밋·push

## 6단계 — 타입 설명 칸 (HTML, 2026-11-05)
- [x] 사용자 선택: 뜻 + 판정 근거 + 실제 환경, 원문 + 한국어 풀이, HTML 에만
- [x] 타이퍼가 판정 정보를 내보냄 — GAFF2 규칙 번호·조정 전 타입·고리/결합 형식, foyer whitelist/blacklist·고리
- [x] web/src/explain.js — DREIDING 라벨 풀이·갈림길, GAFF2 DEF 규칙 풀이·공액 짝, foyer 패턴 나무·우선 규칙·OPLS/2020 짝
- [x] 화면: 타입 카드·3D 원자·원자 목록 줄을 누르면 설명 칸(넓은 화면은 3D 보기 옆에 붙어 따라옴), 이웃으로 이동,
      원자 목록 이웃(타입) 칸, 마우스 정보에 타입 설명, 같은 SMILES 로 역장을 바꾸면 같은 원자를 이어서 설명
- [x] 검증: 세 역장 분자 500여 개에서 설명 생성·원자 누락 없음, DREIDING 이유 = 타이퍼 갈림길(6820 원자),
      GAFF2 근거 = 실제 첫 규칙(5111 원자), 고리 크기 = RDKit MinAtomRingSize(401 분자)
- [x] 테스트: Python 97개, Node 16개, 브라우저 11개 통과 → 커밋·push
