# 체크리스트

## 1단계 — 저장소 골격 (2026-11-04)
- [x] 라이선스 GPL-3.0 으로 변경 (antechamber 이식 때문, 사용자 승인)
- [x] DREIDING 타이퍼를 `fftyper` 패키지로 옮김, `fftyper --ff` 명령
- [x] 기존 DREIDING 테스트(Python 55개, 웹 7개) 통과

## 2단계 — GAFF2 (Python)
- [ ] AmberTools 26 설치(정답지), antechamber 대조 스크립트
- [ ] bondtype(-j part) 처리, 고리·방향족(AR1~AR5) 판정 이식
- [ ] ATOMTYPE_GFF2.DEF 해석기(원자 성질·화학 환경), cc/cd 등 공액 타입 조정
- [ ] 넓은 분자 집합에서 antechamber 와 원자 단위 일치
- [ ] gaff2.lt 에 없는 타입·결합 파라미터 확인, moltemplate 통합 테스트

## 3단계 — GAFF2 (JS) + HTML
- [ ] JS 이식, Python 대조
- [ ] 역장 선택 메뉴, 재빌드, 브라우저 테스트

## 4단계 — OPLS-AA
- [ ] foyer SMARTS 문법 엔진 (Python·JS 같은 결과)
- [ ] foyer 규칙 → oplsaa2024 번호 대응표 (설명 대조)
- [ ] 타입 전하, 타입 없는 원자 오류
- [ ] foyer·참조 분자 대조, moltemplate 통합 테스트

## 5단계 — 문서·배포
- [ ] README·README_kor
- [ ] 커밋·push
