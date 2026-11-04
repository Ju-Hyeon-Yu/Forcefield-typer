# Forcefield-typer

[English](README.md) | **한국어**

Forcefield-typer는 SMILES(필수)와 mol2(선택)로부터 **DREIDING**, **GAFF2**, **OPLS-AA** 원자 타입을 붙인 moltemplate
단일 분자 `.lt` 파일을 만듭니다. 출력은 `import "<역장>.lt"` + `MOL inherits <역장> { ... }` 형태라서 `system.lt`에
바로 불러올 수 있습니다. moltemplate에 들어 있는 역장 파일은 읽기만 하고 바꾸지 않습니다.

| 역장 | moltemplate 파일 | 타입을 정하는 기준 | 전하 |
|---|---|---|---|
| DREIDING | `dreiding.lt` | moltemplate DREIDING 라벨 매뉴얼 본문의 규칙 | mol2의 값 (RESP 등) |
| GAFF2 | `gaff2.lt` (GAFF2 2.11) | AmberTools `antechamber -at gaff2`와 같은 결과 (규칙 파일과 처리 과정을 옮김) | mol2의 값 |
| OPLS-AA | `oplsaa2024.lt` | foyer의 OPLS-AA SMARTS 규칙을 2024 타입 번호에 맞게 씀 | OPLS-AA 타입 전하 |

같은 타이퍼를 오프라인 HTML 파일 하나(`forcefield_typer.html`)로도 쓸 수 있습니다. DREIDING은
[Dreiding_typer](https://github.com/Ju-Hyeon-Yu/Dreiding_typer)에서 옮겨 왔습니다.

## 설치

```bash
conda activate myenv
git clone https://github.com/Ju-Hyeon-Yu/Forcefield-typer.git ~/Forcefield-typer
cd ~/Forcefield-typer
pip install -e .                 # rdkit, moltemplate 가 없으면 함께 설치합니다
```

그 환경에서만 `fftyper` 명령이 생깁니다. Python 3.9 이상이 필요합니다. 역장 파일은 실행할 때마다 설치된
moltemplate(2.22.5에서 확인)의 것을 읽습니다.

## 사용법

```bash
# GAFF2 + mol2 의 RESP 전하 — 원자 순서·좌표·전하는 mol2 에서, 결합 차수·형식 전하는 SMILES 에서 가져옵니다
fftyper --ff gaff2 --smiles "CNC(C)=O" --mol2 nma.mol2 -o nma.lt --name NMA

# OPLS-AA — 전하는 OPLS-AA 타입 전하라서 mol2 는 좌표가 필요할 때만 넣습니다
fftyper --ff oplsaa --smiles "CC(C)(C)O" -o tbuoh.lt

# DREIDING 이온쌍
fftyper --ff dreiding --smiles "Cn1cc[n+](C)c1.[O-]S(=O)(=O)C(F)(F)F" --mol2 ip.mol2 -o ip.lt
```

명령줄에서는 SMILES를 따옴표로 감싸야 합니다. `(`, `[`, `#`를 셸이 먼저 해석하기 때문입니다.

| 옵션 | 설명 |
|---|---|
| `--ff` | 필수. `dreiding`, `gaff2`, `oplsaa` |
| `--smiles` | 필수. 이온쌍은 `.`으로 연결합니다 (예: `A.B`) |
| `--mol2` | 좌표(DREIDING·GAFF2는 전하도)를 가져올 mol2. 없으면 RDKit으로 3D 구조를 만듭니다 |
| `-o`, `--out` | 출력 `.lt` 경로 |
| `--name` | moltemplate 객체 이름. 기본값은 출력 파일 이름입니다 |
| `--ff-file` | 타입과 규칙을 읽을 역장 `.lt`. 기본값은 설치된 moltemplate의 파일입니다 |
| `--seed` | mol2가 없을 때 3D 구조 생성 난수 시드 |

원자별 타입과 전하를 화면에 출력합니다. 역장 파일로 나타낼 수 없는 부분(파라미터가 없는 결합 항 등)은 표준
오류에 경고(`경고:`)로 알립니다. 타입을 정할 수 없으면 `.lt`를 쓰지 않고, 문제 원자를 알려 주고 멈춥니다.

결과는 `system.lt`에서 불러온 뒤 `moltemplate.sh system.lt`로 LAMMPS 입력을 만듭니다.

```
import "nma.lt"
m = new NMA
```

| 역장 | moltemplate 파일이 쓰는 LAMMPS 패키지 | 참고 |
|---|---|---|
| DREIDING | MOLECULE (`hbond/dreiding/lj`, `dihedral charmm`, `improper umbrella`), KSPACE | |
| GAFF2 | EXTRA-MOLECULE (`dihedral fourier`), MOLECULE (`improper cvff`), KSPACE | |
| OPLS-AA | MOLECULE (`dihedral opls`, `improper cvff`), KSPACE | `read_data` 다음에 `system.in.charges`를 읽어야 합니다 (moltemplate가 타입별 전하를 여기서 정함) |

## 오프라인 HTML 도구 (`forcefield_typer.html`)

`forcefield_typer.html`을 내려받아 Edge나 Chrome으로 엽니다. RDKit.js와 모든 규칙표가 들어 있어(약 4.5 MB)
인터넷 연결 없이 동작합니다.

1. 역장(DREIDING, GAFF2, OPLS-AA)을 고릅니다. 고른 역장은 브라우저에 기억됩니다.
2. SMILES를 넣고, 필요하면 좌표 파일(mol2, pdb, xyz)을 끌어다 놓습니다. 원자 순서는 자동으로 맞춥니다.
   mol2는 좌표와 전하를, pdb·xyz는 좌표만(전하 0) 씁니다. OPLS-AA는 항상 타입 전하를 씁니다. 파일이 없으면
   브라우저에서 3D 시작 구조를 만듭니다. MD 전에 에너지 최소화를 하십시오.
3. **타입 할당**을 누르면 3D 보기에 원자마다 타입이 표시되고, 타입 카드에 타입별 개수·원자 번호·설명이 나옵니다
   (누르면 해당 원자만 강조). `.lt` 내용은 복사하거나 내려받을 수 있습니다.
4. 타입 카드, 3D 보기의 원자, 원자 목록의 줄을 누르면 **타입 설명** 칸이 열립니다. 할당된 타입이 맞는지 확인할 때
   씁니다.
   - **뜻**: 역장 파일의 원문 설명(gaff2.lt, oplsaa2024.lt)과 한국어 풀이입니다. DREIDING은 라벨을 부분별로
     풀어 줍니다(혼성, `_b1`·`_b2` 결합 차수 플래그, `_dN`, `_hd`·`_ha`, 평형 결합각).
   - **판정 근거**: 그 타입을 고른 규칙입니다. GAFF2는 처음 맞은 ATOMTYPE_GFF2.DEF 규칙 줄과 그 풀이, 앞에서 맞지
     않은 같은 원소 규칙, 공액 짝(cc/cd 등)·c5/c6 조정을 보여 줍니다. OPLS-AA는 foyer SMARTS 패턴을 나무 모양으로
     풀어 주고, 우선 규칙과 oplsaa2024.lt 번호 대응, OPLS/2020 짝을 보여 줍니다. DREIDING은 타이퍼가 지난
     갈림길(이웃 수, 방향족, 공액, 수소결합 플래그)을 문장으로 보여 줍니다.
   - **원자별 실제 환경**: 원자마다 이웃(결합 차수·이웃의 타입), 붙은 H 수, 형식 전하, 가장 작은 고리 크기를
     보여 줍니다. 역장별로 DREIDING은 결합마다 SMILES 차수와 moltemplate 라벨 쌍 차수를 ✓/✗로 맞대 보고,
     GAFF2는 antechamber의 고리(AR1~AR5)·결합 형식 판정을, OPLS-AA는 맞은 foyer 규칙과 밀린 규칙, foyer 고리
     판정을 보여 줍니다. 이웃을 누르면 그 원자의 설명으로 넘어갑니다.
   - 같은 SMILES로 역장만 바꾸면 보던 원자의 새 타입 설명을 이어서 보여 줍니다. 원자 목록에는 이웃(타입) 칸이
     있습니다.

JS 타이퍼는 Python 타이퍼를 그대로 옮긴 것입니다. 아래의 기준 분자 집합에서 타입·경고·전하·오류가 같고, 원자
순서를 섞은 mol2 입력에서 `.lt` 전체가 Python CLI와 같은지 테스트로 확인합니다.

`web/src`를 고쳤거나 moltemplate의 역장 파일이 바뀌었을 때만 다시 빌드합니다.

```bash
cd web && npm install && cd ..        # RDKit.js (@rdkit/rdkit 2026.9.1)
python3 web/build.py                  # → forcefield_typer.html
cd web && npm test                    # Python 타이퍼와 대조 (Node)
node tests/test_browser.js            # 빌드한 HTML 을 브라우저로 확인 (playwright 필요)
```

## DREIDING

moltemplate DREIDING 라벨 매뉴얼([m-bone/moltemplate-DREIDING](https://github.com/m-bone/moltemplate-DREIDING))
본문의 규칙을 따릅니다. `R`은 방향족 고리 원자에만, `_b1`/`_b2`는 반드시 짝으로, `_d1`/`_d2`는 이웃 수 − 1,
수소결합은 `H_HB`와 `_hd`/`_ha`입니다. 결합 차수 플래그는 `dreiding.lt`의 규칙을 직접 읽어서, moltemplate가
붙이는 결합 차수가 SMILES와 같고 모든 결합에 이면각 규칙이 있도록 고릅니다. 전체 규칙표, 정확한 라벨이 없는
경우, 매뉴얼 그림과의 비교는 [Dreiding_typer 설명](https://github.com/Ju-Hyeon-Yu/Dreiding_typer/blob/main/README_kor.md)에
있습니다.

## GAFF2

- **antechamber와 같은 타입.** GAFF2 타이퍼는 AmberTools antechamber의 처리(`bondtype -j part` 다음
  `atomtype -p gaff2`, 규칙 파일은 AmberTools 26의 `ATOMTYPE_GFF2.DEF`)를 옮긴 것입니다. `antechamber -j 5`처럼
  SMILES의 결합 차수(케쿨레 구조)를 유지합니다. 391개 분자(원자 5279개)에서 antechamber와 타입이 모두 같습니다.
- **GAFF2 2.11.** moltemplate의 `gaff2.lt`는 `gaff211.dat`(GAFF2 2.11)를 amber2lt로 바꾼 것입니다. 2.2에서 생긴
  `c5`/`c6`(5·6원 고리 sp3 C)이 없으므로 이 원자는 2.11 규칙대로 `c3`로 둡니다. `gaff2.lt`에 아예 없는 타입(B, Si,
  금속 등)이 나오면 멈춥니다.
- **빠진 결합 파라미터.** 2.11판에서 `ns`, `nt`, `n7`–`n9`, `nx`–`nz`, `n+`, `nu`, `nv`, `n5`, `n6`, `cs`는 질량과 LJ만
  있고 결합 파라미터가 없습니다. AmberTools에서는 parmchk2가 "같은 타입"(PARMCHK.DAT의 EQUA, 예: `ns`→`n`)의 값을
  복사해 채웁니다. fftyper도 같은 방식으로 찾은 값을 생성 `.lt`의 분자 정의 앞에 `GAFF2` 객체에 덧붙여 씁니다.

  ```
  GAFF2 {
    write_once("In Settings") {
      bond_coeff @bond:c-nt harmonic 356.2 1.3790   # = c-n
      ...
    }
    write_once("Data Bonds By Type") {
      @bond:c-nt @atom:c @atom:nt
    }
    ...
  }
  ```

  찾는 순서는 parmchk2를 단계별로 따릅니다. 378개 분자에서 채운 항이 parmchk2가 쓰는 737개 항과 같습니다. 여러
  `.lt`가 같은 정의를 덧붙여도 됩니다. 근사값(parmchk2의 벌점이 붙은 CORR 타입, 경험식 값)은 넣지 않고, 남은
  누락은 경고로 알립니다. **결합**이 빠지면 moltemplate가 `No bond types defined`로 멈추고(예: 다이사이안아마이드의
  `cg-n2`), 각·이면각이 빠지면 moltemplate가 경고 없이 그 항을 뺍니다.
- **AMBER와 대조.** 같은 좌표로 374개 분자의 모든 결합 항을 tleap(`gaff211.dat` + parmchk2 frcmod)과
  비교했습니다. 결합 4950개는 모두 같습니다. 각·이면각은 fftyper가 누락을 경고하는 항, 그리고 `gaff2.dat` 2.11에
  방향만 반대로 값이 다르게 두 번 들어 있는 이면각 4쌍(35개 항)에서만 다릅니다(tleap·parmchk2는 앞의 것,
  moltemplate는 뒤의 것을 씀). parmchk2가 CORR 타입으로 더하는 improper는 `gaff2.lt`에 없어 넣지 않습니다.

## OPLS-AA

- **규칙: foyer.** 원자 타입은 [foyer](https://github.com/mosdef-hub/foyer) 1.2.0의 OPLS-AA SMARTS 규칙 228개로
  정합니다. foyer의 처리 방식(규칙 순서, 덮어쓰기, 유도 부분그래프 맞춤, 고리 찾기)을 그대로 다시 구현했고,
  foyer 검증 분자 182개와 다른 분자 391개에서 실패를 포함한 결과가 foyer를 직접 돌린 것과 같습니다.
- **번호: oplsaa2024.lt.** foyer는 GROMACS의 `opls_NNN` 이름을, moltemplate의 `oplsaa2024.lt`(Jorgensen 2024)는
  BOSS 번호를 씁니다. 대부분은 같은 번호가 같은 원자를 뜻하므로 `opls_135`는 `@atom:135`가 됩니다. 같은 번호가
  다른 원자를 뜻하거나(`732`, `927`, `930`) 2024 파일에 그 타입이 없으면(foyer의 `151a`–`153c`, `722a`) 문제 원자를
  알려 주고 멈춥니다. 그래서 알케인은 OPLS/2020 타입(54–60)이 아니라 고전 타입(135–140)을 씁니다.
- **OPLS/2020 짝.** 2024 파일은 159(3차 알코올 C)와 189(OCH2O 아세탈 C)의 전하를 바꾸면서 그 산소의 짝 타입을
  새로 두었습니다. 159 옆의 154 O는 160으로, 189 옆의 180 O는 186으로 바꿔서 이 분자들이 중성이 되게 합니다.
- **전하.** `oplsaa2024.lt`의 타입 전하를 씁니다(moltemplate도 실행할 때 타입별로 정함). foyer 규칙은 OPLS-AA의
  분자별 전하 조정(바이페닐·큐멘의 ipso C, 아세테이트의 CH3 등)을 모르므로, 이런 분자는 전하 합이 형식 전하와 달라
  경고가 나옵니다.
- **안 되는 것.** 이미다졸륨, 트리플레이트, TFSI, 술포네이트 같은 이온은 foyer 규칙도 2024 파일의 타입도 없어서
  OPLS-AA로는 멈춥니다.

## mol2 조건

- 원소는 원자 타입 열(SYBYL `C.ar`, GAFF `ca` 등)에서 읽고, 안 되면 원자 이름에서 읽습니다.
- 전하는 9번째 열에서 읽습니다(없으면 0으로 두고 경고).
- **원자 순서가 같지 않아도 됩니다.** SMILES와 mol2 원자를 원소와 연결 관계로 맞추고, 출력은 mol2 순서를
  따릅니다. BOND 절은 연결 관계만 쓰고 결합 차수 열은 쓰지 않습니다.
- DREIDING·GAFF2는 mol2 전하 합이 SMILES 형식 전하와 다르면 경고합니다.

SMILES가 없으면 좌표에서 결합 차수를 먼저 추정해(예: RDKit `rdDetermineBonds.DetermineBonds`) SMILES를 얻고, 눈으로
확인한 뒤 `--smiles`로 넣으십시오. 공명 구조가 여럿인 이온(예: TFSI⁻)은 다른 형태로 나올 수 있습니다.

## 파일

| 파일 | 내용 |
|---|---|
| `fftyper/cli.py` | `fftyper` 명령 |
| `fftyper/molio.py` | mol2 읽기, SMILES↔mol2 원자 대응, 3D 구조, `.lt` 쓰기 |
| `fftyper/dreiding.py` | DREIDING 타이퍼 |
| `fftyper/gaff2.py` | GAFF2 타이퍼(antechamber 이식)와 parmchk2 EQUA 보충 |
| `fftyper/foyer_smarts.py`, `fftyper/oplsaa.py` | foyer 규칙 판정 엔진과 OPLS-AA 타이퍼 |
| `fftyper/ltrules.py` | moltemplate 역장 파일의 타입·By Type 규칙·계수·전하 읽기 |
| `fftyper/data/` | `ATOMTYPE_GFF2.DEF`(AmberTools), `foyer_oplsaa.json`(foyer 규칙)과 foyer 라이선스 |
| `web/` | HTML 도구 소스(`src/`), 빌드 스크립트, 테스트 |
| `tests/` | 테스트, 정답 데이터와 그것을 만드는 스크립트(`antechamber_ref.py`, `parmchk2_ref.py`, `foyer_ref.py`, `foyer_rules.py`) |
| `plan.md`, `checklist.md`, `context-notes.md` | 계획, 할 일 목록, 결정과 그 근거 기록 |

## 테스트

```bash
pip install -e ".[test]"
python3 -m pytest -q tests          # 약 1분 30초, moltemplate·LAMMPS 실행 포함
```

테스트는 저장해 둔 정답 데이터와 비교하므로 AmberTools와 foyer가 없어도 됩니다. 정답 데이터는 AmberTools 26
(`tests/antechamber_ref.py`, `tests/parmchk2_ref.py`)과 foyer 1.2.0(`tests/foyer_ref.py`)으로 만들었습니다. 통합
테스트는 역장마다 moltemplate로 계를 만들어 모든 결합·각·이면각에 타입이 붙는지 보고, LAMMPS 단일점 계산을
돌립니다(`lmp`가 없으면 건너뜀).

## 라이선스

GPL-3.0-or-later (`LICENSE`). GAFF2 타이퍼가 AmberTools antechamber(GPL-2.0-or-later)를 옮긴 것이기 때문입니다.
포함한 외부 자료는 다음과 같습니다.

- `fftyper/data/ATOMTYPE_GFF2.DEF` — AmberTools 26, GPL-2.0-or-later
- `fftyper/data/foyer_oplsaa.json` — foyer `oplsaa.xml`의 규칙, MIT, © 2015 Vanderbilt University
  (`fftyper/data/LICENSE.foyer.txt`)
- `forcefield_typer.html` 안의 RDKit.js — BSD-3-Clause, © RDKit 기여자

## 참고 문헌

- Mayo, Olafson, Goddard, *J. Phys. Chem.* **94**, 8897 (1990) — DREIDING
- Wang, Wang, Kollman, Case, *J. Mol. Graph. Model.* **25**, 247 (2006) — antechamber 원자 타입 판정
  (GAFF2는 AmberTools에 들어 있음)
- Jorgensen et al., *J. Phys. Chem. B* (2024), doi:10.1021/acs.jpcb.3c06602 — `oplsaa2024.lt`의 OPLS-AA 파라미터
- Klein et al., *Comput. Mater. Sci.* **167** (2019), doi:10.1016/j.commatsci.2019.05.026 — foyer
- [jewettaij/moltemplate](https://github.com/jewettaij/moltemplate) — 역장 파일
- [RDKit.js](https://github.com/rdkit/rdkit-js) — HTML 도구의 화학 엔진
