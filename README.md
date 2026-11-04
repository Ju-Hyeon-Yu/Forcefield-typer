# Forcefield-typer

**English** | [한국어](README_kor.md)

Forcefield-typer builds a single-molecule moltemplate `.lt` file with **DREIDING**, **GAFF2** or **OPLS-AA** atom
types from a SMILES string, optionally together with a mol2 file that carries coordinates and charges. The output
has the form `import "<force field>.lt"` + `MOL inherits <FF> { ... }`, so it can be pulled straight into a
`system.lt`. The force-field files shipped with moltemplate are only read, never changed.

| Force field | moltemplate file | How atom types are chosen | Charges |
|---|---|---|---|
| DREIDING | `dreiding.lt` | Rules in the text of the moltemplate DREIDING label manual | from the mol2 (e.g. RESP) |
| GAFF2 | `gaff2.lt` (GAFF2 2.11) | Same result as AmberTools `antechamber -at gaff2` (its rule file and processing are ported) | from the mol2 |
| OPLS-AA | `oplsaa2024.lt` | foyer's OPLS-AA SMARTS rules, mapped to the 2024 type numbers | OPLS-AA type charges |

The same typers are also available as one offline HTML file, `forcefield_typer.html`. DREIDING was moved here from
[Dreiding_typer](https://github.com/Ju-Hyeon-Yu/Dreiding_typer). Program messages and the HTML interface are in
Korean.

## Installation

```bash
conda activate myenv
git clone https://github.com/Ju-Hyeon-Yu/Forcefield-typer.git ~/Forcefield-typer
cd ~/Forcefield-typer
pip install -e .                 # also installs rdkit and moltemplate if missing
```

This creates the `fftyper` command in that environment only. Python 3.9 or later is required. The force-field
files are read from the installed moltemplate (tested with 2.22.5) on every run.

## Usage

```bash
# GAFF2 with RESP charges from a mol2 — atom order, coordinates and charges come from the mol2,
# bond orders and formal charges come from the SMILES
fftyper --ff gaff2 --smiles "CNC(C)=O" --mol2 nma.mol2 -o nma.lt --name NMA

# OPLS-AA — charges are the OPLS-AA type charges, so a mol2 is only needed for coordinates
fftyper --ff oplsaa --smiles "CC(C)(C)O" -o tbuoh.lt

# DREIDING ion pair
fftyper --ff dreiding --smiles "Cn1cc[n+](C)c1.[O-]S(=O)(=O)C(F)(F)F" --mol2 ip.mol2 -o ip.lt
```

Quote the SMILES on the command line, because the shell interprets `(`, `[` and `#`.

| Option | Description |
|---|---|
| `--ff` | Required. `dreiding`, `gaff2` or `oplsaa` |
| `--smiles` | Required. Join the ions of an ion pair with `.` (e.g. `A.B`) |
| `--mol2` | mol2 to take coordinates (and, for DREIDING/GAFF2, charges) from. Without it, RDKit builds a 3D structure |
| `-o`, `--out` | Output `.lt` path |
| `--name` | moltemplate object name. Defaults to the output file name |
| `--ff-file` | Force-field `.lt` to read types and rules from. Defaults to the installed moltemplate's file |
| `--seed` | Random seed for 3D structure generation when there is no mol2 |

The per-atom types and charges are printed. Anything the force field file cannot represent is reported on standard
error as a warning (`경고:`), for example bonded terms that have no parameters. If a type cannot be determined, no
`.lt` is written and the program stops with the atoms concerned.

Load the result from `system.lt` and run `moltemplate.sh system.lt`:

```
import "nma.lt"
m = new NMA
```

| Force field | LAMMPS packages needed by the moltemplate file | Note |
|---|---|---|
| DREIDING | MOLECULE (`hbond/dreiding/lj`, `dihedral charmm`, `improper umbrella`), KSPACE | |
| GAFF2 | EXTRA-MOLECULE (`dihedral fourier`), MOLECULE (`improper cvff`), KSPACE | |
| OPLS-AA | MOLECULE (`dihedral opls`, `improper cvff`), KSPACE | include `system.in.charges` after `read_data` (moltemplate sets charges by type there) |

## Offline HTML tool (`forcefield_typer.html`)

Download `forcefield_typer.html` and open it in Edge or Chrome. It works without an internet connection (RDKit.js
and all rule tables are embedded, about 4.5 MB).

1. Choose the force field (DREIDING, GAFF2, OPLS-AA). The choice is remembered in the browser.
2. Enter a SMILES and optionally drop a coordinate file (mol2, pdb, xyz). Atom order is matched automatically.
   mol2 supplies coordinates and charges; pdb and xyz supply coordinates only (charges 0). OPLS-AA always uses type
   charges. Without a file, a 3D starting structure is built in the browser — minimize it before MD.
3. Press **타입 할당**. The 3D viewer labels every atom, the type cards show the count, atoms and description of each
   type (click to highlight), and the `.lt` text can be copied or downloaded.

The JavaScript typers are ports of the Python ones; tests confirm identical types, warnings, charges and errors for
the reference molecule sets below, and identical `.lt` files for mol2 inputs with shuffled atom order.

Rebuild after changing `web/src` or when moltemplate's force-field files change:

```bash
cd web && npm install && cd ..        # RDKit.js (@rdkit/rdkit 2026.9.1)
python3 web/build.py                  # → forcefield_typer.html
cd web && npm test                    # compare against the Python typers (Node)
node tests/test_browser.js            # check the built HTML in a browser (needs playwright)
```

## DREIDING

The rules follow the text of the moltemplate DREIDING label manual
([m-bone/moltemplate-DREIDING](https://github.com/m-bone/moltemplate-DREIDING)): `R` only for atoms in aromatic
rings, `_b1`/`_b2` only in pairs, `_d1`/`_d2` for the number of neighbors minus one, `H_HB` with `_hd`/`_ha` for
hydrogen bonds. The bond-order flags are chosen by reading the rules in `dreiding.lt`, so that the bond order
moltemplate assigns matches the SMILES and every bond has a dihedral rule. The full rule table, the cases without an
exact label, and the comparison with the manual's figures are in the
[Dreiding_typer README](https://github.com/Ju-Hyeon-Yu/Dreiding_typer#typing-rules).

## GAFF2

- **Same types as antechamber.** The GAFF2 typer is a port of AmberTools antechamber (`bondtype -j part` followed by
  `atomtype -p gaff2`, with `ATOMTYPE_GFF2.DEF` from AmberTools 26). The bond orders of the SMILES (Kekulé form) are
  kept, as with `antechamber -j 5`. For 391 molecules (5279 atoms) the types are identical to antechamber's.
- **GAFF2 2.11.** moltemplate's `gaff2.lt` is `gaff211.dat` (GAFF2 2.11) converted with amber2lt. It has no `c5`/`c6`
  (sp3 C in 5/6-membered rings, added in 2.2), so these atoms become `c3`, which is what the 2.11 rules give. A type
  that is not in `gaff2.lt` at all (B, Si, metals, ...) stops the program.
- **Missing bonded parameters.** In 2.11, `ns`, `nt`, `n7`–`n9`, `nx`–`nz`, `n+`, `nu`, `nv`, `n5`, `n6` and `cs` have
  masses and LJ parameters but no bonded parameters; AmberTools fills them with parmchk2 using "equivalent" types
  (EQUA in PARMCHK.DAT, e.g. `ns`→`n`). fftyper does the same and writes the copied parameters into the generated
  `.lt`, before the molecule, by adding to the `GAFF2` object:

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

  The search follows parmchk2 step by step; for 378 molecules the filled terms are the same 737 terms parmchk2
  writes. Several `.lt` files may add the same definitions. Approximate parameters (parmchk2's CORR types with a
  penalty, empirical values) are not added; the remaining missing terms are reported. A missing **bond** makes
  moltemplate stop (`No bond types defined`, e.g. `cg-n2` in dicyanamide); missing angles and dihedrals are silently
  dropped by moltemplate.
- **Check against AMBER.** With the same coordinates, every bonded term was compared with tleap (`gaff211.dat` +
  parmchk2 frcmod) for 374 molecules: all 4950 bonds are identical; angles and dihedrals differ only where fftyper
  warns about a missing term, and for 35 dihedrals of four pairs that `gaff2.dat` 2.11 defines twice in opposite
  directions with different values (tleap and parmchk2 use the first, moltemplate the later one). Impropers that
  parmchk2 adds through CORR types are not in `gaff2.lt` and are not added.

## OPLS-AA

- **Rules: foyer.** Atom types are decided with the 228 OPLS-AA SMARTS rules of
  [foyer](https://github.com/mosdef-hub/foyer) 1.2.0, by a re-implementation of foyer's procedure (rule order,
  overrides, induced-subgraph matching, its ring finder). For foyer's 182 validation molecules and 391 other
  molecules, the result — including failures — is identical to running foyer itself.
- **Numbers: oplsaa2024.lt.** foyer uses GROMACS `opls_NNN` names; moltemplate's `oplsaa2024.lt` (Jorgensen 2024)
  uses BOSS numbers. Most numbers mean the same atom, so `opls_135` becomes `@atom:135`. Where the same number means a
  different atom (`732`, `927`, `930`) or the 2024 file has no such type (foyer's `151a`–`153c`, `722a`), the program
  stops with the atoms concerned. Alkanes therefore use the classic types 135–140, not the OPLS/2020 types 54–60.
- **OPLS/2020 partners.** The 2024 file changed the charges of 159 (C of tertiary alcohols) and 189 (C of OCH2O
  acetals) together with new partner types for their oxygens; an O typed 154 next to 159 becomes 160, and an O typed
  180 next to 189 becomes 186, so that these molecules stay neutral.
- **Charges.** The type charges of `oplsaa2024.lt` are written (moltemplate also sets them by type). foyer's rules do
  not know OPLS-AA's per-molecule charge adjustments (for example the ipso carbons of biphenyl or cumene, or the CH3 of
  acetate), so for such molecules the charge sum differs from the formal charge and a warning is printed.
- **Not covered.** Ions such as imidazolium, triflate, TFSI and sulfonates have no foyer rules and no types in the 2024
  file, so OPLS-AA typing stops for them.

## mol2 requirements

- Elements are read from the atom-type column (SYBYL `C.ar`, GAFF `ca`, ...), falling back to the atom name.
- Charges are read from the 9th column (missing → 0 with a warning).
- **Atom order does not need to match.** SMILES and mol2 atoms are matched by element and connectivity, and the output
  keeps the mol2 order. Only the connectivity of the BOND section is used, not its bond orders.
- For DREIDING and GAFF2, a warning is printed if the mol2 charges do not add up to the formal charge of the SMILES.

Without a SMILES, infer bond orders from the coordinates first, for example with RDKit
`rdDetermineBonds.DetermineBonds`, check the SMILES by eye, and pass it with `--smiles`. Ions with several resonance
structures (e.g. TFSI⁻) may come out in another form.

## Files

| File | Contents |
|---|---|
| `fftyper/cli.py` | `fftyper` command |
| `fftyper/molio.py` | mol2 reading, SMILES↔mol2 atom matching, 3D structure, `.lt` writing |
| `fftyper/dreiding.py` | DREIDING typer |
| `fftyper/gaff2.py` | GAFF2 typer (antechamber port) and parmchk2 EQUA filling |
| `fftyper/foyer_smarts.py`, `fftyper/oplsaa.py` | foyer rule engine and OPLS-AA typer |
| `fftyper/ltrules.py` | Reading types, By-Type rules, coefficients and charges of moltemplate force-field files |
| `fftyper/data/` | `ATOMTYPE_GFF2.DEF` (AmberTools), `foyer_oplsaa.json` (foyer rules) and foyer's license |
| `web/` | HTML tool sources (`src/`), build scripts and tests |
| `tests/` | Tests, reference data and the scripts that make it (`antechamber_ref.py`, `parmchk2_ref.py`, `foyer_ref.py`, `foyer_rules.py`) |
| `plan.md`, `checklist.md`, `context-notes.md` | Plan, task list and a log of decisions with their reasons (Korean) |

## Tests

```bash
pip install -e ".[test]"
python3 -m pytest -q tests          # about 1.5 minutes, includes moltemplate and LAMMPS runs
```

The tests compare against stored reference data, so AmberTools and foyer are not needed to run them. The reference
data were made with AmberTools 26 (`tests/antechamber_ref.py`, `tests/parmchk2_ref.py`) and foyer 1.2.0
(`tests/foyer_ref.py`). Integration tests build systems with moltemplate for each force field, check that every bond,
angle and dihedral gets a type, and run a LAMMPS single point (skipped if `lmp` is missing).

## License

GPL-3.0-or-later (see `LICENSE`), because the GAFF2 typer ports AmberTools antechamber (GPL-2.0-or-later). Included
third-party material:

- `fftyper/data/ATOMTYPE_GFF2.DEF` — AmberTools 26, GPL-2.0-or-later
- `fftyper/data/foyer_oplsaa.json` — rules from foyer's `oplsaa.xml`, MIT, © 2015 Vanderbilt University
  (`fftyper/data/LICENSE.foyer.txt`)
- RDKit.js inside `forcefield_typer.html` — BSD-3-Clause, © RDKit contributors

## References

- Mayo, Olafson, Goddard, *J. Phys. Chem.* **94**, 8897 (1990) — DREIDING
- Wang, Wang, Kollman, Case, *J. Mol. Graph. Model.* **25**, 247 (2006) — antechamber atom-type perception
  (GAFF2 itself is distributed with AmberTools)
- Jorgensen et al., *J. Phys. Chem. B* (2024), doi:10.1021/acs.jpcb.3c06602 — OPLS-AA parameters in `oplsaa2024.lt`
- Klein et al., *Comput. Mater. Sci.* **167** (2019), doi:10.1016/j.commatsci.2019.05.026 — foyer
- [jewettaij/moltemplate](https://github.com/jewettaij/moltemplate) — force-field files
- [RDKit.js](https://github.com/rdkit/rdkit-js) — chemistry engine of the HTML tool
