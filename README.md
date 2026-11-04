# Forcefield-typer

Builds a single-molecule moltemplate `.lt` file with DREIDING, GAFF2 or OPLS-AA atom types from a SMILES string,
optionally together with a mol2 file that carries coordinates and charges.

Work in progress: DREIDING is available (moved from
[Dreiding_typer](https://github.com/Ju-Hyeon-Yu/Dreiding_typer)); GAFF2 and OPLS-AA are being added.

```bash
pip install -e .
fftyper --ff dreiding --smiles "CCO" -o ethanol.lt
```

License: GPL-3.0-or-later (the GAFF2 typer reproduces AmberTools antechamber rules, which are GPL).
