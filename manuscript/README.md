# Permutation norm dependence

Haobo Ma (ChronoAI Pte Ltd, Singapore) and Wenlin Zhang (National University of Singapore; ChronoAI Pte Ltd, Singapore).

## Contents

- `note.tex`: the main source, with bibliography included.
- `note.pdf`: the typeset note, distributed separately.
- `metadata.json`: title, creators, problem identifier, license and mathematical source link.
- `verify.py`: dependency-free exact checks using Python 3.10 or later.
- `checks.json`: scientific output of the exact checks, distributed separately.

The note extends Clarisse and Wocjan's permutation classification: for r >= 1 and d >= 2 it proves the four-symmetry norm classification and makes their orbit upper bound exact. The same classification holds for rank and Schatten functionals for 0 < p <= infinity, except p = 2; those with p < 1 are quasi-norms. Adaptive identification uses at most 2r - 1 states. Equality on positive-definite density matrices suffices for Schatten equality on all density matrices. Bell states also prove the known pairwise incomparability of the 2^(r-1) - 1 nontrivial partial-transpose cuts. General detection-set incomparability is distinct from norm equality.

## Rebuild

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, use a standard TeX installation:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is embedded; no BibTeX step, external figures, or additional TeX input files are required.

## Exact checks

```sh
python3 verify.py --output checks.json
```

Expected output has `all_checks_satisfied` equal to `true`. Integer Gram identities determine the complete singular-value lists for all 4 mixed-state label placements and all 16 Bell-state label placements. The program checks 5,461 arbitrary paired-label families through 6 parties, every balanced subset through 8 parties, and the row-set convention for every permutation through 3 parties. It uses no floating-point arithmetic.

Adaptive and rank signatures are checked on all balanced subsets through 8 parties. Exact base-two exponents compare the witness Schatten functionals at p = 1/2, 1, 3/2, 2, 4 and infinity through 6 parties; p = 2 gives a single class. Every pair of nontrivial partial-transpose cuts through 8 parties has Bell witnesses in both directions. These finite ranges accompany the all-parameter written proofs; they are not an extrapolation.

For an uncorrelated spectator left as a separate tensor factor, each Schatten functional is multiplied by that of the spectator. In particular, its trace norm is unchanged for any density-matrix spectator; every Schatten functional is unchanged for a pure spectator.

The class counts including the identity for 1 through 8 parties are 1, 3, 7, 23, 71, 252, 890, 3299. These are finite diagnostics; the general proofs are in the note.

The deterministic `checks.json` output has SHA256 `e9977cc777865721b1c955f438d11e59547d1f44cd5eaf56f011891692b7f37a`.

## Mathematical source

The declaration `ClarisseWocjanClassification.result` at the following pinned source states the trace-norm equivalence for positive trace-one states, r >= 1 and d >= 2:

https://github.com/the-omega-institute/trureturing/blob/3464bc4487b2236dab97aaf7bbed237f0db28b71/D5/S3/Quantum/Entanglement/PermutationCriteria/ClarisseWocjanClassification.lean

The class count and Schatten extension have written proofs in the note. The exact checks cover the finite ranges stated above.

## License

CC BY 4.0 (Creative Commons Attribution 4.0 International).
https://creativecommons.org/licenses/by/4.0/
