# uniprot-ptm-lookup

A small reusable Python package for annotating in-line peptide modifications against the UniProt PTM vocabulary.

It supports modified peptide strings such as:

- `A(+27.99)EDPETQVVL`
- `A(+114.04)IPRSPFEVQVSPE`
- `AAAAVRQM(+15.99)NPHIRVT`
- `M(+15.99)PEPTIDEK(+114.04)`

The package ships with a bundled JSON export of the UniProt PTM list so it can be used immediately after installation, and it can also download the current UniProt flatfile when needed.

## Install

Editable install while developing locally:

```bash
python -m pip install -e .
```

Regular install into an environment:

```bash
python -m pip install .
```

## Quick start

```python
from uniprot_ptm_lookup import load_default_lookup

lib = load_default_lookup()

result = lib.annotate_modified_peptide("AAAAVRQM(+15.99)NPHIRVT")
site = result.modifications[0]

print(result.sequence)
print(site.position)
print(site.residue)
print(site.observed_delta)
print(site.assigned_modification_type)
print(site.best_candidate.display_name if site.best_candidate else None)
```

Expected broad assignments for the examples above:

- `A(+27.99)EDPETQVVL` -> `Methylation`
- `A(+114.04)IPRSPFEVQVSPE` -> `Ubiquitination`
- `AAAAVRQM(+15.99)NPHIRVT` -> `Oxidation`

## Batch annotation and CSV export

```python
from uniprot_ptm_lookup import load_default_lookup

peptides = [
    "A(+27.99)EDPETQVVL",
    "A(+114.04)IPRSPFEVQVSPE",
    "AAAAVRQM(+15.99)NPHIRVT",
]

lib = load_default_lookup()
results = lib.annotate_modified_peptides_to_csv(
    peptides,
    "annotated_ptms.csv",
    tolerance=0.05,
    max_candidates=5,
)
```

## Command line usage

Use the bundled PTM library that ships with the package:

```bash
python -m uniprot_ptm_lookup \
  --use-packaged-library \
  --modified-peptide "A(+27.99)EDPETQVVL" \
  --modified-peptide "AAAAVRQM(+15.99)NPHIRVT"
```

Or export a CSV from a text file containing one modified peptide per line:

```bash
uniprot-ptm-lookup \
  --use-packaged-library \
  --modified-peptide-file peptides.txt \
  --export-modified-csv annotated_ptms.csv
```

## Installing from Git in other projects

Once this directory is committed to a Git repository, another project can install it directly:

```bash
python -m pip install "git+https://YOUR_GIT_SERVER/YOUR_ORG/uniprot-ptm-lookup.git"
```

For an editable Git checkout during development:

```bash
git clone https://YOUR_GIT_SERVER/YOUR_ORG/uniprot-ptm-lookup.git
cd uniprot-ptm-lookup
python -m pip install -e .
```

## Notes

- The `+114.0429` ubiquitination / digly remnant mapping is implemented as a curated fallback alias, because this shorthand is commonly used in peptide search outputs but is not exposed as a direct site entry in the UniProt PTM list.
- For strict UniProt-only behavior, pass `allow_fallback_aliases=False` when annotating modified peptides.
- The bundled PTM JSON library is based on UniProt release `2026_01`.


Compatibility
-------------
- Python 3.7+
