# uniprot-ptm-lookup

A Python package for annotating inline peptide modifications against the UniProt PTM vocabulary.

It supports modified peptide strings such as:

- `A(+27.99)EDPETQVVL`
- `A(+114.04)IPRSPFEVQVSPE`
- `AAAAVRQM(+15.99)NPHIRVT`
- `M(+15.99)PEPTIDEK(+114.04)`

The package ships with a bundled JSON export of the UniProt PTM list so it can be used immediately after installation. It can also download the current UniProt flatfile when needed.

## Installation

Install from a published release asset:

```bash
python -m pip install "uniprot-ptm-lookup @ https://github.com/JasonPBennett/UniProt-PTM-Lookup/releases/download/v0.1.3/uniprot_ptm_lookup-0.1.3-py3-none-any.whl"
```

Install from a local checkout:

```bash
python -m pip install .
```

Install in editable mode while developing locally:

```bash
python -m pip install -e .
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

## CLI usage

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

## Using in other projects

Published release artifacts are attached to GitHub Releases. To pin this package in another project, add the release wheel URL to `requirements.txt`:

```text
uniprot-ptm-lookup @ https://github.com/JasonPBennett/UniProt-PTM-Lookup/releases/download/v0.1.3/uniprot_ptm_lookup-0.1.3-py3-none-any.whl
```

## Developer checkout from GitLab

For institutional development, clone the project from your GitLab namespace and install it in editable mode.

In GitLab, open the repository, select **Code**, copy either the HTTPS or SSH clone URL, and substitute it into the commands below.

HTTPS:

```bash
git clone https://<gitlab-host>/<group>/uniprot-ptm-lookup.git
cd uniprot-ptm-lookup
python -m pip install -e .
```

SSH:

```bash
git clone git@<gitlab-host>:<group>/uniprot-ptm-lookup.git
cd uniprot-ptm-lookup
python -m pip install -e .
```

## Notes

- The `+114.0429` ubiquitination / digly remnant mapping is implemented as a curated fallback alias, because this shorthand is commonly used in peptide search outputs but is not exposed as a direct site entry in the UniProt PTM list.
- For strict UniProt-only behavior, pass `allow_fallback_aliases=False` when annotating modified peptides.
- The bundled PTM JSON library is based on UniProt release `2026_01`.

## Compatibility

- Python 3.7+
