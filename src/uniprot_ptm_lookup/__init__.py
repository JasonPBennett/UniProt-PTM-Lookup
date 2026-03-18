"""Reusable UniProt PTM lookup package.

Typical usage:

    from uniprot_ptm_lookup import load_default_lookup

    lib = load_default_lookup()
    result = lib.annotate_modified_peptide("AAAAVRQM(+15.99)NPHIRVT")
"""

from .core import (
    DEFAULT_MODIFIED_PEPTIDE_TOLERANCE,
    DEFAULT_PTM_URL,
    MassAlias,
    ModifiedPeptideParseResult,
    ModifiedPeptideResult,
    ModifiedSiteResult,
    PTMCandidate,
    PTMEntry,
    PTMHit,
    SiteAnnotation,
    SiteObservation,
    UniProtPTMLookup,
    download_ptm_list,
    infer_modification_type,
    load_default_lookup,
    parse_modified_peptide_string,
)

__all__ = [
    "DEFAULT_MODIFIED_PEPTIDE_TOLERANCE",
    "DEFAULT_PTM_URL",
    "MassAlias",
    "ModifiedPeptideParseResult",
    "ModifiedPeptideResult",
    "ModifiedSiteResult",
    "PTMCandidate",
    "PTMEntry",
    "PTMHit",
    "SiteAnnotation",
    "SiteObservation",
    "UniProtPTMLookup",
    "download_ptm_list",
    "infer_modification_type",
    "load_default_lookup",
    "parse_modified_peptide_string",
]

__version__ = "0.1.0"
