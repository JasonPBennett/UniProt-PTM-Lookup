from __future__ import annotations

"""
UniProt PTM lookup utility.

This module downloads and parses the UniProt PTM flatfile and builds in-memory
indexes that are useful for proteomics-style annotation workflows.

What it gives you
-----------------
- Parsed PTM entries with names, accessions, target residue(s), mass deltas,
  keywords, and external cross-references.
- Fast exact lookup keyed by (residue, quantized mass delta).
- Fast tolerant lookup keyed by residue plus a bisected sorted mass index.
- Separate crosslink lookup keyed by the participating residue signature.
- High-level parsing of modified peptide strings such as:
      A(+27.99)EDPETQVVL
      A(+114.04)IPRSPFEVQVSPE
- Batch annotation helpers that can be called one peptide at a time or on a
  whole list of peptide strings.

No third-party packages are required.

Low-level example
-----------------
>>> lib = UniProtPTMLookup.from_uniprot(cache_path="ptmlist.txt")
>>> hits = lib.lookup_site("S", 79.966331, tolerance=0.01)
>>> [hit.entry.name for hit in hits][:3]
['Phosphoserine']

High-level example
------------------
>>> lib = UniProtPTMLookup.from_uniprot(cache_path="ptmlist.txt")
>>> result = lib.annotate_modified_peptide("A(+27.99)EDPETQVVL", tolerance=0.05)
>>> result.modifications[0].best_candidate.modification_type
'Methylation'
"""

import argparse
import bisect
import csv
import json
import re
import importlib.resources as importlib_resources
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import DefaultDict, Dict, FrozenSet, Iterable, List, Optional, Sequence, Tuple, Union
from urllib.request import urlopen

DEFAULT_PTM_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/complete/docs/ptmlist"
)

PACKAGED_LIBRARY_FILENAME = "uniprot_ptm_entries_2026_01.json"

LINE_CODES = frozenset({"ID", "AC", "FT", "TG", "PA", "PP", "CF", "MM", "MA", "LC", "TR", "KW", "DR"})
FIELD_RE = re.compile(r"^(ID|AC|FT|TG|PA|PP|CF|MM|MA|LC|TR|KW|DR)\s{2,}(.*)$")
ENTRY_SEPARATOR_PREFIX = "_____"

AA_NAME_TO_CODE: Dict[str, str] = {
    "Alanine": "A",
    "Arginine": "R",
    "Asparagine": "N",
    "Aspartate": "D",
    "Cysteine": "C",
    "Glutamine": "Q",
    "Glutamate": "E",
    "Glycine": "G",
    "Histidine": "H",
    "Isoleucine": "I",
    "Leucine": "L",
    "Lysine": "K",
    "Methionine": "M",
    "Phenylalanine": "F",
    "Proline": "P",
    "Serine": "S",
    "Threonine": "T",
    "Tryptophan": "W",
    "Tyrosine": "Y",
    "Valine": "V",
    "Selenocysteine": "U",
    "Pyrrolysine": "O",
}

AA_THREE_TO_CODE: Dict[str, str] = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Sec": "U",
    "Pyl": "O",
}

AA_CODE_TO_NAME: Dict[str, str] = {code: name for name, code in AA_NAME_TO_CODE.items()}
RESIDUE_ALIASES: Dict[str, str] = {}
for _name, _code in AA_NAME_TO_CODE.items():
    RESIDUE_ALIASES[_name.lower()] = _code
for _name, _code in AA_THREE_TO_CODE.items():
    RESIDUE_ALIASES[_name.lower()] = _code
for _code in AA_CODE_TO_NAME:
    RESIDUE_ALIASES[_code.lower()] = _code

DEFAULT_SITE_FEATURE_KEYS = frozenset({"MOD_RES", "CARBOHYD", "LIPID"})
DEFAULT_CROSSLINK_FEATURE_KEYS = frozenset({"CROSSLNK"})
DEFAULT_MODIFIED_PEPTIDE_TOLERANCE = 0.05

# Parsed peptide strings frequently use rounded masses. This small alias table
# is intentionally conservative and only covers common search-engine shorthand
# that is not represented directly in the UniProt PTM list in a convenient form.
# The canonical UniProt entries remain the primary source of truth whenever
# a match exists there.
DEFAULT_FALLBACK_MASS_ALIASES: Tuple["MassAlias", ...] = ()


@dataclass
class PTMEntry:
    accession: str
    name: str
    feature_key: str
    raw_target: str
    target_groups: Tuple[Tuple[str, ...], ...]
    position_amino_acid: Optional[str] = None
    position_polypeptide: Optional[str] = None
    correction_formula: Optional[str] = None
    mono_mass_delta: Optional[float] = None
    avg_mass_delta: Optional[float] = None
    cellular_location: Optional[str] = None
    taxonomic_range: Tuple[str, ...] = ()
    keywords: Tuple[str, ...] = ()
    xrefs: Dict[str, Tuple[str, ...]] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, payload: Dict[str, object]) -> "PTMEntry":
        return cls(
            accession=str(payload["accession"]),
            name=str(payload["name"]),
            feature_key=str(payload["feature_key"]),
            raw_target=str(payload["raw_target"]),
            target_groups=tuple(tuple(group) for group in payload["target_groups"]),
            position_amino_acid=payload.get("position_amino_acid"),
            position_polypeptide=payload.get("position_polypeptide"),
            correction_formula=payload.get("correction_formula"),
            mono_mass_delta=payload.get("mono_mass_delta"),
            avg_mass_delta=payload.get("avg_mass_delta"),
            cellular_location=payload.get("cellular_location"),
            taxonomic_range=tuple(payload.get("taxonomic_range", ())),
            keywords=tuple(payload.get("keywords", ())),
            xrefs={db: tuple(ids) for db, ids in payload.get("xrefs", {}).items()},
        )

    @property
    def target_residues(self) -> Tuple[str, ...]:
        seen: List[str] = []
        for group in self.target_groups:
            for residue in group:
                if residue not in seen:
                    seen.append(residue)
        return tuple(seen)

    @property
    def is_crosslink(self) -> bool:
        return self.feature_key == "CROSSLNK"

    def to_dict(self) -> Dict[str, object]:
        return {
            "accession": self.accession,
            "name": self.name,
            "feature_key": self.feature_key,
            "raw_target": self.raw_target,
            "target_groups": [list(group) for group in self.target_groups],
            "position_amino_acid": self.position_amino_acid,
            "position_polypeptide": self.position_polypeptide,
            "correction_formula": self.correction_formula,
            "mono_mass_delta": self.mono_mass_delta,
            "avg_mass_delta": self.avg_mass_delta,
            "cellular_location": self.cellular_location,
            "taxonomic_range": list(self.taxonomic_range),
            "keywords": list(self.keywords),
            "xrefs": {db: list(ids) for db, ids in self.xrefs.items()},
        }


@dataclass
class PTMHit:
    entry: PTMEntry
    matched_target: Tuple[str, ...]
    observed_delta: float
    error_da: float

    def to_dict(self) -> Dict[str, object]:
        return {
            "entry": self.entry.to_dict(),
            "matched_target": list(self.matched_target),
            "observed_delta": self.observed_delta,
            "error_da": self.error_da,
        }


@dataclass
class SiteObservation:
    position: int  # 1-based peptide position
    delta_mass: float
    polypeptide_position: Optional[str] = None
    amino_acid_position: Optional[str] = None


@dataclass
class SiteAnnotation:
    position: int
    residue: str
    observed_delta: float
    hits: Tuple[PTMHit, ...]

    def to_dict(self) -> Dict[str, object]:
        return {
            "position": self.position,
            "residue": self.residue,
            "observed_delta": self.observed_delta,
            "hits": [hit.to_dict() for hit in self.hits],
        }


@dataclass
class ModifiedPeptideParseResult:
    raw_sequence: str
    sequence: str
    observations: Tuple[SiteObservation, ...]

    def to_dict(self) -> Dict[str, object]:
        return {
            "raw_sequence": self.raw_sequence,
            "sequence": self.sequence,
            "observations": [
                {
                    "position": obs.position,
                    "delta_mass": obs.delta_mass,
                    "polypeptide_position": obs.polypeptide_position,
                    "amino_acid_position": obs.amino_acid_position,
                }
                for obs in self.observations
            ],
        }


@dataclass
class MassAlias:
    name: str
    modification_type: str
    mono_mass_delta: float
    allowed_residues: Tuple[str, ...] = ()
    allow_peptide_n_term_any_residue: bool = False
    allow_peptide_c_term_any_residue: bool = False
    note: Optional[str] = None

    def matches(self, residue: str, position: int, peptide_length: int) -> bool:
        residue_code = residue.upper()
        if self.allowed_residues and residue_code in self.allowed_residues:
            return True
        if self.allow_peptide_n_term_any_residue and position == 1:
            return True
        if self.allow_peptide_c_term_any_residue and position == peptide_length:
            return True
        return False


@dataclass
class PTMCandidate:
    display_name: str
    modification_type: Optional[str]
    source: str
    accession: Optional[str]
    observed_delta: float
    reference_delta: Optional[float]
    error_da: Optional[float]
    keywords: Tuple[str, ...] = ()
    note: Optional[str] = None

    def to_dict(self) -> Dict[str, object]:
        return {
            "display_name": self.display_name,
            "modification_type": self.modification_type,
            "source": self.source,
            "accession": self.accession,
            "observed_delta": self.observed_delta,
            "reference_delta": self.reference_delta,
            "error_da": self.error_da,
            "keywords": list(self.keywords),
            "note": self.note,
        }


@dataclass
class ModifiedSiteResult:
    position: int
    residue: str
    observed_delta: float
    candidates: Tuple[PTMCandidate, ...]

    @property
    def best_candidate(self) -> Optional[PTMCandidate]:
        return self.candidates[0] if self.candidates else None

    @property
    def assigned_modification_type(self) -> Optional[str]:
        best = self.best_candidate
        if best is None:
            return None
        return best.modification_type or best.display_name

    @property
    def assigned_label(self) -> Optional[str]:
        best = self.best_candidate
        if best is None:
            return None
        return self.assigned_modification_type or best.display_name

    def to_dict(self) -> Dict[str, object]:
        return {
            "position": self.position,
            "residue": self.residue,
            "observed_delta": self.observed_delta,
            "assigned_modification_type": self.assigned_modification_type,
            "assigned_label": self.assigned_label,
            "best_candidate": self.best_candidate.to_dict() if self.best_candidate is not None else None,
            "candidates": [candidate.to_dict() for candidate in self.candidates],
        }


@dataclass
class ModifiedPeptideResult:
    raw_sequence: str
    sequence: str
    modifications: Tuple[ModifiedSiteResult, ...]

    def to_dict(self) -> Dict[str, object]:
        return {
            "raw_sequence": self.raw_sequence,
            "sequence": self.sequence,
            "modifications": [site.to_dict() for site in self.modifications],
        }


@dataclass
class _MassIndex:
    masses: Tuple[float, ...]
    targets: Tuple[Tuple[str, ...], ...]
    entries: Tuple[PTMEntry, ...]

    @classmethod
    def build(cls, rows: Iterable[Tuple[float, Tuple[str, ...], PTMEntry]]) -> "_MassIndex":
        ordered = sorted(rows, key=lambda row: row[0])
        return cls(
            masses=tuple(row[0] for row in ordered),
            targets=tuple(row[1] for row in ordered),
            entries=tuple(row[2] for row in ordered),
        )

    def search(self, observed_delta: float, tolerance: float) -> List[PTMHit]:
        left = bisect.bisect_left(self.masses, observed_delta - tolerance)
        right = bisect.bisect_right(self.masses, observed_delta + tolerance)
        hits: List[PTMHit] = []
        for i in range(left, right):
            hits.append(
                PTMHit(
                    entry=self.entries[i],
                    matched_target=self.targets[i],
                    observed_delta=observed_delta,
                    error_da=observed_delta - self.masses[i],
                )
            )
        hits.sort(key=lambda hit: (abs(hit.error_da), hit.entry.name, hit.entry.accession))
        return hits


def _strip_terminal_period(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    value = value.strip()
    if value.endswith("."):
        return value[:-1]
    return value


def _normalize_residue(residue: str) -> str:
    key = residue.strip()
    if not key:
        raise ValueError("Residue cannot be empty.")
    alias = RESIDUE_ALIASES.get(key.lower())
    if alias is None:
        raise ValueError(f"Unknown residue: {residue!r}")
    return alias


def _normalize_residue_group(residues: Sequence[str]) -> Tuple[str, ...]:
    return tuple(_normalize_residue(residue) for residue in residues)


def _target_sort_key(group: Tuple[str, ...]) -> Tuple[str, ...]:
    return tuple(sorted(group))


def _mass_key(value: float, scale: int) -> int:
    return int(round(value * scale))


def _parse_target_groups(raw_target: str) -> Tuple[Tuple[str, ...], ...]:
    text = _strip_terminal_period(raw_target) or ""
    if not text or text == "Undefined":
        return tuple()
    groups: List[Tuple[str, ...]] = []
    for alternative in text.split(" or "):
        parts = [part.strip() for part in alternative.split("-")]
        groups.append(tuple(_normalize_residue(part) for part in parts))
    return tuple(groups)


def _parse_xrefs(values: Sequence[str]) -> Dict[str, Tuple[str, ...]]:
    grouped: DefaultDict[str, List[str]] = defaultdict(list)
    for value in values:
        db, _, ident = value.partition(";")
        db = db.strip()
        ident = _strip_terminal_period(ident.strip()) or ""
        grouped[db].append(ident)
    return {db: tuple(ids) for db, ids in grouped.items()}


def _extract_header_value(text: str, field_name: str) -> Optional[str]:
    prefix = f"{field_name}:"
    for line in text.splitlines()[:30]:
        if line.startswith(prefix):
            return line.split(":", 1)[1].strip()
    return None


def download_ptm_list(url: str = DEFAULT_PTM_URL, timeout: int = 60) -> str:
    with urlopen(url, timeout=timeout) as response:
        return response.read().decode("utf-8")


def parse_ptm_flatfile(text: str) -> List[PTMEntry]:
    lines = text.splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if line.startswith(ENTRY_SEPARATOR_PREFIX))
    except StopIteration as exc:
        raise ValueError("Could not locate the start of the PTM entry block.") from exc

    raw_entries: List[Dict[str, List[str]]] = []
    current: Optional[DefaultDict[str, List[str]]] = None
    current_code: Optional[str] = None

    for line in lines[start + 1 :]:
        if line == "//":
            if current is not None:
                raw_entries.append(dict(current))
                current = None
                current_code = None
            continue

        if not line.strip():
            continue

        match = FIELD_RE.match(line)
        if match:
            code, value = match.groups()
            if code not in LINE_CODES:
                continue
            if current is None:
                current = defaultdict(list)
            current[code].append(value.rstrip())
            current_code = code
            continue

        if current is not None and current_code is not None:
            current[current_code][-1] = f"{current[current_code][-1]} {line.strip()}".strip()
            continue

        if line.startswith("-") or line.startswith("Copyrighted") or line.startswith("Distributed"):
            continue

    entries: List[PTMEntry] = []
    for record in raw_entries:
        if "ID" not in record or "AC" not in record or "FT" not in record or "TG" not in record:
            continue
        entry = PTMEntry(
            accession=record["AC"][0].strip(),
            name=record["ID"][0].strip(),
            feature_key=record["FT"][0].strip(),
            raw_target=record["TG"][0].strip(),
            target_groups=_parse_target_groups(record["TG"][0].strip()),
            position_amino_acid=_strip_terminal_period(record.get("PA", [None])[0]),
            position_polypeptide=_strip_terminal_period(record.get("PP", [None])[0]),
            correction_formula=record.get("CF", [None])[0],
            mono_mass_delta=float(record["MM"][0]) if "MM" in record else None,
            avg_mass_delta=float(record["MA"][0]) if "MA" in record else None,
            cellular_location=_strip_terminal_period(record.get("LC", [None])[0]),
            taxonomic_range=tuple(_strip_terminal_period(value) or value for value in record.get("TR", [])),
            keywords=tuple(_strip_terminal_period(value) or value for value in record.get("KW", [])),
            xrefs=_parse_xrefs(record.get("DR", [])),
        )
        entries.append(entry)

    return entries


def _consume_bracketed_mass(text: str, start: int) -> Tuple[str, int]:
    opener = text[start]
    closer = ")" if opener == "(" else "]"
    end = text.find(closer, start + 1)
    if end == -1:
        raise ValueError(f"Unclosed modification mass starting at index {start}: {text[start:]!r}")
    mass_text = text[start + 1 : end].strip()
    if not mass_text:
        raise ValueError("Encountered an empty modification mass annotation.")
    return mass_text, end + 1


def parse_modified_peptide_string(modified_sequence: str) -> ModifiedPeptideParseResult:
    """
    Parse peptide strings where mass shifts are written inline immediately after
    the modified residue, for example:
        A(+27.99)EDPETQVVL
        A(+114.04)IPRSPFEVQVSPE
        AAAAVRQM(+15.99)NPHIRVT
        M(+15.99)PEPTIDEK(+114.04)

    Multiple inline modifications within the same peptide are supported. The
    function also tolerates square brackets and a leading N-terminal mass
    annotation such as:
        (+42.01)ACDEFGHIK
    """
    raw = modified_sequence.strip()
    if not raw:
        raise ValueError("Modified peptide string cannot be empty.")

    text = raw.replace(" ", "")
    residues: List[str] = []
    observations: List[SiteObservation] = []
    pending_prefix_masses: List[float] = []

    position = 0
    i = 0
    while i < len(text):
        ch = text[i]

        if ch in ".-_":
            i += 1
            continue

        if ch in "([":
            mass_text, next_index = _consume_bracketed_mass(text, i)
            try:
                delta_mass = float(mass_text)
            except ValueError as exc:
                raise ValueError(f"Could not parse mass shift {mass_text!r} in {raw!r}.") from exc
            if position == 0:
                pending_prefix_masses.append(delta_mass)
            else:
                observations.append(SiteObservation(position=position, delta_mass=delta_mass))
            i = next_index
            continue

        if ch.isalpha():
            residue = ch.upper()
            residues.append(residue)
            position += 1
            if pending_prefix_masses:
                for delta_mass in pending_prefix_masses:
                    observations.append(SiteObservation(position=position, delta_mass=delta_mass))
                pending_prefix_masses.clear()
            i += 1
            continue

        raise ValueError(f"Unexpected character {ch!r} in modified peptide string {raw!r}.")

    if pending_prefix_masses:
        raise ValueError(f"Found leading modification mass with no residue to attach it to in {raw!r}.")

    return ModifiedPeptideParseResult(raw_sequence=raw, sequence="".join(residues), observations=tuple(observations))


def _normalize_keyword(keyword: str) -> Optional[str]:
    keyword = keyword.strip()
    if not keyword:
        return None

    keyword_map = {
        "Phosphoprotein": "Phosphorylation",
        "Glycoprotein": "Glycosylation",
        "Proteoglycan": "Glycosylation",
        "Lipoprotein": "Lipidation",
        "GPI-anchor": "Lipidation",
        "Palmitate": "Palmitoylation",
        "Myristate": "Myristoylation",
        "Thioether bond": "Cross-link",
        "Thioester bond": "Cross-link",
        "Isopeptide bond": "Cross-link",
    }
    return keyword_map.get(keyword, keyword)


def _keyword_priority(keyword: str) -> Tuple[int, str]:
    preferred = {
        "Ubiquitination": 0,
        "Phosphorylation": 1,
        "Methylation": 2,
        "Acetylation": 3,
        "Hydroxylation": 4,
        "Oxidation": 5,
        "Sulfation": 6,
        "Glycosylation": 7,
        "Lipidation": 8,
        "Cross-link": 9,
    }
    return preferred.get(keyword, 100), keyword


def infer_modification_type(entry: PTMEntry) -> Optional[str]:
    """
    Derive a broad PTM class such as 'Methylation' or 'Phosphorylation' from a
    UniProt PTM entry. Keywords are preferred when present, with fallbacks based
    on the entry name and feature key.
    """
    normalized_keywords: List[str] = []
    for raw_keyword in entry.keywords:
        for piece in raw_keyword.split(";"):
            normalized = _normalize_keyword(piece)
            if normalized is not None:
                normalized_keywords.append(normalized)

    if normalized_keywords:
        keyword = sorted(normalized_keywords, key=_keyword_priority)[0]
        return keyword

    name = entry.name.lower()

    pattern_map = (
        ("ubiquit", "Ubiquitination"),
        ("digly", "Ubiquitination"),
        ("glygly", "Ubiquitination"),
        ("methyl", "Methylation"),
        ("acetyl", "Acetylation"),
        ("formyl", "Formylation"),
        ("phospho", "Phosphorylation"),
        ("hydroxy", "Hydroxylation"),
        ("oxid", "Oxidation"),
        ("sulf", "Sulfation"),
        ("nitro", "Nitration"),
        ("iodo", "Iodination"),
        ("biotinyl", "Biotinylation"),
        ("malonyl", "Malonylation"),
        ("succinyl", "Succinylation"),
        ("glutaryl", "Glutarylation"),
        ("crotonyl", "Crotonylation"),
        ("butyryl", "Butyrylation"),
        ("propionyl", "Propionylation"),
        ("lactyl", "Lactylation"),
        ("myristoyl", "Myristoylation"),
        ("palmitoyl", "Palmitoylation"),
        ("prenyl", "Prenylation"),
        ("farnesyl", "Prenylation"),
        ("geranylgeranyl", "Prenylation"),
        ("glycyl", "Glycylation"),
        ("glutamyl", "Glutamylation"),
    )
    for pattern, label in pattern_map:
        if pattern in name:
            return label

    if entry.feature_key == "CARBOHYD":
        return "Glycosylation"
    if entry.feature_key == "LIPID":
        return "Lipidation"
    if entry.feature_key == "CROSSLNK":
        return "Cross-link"

    return None


# Conservative alias list. It is outside UniProt curation and only used if no
# direct UniProt candidate exists for the modified residue + mass pair.
DEFAULT_FALLBACK_MASS_ALIASES = (
    MassAlias(
        name="Ubiquitination (digly remnant)",
        modification_type="Ubiquitination",
        mono_mass_delta=114.042927,
        allowed_residues=("K",),
        allow_peptide_n_term_any_residue=True,
        note=(
            "Curated fallback alias for the common +114.0429 digly remnant used "
            "in peptide search outputs. This shorthand is not listed as a direct "
            "site entry in the UniProt PTM list."
        ),
    ),
)


class UniProtPTMLookup:
    def __init__(
        self,
        entries: Sequence[PTMEntry],
        *,
        release: Optional[str] = None,
        source_url: Optional[str] = None,
        exact_mass_scale: int = 1_000_000,
        fallback_mass_aliases: Sequence[MassAlias] = DEFAULT_FALLBACK_MASS_ALIASES,
    ) -> None:
        self.entries: Tuple[PTMEntry, ...] = tuple(entries)
        self.release = release
        self.source_url = source_url
        self.exact_mass_scale = exact_mass_scale
        self.fallback_mass_aliases: Tuple[MassAlias, ...] = tuple(fallback_mass_aliases)

        self.by_accession: Dict[str, PTMEntry] = {entry.accession: entry for entry in self.entries}
        self.by_name: DefaultDict[str, List[PTMEntry]] = defaultdict(list)
        self.by_xref: DefaultDict[Tuple[str, str], List[PTMEntry]] = defaultdict(list)

        site_mono_rows: DefaultDict[str, List[Tuple[float, Tuple[str, ...], PTMEntry]]] = defaultdict(list)
        site_avg_rows: DefaultDict[str, List[Tuple[float, Tuple[str, ...], PTMEntry]]] = defaultdict(list)
        crosslink_mono_rows_ordered: DefaultDict[Tuple[str, ...], List[Tuple[float, Tuple[str, ...], PTMEntry]]] = defaultdict(list)
        crosslink_avg_rows_ordered: DefaultDict[Tuple[str, ...], List[Tuple[float, Tuple[str, ...], PTMEntry]]] = defaultdict(list)
        crosslink_mono_rows_unordered: DefaultDict[Tuple[str, ...], List[Tuple[float, Tuple[str, ...], PTMEntry]]] = defaultdict(list)
        crosslink_avg_rows_unordered: DefaultDict[Tuple[str, ...], List[Tuple[float, Tuple[str, ...], PTMEntry]]] = defaultdict(list)

        self._site_exact_mono: DefaultDict[Tuple[str, int], List[Tuple[Tuple[str, ...], PTMEntry]]] = defaultdict(list)
        self._site_exact_avg: DefaultDict[Tuple[str, int], List[Tuple[Tuple[str, ...], PTMEntry]]] = defaultdict(list)
        self._crosslink_exact_mono_ordered: DefaultDict[Tuple[Tuple[str, ...], int], List[Tuple[Tuple[str, ...], PTMEntry]]] = defaultdict(list)
        self._crosslink_exact_avg_ordered: DefaultDict[Tuple[Tuple[str, ...], int], List[Tuple[Tuple[str, ...], PTMEntry]]] = defaultdict(list)
        self._crosslink_exact_mono_unordered: DefaultDict[Tuple[Tuple[str, ...], int], List[Tuple[Tuple[str, ...], PTMEntry]]] = defaultdict(list)
        self._crosslink_exact_avg_unordered: DefaultDict[Tuple[Tuple[str, ...], int], List[Tuple[Tuple[str, ...], PTMEntry]]] = defaultdict(list)

        for entry in self.entries:
            self.by_name[entry.name.lower()].append(entry)
            for db, identifiers in entry.xrefs.items():
                for identifier in identifiers:
                    self.by_xref[(db.upper(), identifier)].append(entry)

            for target_group in entry.target_groups:
                if len(target_group) == 1 and not entry.is_crosslink:
                    residue = target_group[0]
                    if entry.mono_mass_delta is not None:
                        site_mono_rows[residue].append((entry.mono_mass_delta, target_group, entry))
                        key = (residue, _mass_key(entry.mono_mass_delta, self.exact_mass_scale))
                        self._site_exact_mono[key].append((target_group, entry))
                    if entry.avg_mass_delta is not None:
                        site_avg_rows[residue].append((entry.avg_mass_delta, target_group, entry))
                        key = (residue, _mass_key(entry.avg_mass_delta, self.exact_mass_scale))
                        self._site_exact_avg[key].append((target_group, entry))
                else:
                    ordered_key = target_group
                    unordered_key = _target_sort_key(target_group)
                    if entry.mono_mass_delta is not None:
                        crosslink_mono_rows_ordered[ordered_key].append((entry.mono_mass_delta, target_group, entry))
                        crosslink_mono_rows_unordered[unordered_key].append((entry.mono_mass_delta, target_group, entry))
                        ordered_index_key = (ordered_key, _mass_key(entry.mono_mass_delta, self.exact_mass_scale))
                        unordered_index_key = (unordered_key, _mass_key(entry.mono_mass_delta, self.exact_mass_scale))
                        self._crosslink_exact_mono_ordered[ordered_index_key].append((target_group, entry))
                        self._crosslink_exact_mono_unordered[unordered_index_key].append((target_group, entry))
                    if entry.avg_mass_delta is not None:
                        crosslink_avg_rows_ordered[ordered_key].append((entry.avg_mass_delta, target_group, entry))
                        crosslink_avg_rows_unordered[unordered_key].append((entry.avg_mass_delta, target_group, entry))
                        ordered_index_key = (ordered_key, _mass_key(entry.avg_mass_delta, self.exact_mass_scale))
                        unordered_index_key = (unordered_key, _mass_key(entry.avg_mass_delta, self.exact_mass_scale))
                        self._crosslink_exact_avg_ordered[ordered_index_key].append((target_group, entry))
                        self._crosslink_exact_avg_unordered[unordered_index_key].append((target_group, entry))

        self._site_mono_index = {key: _MassIndex.build(rows) for key, rows in site_mono_rows.items()}
        self._site_avg_index = {key: _MassIndex.build(rows) for key, rows in site_avg_rows.items()}
        self._crosslink_mono_index_ordered = {key: _MassIndex.build(rows) for key, rows in crosslink_mono_rows_ordered.items()}
        self._crosslink_avg_index_ordered = {key: _MassIndex.build(rows) for key, rows in crosslink_avg_rows_ordered.items()}
        self._crosslink_mono_index_unordered = {key: _MassIndex.build(rows) for key, rows in crosslink_mono_rows_unordered.items()}
        self._crosslink_avg_index_unordered = {key: _MassIndex.build(rows) for key, rows in crosslink_avg_rows_unordered.items()}

    @classmethod
    def from_text(cls, text: str, *, source_url: Optional[str] = None) -> "UniProtPTMLookup":
        return cls(
            parse_ptm_flatfile(text),
            release=_extract_header_value(text, "Release"),
            source_url=source_url,
        )

    @classmethod
    def from_path(cls, path: Union[str, Path], *, source_url: Optional[str] = None) -> "UniProtPTMLookup":
        text = Path(path).read_text(encoding="utf-8")
        return cls.from_text(text, source_url=source_url)

    @classmethod
    def from_json_export(cls, path: Union[str, Path]) -> "UniProtPTMLookup":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        entries = [PTMEntry.from_dict(entry) for entry in payload["entries"]]
        return cls(
            entries,
            release=payload.get("release"),
            source_url=payload.get("source_url"),
        )

    @classmethod
    def from_packaged_library(cls) -> "UniProtPTMLookup":
        if hasattr(importlib_resources, "files"):
            package_root = importlib_resources.files("uniprot_ptm_lookup")
            payload_text = package_root.joinpath(PACKAGED_LIBRARY_FILENAME).read_text(encoding="utf-8")
        else:
            payload_text = importlib_resources.read_text(
                "uniprot_ptm_lookup",
                PACKAGED_LIBRARY_FILENAME,
                encoding="utf-8",
            )
        payload = json.loads(payload_text)
        entries = [PTMEntry.from_dict(entry) for entry in payload["entries"]]
        return cls(
            entries,
            release=payload.get("release"),
            source_url=payload.get("source_url"),
        )

    @classmethod
    def from_uniprot(
        cls,
        *,
        url: str = DEFAULT_PTM_URL,
        cache_path: Optional[Union[str, Path]] = None,
        refresh: bool = False,
        timeout: int = 60,
    ) -> "UniProtPTMLookup":
        text: str
        cache_file = Path(cache_path) if cache_path is not None else None
        if cache_file is not None and cache_file.exists() and not refresh:
            text = cache_file.read_text(encoding="utf-8")
        else:
            text = download_ptm_list(url=url, timeout=timeout)
            if cache_file is not None:
                cache_file.write_text(text, encoding="utf-8")
        return cls.from_text(text, source_url=url)

    def __len__(self) -> int:
        return len(self.entries)

    def get(self, accession: str) -> PTMEntry:
        return self.by_accession[accession]

    def find_by_name(self, query: str) -> List[PTMEntry]:
        query_lc = query.strip().lower()
        return [entry for entry in self.entries if query_lc in entry.name.lower()]

    def find_by_xref(self, db: str, identifier: str) -> List[PTMEntry]:
        return list(self.by_xref.get((db.upper(), identifier), []))

    def stats(self) -> Dict[str, object]:
        feature_counts: DefaultDict[str, int] = defaultdict(int)
        with_mono = 0
        with_avg = 0
        for entry in self.entries:
            feature_counts[entry.feature_key] += 1
            if entry.mono_mass_delta is not None:
                with_mono += 1
            if entry.avg_mass_delta is not None:
                with_avg += 1
        return {
            "release": self.release,
            "source_url": self.source_url,
            "n_entries": len(self.entries),
            "n_with_mono_mass": with_mono,
            "n_with_avg_mass": with_avg,
            "feature_counts": dict(sorted(feature_counts.items())),
        }

    def lookup_site(
        self,
        residue: str,
        observed_delta: float,
        *,
        tolerance: float = 0.01,
        feature_keys: Optional[Iterable[str]] = None,
        polypeptide_position: Optional[str] = None,
        amino_acid_position: Optional[str] = None,
        use_average_mass: bool = False,
        exact: bool = False,
    ) -> List[PTMHit]:
        residue_code = _normalize_residue(residue)
        allowed_feature_keys = frozenset(feature_keys) if feature_keys is not None else DEFAULT_SITE_FEATURE_KEYS
        polypeptide_position = _strip_terminal_period(polypeptide_position)
        amino_acid_position = _strip_terminal_period(amino_acid_position)

        if exact or tolerance == 0:
            exact_index = self._site_exact_avg if use_average_mass else self._site_exact_mono
            exact_key = (residue_code, _mass_key(observed_delta, self.exact_mass_scale))
            hits = [
                PTMHit(entry=entry, matched_target=target_group, observed_delta=observed_delta, error_da=0.0)
                for target_group, entry in exact_index.get(exact_key, [])
            ]
        else:
            mass_index = self._site_avg_index if use_average_mass else self._site_mono_index
            index = mass_index.get(residue_code)
            if index is None:
                return []
            hits = index.search(observed_delta, tolerance)

        return [
            hit
            for hit in hits
            if self._site_hit_matches(
                hit,
                allowed_feature_keys=allowed_feature_keys,
                polypeptide_position=polypeptide_position,
                amino_acid_position=amino_acid_position,
            )
        ]

    def lookup_crosslink(
        self,
        residues: Sequence[str],
        observed_delta: float,
        *,
        tolerance: float = 0.01,
        feature_keys: Optional[Iterable[str]] = None,
        use_average_mass: bool = False,
        exact: bool = False,
        ordered: bool = False,
    ) -> List[PTMHit]:
        residue_group = _normalize_residue_group(residues)
        lookup_key = residue_group if ordered else _target_sort_key(residue_group)
        allowed_feature_keys = frozenset(feature_keys) if feature_keys is not None else DEFAULT_CROSSLINK_FEATURE_KEYS

        if exact or tolerance == 0:
            if use_average_mass:
                exact_index = self._crosslink_exact_avg_ordered if ordered else self._crosslink_exact_avg_unordered
            else:
                exact_index = self._crosslink_exact_mono_ordered if ordered else self._crosslink_exact_mono_unordered
            exact_key = (lookup_key, _mass_key(observed_delta, self.exact_mass_scale))
            hits = [
                PTMHit(entry=entry, matched_target=target_group, observed_delta=observed_delta, error_da=0.0)
                for target_group, entry in exact_index.get(exact_key, [])
            ]
        else:
            if use_average_mass:
                mass_index = self._crosslink_avg_index_ordered if ordered else self._crosslink_avg_index_unordered
            else:
                mass_index = self._crosslink_mono_index_ordered if ordered else self._crosslink_mono_index_unordered
            index = mass_index.get(lookup_key)
            if index is None:
                return []
            hits = index.search(observed_delta, tolerance)

        return [hit for hit in hits if hit.entry.feature_key in allowed_feature_keys]

    def annotate_peptide(
        self,
        sequence: str,
        observations: Sequence[SiteObservation],
        *,
        tolerance: float = 0.01,
        feature_keys: Optional[Iterable[str]] = None,
        use_average_mass: bool = False,
        exact: bool = False,
    ) -> List[SiteAnnotation]:
        peptide = sequence.strip().upper()
        annotations: List[SiteAnnotation] = []
        for obs in observations:
            if obs.position < 1 or obs.position > len(peptide):
                raise IndexError(
                    f"Observation position {obs.position} is outside peptide length {len(peptide)}."
                )
            residue = peptide[obs.position - 1]
            hits = self.lookup_site(
                residue,
                obs.delta_mass,
                tolerance=tolerance,
                feature_keys=feature_keys,
                polypeptide_position=obs.polypeptide_position,
                amino_acid_position=obs.amino_acid_position,
                use_average_mass=use_average_mass,
                exact=exact,
            )
            annotations.append(
                SiteAnnotation(
                    position=obs.position,
                    residue=residue,
                    observed_delta=obs.delta_mass,
                    hits=tuple(hits),
                )
            )
        return annotations

    @staticmethod
    def parse_modified_peptide(modified_sequence: str) -> ModifiedPeptideParseResult:
        return parse_modified_peptide_string(modified_sequence)

    def annotate_modified_peptide(
        self,
        modified_sequence: str,
        *,
        tolerance: float = DEFAULT_MODIFIED_PEPTIDE_TOLERANCE,
        feature_keys: Optional[Iterable[str]] = None,
        use_average_mass: bool = False,
        exact: bool = False,
        allow_fallback_aliases: bool = True,
        max_candidates: int = 10,
    ) -> ModifiedPeptideResult:
        parsed = parse_modified_peptide_string(modified_sequence)
        peptide = parsed.sequence
        modifications: List[ModifiedSiteResult] = []

        for obs in parsed.observations:
            residue = peptide[obs.position - 1]
            try:
                site_hits = self.lookup_site(
                    residue,
                    obs.delta_mass,
                    tolerance=tolerance,
                    feature_keys=feature_keys,
                    polypeptide_position=obs.polypeptide_position,
                    amino_acid_position=obs.amino_acid_position,
                    use_average_mass=use_average_mass,
                    exact=exact,
                )
            except ValueError:
                site_hits = []

            candidates = self._site_hits_to_candidates(
                site_hits,
                use_average_mass=use_average_mass,
            )

            if allow_fallback_aliases and not use_average_mass:
                candidates.extend(
                    self._lookup_fallback_alias_candidates(
                        residue=residue,
                        position=obs.position,
                        peptide_length=len(peptide),
                        observed_delta=obs.delta_mass,
                        tolerance=tolerance,
                    )
                )

            self._sort_candidates_in_place(candidates)

            if max_candidates > 0:
                candidates = candidates[:max_candidates]

            modifications.append(
                ModifiedSiteResult(
                    position=obs.position,
                    residue=residue,
                    observed_delta=obs.delta_mass,
                    candidates=tuple(candidates),
                )
            )

        return ModifiedPeptideResult(
            raw_sequence=parsed.raw_sequence,
            sequence=parsed.sequence,
            modifications=tuple(modifications),
        )

    def annotate_modified_peptides(
        self,
        modified_sequences: Sequence[str],
        *,
        tolerance: float = DEFAULT_MODIFIED_PEPTIDE_TOLERANCE,
        feature_keys: Optional[Iterable[str]] = None,
        use_average_mass: bool = False,
        exact: bool = False,
        allow_fallback_aliases: bool = True,
        max_candidates: int = 10,
    ) -> List[ModifiedPeptideResult]:
        return [
            self.annotate_modified_peptide(
                modified_sequence,
                tolerance=tolerance,
                feature_keys=feature_keys,
                use_average_mass=use_average_mass,
                exact=exact,
                allow_fallback_aliases=allow_fallback_aliases,
                max_candidates=max_candidates,
            )
            for modified_sequence in modified_sequences
        ]

    def export_entries_json(self, path: Union[str, Path]) -> None:
        payload = {
            "release": self.release,
            "source_url": self.source_url,
            "entries": [entry.to_dict() for entry in self.entries],
        }
        Path(path).write_text(json.dumps(payload, indent=2), encoding="utf-8")

    @staticmethod
    def results_to_flat_rows(
        results: Sequence[ModifiedPeptideResult],
        *,
        include_all_candidates: bool = True,
        candidate_separator: str = " | ",
    ) -> List[Dict[str, object]]:
        rows: List[Dict[str, object]] = []
        for result in results:
            if not result.modifications:
                rows.append(
                    {
                        "raw_sequence": result.raw_sequence,
                        "sequence": result.sequence,
                        "status": "no_inline_modification_found",
                        "position": "",
                        "residue": "",
                        "observed_delta": "",
                        "assigned_modification_type": "",
                        "assigned_label": "",
                        "best_candidate_name": "",
                        "best_candidate_source": "",
                        "best_candidate_accession": "",
                        "best_reference_delta": "",
                        "best_error_da": "",
                        "n_candidates": 0,
                    }
                )
                continue

            for site in result.modifications:
                best = site.best_candidate
                row: Dict[str, object] = {
                    "raw_sequence": result.raw_sequence,
                    "sequence": result.sequence,
                    "status": "ok",
                    "position": site.position,
                    "residue": site.residue,
                    "observed_delta": site.observed_delta,
                    "assigned_modification_type": site.assigned_modification_type or "",
                    "assigned_label": site.assigned_label or "",
                    "best_candidate_name": best.display_name if best is not None else "",
                    "best_candidate_source": best.source if best is not None else "",
                    "best_candidate_accession": best.accession if best is not None else "",
                    "best_reference_delta": best.reference_delta if best is not None else "",
                    "best_error_da": best.error_da if best is not None else "",
                    "n_candidates": len(site.candidates),
                }
                if include_all_candidates:
                    row["all_candidate_names"] = candidate_separator.join(candidate.display_name for candidate in site.candidates)
                    row["all_candidate_types"] = candidate_separator.join(
                        candidate.modification_type or "" for candidate in site.candidates
                    )
                    row["all_candidate_sources"] = candidate_separator.join(candidate.source for candidate in site.candidates)
                    row["all_candidate_accessions"] = candidate_separator.join(
                        candidate.accession or "" for candidate in site.candidates
                    )
                rows.append(row)
        return rows

    @classmethod
    def export_modified_peptide_results_csv(
        cls,
        results: Sequence[ModifiedPeptideResult],
        path: Union[str, Path],
        *,
        include_all_candidates: bool = True,
        candidate_separator: str = " | ",
    ) -> None:
        rows = cls.results_to_flat_rows(
            results,
            include_all_candidates=include_all_candidates,
            candidate_separator=candidate_separator,
        )
        fieldnames = [
            "raw_sequence",
            "sequence",
            "status",
            "position",
            "residue",
            "observed_delta",
            "assigned_modification_type",
            "assigned_label",
            "best_candidate_name",
            "best_candidate_source",
            "best_candidate_accession",
            "best_reference_delta",
            "best_error_da",
            "n_candidates",
        ]
        if include_all_candidates:
            fieldnames.extend(
                [
                    "all_candidate_names",
                    "all_candidate_types",
                    "all_candidate_sources",
                    "all_candidate_accessions",
                ]
            )
        output_path = Path(path)
        with output_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

    def annotate_modified_peptides_to_csv(
        self,
        modified_sequences: Sequence[str],
        path: Union[str, Path],
        *,
        tolerance: float = DEFAULT_MODIFIED_PEPTIDE_TOLERANCE,
        feature_keys: Optional[Iterable[str]] = None,
        use_average_mass: bool = False,
        exact: bool = False,
        allow_fallback_aliases: bool = True,
        max_candidates: int = 10,
        include_all_candidates: bool = True,
        candidate_separator: str = " | ",
    ) -> List[ModifiedPeptideResult]:
        results = self.annotate_modified_peptides(
            modified_sequences,
            tolerance=tolerance,
            feature_keys=feature_keys,
            use_average_mass=use_average_mass,
            exact=exact,
            allow_fallback_aliases=allow_fallback_aliases,
            max_candidates=max_candidates,
        )
        self.export_modified_peptide_results_csv(
            results,
            path,
            include_all_candidates=include_all_candidates,
            candidate_separator=candidate_separator,
        )
        return results

    @staticmethod
    def _site_hit_matches(
        hit: PTMHit,
        *,
        allowed_feature_keys: FrozenSet[str],
        polypeptide_position: Optional[str],
        amino_acid_position: Optional[str],
    ) -> bool:
        entry = hit.entry
        if entry.feature_key not in allowed_feature_keys:
            return False

        if amino_acid_position is not None and entry.position_amino_acid not in (None, amino_acid_position):
            return False

        if polypeptide_position is not None:
            if entry.position_polypeptide is None:
                return True
            if entry.position_polypeptide == "Anywhere":
                return True
            if entry.position_polypeptide != polypeptide_position:
                return False

        return True

    @staticmethod
    def _site_hits_to_candidates(
        hits: Sequence[PTMHit],
        *,
        use_average_mass: bool,
    ) -> List[PTMCandidate]:
        candidates: List[PTMCandidate] = []
        for hit in hits:
            reference_delta = hit.entry.avg_mass_delta if use_average_mass else hit.entry.mono_mass_delta
            candidates.append(
                PTMCandidate(
                    display_name=hit.entry.name,
                    modification_type=infer_modification_type(hit.entry),
                    source="UniProt",
                    accession=hit.entry.accession,
                    observed_delta=hit.observed_delta,
                    reference_delta=reference_delta,
                    error_da=hit.error_da,
                    keywords=hit.entry.keywords,
                    note=None,
                )
            )
        return candidates

    @staticmethod
    def _candidate_source_priority(source: str) -> int:
        return 0 if source == "UniProt" else 1

    @staticmethod
    def _candidate_specificity_penalty(display_name: str) -> Tuple[int, int, str]:
        name = display_name.lower()
        stereochemistry_penalty = 1 if any(tag in name for tag in ("(r)", "(s)", "(e)", "(z)")) else 0
        parenthetical_penalty = name.count("(")
        return stereochemistry_penalty, parenthetical_penalty, display_name

    @classmethod
    def _candidate_sort_key(cls, candidate: PTMCandidate) -> Tuple[float, int, Tuple[int, int, str], str]:
        error = abs(candidate.error_da) if candidate.error_da is not None else float("inf")
        return (
            error,
            cls._candidate_source_priority(candidate.source),
            cls._candidate_specificity_penalty(candidate.display_name),
            candidate.display_name,
        )

    @classmethod
    def _sort_candidates_in_place(cls, candidates: List[PTMCandidate]) -> None:
        candidates.sort(key=cls._candidate_sort_key)

    def _lookup_fallback_alias_candidates(
        self,
        *,
        residue: str,
        position: int,
        peptide_length: int,
        observed_delta: float,
        tolerance: float,
    ) -> List[PTMCandidate]:
        candidates: List[PTMCandidate] = []
        for alias in self.fallback_mass_aliases:
            if abs(observed_delta - alias.mono_mass_delta) <= tolerance and alias.matches(
                residue=residue,
                position=position,
                peptide_length=peptide_length,
            ):
                candidates.append(
                    PTMCandidate(
                        display_name=alias.name,
                        modification_type=alias.modification_type,
                        source="Curated mass alias",
                        accession=None,
                        observed_delta=observed_delta,
                        reference_delta=alias.mono_mass_delta,
                        error_da=observed_delta - alias.mono_mass_delta,
                        keywords=(),
                        note=alias.note,
                    )
                )
        self._sort_candidates_in_place(candidates)
        return candidates


_DEFAULT_LOOKUP: Optional[UniProtPTMLookup] = None


def load_default_lookup(*, cache: bool = True) -> UniProtPTMLookup:
    """Load the packaged PTM library that ships with this distribution.

    When ``cache`` is True, the parsed lookup object is reused within the
    current Python process so repeated calls are cheap.
    """
    global _DEFAULT_LOOKUP
    if not cache:
        return UniProtPTMLookup.from_packaged_library()
    if _DEFAULT_LOOKUP is None:
        _DEFAULT_LOOKUP = UniProtPTMLookup.from_packaged_library()
    return _DEFAULT_LOOKUP


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build and query a UniProt PTM lookup index.")
    parser.add_argument(
        "--use-packaged-library",
        action="store_true",
        help="Use the PTM JSON library bundled with this installed package instead of downloading UniProt.",
    )
    parser.add_argument("--cache", type=Path, default=Path("ptmlist.txt"), help="Path for raw UniProt PTM flatfile cache.")
    parser.add_argument("--library-json", type=Path, help="Path to a previously exported PTM JSON library file.")
    parser.add_argument("--refresh", action="store_true", help="Re-download the flatfile even if the cache exists.")
    parser.add_argument("--stats", action="store_true", help="Print a JSON summary of the parsed PTM library.")
    parser.add_argument("--residue", help="Residue for site lookup, e.g. S or Serine.")
    parser.add_argument("--delta", type=float, help="Observed mass delta for lookup.")
    parser.add_argument("--tolerance", type=float, default=0.01, help="Mass tolerance in Da for low-level lookup.")
    parser.add_argument("--sequence", help="Peptide sequence for low-level annotation.")
    parser.add_argument(
        "--obs",
        action="append",
        default=[],
        help=(
            "Observed site modification in the form position:delta_mass, e.g. 2:79.966331. "
            "Can be provided more than once."
        ),
    )
    parser.add_argument(
        "--modified-peptide",
        action="append",
        default=[],
        help=(
            "Modified peptide string such as A(+27.99)EDPETQVVL. "
            "Can be provided more than once for batch annotation."
        ),
    )
    parser.add_argument(
        "--modified-peptide-file",
        type=Path,
        help="Text file containing one modified peptide string per line.",
    )
    parser.add_argument(
        "--modified-tolerance",
        type=float,
        default=DEFAULT_MODIFIED_PEPTIDE_TOLERANCE,
        help="Mass tolerance in Da for modified peptide string annotation.",
    )
    parser.add_argument(
        "--max-candidates",
        type=int,
        default=10,
        help="Maximum number of candidates to retain per modified site.",
    )
    parser.add_argument(
        "--export-modified-csv",
        type=Path,
        help="Write flattened modified peptide annotation results to a CSV file.",
    )
    parser.add_argument(
        "--compact-csv",
        action="store_true",
        help="When exporting CSV, omit the semicolon-separated all-candidate columns.",
    )
    parser.add_argument("--export-json", type=Path, help="Write parsed entries to a JSON file.")
    return parser


def _parse_cli_observations(items: Sequence[str]) -> List[SiteObservation]:
    observations: List[SiteObservation] = []
    for item in items:
        position_text, sep, delta_text = item.partition(":")
        if not sep:
            raise ValueError(f"Invalid observation {item!r}; expected position:delta_mass")
        observations.append(SiteObservation(position=int(position_text), delta_mass=float(delta_text)))
    return observations


def _read_modified_peptide_file(path: Path) -> List[str]:
    items: List[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        items.append(text)
    return items


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()

    if args.library_json is not None:
        lib = UniProtPTMLookup.from_json_export(args.library_json)
    elif args.use_packaged_library:
        lib = UniProtPTMLookup.from_packaged_library()
    else:
        lib = UniProtPTMLookup.from_uniprot(cache_path=args.cache, refresh=args.refresh)

    if args.export_json is not None:
        lib.export_entries_json(args.export_json)

    if args.stats:
        print(json.dumps(lib.stats(), indent=2))

    if args.residue is not None and args.delta is not None:
        hits = lib.lookup_site(args.residue, args.delta, tolerance=args.tolerance)
        print(json.dumps([hit.to_dict() for hit in hits], indent=2))

    if args.sequence and args.obs:
        observations = _parse_cli_observations(args.obs)
        annotations = lib.annotate_peptide(args.sequence, observations, tolerance=args.tolerance)
        print(json.dumps([annotation.to_dict() for annotation in annotations], indent=2))

    modified_sequences = list(args.modified_peptide)
    if args.modified_peptide_file is not None:
        modified_sequences.extend(_read_modified_peptide_file(args.modified_peptide_file))

    if modified_sequences:
        results = lib.annotate_modified_peptides(
            modified_sequences,
            tolerance=args.modified_tolerance,
            max_candidates=args.max_candidates,
        )
        print(json.dumps([result.to_dict() for result in results], indent=2))

        if args.export_modified_csv is not None:
            lib.export_modified_peptide_results_csv(
                results,
                args.export_modified_csv,
                include_all_candidates=not args.compact_csv,
            )
    elif args.export_modified_csv is not None:
        raise ValueError("--export-modified-csv was provided, but no modified peptide inputs were supplied.")


if __name__ == "__main__":
    main()
