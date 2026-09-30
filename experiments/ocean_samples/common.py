
"""Utilities for the 159-sample ocean-community experiment."""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import pandas as pd
import yaml


MODEL_SUFFIXES = (".xml", ".sbml")


def safe_name(value: str) -> str:
    """Return a filesystem/solver-friendly name fragment."""
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._-")
    return cleaned or "unnamed"


def normalized_id(value: str) -> str:
    """Normalize an identifier for conservative alias matching."""
    text = str(value).strip().lower()

    changed = True
    while changed:
        changed = False
        for suffix in (
            ".xml", ".sbml", ".gz", "_model", "-model", ".model", "_gapseq", "-gapseq",
        ):
            if text.endswith(suffix):
                text = text[: -len(suffix)]
                changed = True

    for prefix in ("model_", "model-", "gapseq_", "gapseq-"):
        if text.startswith(prefix):
            text = text[len(prefix) :]

    return re.sub(r"[^a-z0-9]+", "", text)


def _choose_unique_id_column(
    raw: pd.DataFrame,
    candidate_columns: list[str],
    expected_count: int,
    entity_label: str,
) -> str:
    """Choose a unique identifier column from non-abundance metadata."""

    unique_candidates = []

    for column in candidate_columns:
        values = raw[column].map(
            lambda value: str(value).strip()
        )

        if (
            len(values) == expected_count
            and not values.eq("").any()
            and values.nunique(dropna=False) == expected_count
        ):
            lower = column.lower()

            score = 0

            if lower in {
                "mag_id",
                "magid",
                "genome_id",
                "genomeid",
                "bin_id",
                "binid",
                "sample_id",
                "sampleid",
            }:
                score += 100

            if "mag" in lower:
                score += 30

            if "genome" in lower:
                score += 25

            if "bin" in lower:
                score += 20

            if "sample" in lower:
                score += 20

            if lower == "id" or lower.endswith("_id"):
                score += 15

            if "tax" in lower:
                score -= 50

            unique_candidates.append(
                (score, column)
            )

    if not unique_candidates:
        raise ValueError(
            f"Could not identify a unique {entity_label} ID column. "
            "Candidate metadata columns were: "
            + ", ".join(candidate_columns)
        )

    unique_candidates.sort(
        key=lambda item: (
            -item[0],
            candidate_columns.index(item[1]),
        )
    )

    best_score = unique_candidates[0][0]
    best = [
        column
        for score, column in unique_candidates
        if score == best_score
    ]

    if len(best) > 1 and best_score <= 0:
        raise ValueError(
            f"Multiple possible unique {entity_label} ID columns "
            f"were found: {best}. Rename the intended column to "
            f"'{entity_label}_id' or provide a simplified matrix."
        )

    return best[0]


def read_abundance_matrix(
    path: Path,
    expected_samples: int = 159,
    expected_mags: int = 1375,
) -> pd.DataFrame:
    """Read and validate a sample x MAG abundance matrix.

    The CEODOS table is allowed to contain metadata columns in addition to
    abundances. In particular, a MAG x sample table may contain a repeated
    taxonomy column plus a separate unique MAG identifier column.

    The function detects numeric abundance columns first, finds a unique ID
    column among the remaining metadata, and always returns samples x MAGs.
    """

    path = Path(path).expanduser().resolve()

    if not path.exists():
        raise FileNotFoundError(
            f"Abundance matrix not found: {path}"
        )

    if path.suffix.lower() in {".tsv", ".tab"}:
        sep, engine = "\t", "c"
    elif path.suffix.lower() == ".csv":
        sep, engine = ",", "c"
    else:
        sep, engine = None, "python"

    raw = pd.read_csv(
        path,
        sep=sep,
        engine=engine,
        header=0,
    )

    if raw.empty:
        raise ValueError(
            f"Abundance matrix is empty: {path}"
        )

    raw.columns = [
        str(value).strip()
        for value in raw.columns
    ]

    if pd.Index(raw.columns).has_duplicates:
        duplicates = (
            pd.Index(raw.columns)[
                pd.Index(raw.columns).duplicated()
            ]
            .unique()
            .tolist()
        )

        raise ValueError(
            "Duplicate column names in abundance matrix: "
            + ", ".join(
                map(
                    str,
                    duplicates[:20],
                )
            )
        )

    numeric_columns = []
    converted_columns = {}

    for column in raw.columns:
        converted = pd.to_numeric(
            raw[column],
            errors="coerce",
        )

        original_nonempty = (
            raw[column]
            .map(
                lambda value: (
                    ""
                    if pd.isna(value)
                    else str(value).strip()
                )
            )
            .ne("")
        )

        invalid = (
            original_nonempty
            & converted.isna()
        )

        if not invalid.any():
            numeric_columns.append(
                column
            )
            converted_columns[
                column
            ] = converted.fillna(0.0)

    metadata_columns = [
        column
        for column in raw.columns
        if column not in numeric_columns
    ]

    # --------------------------------------------------------------
    # CEODOS/native MAG x sample orientation.
    #
    # Expected structure:
    #
    #     1375 MAG rows
    #     159 numeric sample columns
    #     one or more metadata columns (e.g. taxonomy + MAG ID)
    #
    # Taxonomy is not required to be unique; the MAG ID column is.
    # --------------------------------------------------------------
    if (
        len(raw) == expected_mags
        and len(numeric_columns) == expected_samples
    ):
        mag_id_column = _choose_unique_id_column(
            raw=raw,
            candidate_columns=metadata_columns,
            expected_count=expected_mags,
            entity_label="mag",
        )

        mag_ids = (
            raw[mag_id_column]
            .map(
                lambda value: str(value).strip()
            )
            .tolist()
        )

        abundance = pd.DataFrame(
            {
                column: converted_columns[column]
                for column in numeric_columns
            }
        )

        abundance.index = mag_ids
        converted = abundance.T

        converted.index = [
            str(value).strip()
            for value in converted.index
        ]

        converted.columns = mag_ids

        converted.attrs[
            "detected_orientation"
        ] = "MAGs_x_samples_transposed_to_samples_x_MAGs"

        converted.attrs[
            "mag_id_column"
        ] = mag_id_column

    # --------------------------------------------------------------
    # Conventional sample x MAG orientation.
    #
    # Expected structure:
    #
    #     159 sample rows
    #     1375 numeric MAG columns
    #     one metadata column holding unique sample IDs
    # --------------------------------------------------------------
    elif (
        len(raw) == expected_samples
        and len(numeric_columns) == expected_mags
    ):
        sample_id_column = _choose_unique_id_column(
            raw=raw,
            candidate_columns=metadata_columns,
            expected_count=expected_samples,
            entity_label="sample",
        )

        sample_ids = (
            raw[sample_id_column]
            .map(
                lambda value: str(value).strip()
            )
            .tolist()
        )

        converted = pd.DataFrame(
            {
                column: converted_columns[column]
                for column in numeric_columns
            }
        )

        converted.index = sample_ids
        converted.columns = [
            str(value).strip()
            for value in numeric_columns
        ]

        converted.attrs[
            "detected_orientation"
        ] = "samples_x_MAGs"

        converted.attrs[
            "sample_id_column"
        ] = sample_id_column

    else:
        raise ValueError(
            "Unable to identify the abundance block automatically. "
            f"Observed {len(raw)} data rows, "
            f"{len(raw.columns)} total columns, "
            f"{len(numeric_columns)} fully numeric columns, and "
            f"{len(metadata_columns)} metadata columns. "
            f"Expected either {expected_mags} MAG rows with "
            f"{expected_samples} numeric sample columns, or "
            f"{expected_samples} sample rows with "
            f"{expected_mags} numeric MAG columns. "
            "Metadata columns such as taxonomy are allowed."
        )

    if converted.index.has_duplicates:
        duplicates = (
            converted.index[
                converted.index.duplicated()
            ]
            .unique()
            .tolist()
        )

        raise ValueError(
            "Duplicate sample IDs after matrix parsing: "
            + ", ".join(
                map(
                    str,
                    duplicates[:20],
                )
            )
        )

    if pd.Index(
        converted.columns
    ).has_duplicates:
        duplicates = (
            pd.Index(converted.columns)[
                pd.Index(
                    converted.columns
                ).duplicated()
            ]
            .unique()
            .tolist()
        )

        raise ValueError(
            "Duplicate MAG IDs after matrix parsing: "
            + ", ".join(
                map(
                    str,
                    duplicates[:20],
                )
            )
        )

    values = converted.to_numpy()

    if not all(
        math.isfinite(
            float(value)
        )
        for value in values.ravel()
    ):
        raise ValueError(
            "Abundance matrix contains non-finite values."
        )

    if (
        converted < 0
    ).any().any():
        raise ValueError(
            "Abundance matrix contains negative values."
        )

    return converted


def discover_model_files(models_dir: Path) -> list[Path]:
    models_dir = Path(models_dir).expanduser().resolve()
    if not models_dir.exists():
        raise FileNotFoundError(f"Models directory not found: {models_dir}")
    paths = []
    for suffix in MODEL_SUFFIXES:
        paths.extend(models_dir.rglob(f"*{suffix}"))
    paths = sorted({path.resolve() for path in paths})
    if not paths:
        raise FileNotFoundError(f"No .xml or .sbml models found below {models_dir}")
    return paths


def _model_aliases(path: Path) -> set[str]:
    name, stem = path.name.strip(), path.stem.strip()
    aliases = {name.lower(), stem.lower(), normalized_id(name), normalized_id(stem)}
    for value in (name, stem):
        lower = value.lower()
        for suffix in ("_model", "-model", ".model", "_gapseq", "-gapseq"):
            if lower.endswith(suffix):
                stripped = value[: -len(suffix)]
                aliases.add(stripped.lower())
                aliases.add(normalized_id(stripped))
    return {alias for alias in aliases if alias}


def build_model_index(model_paths: list[Path]) -> tuple[dict[str, set[Path]], dict[str, set[Path]]]:
    exact: dict[str, set[Path]] = {}
    normalized: dict[str, set[Path]] = {}
    for path in model_paths:
        for alias in _model_aliases(path):
            target = normalized if re.fullmatch(r"[a-z0-9]+", alias) else exact
            target.setdefault(alias, set()).add(path)
        exact.setdefault(path.name.lower(), set()).add(path)
        exact.setdefault(path.stem.lower(), set()).add(path)
        normalized.setdefault(normalized_id(path.name), set()).add(path)
        normalized.setdefault(normalized_id(path.stem), set()).add(path)
    return exact, normalized


def read_explicit_model_map(path: Path | None) -> dict[str, Path]:
    if path is None:
        return {}
    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Model map not found: {path}")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel_tab
        raw_rows = [row for row in csv.reader(handle, dialect) if row]

    if not raw_rows:
        return {}
    header = [cell.strip().lower() for cell in raw_rows[0]]
    if "mag_id" in header and "model_path" in header:
        mag_col, path_col, data_rows = header.index("mag_id"), header.index("model_path"), raw_rows[1:]
    else:
        mag_col, path_col, data_rows = 0, 1, raw_rows

    mapping: dict[str, Path] = {}
    for row in data_rows:
        if len(row) <= max(mag_col, path_col):
            continue
        mag_id, model_path = row[mag_col].strip(), row[path_col].strip()
        if not mag_id or not model_path:
            continue
        candidate = Path(model_path).expanduser()
        candidate = (path.parent / candidate).resolve() if not candidate.is_absolute() else candidate.resolve()
        if not candidate.exists():
            raise FileNotFoundError(f"Model map points to missing file for {mag_id}: {candidate}")
        mapping[mag_id] = candidate
    return mapping


def resolve_mag_model(mag_id: str, exact_index: dict[str, set[Path]], normalized_index: dict[str, set[Path]], explicit_map: dict[str, Path] | None = None) -> Path:
    explicit_map = explicit_map or {}
    if mag_id in explicit_map:
        return explicit_map[mag_id]

    exact_candidates: set[Path] = set()
    for key in {mag_id.lower(), f"{mag_id}.xml".lower(), f"{mag_id}.sbml".lower()}:
        exact_candidates.update(exact_index.get(key, set()))
    if len(exact_candidates) == 1:
        return next(iter(exact_candidates))
    if len(exact_candidates) > 1:
        raise ValueError(f"Ambiguous exact model match for MAG {mag_id}: " + ", ".join(str(path) for path in sorted(exact_candidates)))

    normalized_candidates = normalized_index.get(normalized_id(mag_id), set())
    if len(normalized_candidates) == 1:
        return next(iter(normalized_candidates))
    if len(normalized_candidates) > 1:
        raise ValueError(f"Ambiguous normalized model match for MAG {mag_id}: " + ", ".join(str(path) for path in sorted(normalized_candidates)))
    raise FileNotFoundError(f"No metabolic model could be matched to MAG {mag_id}")


def read_medium_spec(
    path: Path,
    default_uptake_bound: float = -1000.0,
) -> tuple[list[dict], dict]:
    """Read either the temporary gapseq CSV medium or future mathomics.txt.

    Supported formats:
    1. CSV with columns compound, name, maxFlux. A positive maxFlux is
       interpreted as uptake magnitude, so the lower bound is -abs(maxFlux).
    2. Headerless text with one metabolite identifier per row. Those entries
       use default_uptake_bound.
    """

    path = Path(path).expanduser().resolve()

    if not path.exists():
        raise FileNotFoundError(f"Medium file not found: {path}")

    if default_uptake_bound >= 0:
        raise ValueError("default_uptake_bound must be negative.")

    entries: list[dict] = []
    source_format = "headerless_one_metabolite_per_row"

    if path.suffix.lower() == ".csv":
        with path.open(
            "r",
            encoding="utf-8-sig",
            newline="",
        ) as handle:
            reader = csv.DictReader(handle)
            fieldnames = [
                str(value).strip()
                for value in (reader.fieldnames or [])
            ]
            normalized_fields = {
                value.lower(): value
                for value in fieldnames
            }

            if "compound" in normalized_fields:
                source_format = "csv_compound_name_maxFlux"

                compound_field = normalized_fields["compound"]
                name_field = normalized_fields.get("name")
                flux_field = (
                    normalized_fields.get("maxflux")
                    or normalized_fields.get("max_flux")
                )

                for line_number, row in enumerate(
                    reader,
                    start=2,
                ):
                    token = str(
                        row.get(compound_field, "")
                    ).strip()

                    if not token:
                        continue

                    name = (
                        str(row.get(name_field, "")).strip()
                        if name_field is not None
                        else ""
                    )

                    max_flux = None
                    lower_bound = float(
                        default_uptake_bound
                    )

                    if flux_field is not None:
                        raw_flux = str(
                            row.get(flux_field, "")
                        ).strip()

                        if raw_flux:
                            try:
                                max_flux = float(raw_flux)
                            except ValueError as error:
                                raise ValueError(
                                    "Invalid maxFlux value "
                                    f"{raw_flux!r} at line "
                                    f"{line_number} in {path}."
                                ) from error

                            if (
                                not math.isfinite(max_flux)
                                or max_flux <= 0
                            ):
                                raise ValueError(
                                    "maxFlux must be finite "
                                    "and > 0 at line "
                                    f"{line_number} in {path}; "
                                    f"observed {max_flux}."
                                )

                            lower_bound = -abs(max_flux)

                    entries.append(
                        {
                            "token": token,
                            "name": name or None,
                            "max_flux": max_flux,
                            "lower_bound": float(
                                lower_bound
                            ),
                        }
                    )

    if not entries:
        source_format = "headerless_one_metabolite_per_row"

        with path.open(
            "r",
            encoding="utf-8-sig",
        ) as handle:
            for line in handle:
                stripped = line.strip()

                if (
                    not stripped
                    or stripped.startswith("#")
                ):
                    continue

                if "\t" in stripped:
                    stripped = stripped.split(
                        "\t",
                        1,
                    )[0].strip()
                elif "," in stripped:
                    stripped = stripped.split(
                        ",",
                        1,
                    )[0].strip()

                if stripped:
                    entries.append(
                        {
                            "token": stripped,
                            "name": None,
                            "max_flux": None,
                            "lower_bound": float(
                                default_uptake_bound
                            ),
                        }
                    )

    if not entries:
        raise ValueError(
            f"No metabolites found in medium file: {path}"
        )

    tokens = [
        entry["token"]
        for entry in entries
    ]

    if len(tokens) != len(set(tokens)):
        duplicates = sorted(
            {
                token
                for token in tokens
                if tokens.count(token) > 1
            }
        )
        raise ValueError(
            "Duplicate compounds in medium file: "
            + ", ".join(duplicates[:20])
        )

    audit = {
        "file": str(path),
        "source_format": source_format,
        "number_entries": len(entries),
        "uses_per_compound_max_flux": any(
            entry["max_flux"] is not None
            for entry in entries
        ),
        "default_uptake_bound": float(
            default_uptake_bound
        ),
    }

    return entries, audit


def read_medium_tokens(path: Path) -> list[str]:
    """Backward-compatible token-only medium reader."""
    entries, _ = read_medium_spec(path)
    return [
        entry["token"]
        for entry in entries
    ]


def _medium_aliases_for_exchange(community, rid: str) -> set[str]:
    aliases = {rid, rid.lower()}
    if rid.startswith("R_EX_"):
        aliases.add(rid[len("R_EX_"):])
        aliases.add(rid[len("R_EX_"):].lower())
    reaction = community.merged_model.reactions[rid]
    if len(reaction.stoichiometry) == 1:
        metabolite_id = next(iter(reaction.stoichiometry))
        aliases.update({metabolite_id, metabolite_id.lower()})
        if metabolite_id.startswith("M_"):
            aliases.update({metabolite_id[2:], metabolite_id[2:].lower()})
        metabolite = community.merged_model.metabolites.get(metabolite_id)
        if metabolite is not None and getattr(metabolite, "name", None):
            aliases.update({str(metabolite.name), str(metabolite.name).lower()})
    aliases.update({normalized_id(alias) for alias in aliases if alias})
    return {alias for alias in aliases if alias}


def resolve_medium_for_community(
    community,
    entries,
    uptake_bound: float = -1000.0,
) -> tuple[dict[str, float], dict]:
    """Resolve medium entries against global community exchanges."""

    if uptake_bound >= 0:
        raise ValueError(
            "uptake_bound must be negative."
        )

    global_exchanges = [
        rid
        for rid in community.merged_model.reactions
        if (
            rid.startswith("R_EX_")
            and not rid.endswith("_i")
        )
    ]

    exact_index: dict[str, set[str]] = {}
    normalized_index: dict[str, set[str]] = {}

    for rid in global_exchanges:
        for alias in _medium_aliases_for_exchange(
            community,
            rid,
        ):
            exact_index.setdefault(
                alias.lower(),
                set(),
            ).add(rid)

            normalized_index.setdefault(
                normalized_id(alias),
                set(),
            ).add(rid)

    medium: dict[str, float] = {}
    matched: dict[str, dict] = {}
    ambiguous: dict[str, list[str]] = {}
    unresolved: list[str] = []

    normalized_entries = []

    for entry in entries:
        if isinstance(entry, str):
            normalized_entries.append(
                {
                    "token": entry,
                    "name": None,
                    "max_flux": None,
                    "lower_bound": float(
                        uptake_bound
                    ),
                }
            )
        else:
            normalized_entries.append(
                {
                    "token": str(
                        entry["token"]
                    ).strip(),
                    "name": entry.get("name"),
                    "max_flux": entry.get(
                        "max_flux"
                    ),
                    "lower_bound": float(
                        entry.get(
                            "lower_bound",
                            uptake_bound,
                        )
                    ),
                }
            )

    for entry in normalized_entries:
        token = entry["token"]

        if not token:
            continue

        lower_bound = float(
            entry["lower_bound"]
        )

        if lower_bound >= 0:
            raise ValueError(
                f"Medium entry {token} has "
                "non-negative uptake lower bound "
                f"{lower_bound}."
            )

        candidates: set[str] = set()

        direct_forms = {
            token,
            f"R_EX_{token}",
            f"R_EX_{token}_e",
            f"R_EX_{token}_e0",
        }

        if token.startswith("EX_"):
            direct_forms.add(
                f"R_{token}"
            )

        if token.startswith("M_"):
            direct_forms.add(
                f"R_EX_{token[2:]}"
            )

        for form in direct_forms:
            if (
                form
                in community.merged_model.reactions
                and form.startswith("R_EX_")
                and not form.endswith("_i")
            ):
                candidates.add(form)

            candidates.update(
                exact_index.get(
                    form.lower(),
                    set(),
                )
            )

        if not candidates:
            candidates.update(
                normalized_index.get(
                    normalized_id(token),
                    set(),
                )
            )

        if len(candidates) == 1:
            rid = next(iter(candidates))

            previous_bound = medium.get(rid)
            if (
                previous_bound is not None
                and abs(
                    previous_bound
                    - lower_bound
                ) > 1e-12
            ):
                raise ValueError(
                    "Two medium entries resolve "
                    f"to {rid} with different "
                    "bounds: "
                    f"{previous_bound} and "
                    f"{lower_bound}."
                )

            medium[rid] = lower_bound

            matched[token] = {
                "exchange_reaction": rid,
                "lower_bound": lower_bound,
                "name": entry.get("name"),
                "max_flux": entry.get(
                    "max_flux"
                ),
            }

        elif len(candidates) > 1:
            ambiguous[token] = sorted(
                candidates
            )

        else:
            unresolved.append(token)

    if ambiguous:
        formatted = "; ".join(
            f"{token} -> {values}"
            for token, values
            in list(
                ambiguous.items()
            )[:20]
        )

        raise ValueError(
            "Ambiguous medium-metabolite "
            "mappings detected: "
            + formatted
        )

    if not medium:
        raise ValueError(
            "None of the medium metabolites "
            "matched a global exchange reaction "
            "in this sample community."
        )

    audit = {
        "number_tokens": len(
            normalized_entries
        ),
        "number_matched_tokens": len(
            matched
        ),
        "number_unresolved_tokens": len(
            unresolved
        ),
        "matched": matched,
        "unresolved": unresolved,
        "fallback_uptake_bound": float(
            uptake_bound
        ),
    }

    return medium, audit


def yaml_dump_atomic(data: dict, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False)
    temporary.replace(path)


def yaml_load(path: Path) -> dict:
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)
