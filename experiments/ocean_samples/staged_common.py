"""Shared helpers for staged ocean-sample execution."""

from __future__ import annotations

import json
import traceback
from pathlib import Path
from time import perf_counter

from reframed.solvers.solver import Parameter

from common import (
    read_medium_spec,
    resolve_medium_for_community,
    yaml_dump_atomic,
    yaml_load,
)
from misosoup.library.readwrite import load_models


def solver_params(prod_tol: float, int_tol: float) -> dict:
    return {
        Parameter.OPTIMALITY_TOL: prod_tol,
        Parameter.FEASIBILITY_TOL: prod_tol,
        Parameter.INT_FEASIBILITY_TOL: int_tol,
    }


def get_sample(manifest_path: Path, sample_index: int) -> tuple[dict, dict]:
    manifest = yaml_load(manifest_path.expanduser().resolve())
    samples = manifest["samples"]
    if not 0 <= sample_index < len(samples):
        raise IndexError(
            f"sample index {sample_index} outside 0..{len(samples) - 1}"
        )
    return manifest, samples[sample_index]


def sample_dir(output_root: Path, sample_index: int) -> Path:
    path = output_root.expanduser().resolve() / "samples" / f"{sample_index:03d}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_sample_models(sample: dict):
    paths = [
        Path(entry["model_path"]).expanduser().resolve()
        for entry in sample["mags"]
    ]
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing model files: " + ", ".join(missing[:20])
        )

    start = perf_counter()
    models = load_models([str(path) for path in paths])
    load_seconds = perf_counter() - start

    if len(models) != len(sample["mags"]):
        raise RuntimeError("Loaded-model count differs from manifest")

    solver_to_mag = {}
    mag_to_solver = {}

    for position, (model, entry) in enumerate(zip(models, sample["mags"])):
        solver_id = f"org_{position:05d}"
        mag_id = str(entry["mag_id"])
        model.id = solver_id
        solver_to_mag[solver_id] = mag_id
        mag_to_solver[mag_id] = solver_id

    return models, solver_to_mag, mag_to_solver, load_seconds


def resolve_medium(
    community,
    medium_file: Path,
    uptake_bound: float,
):
    entries, file_audit = read_medium_spec(
        medium_file.expanduser().resolve(),
        default_uptake_bound=uptake_bound,
    )
    medium, audit = resolve_medium_for_community(
        community,
        entries,
        uptake_bound,
    )
    return medium, audit, file_audit


def write_status(path: Path, **payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    temporary.replace(path)


def mark_complete(stage_dir: Path) -> None:
    (stage_dir / "COMPLETE").write_text("ok\n", encoding="utf-8")


def write_stage_failure(
    *,
    stage_name: str,
    stage_dir: Path,
    sample_index: int,
    sample_id: str,
    error: Exception,
    elapsed_seconds: float,
    sample_root: Path | None = None,
) -> None:
    failure = {
        "status": "failed",
        "stage": stage_name,
        "sample": {
            "index": int(sample_index),
            "sample_id": str(sample_id),
        },
        "error_type": type(error).__name__,
        "error": str(error),
        "timing_seconds": {
            "total": float(elapsed_seconds),
        },
        "traceback": traceback.format_exc(),
    }
    yaml_dump_atomic(failure, stage_dir / "failure.yaml")
    if sample_root is not None:
        yaml_dump_atomic(failure, sample_root / "failure.yaml")
