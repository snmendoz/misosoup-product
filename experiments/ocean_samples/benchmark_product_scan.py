"""Benchmark the pure-LP product scan on a small deterministic exchange subset."""

from __future__ import annotations

import argparse
from pathlib import Path
from time import perf_counter

from reframed.solvers.solution import Status
from reframed.solvers.solver import Parameter

from common import (
    read_medium_spec,
    resolve_medium_for_community,
    yaml_load,
)
from misosoup.library.product_reference import (
    constrain_full_community_lp,
    find_producible_exchanges,
)
from misosoup.library.readwrite import load_models
from misosoup.reframed.layered_community import LayeredCommunity


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--sample-index", type=int, default=0)
    parser.add_argument("--medium-file", required=True, type=Path)
    parser.add_argument("--exchanges", type=int, default=10)
    parser.add_argument("--minimal-growth", type=float, default=0.01)
    parser.add_argument("--production-tolerance", type=float, default=1e-6)
    parser.add_argument("--integer-tolerance", type=float, default=1e-9)
    parser.add_argument("--uptake-bound", type=float, default=-1000.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.exchanges < 1:
        raise ValueError("--exchanges must be >= 1")
    if args.minimal_growth <= 0:
        raise ValueError("--minimal-growth must be > 0")
    if args.uptake_bound >= 0:
        raise ValueError("--uptake-bound must be negative")

    manifest = yaml_load(args.manifest.expanduser().resolve())
    samples = manifest["samples"]

    if not 0 <= args.sample_index < len(samples):
        raise IndexError(
            f"sample index {args.sample_index} outside 0..{len(samples) - 1}"
        )

    sample = samples[args.sample_index]
    model_paths = [
        Path(entry["model_path"]).expanduser().resolve()
        for entry in sample["mags"]
    ]

    missing = [str(path) for path in model_paths if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing model files: " + ", ".join(missing[:20])
        )

    print("=" * 78, flush=True)
    print(
        f"PURE-LP PRODUCT-SCAN BENCHMARK: "
        f"sample={sample['sample_id']} "
        f"models={len(model_paths)} "
        f"exchanges={args.exchanges}",
        flush=True,
    )
    print("=" * 78, flush=True)

    total_start = perf_counter()

    load_start = perf_counter()
    models = load_models([str(path) for path in model_paths])
    print(
        f"Loaded {len(models)} models in "
        f"{perf_counter() - load_start:.3f}s.",
        flush=True,
    )

    for position, model in enumerate(models):
        model.id = f"org_{position:05d}"

    params = {
        Parameter.OPTIMALITY_TOL: args.production_tolerance,
        Parameter.FEASIBILITY_TOL: args.production_tolerance,
        Parameter.INT_FEASIBILITY_TOL: args.integer_tolerance,
    }

    build_start = perf_counter()
    community = LayeredCommunity(
        f"benchmark_product_scan_{args.sample_index:03d}",
        models,
        copy_models=False,
        params=params,
    )
    print(
        f"Built merged community in "
        f"{perf_counter() - build_start:.3f}s.",
        flush=True,
    )

    constrain_full_community_lp(
        community,
        minimal_growth=args.minimal_growth,
    )

    medium_entries, _ = read_medium_spec(
        args.medium_file.expanduser().resolve(),
        default_uptake_bound=args.uptake_bound,
    )
    medium, medium_audit = resolve_medium_for_community(
        community,
        medium_entries,
        args.uptake_bound,
    )
    community.setup_medium(medium)

    print(
        f"Medium matched {medium_audit['number_matched_tokens']} / "
        f"{medium_audit['number_tokens']} tokens.",
        flush=True,
    )

    feasibility_start = perf_counter()
    feasibility = community.check_feasibility(["community_growth"])
    feasibility_time = perf_counter() - feasibility_start

    print(
        f"Full-community LP feasibility: {feasibility.status} "
        f"in {feasibility_time:.3f}s.",
        flush=True,
    )

    if feasibility.status != Status.OPTIMAL:
        raise RuntimeError(
            "Benchmark community is not feasible under the selected medium."
        )

    scan_start = perf_counter()
    producible = find_producible_exchanges(
        community,
        tolerance=args.production_tolerance,
        max_exchanges=args.exchanges,
    )
    scan_time = perf_counter() - scan_start

    print("=" * 78, flush=True)
    print(
        f"BENCHMARK COMPLETE: {args.exchanges} exchanges in "
        f"{scan_time:.3f}s; average={scan_time / args.exchanges:.3f}s/exchange; "
        f"producible={len(producible)}; "
        f"total_wall={perf_counter() - total_start:.3f}s.",
        flush=True,
    )
    print("=" * 78, flush=True)


if __name__ == "__main__":
    main()
