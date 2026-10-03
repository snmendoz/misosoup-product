"""Profile the serial pre-solve path of one ocean Product Scan sample.

This diagnostic does not write Product Scan checkpoints and does not use Gurobi.
It times model loading one file at a time, community merge, medium resolution,
LP construction, and the initial HiGHS feasibility solve.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from time import perf_counter

from reframed.io.sbml import load_cbmodel

from common import yaml_load
from misosoup.library.highs_product_scan import (
    build_product_scan_lp,
    check_product_scan_feasibility,
)
from misosoup.library.product_filter import filter_exchange_candidates
from misosoup.library.product_reference import get_community_exchanges
from misosoup.reframed.layered_community import LayeredCommunity
from staged_common import get_sample, resolve_medium


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--sample-index", required=True, type=int)
    parser.add_argument("--medium-file", required=True, type=Path)
    parser.add_argument("--minimal-growth", type=float, default=0.01)
    parser.add_argument("--production-tolerance", type=float, default=1e-6)
    parser.add_argument("--uptake-bound", type=float, default=-1000.0)
    parser.add_argument(
        "--stop-after",
        choices=("load", "merge", "lp", "feasibility"),
        default="feasibility",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    _, sample = get_sample(args.manifest, args.sample_index)
    print(
        f"DIAGNOSTIC sample={sample['sample_id']} "
        f"index={args.sample_index} MAGs={sample['number_modeled_mags']}",
        flush=True,
    )

    models = []
    load_total = perf_counter()
    slowest = []
    for position, entry in enumerate(sample["mags"], start=1):
        path = Path(entry["model_path"]).expanduser().resolve()
        start = perf_counter()
        model = load_cbmodel(str(path), flavor="fbc2")
        elapsed = perf_counter() - start
        model.id = f"org_{position - 1:05d}"
        models.append(model)
        slowest.append((elapsed, str(entry["mag_id"]), str(path)))
        print(
            f"LOAD {position}/{len(sample['mags'])} "
            f"{entry['mag_id']} {elapsed:.3f}s "
            f"rxns={len(model.reactions)} mets={len(model.metabolites)} "
            f"genes={len(model.genes)}",
            flush=True,
        )
    total_load = perf_counter() - load_total
    print(f"LOAD COMPLETE total={total_load:.3f}s", flush=True)
    for elapsed, mag_id, path in sorted(slowest, reverse=True)[:10]:
        print(
            f"SLOW LOAD {mag_id} {elapsed:.3f}s {path}",
            flush=True,
        )
    if args.stop_after == "load":
        return

    start = perf_counter()
    community = LayeredCommunity(
        f"diagnostic_sample_{args.sample_index:03d}",
        models,
        copy_models=False,
        create_solver=False,
    )
    elapsed = perf_counter() - start
    print(
        f"MERGE COMPLETE {elapsed:.3f}s "
        f"rxns={len(community.merged_model.reactions)} "
        f"mets={len(community.merged_model.metabolites)} "
        f"genes={len(community.merged_model.genes)}",
        flush=True,
    )
    if args.stop_after == "merge":
        return

    start = perf_counter()
    medium, _, _ = resolve_medium(
        community,
        args.medium_file,
        args.uptake_bound,
    )
    print(
        f"MEDIUM COMPLETE {perf_counter() - start:.3f}s "
        f"mapped={len(medium)}",
        flush=True,
    )

    exchanges = sorted(get_community_exchanges(community))
    candidates, _ = filter_exchange_candidates(
        community=community,
        exchange_reactions=exchanges,
        keep_oxygen=True,
    )

    start = perf_counter()
    lp = build_product_scan_lp(
        community,
        minimal_growth=args.minimal_growth,
        medium=medium,
    )
    elapsed = perf_counter() - start
    print(
        f"LP COMPLETE {elapsed:.3f}s "
        f"vars={len(lp.reaction_ids)} rows={lp.a_eq.shape[0]} "
        f"nnz={lp.a_eq.nnz} candidates={len(candidates)}",
        flush=True,
    )
    if args.stop_after == "lp":
        return

    start = perf_counter()
    result = check_product_scan_feasibility(
        lp,
        tolerance=args.production_tolerance,
    )
    print(
        f"FEASIBILITY COMPLETE {perf_counter() - start:.3f}s "
        f"optimal={result['optimal']} status={result['status']} "
        f"message={result['message']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
