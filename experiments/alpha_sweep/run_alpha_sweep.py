"""Sensitivity analysis for the simultaneous-product workflow over alpha.

The individual product maxima M_i and the C/N/S/P filter are computed once,
because they do not depend on alpha. For each alpha, Stages A-D are rebuilt
with fresh LayeredCommunity solver objects.

Outputs:
  alpha_sweep_results.yaml
  alpha_sweep_summary.csv

If one alpha fails, the failure is recorded and the sweep continues. The
script exits non-zero at the end if any alpha failed, while preserving partial
results and any IIS files produced by the stage solvers.
"""

from __future__ import annotations

import argparse
import csv
import gc
import os
import traceback
from pathlib import Path
from time import perf_counter

import yaml

from misosoup.cobra.solver import Status, Parameter

from misosoup.library.product_filter import filter_product_candidates
from misosoup.library.product_reference import (
    constrain_full_community,
    find_producible_exchanges,
)
from misosoup.library.minimal_product_communities import (
    findMinimalProductCommunities,
)
from misosoup.library.product_selection import (
    getMaxProduct,
    getSelectedProductsFromProductMaximization,
    maximizeProducts,
    solvepFBAUsingFixProducts,
)
from misosoup.library.readwrite import load_models, read_compounds
from misosoup.cobra.layered_community import LayeredCommunity


ROOT = Path(__file__).resolve().parents[2]

MODEL_DIR = ROOT / "examples" / "marine" / "strains"
MEDIA_FILE = ROOT / "examples" / "marine" / "media.yaml"

MODEL_FILES = [
    MODEL_DIR / "A1R12.xml",
    MODEL_DIR / "I2R16.xml",
    MODEL_DIR / "I3M07.xml",
]

MEDIUM_NAME = "ac"

DEFAULT_ALPHAS = [
    0.01,
    0.025,
    0.05,
    0.10,
    0.15,
    0.20,
    0.30,
    0.40,
    0.50,
]

MINIMAL_GROWTH = 0.01
PRODUCTION_TOLERANCE = 1e-6
INTEGER_TOLERANCE = 1e-9
LEXICOGRAPHIC_TOLERANCE = 1e-5
KEEP_OXYGEN = True
MIN_COMMUNITY_PRODUCT_RETENTION = 0.90
MAX_MINIMAL_COMMUNITIES = 100


def elapsed_seconds(start: float) -> float:
    return perf_counter() - start


def alpha_tag(alpha: float) -> str:
    return f"{alpha:g}".replace(".", "p")


def solver_params() -> dict:
    return {
        Parameter.OPTIMALITY_TOL: PRODUCTION_TOLERANCE,
        Parameter.FEASIBILITY_TOL: PRODUCTION_TOLERANCE,
        Parameter.INT_FEASIBILITY_TOL: INTEGER_TOLERANCE,
    }


def load_medium() -> dict:
    media_data = read_compounds(str(MEDIA_FILE))
    base_medium = media_data.get("base_medium", {})
    selected_medium = media_data.get(MEDIUM_NAME)

    if selected_medium is None:
        raise KeyError(
            f"Medium '{MEDIUM_NAME}' not found in {MEDIA_FILE}"
        )

    return {
        **base_medium,
        **selected_medium,
    }


def serializable_stage_b_info(product_selection: dict) -> dict:
    return {
        rid: {
            key: float(value)
            for key, value in info.items()
        }
        for rid, info in product_selection[
            "selected_product_info"
        ].items()
    }


def serializable_reference_products(stage_c: dict) -> dict:
    return {
        rid: {
            key: float(value)
            for key, value in info.items()
        }
        for rid, info in stage_c["reference_products"].items()
    }


def write_outputs(output_dir: Path, output: dict) -> None:
    yaml_path = output_dir / "alpha_sweep_results.yaml"
    csv_path = output_dir / "alpha_sweep_summary.csv"

    with yaml_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(
            output,
            handle,
            sort_keys=False,
        )

    fieldnames = [
        "alpha",
        "status",
        "k_star",
        "q_star",
        "q_star_from_values",
        "q_star_difference",
        "number_selected_products",
        "selected_products",
        "pfba_objective",
        "stage_c_q_sum",
        "minimum_size",
        "number_communities",
        "enumeration_complete",
        "communities",
        "stage_a_seconds",
        "stage_b1_seconds",
        "stage_b2_seconds",
        "stage_c_seconds",
        "stage_d_seconds",
        "alpha_total_seconds",
        "error_type",
        "error",
    ]

    with csv_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()

        for result in output["results"]:
            row = {
                field: ""
                for field in fieldnames
            }
            row["alpha"] = result["alpha"]
            row["status"] = result["status"]

            if result["status"] == "ok":
                stage_a = result["stage_a"]
                stage_b = result["stage_b"]
                stage_c = result["stage_c"]
                stage_d = result["stage_d"]
                timing = result["timing_seconds"]

                row.update(
                    {
                        "k_star": stage_a["k_star"],
                        "q_star": stage_b["q_star"],
                        "q_star_from_values": stage_b[
                            "q_star_from_values"
                        ],
                        "q_star_difference": stage_b[
                            "q_star_difference"
                        ],
                        "number_selected_products": len(
                            stage_b["selected_products"]
                        ),
                        "selected_products": ";".join(
                            stage_b["selected_products"]
                        ),
                        "pfba_objective": stage_c[
                            "pfba_objective"
                        ],
                        "stage_c_q_sum": stage_c[
                            "stage_c_q_sum"
                        ],
                        "minimum_size": stage_d[
                            "minimum_size"
                        ],
                        "number_communities": stage_d[
                            "number_communities"
                        ],
                        "enumeration_complete": stage_d[
                            "enumeration_complete"
                        ],
                        "communities": ";".join(
                            "+".join(item["organisms"])
                            for item in stage_d[
                                "communities"
                            ]
                        ),
                        "stage_a_seconds": timing[
                            "stage_a"
                        ],
                        "stage_b1_seconds": timing[
                            "stage_b1"
                        ],
                        "stage_b2_seconds": timing[
                            "stage_b2"
                        ],
                        "stage_c_seconds": timing[
                            "stage_c"
                        ],
                        "stage_d_seconds": timing[
                            "stage_d"
                        ],
                        "alpha_total_seconds": timing[
                            "alpha_total"
                        ],
                    }
                )
            else:
                row["error_type"] = result.get(
                    "error_type",
                    "",
                )
                row["error"] = result.get(
                    "error",
                    "",
                )

            writer.writerow(row)


def run_one_alpha(
    alpha: float,
    models: list,
    medium: dict,
    max_secretion: dict,
    alpha_dir: Path,
) -> dict:
    run_start = perf_counter()
    old_cwd = Path.cwd()

    alpha_dir.mkdir(
        parents=True,
        exist_ok=True,
    )
    os.chdir(alpha_dir)

    try:
        # A/B: fresh solver for this alpha.
        ab_community = LayeredCommunity(
            f"alpha_{alpha_tag(alpha)}_ab",
            models,
            copy_models=False,
            params=solver_params(),
        )

        constrain_full_community(
            ab_community,
            minimal_growth=MINIMAL_GROWTH,
        )
        ab_community.setup_medium(medium)

        feasibility = ab_community.check_feasibility(
            ["community_growth"]
        )

        if feasibility.status != Status.OPTIMAL:
            raise RuntimeError(
                "Full community is not feasible at "
                f"alpha={alpha}. Solver status: "
                f"{feasibility.status}"
            )

        stage_a_start = perf_counter()
        stage_a = getMaxProduct(
            community=ab_community,
            max_secretion=max_secretion,
            fraction=alpha,
            tolerance=PRODUCTION_TOLERANCE,
            integer_tolerance=INTEGER_TOLERANCE,
        )
        stage_a_time = elapsed_seconds(stage_a_start)

        stage_b_start = perf_counter()
        stage_b = maximizeProducts(
            community=ab_community,
            stage_a=stage_a,
        )
        stage_b_time = elapsed_seconds(stage_b_start)

        stage_b2_start = perf_counter()
        product_selection = (
            getSelectedProductsFromProductMaximization(
                stage_b=stage_b,
            )
        )
        stage_b2_time = elapsed_seconds(stage_b2_start)

        # C: fresh solver.
        stage_c_start = perf_counter()
        stage_c_community = LayeredCommunity(
            f"alpha_{alpha_tag(alpha)}_stage_c",
            models,
            copy_models=False,
            params=solver_params(),
        )

        stage_c = solvepFBAUsingFixProducts(
            community=stage_c_community,
            stage_b=stage_b,
            product_selection=product_selection,
            medium=medium,
            minimal_growth=MINIMAL_GROWTH,
            lexicographic_tolerance=(
                LEXICOGRAPHIC_TOLERANCE
            ),
            tolerance=PRODUCTION_TOLERANCE,
            check_feasibility=True,
        )
        stage_c_time = elapsed_seconds(stage_c_start)

        if set(stage_c["selected_products"]) != set(
            product_selection["selected_products"]
        ):
            raise RuntimeError(
                "Stage C product set differs from the "
                "Stage B selected set."
            )

        # D: fresh solver.
        stage_d_start = perf_counter()
        minimal_community = LayeredCommunity(
            f"alpha_{alpha_tag(alpha)}_stage_d",
            models,
            copy_models=False,
            params=solver_params(),
        )

        stage_d = findMinimalProductCommunities(
            community=minimal_community,
            medium=medium,
            reference_products=stage_c[
                "reference_products"
            ],
            minimal_growth=MINIMAL_GROWTH,
            product_retention=(
                MIN_COMMUNITY_PRODUCT_RETENTION
            ),
            tolerance=PRODUCTION_TOLERANCE,
            max_communities=MAX_MINIMAL_COMMUNITIES,
        )
        stage_d_time = elapsed_seconds(stage_d_start)

        return {
            "alpha": float(alpha),
            "status": "ok",
            "stage_a": {
                "k_star": int(stage_a["k_star"]),
                "selected_products": list(
                    stage_a[
                        "stage_a_selected_products"
                    ]
                ),
            },
            "stage_b": {
                "q_star": float(stage_b["q_star"]),
                "q_star_from_values": float(
                    stage_b["q_star_from_values"]
                ),
                "q_star_difference": float(
                    stage_b["q_star_difference"]
                ),
                "selected_products": list(
                    product_selection[
                        "selected_products"
                    ]
                ),
                "selected_product_info": (
                    serializable_stage_b_info(
                        product_selection
                    )
                ),
            },
            "stage_c": {
                "selected_products": list(
                    stage_c["selected_products"]
                ),
                "stage_b_objective_floor": float(
                    stage_c[
                        "stage_b_objective_floor"
                    ]
                ),
                "stage_c_q_sum": float(
                    stage_c["stage_c_q_sum"]
                ),
                "lexicographic_tolerance": float(
                    stage_c[
                        "lexicographic_tolerance"
                    ]
                ),
                "pfba_objective": float(
                    stage_c["pfba_objective"]
                ),
                "reference_products": (
                    serializable_reference_products(
                        stage_c
                    )
                ),
            },
            "stage_d": {
                "minimum_size": int(
                    stage_d["minimum_size"]
                ),
                "product_retention": float(
                    stage_d["product_retention"]
                ),
                "product_requirements": {
                    rid: float(value)
                    for rid, value in stage_d[
                        "product_requirements"
                    ].items()
                },
                "number_communities": int(
                    stage_d["number_communities"]
                ),
                "enumeration_complete": bool(
                    stage_d[
                        "enumeration_complete"
                    ]
                ),
                "communities": stage_d[
                    "communities"
                ],
            },
            "timing_seconds": {
                "stage_a": stage_a_time,
                "stage_b1": stage_b_time,
                "stage_b2": stage_b2_time,
                "stage_c": stage_c_time,
                "stage_d": stage_d_time,
                "alpha_total": elapsed_seconds(
                    run_start
                ),
            },
        }

    finally:
        os.chdir(old_cwd)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep the simultaneous-product workflow "
            "over Stage-A alpha values."
        )
    )
    parser.add_argument(
        "--alphas",
        nargs="+",
        type=float,
        default=DEFAULT_ALPHAS,
        help=(
            "Alpha values to test. Default: "
            + " ".join(
                str(value)
                for value in DEFAULT_ALPHAS
            )
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path.cwd() / "alpha_sweep_results",
        help=(
            "Directory for YAML, CSV, and per-alpha "
            "diagnostic files."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    alphas = [
        float(value)
        for value in args.alphas
    ]

    if not alphas:
        raise ValueError(
            "At least one alpha value is required."
        )

    for alpha in alphas:
        if not 0 < alpha <= 1:
            raise ValueError(
                f"Invalid alpha={alpha}; alpha must "
                "be >0 and <=1."
            )

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    total_start = perf_counter()

    print("=" * 72)
    print("MiSoSoup alpha sensitivity sweep")
    print(f"Repository root : {ROOT}")
    print(f"Output directory: {output_dir}")
    print(f"Alphas          : {alphas}")
    print("=" * 72)
    print(flush=True)

    models = load_models(
        [str(path) for path in MODEL_FILES]
    )
    medium = load_medium()

    # M_i and the filter are alpha-independent, so calculate them once.
    reference_start = perf_counter()

    reference_community = LayeredCommunity(
        "alpha_sweep_reference",
        models,
        copy_models=False,
        params=solver_params(),
    )
    constrain_full_community(
        reference_community,
        minimal_growth=MINIMAL_GROWTH,
    )
    reference_community.setup_medium(medium)

    feasibility = reference_community.check_feasibility(
        ["community_growth"]
    )

    if feasibility.status != Status.OPTIMAL:
        raise RuntimeError(
            "Reference full community is not feasible. "
            f"Solver status: {feasibility.status}"
        )

    unfiltered_max_secretion = (
        find_producible_exchanges(
            reference_community,
            tolerance=PRODUCTION_TOLERANCE,
        )
    )

    (
        max_secretion,
        product_filter_audit,
    ) = filter_product_candidates(
        community=reference_community,
        max_secretion=unfiltered_max_secretion,
        keep_oxygen=KEEP_OXYGEN,
    )

    reference_time = elapsed_seconds(
        reference_start
    )

    output = {
        "configuration": {
            "models": [
                path.name
                for path in MODEL_FILES
            ],
            "medium": MEDIUM_NAME,
            "minimal_growth": MINIMAL_GROWTH,
            "production_tolerance": (
                PRODUCTION_TOLERANCE
            ),
            "integer_tolerance": (
                INTEGER_TOLERANCE
            ),
            "lexicographic_tolerance": (
                LEXICOGRAPHIC_TOLERANCE
            ),
            "keep_oxygen": KEEP_OXYGEN,
            "min_community_product_retention": (
                MIN_COMMUNITY_PRODUCT_RETENTION
            ),
            "max_minimal_communities": (
                MAX_MINIMAL_COMMUNITIES
            ),
            "alphas": alphas,
        },
        "reference_product_scan": {
            "number_unfiltered_producible_exchanges": len(
                unfiltered_max_secretion
            ),
            "number_filtered_products": len(
                max_secretion
            ),
            "max_secretion": {
                rid: float(value)
                for rid, value in max_secretion.items()
            },
            "filter_audit": product_filter_audit,
            "timing_seconds": reference_time,
        },
        "results": [],
    }

    write_outputs(
        output_dir,
        output,
    )

    print(
        "Reference product scan: "
        f"{len(unfiltered_max_secretion)} producible -> "
        f"{len(max_secretion)} filtered candidates "
        f"in {reference_time:.3f} s",
        flush=True,
    )

    failures = []

    del reference_community
    gc.collect()

    for index, alpha in enumerate(
        alphas,
        start=1,
    ):
        print()
        print("=" * 72)
        print(
            f"[{index}/{len(alphas)}] alpha = {alpha:g}"
        )
        print("=" * 72)
        print(flush=True)

        alpha_dir = (
            output_dir
            / f"alpha_{alpha_tag(alpha)}"
        )

        try:
            result = run_one_alpha(
                alpha=alpha,
                models=models,
                medium=medium,
                max_secretion=max_secretion,
                alpha_dir=alpha_dir,
            )

            print(
                "Result: "
                f"K*={result['stage_a']['k_star']}, "
                f"Q*={result['stage_b']['q_star']:.12g}, "
                f"pFBA="
                f"{result['stage_c']['pfba_objective']:.12g}, "
                f"N_min="
                f"{result['stage_d']['minimum_size']}, "
                f"communities="
                f"{result['stage_d']['number_communities']}",
                flush=True,
            )

        except Exception as error:
            traceback.print_exc()

            result = {
                "alpha": float(alpha),
                "status": "failed",
                "error_type": type(error).__name__,
                "error": str(error),
            }

            failures.append(
                {
                    "alpha": float(alpha),
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
            )

        output["results"].append(result)

        # Save a checkpoint after every alpha.
        write_outputs(
            output_dir,
            output,
        )
        gc.collect()

    output["timing_seconds"] = {
        "reference_product_scan": reference_time,
        "total": elapsed_seconds(total_start),
    }
    output["number_successful"] = sum(
        result["status"] == "ok"
        for result in output["results"]
    )
    output["number_failed"] = len(failures)

    write_outputs(
        output_dir,
        output,
    )

    print()
    print("=" * 72)
    print("ALPHA SWEEP COMPLETE")
    print("=" * 72)
    print(
        f"Successful: "
        f"{output['number_successful']} / {len(alphas)}"
    )
    print(
        f"Failed    : "
        f"{output['number_failed']} / {len(alphas)}"
    )
    print(
        f"YAML      : "
        f"{output_dir / 'alpha_sweep_results.yaml'}"
    )
    print(
        f"CSV       : "
        f"{output_dir / 'alpha_sweep_summary.csv'}"
    )
    print(
        "Total time: "
        f"{output['timing_seconds']['total']:.3f} s"
    )
    print(flush=True)

    if failures:
        details = "; ".join(
            (
                f"alpha={item['alpha']}: "
                f"{item['error_type']}: "
                f"{item['error']}"
            )
            for item in failures
        )
        raise RuntimeError(
            "One or more alpha values failed. "
            f"Partial results were saved. {details}"
        )


if __name__ == "__main__":
    main()
