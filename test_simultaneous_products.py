"""Integration test for the A -> B.1 -> B.2 -> C product workflow.

Pipeline
========

1. Build the full microbial community.
2. Require every organism to grow at least MINIMAL_GROWTH.
3. Determine M_i for each individually producible community exchange.
4. Stage A:
       K* = max sum_i z_i
5. Stage B.1:
       Q* = max sum_i q_i, subject to sum_i z_i = K*
6. Stage B.2:
       Extract the selected product set S_B.
7. Stage C:
       Fix S_B, preserve the Stage-B objective:

           sum_i q_i >= Q* - epsilon

       and solve classical split-flux pFBA:

           v_r = v_r^+ - v_r^-

           min sum_r(v_r^+ + v_r^-)

8. Stage D:
       Build a fresh community and solve:

           N_min = min sum_j y_j

       subject to preserving every selected Stage-C product:

           v_i >= beta * v_i^ref

       Then enumerate alternative communities with exactly N_min organisms.

The final output is written to YAML without serializing solver objects.
"""

from pathlib import Path
from time import perf_counter

import yaml

from reframed.solvers.solution import Status
from reframed.solvers.solver import Parameter

from misosoup.library.product_reference import (
    constrain_full_community,
    find_producible_exchanges,
)
from misosoup.library.product_filter import (
    filter_product_candidates,
)
from misosoup.library.minimal_product_communities import (
    findMinimalProductCommunities,
)
from misosoup.library.product_selection import (
    getMaxProduct,
    maximizeProducts,
    getSelectedProductsFromProductMaximization,
    solvepFBAUsingFixProducts,
)
from misosoup.library.readwrite import (
    load_models,
    read_compounds,
)
from misosoup.reframed.layered_community import LayeredCommunity


# ======================================================================
# CONFIGURATION
# ======================================================================

ROOT = Path(__file__).resolve().parent

MODEL_DIR = ROOT / "examples" / "marine" / "strains"
MEDIA_FILE = ROOT / "examples" / "marine" / "media.yaml"

MODEL_FILES = [
    MODEL_DIR / "A1R12.xml",
    MODEL_DIR / "I2R16.xml",
    MODEL_DIR / "I3M07.xml",
]

MEDIUM_NAME = "ac"

# Every organism in the full reference community must satisfy:
#
#     growth_j >= MINIMAL_GROWTH
MINIMAL_GROWTH = 0.01

# Numerical tolerance used for feasibility/optimality and for deciding
# whether an exchange has a positive individual maximum M_i.
PRODUCTION_TOLERANCE = 1e-6

# Strict integrality tolerance for z_i.
INTEGER_TOLERANCE = 1e-9

# Stage A threshold:
#
#     z_i = 1  =>  v_i >= alpha M_i
PRODUCTION_FRACTION = 0.20

# Stage C preserves the Stage-B objective lexicographically:
#
#     sum_i q_i >= Q* - epsilon
#
# A small epsilon avoids demanding the numerically exact MIP optimum.
LEXICOGRAPHIC_TOLERANCE = 1e-5

# Keep molecular oxygen as a candidate even though it contains no C/N/S/P.
# Set False later if O2 should be treated as a currency/redox species only.
KEEP_OXYGEN = True

# Stage D preservation threshold:
#
#     v_i >= beta * v_i^ref
#
# beta=0.90 means that every minimum community must preserve at least
# 90% of EACH Stage-C reference product flux.
MIN_COMMUNITY_PRODUCT_RETENTION = 0.90

# Cap enumeration for safety on larger communities.
# Use None for exhaustive enumeration without a cap.
MAX_MINIMAL_COMMUNITIES = 100

OUTPUT_FILE = ROOT / "simultaneous_products_test.yaml"


def elapsed_seconds(start: float) -> float:
    """Return elapsed wall-clock seconds from ``start``."""

    return perf_counter() - start


# ======================================================================
# LOAD MODELS
# ======================================================================

total_start = perf_counter()

print(
    f"Loading {len(MODEL_FILES)} models:",
    flush=True,
)

for model_file in MODEL_FILES:
    print(
        f"  - {model_file.name}",
        flush=True,
    )

models = load_models(
    [str(path) for path in MODEL_FILES]
)


# ======================================================================
# LOAD MEDIUM
# ======================================================================

media_data = read_compounds(
    str(MEDIA_FILE)
)

base_medium = media_data.get(
    "base_medium",
    {},
)

selected_medium = media_data.get(
    MEDIUM_NAME
)

if selected_medium is None:
    raise KeyError(
        f"Medium '{MEDIUM_NAME}' not found in {MEDIA_FILE}"
    )

# Specific-medium entries override base-medium entries if the same exchange
# appears in both dictionaries.
medium = {
    **base_medium,
    **selected_medium,
}


# ======================================================================
# BUILD FULL COMMUNITY
# ======================================================================

community = LayeredCommunity(
    "full_community",
    models,
    copy_models=False,
    params={
        Parameter.OPTIMALITY_TOL: (
            PRODUCTION_TOLERANCE
        ),
        Parameter.FEASIBILITY_TOL: (
            PRODUCTION_TOLERANCE
        ),
        Parameter.INT_FEASIBILITY_TOL: (
            INTEGER_TOLERANCE
        ),
    },
)


# ======================================================================
# REQUIRE ALL ORGANISMS TO BE PRESENT
# ======================================================================
#
# constrain_full_community() creates the organism activity variables y_j
# and imposes:
#
#     sum_j y_j = N
#
# Since every y_j is binary and there are N organisms:
#
#     y_j = 1  for every organism j
#
# MiSoSoup then requires:
#
#     growth_j >= MINIMAL_GROWTH
#
# for every active organism.
# ======================================================================

constrain_full_community(
    community,
    minimal_growth=MINIMAL_GROWTH,
)


# ======================================================================
# APPLY MEDIUM
# ======================================================================

community.setup_medium(
    medium
)


# ======================================================================
# CHECK FULL-COMMUNITY FEASIBILITY
# ======================================================================

feasibility = community.check_feasibility(
    ["community_growth"]
)

if feasibility.status != Status.OPTIMAL:
    raise RuntimeError(
        "Full community is not feasible. "
        f"Solver status: {feasibility.status}"
    )

print()
print(
    "Full community is feasible.",
    flush=True,
)

print(
    f"Minimal growth per organism: {MINIMAL_GROWTH}",
    flush=True,
)


# ======================================================================
# STEP 0: INDIVIDUAL PRODUCT MAXIMA
# ======================================================================
#
# For each community exchange i:
#
#     M_i = max v_i
#
# while:
#
#     - all organisms are present,
#     - every organism satisfies minimum growth,
#     - the medium is fixed.
# ======================================================================

reference_start = perf_counter()

max_secretion = find_producible_exchanges(
    community,
    tolerance=PRODUCTION_TOLERANCE,
)

reference_time = elapsed_seconds(
    reference_start
)

# ======================================================================
# SYSTEMATIC PRODUCT FILTER
# ======================================================================
#
# Keep extracellular products relevant to C/N/S/P metabolism.
#
# A metabolite is retained if its formula contains C, N, S, or P.
# Thus CO2, H2S, elemental S, inorganic N species, sulfate and phosphate
# remain candidates, while H+, H2O and pure metal ions are excluded.
#
# O2 is controlled separately by KEEP_OXYGEN.
# Missing formulas are kept conservatively and flagged in the audit.
# ======================================================================

unfiltered_max_secretion = dict(max_secretion)

max_secretion, product_filter_audit = filter_product_candidates(
    community=community,
    max_secretion=unfiltered_max_secretion,
    keep_oxygen=KEEP_OXYGEN,
)

excluded_products = {
    rid: info
    for rid, info in product_filter_audit.items()
    if not info["keep"]
}

unknown_formula_products = {
    rid: info
    for rid, info in product_filter_audit.items()
    if info["reason"] == "missing_formula_kept_for_review"
}

print()
print(
    "Product filtering: "
    f"{len(unfiltered_max_secretion)} -> {len(max_secretion)} candidates",
    flush=True,
)

print(
    f"Excluded products: {len(excluded_products)}",
    flush=True,
)

for rid, info in sorted(excluded_products.items()):
    print(
        f"  EXCLUDE {rid:20s} "
        f"formula={str(info['formula']):12s} "
        f"reason={info['reason']}",
        flush=True,
    )

if unknown_formula_products:
    print()
    print(
        f"WARNING: {len(unknown_formula_products)} retained products "
        "have no molecular formula and require review:",
        flush=True,
    )
    for rid in sorted(unknown_formula_products):
        print(
            f"  REVIEW {rid}",
            flush=True,
        )

print()
print(
    f"Number of candidate products |L| = {len(max_secretion)}",
    flush=True,
)

print(
    f"Individual product scan time: {reference_time:.3f} s",
    flush=True,
)


# ======================================================================
# STAGE A
# ======================================================================
#
# Solve:
#
#     K* = max sum_i z_i
#
# subject to:
#
#     z_i = 1  =>  v_i >= alpha M_i
#
# using native Gurobi indicator constraints.
# ======================================================================

print()
print("=" * 70)
print("STAGE A: MAXIMUM NUMBER OF SIMULTANEOUS PRODUCTS")
print("=" * 70)

stage_a_start = perf_counter()

stage_a = getMaxProduct(
    community=community,
    max_secretion=max_secretion,
    fraction=PRODUCTION_FRACTION,
    tolerance=PRODUCTION_TOLERANCE,
    integer_tolerance=INTEGER_TOLERANCE,
)

stage_a_time = elapsed_seconds(
    stage_a_start
)

print(
    f"K* = {stage_a['k_star']} / {len(stage_a['product_ids'])}",
    flush=True,
)

print(
    f"Stage A time: {stage_a_time:.3f} s",
    flush=True,
)

print()
print(
    "Stage A selected products:",
    flush=True,
)

for rid in stage_a["stage_a_selected_products"]:
    print(
        f"  {rid}",
        flush=True,
    )


# ======================================================================
# STAGE B.1
# ======================================================================
#
# Fix:
#
#     sum_i z_i = K*
#
# and solve:
#
#     Q* = max sum_i q_i
#
# subject to:
#
#     q_i <= z_i
#
#     z_i = 1  =>  v_i >= M_i q_i
# ======================================================================

print()
print("=" * 70)
print("STAGE B.1: MAXIMIZE NORMALIZED PRODUCTION")
print("=" * 70)

stage_b_start = perf_counter()

stage_b = maximizeProducts(
    community=community,
    stage_a=stage_a,
)

stage_b_time = elapsed_seconds(
    stage_b_start
)

print(
    f"Q* = {stage_b['q_star']:.12g}",
    flush=True,
)

print(
    "Q* reconstructed from q_i = "
    f"{stage_b['q_star_from_values']:.12g}",
    flush=True,
)

print(
    "Q* difference = "
    f"{stage_b['q_star_difference']:.12g}",
    flush=True,
)

print(
    f"Stage B.1 time: {stage_b_time:.3f} s",
    flush=True,
)


# ======================================================================
# STAGE B.2
# ======================================================================
#
# Extract:
#
#     S_B = {i : z_i = 1}
#
# No optimization is performed in this step.
# ======================================================================

print()
print("=" * 70)
print("STAGE B.2: EXTRACT SELECTED PRODUCT SET")
print("=" * 70)

stage_b2_start = perf_counter()

product_selection = (
    getSelectedProductsFromProductMaximization(
        stage_b=stage_b,
    )
)

stage_b2_time = elapsed_seconds(
    stage_b2_start
)

selected_products = product_selection[
    "selected_products"
]

print(
    f"Selected products: {len(selected_products)}",
    flush=True,
)

for rid in selected_products:
    info = product_selection[
        "selected_product_info"
    ][rid]

    print()
    print(
        rid,
        flush=True,
    )

    print(
        "    individual maximum M_i : "
        f"{info['max_individual']:.8g}",
        flush=True,
    )

    print(
        "    required >= alpha M_i  : "
        f"{info['required_flux']:.8g}",
        flush=True,
    )

    print(
        "    Stage B flux           : "
        f"{info['stage_b_flux']:.8g}",
        flush=True,
    )

    print(
        "    Stage B q_i            : "
        f"{info['stage_b_q']:.8g}",
        flush=True,
    )

    print(
        "    Stage B v_i/M_i        : "
        f"{info['stage_b_normalized_flux']:.8g}",
        flush=True,
    )

print()
print(
    f"Stage B.2 time: {stage_b2_time:.3f} s",
    flush=True,
)


# ======================================================================
# COMPARE STAGE A AND STAGE B PRODUCT IDENTITIES
# ======================================================================

stage_a_set = set(
    stage_a["stage_a_selected_products"]
)

stage_b_set = set(
    selected_products
)

removed_from_a = sorted(
    stage_a_set - stage_b_set
)

added_in_b = sorted(
    stage_b_set - stage_a_set
)

print()
print("=" * 70)
print("STAGE A -> STAGE B PRODUCT-SET COMPARISON")
print("=" * 70)

if not removed_from_a and not added_in_b:
    print(
        "No change in product identity.",
        flush=True,
    )
else:
    if removed_from_a:
        print(
            "Removed after Stage A:",
            flush=True,
        )
        for rid in removed_from_a:
            print(
                f"  {rid}",
                flush=True,
            )

    if added_in_b:
        print(
            "Added in Stage B:",
            flush=True,
        )
        for rid in added_in_b:
            print(
                f"  {rid}",
                flush=True,
            )


# ======================================================================
# STAGE C
# ======================================================================
#
# Rebuild the Stage-B product set S_B in a FRESH LayeredCommunity.
#
# Preserve exactly the objective optimized in Stage B:
#
#     sum_i q_i >= Q* - epsilon
#
# Diagnostic sequence:
#
#     C0 = fresh biological model before split variables
#     C1 = split representation added, no pFBA objective
#     C2 = pFBA minimization
#
# Then solve classical split-flux pFBA:
#
#     v_r = v_r^+ - v_r^-
#
#     v_r^+ >= 0
#     v_r^- >= 0
#
#     min sum_r(v_r^+ + v_r^-)
# ======================================================================

print()
print("=" * 70)
print("STAGE C: FIX PRODUCTS AND SOLVE pFBA")
print("=" * 70)

stage_c_start = perf_counter()

# Stage C intentionally uses a NEW community/solver.
#
# The A/B solver contains z_i variables, q_i variables, indicators and
# the K* constraint. None of those solver objects is inherited by Stage C.
stage_c_community = LayeredCommunity(
    "stage_c_reference",
    models,
    copy_models=False,
    params={
        Parameter.OPTIMALITY_TOL: (
            PRODUCTION_TOLERANCE
        ),
        Parameter.FEASIBILITY_TOL: (
            PRODUCTION_TOLERANCE
        ),
        Parameter.INT_FEASIBILITY_TOL: (
            INTEGER_TOLERANCE
        ),
    },
)

stage_c = solvepFBAUsingFixProducts(
    community=stage_c_community,
    stage_b=stage_b,
    product_selection=product_selection,
    medium=medium,
    minimal_growth=MINIMAL_GROWTH,
    lexicographic_tolerance=LEXICOGRAPHIC_TOLERANCE,
    tolerance=PRODUCTION_TOLERANCE,
    check_feasibility=True,
)

stage_c_time = elapsed_seconds(
    stage_c_start
)

print(
    "Stage B Q*                : "
    f"{stage_c['stage_b_q_star']:.12g}",
    flush=True,
)

print(
    "Lexicographic epsilon     : "
    f"{stage_c['lexicographic_tolerance']:.12g}",
    flush=True,
)

print(
    "Required sum(q_i)         : "
    f"{stage_c['stage_b_objective_floor']:.12g}",
    flush=True,
)

print(
    "Final sum(q_i)            : "
    f"{stage_c['stage_c_q_sum']:.12g}",
    flush=True,
)

print(
    "pFBA objective            : "
    f"{stage_c['pfba_objective']:.12g}",
    flush=True,
)

print(
    f"Stage C time: {stage_c_time:.3f} s",
    flush=True,
)

print()
print(
    "Final reference products:",
    flush=True,
)

for rid in stage_c["selected_products"]:
    info = stage_c[
        "reference_products"
    ][rid]

    print()
    print(
        rid,
        flush=True,
    )

    print(
        "    individual maximum M_i : "
        f"{info['max_individual']:.8g}",
        flush=True,
    )

    print(
        "    required fraction alpha: "
        f"{info['required_fraction']:.8g}",
        flush=True,
    )

    print(
        "    required >= alpha M_i  : "
        f"{info['required_flux']:.8g}",
        flush=True,
    )

    print(
        "    reference flux         : "
        f"{info['reference_flux']:.8g}",
        flush=True,
    )

    print(
        "    normalized flux v_i/M_i: "
        f"{info['normalized_flux']:.8g}",
        flush=True,
    )

    print(
        "    Stage C q_i           : "
        f"{info['q_value']:.8g}",
        flush=True,
    )


# ======================================================================
# FINAL CONSISTENCY CHECK
# ======================================================================

final_set = set(
    stage_c["selected_products"]
)

if final_set != stage_b_set:
    raise RuntimeError(
        "Stage C product set differs from the Stage B selected set."
    )

print()
print(
    "Stage C preserves exactly the Stage B product identities.",
    flush=True,
)


# ======================================================================
# STAGE D: MINIMUM PRODUCT-PRESERVING COMMUNITIES
# ======================================================================
#
# IMPORTANT:
# Stage D uses a FRESH LayeredCommunity.
#
# The reference community used by A/B/C contains:
#
#     sum_j y_j = N
#
# because all organisms were forced active while defining the reference
# phenotype. Reusing that solver would make community minimization
# impossible.
#
# Stage D instead solves:
#
#     N_min = min sum_j y_j
#
# subject to, for every selected Stage-C product:
#
#     v_i >= beta * v_i^ref
#
# and then enumerates alternative organism sets of exactly size N_min.
# ======================================================================

print()
print("=" * 70)
print("STAGE D: MINIMUM PRODUCT-PRESERVING COMMUNITIES")
print("=" * 70)

stage_d_start = perf_counter()

minimal_community = LayeredCommunity(
    "minimum_product_community",
    models,
    copy_models=False,
    params={
        Parameter.OPTIMALITY_TOL: (
            PRODUCTION_TOLERANCE
        ),
        Parameter.FEASIBILITY_TOL: (
            PRODUCTION_TOLERANCE
        ),
        Parameter.INT_FEASIBILITY_TOL: (
            INTEGER_TOLERANCE
        ),
    },
)

stage_d = findMinimalProductCommunities(
    community=minimal_community,
    medium=medium,
    reference_products=stage_c["reference_products"],
    minimal_growth=MINIMAL_GROWTH,
    product_retention=MIN_COMMUNITY_PRODUCT_RETENTION,
    tolerance=PRODUCTION_TOLERANCE,
    max_communities=MAX_MINIMAL_COMMUNITIES,
)

stage_d_time = elapsed_seconds(
    stage_d_start
)

print(
    f"Minimum community size N_min: {stage_d['minimum_size']}",
    flush=True,
)

print(
    f"Minimum communities found  : {stage_d['number_communities']}",
    flush=True,
)

print(
    f"Enumeration complete       : {stage_d['enumeration_complete']}",
    flush=True,
)

print(
    "Product retention beta     : "
    f"{stage_d['product_retention']:.6g}",
    flush=True,
)

for index, result in enumerate(
    stage_d["communities"],
    start=1,
):
    print()
    print(
        f"Minimum community {index}: "
        f"{', '.join(result['organisms'])}",
        flush=True,
    )

    print(
        f"    size             : {result['size']}",
        flush=True,
    )

    print(
        "    community growth : "
        f"{result['community_growth']:.8g}",
        flush=True,
    )

    print(
        "    product fluxes:",
        flush=True,
    )

    for rid in stage_c["selected_products"]:
        observed = result[
            "product_fluxes"
        ][rid]

        required = stage_d[
            "product_requirements"
        ][rid]

        print(
            f"      {rid:20s} "
            f"{observed:.8g} "
            f"(required >= {required:.8g})",
            flush=True,
        )

print()
print(
    f"Stage D time: {stage_d_time:.3f} s",
    flush=True,
)


# ======================================================================
# SERIALIZABLE YAML OUTPUT
# ======================================================================
#
# Do not dump stage_a/stage_b/stage_c directly because they contain
# ReFramed Solution objects. Instead, save only plain Python structures.
# ======================================================================

total_time = elapsed_seconds(
    total_start
)

output = {
    "configuration": {
        "community_size": len(models),
        "models": [
            path.name
            for path in MODEL_FILES
        ],
        "medium": MEDIUM_NAME,
        "minimal_growth": MINIMAL_GROWTH,
        "production_tolerance": PRODUCTION_TOLERANCE,
        "integer_tolerance": INTEGER_TOLERANCE,
        "production_fraction": PRODUCTION_FRACTION,
        "lexicographic_tolerance": LEXICOGRAPHIC_TOLERANCE,
        "keep_oxygen": KEEP_OXYGEN,
        "min_community_product_retention": (
            MIN_COMMUNITY_PRODUCT_RETENTION
        ),
        "max_minimal_communities": (
            MAX_MINIMAL_COMMUNITIES
        ),
    },

    "individual_product_maxima": {
        "number_unfiltered_producible_exchanges": (
            len(unfiltered_max_secretion)
        ),
        "number_filtered_products": (
            len(max_secretion)
        ),
        "max_secretion": max_secretion,
        "filter_audit": product_filter_audit,
    },

    "stage_a": {
        "k_star": stage_a["k_star"],
        "selected_products": (
            stage_a["stage_a_selected_products"]
        ),
    },

    "stage_b": {
        "q_star": stage_b["q_star"],
        "q_star_from_values": (
            stage_b["q_star_from_values"]
        ),
        "q_star_difference": (
            stage_b["q_star_difference"]
        ),
        "selected_products": (
            product_selection[
                "selected_products"
            ]
        ),
        "selected_product_info": (
            product_selection[
                "selected_product_info"
            ]
        ),
    },

    "stage_c": {
        "selected_products": (
            stage_c["selected_products"]
        ),
        "reference_products": (
            stage_c["reference_products"]
        ),
        "stage_b_objective_floor": (
            stage_c[
                "stage_b_objective_floor"
            ]
        ),
        "stage_c_q_sum": (
            stage_c[
                "stage_c_q_sum"
            ]
        ),
        "lexicographic_tolerance": (
            stage_c[
                "lexicographic_tolerance"
            ]
        ),
        "pfba_objective": (
            stage_c["pfba_objective"]
        ),
        "exchange_fluxes": (
            stage_c["exchange_fluxes"]
        ),
    },

    "stage_d": {
        "minimum_size": stage_d["minimum_size"],
        "product_retention": stage_d["product_retention"],
        "product_requirements": (
            stage_d["product_requirements"]
        ),
        "number_communities": (
            stage_d["number_communities"]
        ),
        "enumeration_complete": (
            stage_d["enumeration_complete"]
        ),
        "communities": (
            stage_d["communities"]
        ),
    },

    "timing_seconds": {
        "individual_product_scan": (
            reference_time
        ),
        "stage_a": (
            stage_a_time
        ),
        "stage_b1": (
            stage_b_time
        ),
        "stage_b2": (
            stage_b2_time
        ),
        "stage_c": (
            stage_c_time
        ),
        "stage_d": (
            stage_d_time
        ),
        "total": (
            total_time
        ),
    },
}


with open(
    OUTPUT_FILE,
    "w",
    encoding="utf8",
) as file_descriptor:
    yaml.safe_dump(
        output,
        file_descriptor,
        sort_keys=False,
    )


print()
print("=" * 70)
print("TEST COMPLETED SUCCESSFULLY")
print("=" * 70)

print(
    f"Output: {OUTPUT_FILE}",
    flush=True,
)

print(
    f"Total elapsed time: {total_time:.3f} s",
    flush=True,
)
