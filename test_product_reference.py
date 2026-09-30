"""Test the product-reference workflow on the marine MiSoSoup example."""

from pathlib import Path

import yaml
from reframed.solvers.solution import Status
from reframed.solvers.solver import Parameter

from misosoup.library.product_reference import (
    constrain_full_community,
    find_producible_exchanges,
)
from misosoup.library.readwrite import load_models, read_compounds
from misosoup.reframed.layered_community import LayeredCommunity


ROOT = Path(__file__).resolve().parent

STRAIN_DIR = ROOT / "examples" / "marine" / "strains"
MEDIA_FILE = ROOT / "examples" / "marine" / "media.yaml"
OUTPUT_FILE = ROOT / "product_reference_test.yaml"

MINIMAL_GROWTH = 0.01
TOLERANCE = 1e-6


# ---------------------------------------------------------
# 1. Load the complete community
# ---------------------------------------------------------

model_paths = sorted(STRAIN_DIR.glob("*.xml"))

print(f"Loading {len(model_paths)} models:")
for path in model_paths:
    print(f"  - {path.name}")

models = load_models([str(path) for path in model_paths])


# ---------------------------------------------------------
# 2. Load medium
# ---------------------------------------------------------

media = read_compounds(str(MEDIA_FILE))

base_medium = media.get("base_medium", {})
selected_medium = media["ac"]

medium = {
    **selected_medium,
    **base_medium,
}


# ---------------------------------------------------------
# 3. Construct the N-organism community
# ---------------------------------------------------------

community = LayeredCommunity(
    "full_community",
    models,
    copy_models=False,
    params={
        Parameter.OPTIMALITY_TOL: TOLERANCE,
        Parameter.FEASIBILITY_TOL: TOLERANCE,
    },
)


# ---------------------------------------------------------
# 4. Force all organisms to be active
# ---------------------------------------------------------

constrain_full_community(
    community,
    minimal_growth=MINIMAL_GROWTH,
)

community.setup_medium(medium)


# ---------------------------------------------------------
# 5. Verify feasibility of the complete community
# ---------------------------------------------------------

feasibility = community.check_feasibility(
    ["community_growth"]
)

if feasibility.status != Status.OPTIMAL:
    raise RuntimeError(
        f"Full community is not feasible. "
        f"Solver status: {feasibility.status}"
    )

print()
print("Full community is feasible.")
print(
    "Community growth:",
    feasibility.values.get("community_growth"),
)


# ---------------------------------------------------------
# 6. Positive side of FVA:
#    maximize every global exchange independently
# ---------------------------------------------------------

products = find_producible_exchanges(
    community,
    tolerance=TOLERANCE,
)


# ---------------------------------------------------------
# 7. Save results
# ---------------------------------------------------------

output = {
    "community_size": len(models),
    "medium": "ac",
    "minimal_growth": MINIMAL_GROWTH,
    "production_tolerance": TOLERANCE,
    "number_producible_exchanges": len(products),
    "max_secretion": dict(sorted(products.items())),
}

with open(OUTPUT_FILE, "w", encoding="utf8") as handle:
    yaml.safe_dump(
        output,
        handle,
        sort_keys=False,
    )


# ---------------------------------------------------------
# 8. Human-readable summary
# ---------------------------------------------------------

print()
print("=" * 70)
print("PRODUCT REFERENCE TEST")
print("=" * 70)

print(f"Community size: {len(models)}")
print(f"Producible exchanges: {len(products)}")

print()
print("Maximum secretion fluxes:")

for reaction, maximum in sorted(
    products.items(),
    key=lambda item: item[1],
    reverse=True,
):
    print(f"{reaction:35s} {maximum:14.6g}")

print()
print(f"Results written to: {OUTPUT_FILE}")