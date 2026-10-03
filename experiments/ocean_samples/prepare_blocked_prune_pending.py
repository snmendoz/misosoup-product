"""Build a list of source models that still need blocked-reaction pruning."""

from __future__ import annotations

import argparse
from pathlib import Path

from prune_blocked_model_highs import (
    audit_path_for,
    cache_valid,
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--output-list", required=True, type=Path)
    parser.add_argument("--blocked-tolerance", type=float, default=1e-9)
    parser.add_argument("--open-exchange-bound", type=float, default=1000.0)
    return parser.parse_args()


def main():
    args = parse_args()
    source_dir = args.source_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_list = args.output_list.expanduser().resolve()

    models = sorted(source_dir.glob("*.xml"))
    if not models:
        raise FileNotFoundError(
            f"No XML models found in {source_dir}"
        )

    pending = []

    for source in models:
        output_model = output_dir / (
            source.stem + "_without_blocked_reactions.xml"
        )
        audit_path = audit_path_for(output_model)

        if not cache_valid(
            source,
            output_model,
            audit_path,
            args.blocked_tolerance,
            args.open_exchange_bound,
        ):
            pending.append(source)

    output_list.parent.mkdir(parents=True, exist_ok=True)
    output_list.write_text(
        "".join(f"{path}\n" for path in pending),
        encoding="utf-8",
    )

    print(
        f"Source models: {len(models)}; "
        f"already complete: {len(models) - len(pending)}; "
        f"pending: {len(pending)}; "
        f"list: {output_list}",
        flush=True,
    )


if __name__ == "__main__":
    main()
