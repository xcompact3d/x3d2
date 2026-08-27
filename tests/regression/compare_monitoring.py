#!/usr/bin/env python3
"""Regression comparator for x3d2 ``monitoring.csv`` output.

Element-wise comparison of a freshly produced ``monitoring.csv`` against a
committed golden file. Columns are matched by name (so column order or the
addition of new, unrelated columns does not break the test). A value passes
when ``|produced - golden| <= atol + rtol * |golden|`` (numpy ``allclose``
semantics), which handles both O(1) quantities (via ``rtol``) and near-zero
quantities such as divergence (via ``atol``). This is the day-to-day "did any
commit change the numbers" guard. Physics validation against reference data is
performed separately and manually.

Blessing
    Passing ``--bless`` copies the produced file to the golden path and exits
    0. This is the *only* way to write a golden: a ``--golden`` path that does
    not exist is a failure, so a typo or a rename cannot quietly turn the
    regression guard off.
    Use this to (re)generate the reference in the *same* environment that CI
    runs in — golden files are only reproducible for a fixed
    compiler/precision/rank-count. New columns are compared only once they
    exist in the golden, so re-bless after a column is added to monitoring.csv
    (e.g. kinetic_energy from PR #338) to bring it under regression coverage.

Exit status is 0 on pass and 1 on failure, so it plugs straight into CTest.
"""

from __future__ import annotations

import argparse
import math
import shutil
import sys
from pathlib import Path


def read_monitoring(path: Path) -> tuple[list[str], list[list[float]]]:
    """Parse a monitoring CSV into (column_names, rows).

    The header line looks like ``# time, enstrophy, div_u_max, ...`` and data
    lines are comma-separated floats. Leading ``#`` and surrounding whitespace
    are stripped from every field.
    """
    columns: list[str] = []
    rows: list[list[float]] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if not columns:
                header = line.lstrip("#").strip()
                columns = [c.strip() for c in header.split(",")]
                continue
            rows.append([float(v) for v in line.split(",")])
    if not columns:
        raise ValueError(f"{path}: no header found")
    return columns, rows


def column(name: str, columns: list[str], rows: list[list[float]]) -> list[float]:
    idx = columns.index(name)
    return [row[idx] for row in rows]


def compare_golden(produced: Path, golden: Path, rtol: float, atol: float) -> int:
    prod_cols, prod_rows = read_monitoring(produced)
    gold_cols, gold_rows = read_monitoring(golden)

    if len(prod_rows) != len(gold_rows):
        print(
            f"FAIL: row count mismatch: produced {len(prod_rows)} "
            f"vs golden {len(gold_rows)} rows",
            file=sys.stderr,
        )
        return 1

    # Rows are matched by index, so bail early with a clear message if the time
    # columns disagree rather than emitting a cascade of value mismatches.
    if "time" in prod_cols and "time" in gold_cols:
        prod_t = column("time", prod_cols, prod_rows)
        gold_t = column("time", gold_cols, gold_rows)
        for i, (pt, gt) in enumerate(zip(prod_t, gold_t)):
            if abs(pt - gt) > 1e-9 + 1e-9 * abs(gt):
                print(
                    f"FAIL: time column misaligned at row {i}: "
                    f"produced t={pt:.12e} vs golden t={gt:.12e}",
                    file=sys.stderr,
                )
                return 1

    missing = [c for c in gold_cols if c not in prod_cols]
    if missing:
        print(
            "FAIL: produced output is missing golden column(s): "
            + ", ".join(missing),
            file=sys.stderr,
        )
        return 1

    failures = 0
    worst = (0.0, "", -1)  # (excess, column, row)
    for name in gold_cols:
        prod = column(name, prod_cols, prod_rows)
        gold = column(name, gold_cols, gold_rows)
        for i, (p, g) in enumerate(zip(prod, gold)):
            if not math.isfinite(p) or not math.isfinite(g):
                failures += 1
                if worst[2] == -1:
                    worst = (math.inf, name, i)
                if failures <= 10:
                    print(
                        f"FAIL: {name} row {i}: non-finite value "
                        f"(produced={p}, golden={g})",
                        file=sys.stderr,
                    )
                continue
            allowed = atol + rtol * abs(g)
            diff = abs(p - g)
            if diff > allowed:
                failures += 1
                excess = diff - allowed
                if excess > worst[0]:
                    worst = (excess, name, i)
                if failures <= 10:
                    print(
                        f"FAIL: {name} row {i}: produced={p:.12e} "
                        f"golden={g:.12e} |diff|={diff:.3e} > {allowed:.3e}",
                        file=sys.stderr,
                    )

    if failures:
        print(
            f"FAIL: {failures} value(s) exceeded tolerance "
            f"(rtol={rtol:g}, atol={atol:g}); worst: column '{worst[1]}' "
            f"row {worst[2]}",
            file=sys.stderr,
        )
        return 1

    print(
        f"PASS: {len(gold_cols)} column(s) x {len(prod_rows)} row(s) "
        f"within tolerance (rtol={rtol:g}, atol={atol:g})"
    )
    return 0



def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)

    g = sub.add_parser("golden", help="element-wise diff vs a golden CSV")
    g.add_argument("--produced", type=Path, required=True)
    g.add_argument("--golden", type=Path, required=True)
    g.add_argument("--rtol", type=float, default=1e-4)
    g.add_argument("--atol", type=float, default=1e-10)
    g.add_argument("--bless", action="store_true", help="write produced as golden")

    args = parser.parse_args(argv)

    if args.bless:
        args.golden.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(args.produced, args.golden)
        print(f"BLESSED: wrote golden {args.golden}")
        return 0

    # Never bless implicitly: under CTest a missing golden would otherwise be
    # an unconditional pass that also writes a new reference into the source
    # tree, leaving a permanently green test that compares nothing.
    if not args.golden.exists():
        print(
            f"FAIL: golden file not found: {args.golden}\n"
            "      Check the path, or pass --bless to create it deliberately.",
            file=sys.stderr,
        )
        return 1

    return compare_golden(args.produced, args.golden, args.rtol, args.atol)


if __name__ == "__main__":
    sys.exit(main())
