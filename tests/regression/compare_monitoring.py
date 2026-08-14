#!/usr/bin/env python3
"""Regression comparator for x3d2 ``monitoring.csv`` output.

Two modes:

``golden``
    Element-wise comparison of a freshly produced ``monitoring.csv`` against a
    committed golden file. Columns are matched by name (so column order or the
    addition of new, unrelated columns does not break the test). A value passes
    when ``|produced - golden| <= atol + rtol * |golden|`` (numpy ``allclose``
    semantics), which handles both O(1) quantities (via ``rtol``) and
    near-zero quantities such as divergence (via ``atol``). This is the
    day-to-day "did any commit change the numbers" guard.

``peak``
    Extracts the dissipation peak (value and time) from a produced
    ``monitoring.csv`` and checks it against reference values within a relative
    tolerance (default 5%). Dissipation is taken as ``eps = -dEk/dt`` when a
    ``kinetic_energy`` column is present, otherwise it falls back to
    ``eps = (2/Re) * enstrophy``. This is the physics-oriented check used by the
    nightly full-size run.

Blessing
    Passing ``--bless`` (or pointing ``--golden`` at a non-existent file in
    ``golden`` mode) copies the produced file to the golden path and exits 0.
    Use this to (re)generate the reference in the *same* environment that CI
    runs in — golden files are only reproducible for a fixed
    compiler/precision/rank-count. New columns are compared only once they
    exist in the golden, so re-bless after a column is added to monitoring.csv
    (e.g. kinetic_energy from PR #338) to bring it under regression coverage.

Exit status is 0 on pass and 1 on failure, so it plugs straight into CTest.
"""

from __future__ import annotations

import argparse
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

    shared = [c for c in gold_cols if c in prod_cols]
    if not shared:
        print("FAIL: no shared columns between produced and golden", file=sys.stderr)
        return 1

    failures = 0
    worst = (0.0, "", -1)  # (excess, column, row)
    for name in shared:
        prod = column(name, prod_cols, prod_rows)
        gold = column(name, gold_cols, gold_rows)
        for i, (p, g) in enumerate(zip(prod, gold)):
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
        f"PASS: {len(shared)} column(s) x {len(prod_rows)} row(s) "
        f"within tolerance (rtol={rtol:g}, atol={atol:g})"
    )
    return 0


def dissipation(columns: list[str], rows: list[list[float]], re: float):
    """Return (times, eps) for the dissipation rate.

    Prefers eps = -dEk/dt (central differences) when kinetic energy is logged;
    otherwise uses eps = (2/Re) * enstrophy.
    """
    times = column("time", columns, rows)
    if "kinetic_energy" in columns:
        ek = column("kinetic_energy", columns, rows)
        eps = []
        for i in range(len(times)):
            lo = max(i - 1, 0)
            hi = min(i + 1, len(times) - 1)
            if hi == lo:
                eps.append(0.0)
            else:
                eps.append(-(ek[hi] - ek[lo]) / (times[hi] - times[lo]))
        return times, eps
    if "enstrophy" in columns:
        if re is None:
            raise ValueError("peak mode needs --re when kinetic_energy is absent")
        ens = column("enstrophy", columns, rows)
        return times, [2.0 / re * e for e in ens]
    raise ValueError("no 'kinetic_energy' or 'enstrophy' column for dissipation")


def compare_peak(
    produced: Path,
    re: float,
    ref_value: float,
    ref_time: float,
    rtol: float,
) -> int:
    columns, rows = read_monitoring(produced)
    times, eps = dissipation(columns, rows, re)

    peak_i = max(range(len(eps)), key=lambda i: eps[i])
    peak_value = eps[peak_i]
    peak_time = times[peak_i]

    err_value = abs(peak_value - ref_value) / abs(ref_value)
    err_time = abs(peak_time - ref_time) / abs(ref_time)

    print(
        f"dissipation peak: value={peak_value:.6e} (ref {ref_value:.6e}, "
        f"err {err_value * 100:.2f}%), time={peak_time:.4f} "
        f"(ref {ref_time:.4f}, err {err_time * 100:.2f}%)"
    )

    ok = err_value <= rtol and err_time <= rtol
    if not ok:
        print(
            f"FAIL: dissipation peak outside {rtol * 100:.1f}% tolerance",
            file=sys.stderr,
        )
        return 1
    print(f"PASS: dissipation peak within {rtol * 100:.1f}% tolerance")
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

    p = sub.add_parser("peak", help="dissipation peak vs reference values")
    p.add_argument("--produced", type=Path, required=True)
    p.add_argument("--re", type=float, default=None, help="Reynolds number")
    p.add_argument("--ref-value", type=float, required=True)
    p.add_argument("--ref-time", type=float, required=True)
    p.add_argument("--rtol", type=float, default=0.05)

    args = parser.parse_args(argv)

    if args.mode == "golden":
        if args.bless or not args.golden.exists():
            args.golden.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(args.produced, args.golden)
            print(f"BLESSED: wrote golden {args.golden}")
            return 0
        return compare_golden(args.produced, args.golden, args.rtol, args.atol)

    return compare_peak(
        args.produced, args.re, args.ref_value, args.ref_time, args.rtol
    )


if __name__ == "__main__":
    sys.exit(main())
