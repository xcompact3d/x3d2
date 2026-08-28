#!/usr/bin/env python3
"""Analyse neutral-ABL mean profiles produced by x3d2."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Profile:
    y: np.ndarray
    u: np.ndarray
    v: np.ndarray
    w: np.ndarray
    u_log: np.ndarray
    metadata: dict[str, float]


@dataclass(frozen=True)
class NeutralMetrics:
    log_law_relative_l2: float
    phi_outer_mean: float
    phi_surface_max: float
    friction_velocity_relative_error: float


def load_profile(path: Path) -> Profile:
    metadata: dict[str, float] = {}
    column_names: list[str] | None = None
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if not line.startswith("#"):
                break
            item = line[1:].strip()
            if "=" in item:
                name, value = item.split("=", 1)
                metadata[name.strip()] = float(value)
            elif item.lower().startswith("y,"):
                column_names = [name.strip().lower() for name in item.split(",")]

    if column_names is None:
        raise ValueError(f"{path}: missing named profile header")
    data = np.loadtxt(path, delimiter=",", comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != len(column_names):
        raise ValueError(f"{path}: header and data column counts differ")

    required = ["y", "u_mean", "v_mean", "w_mean", "u_log"]
    missing = [name for name in required if name not in column_names]
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(missing)}")
    if data.shape[0] < 3 or not np.all(np.diff(data[:, column_names.index("y")]) > 0):
        raise ValueError(f"{path}: profile requires at least three increasing y values")

    def column(name: str) -> np.ndarray:
        return data[:, column_names.index(name)]

    return Profile(
        column("y"), column("u_mean"), column("v_mean"),
        column("w_mean"), column("u_log"), metadata,
    )


def phi(profile: Profile, kappa: float, imposed_u_star: float) -> np.ndarray:
    if kappa <= 0.0 or imposed_u_star <= 0.0:
        raise ValueError("kappa and imposed u_star must be positive")
    return kappa*profile.y/imposed_u_star*np.gradient(
        profile.u, profile.y, edge_order=2
    )


def _height_mask(y: np.ndarray, height: float, lower: float, upper: float) -> np.ndarray:
    mask = (y/height >= lower) & (y/height <= upper)
    if np.count_nonzero(mask) < 2:
        raise ValueError(f"profile has too few points in [{lower}, {upper}] H")
    return mask


def neutral_metrics(
    profile: Profile,
    height: float,
    kappa: float,
    surface_layer: tuple[float, float],
    outer_layer: tuple[float, float],
) -> NeutralMetrics:
    imposed = profile.metadata["imposed_u_star"]
    diagnosed = profile.metadata["diagnosed_u_star"]
    phi_values = phi(profile, kappa, imposed)
    surface = _height_mask(profile.y, height, *surface_layer)
    outer = _height_mask(profile.y, height, *outer_layer)
    reference = profile.u_log[surface]
    log_error = np.linalg.norm(profile.u[surface] - reference)/np.linalg.norm(reference)
    return NeutralMetrics(
        log_error,
        float(np.mean(phi_values[outer])),
        float(np.max(phi_values[surface])),
        abs(diagnosed - imposed)/imposed,
    )


def ekman_veer_degrees(profile: Profile, height: float) -> float:
    layer = _height_mask(profile.y, height, 0.05, 0.8)
    angle = np.unwrap(np.arctan2(profile.w[layer], profile.u[layer]))
    return float(np.degrees(np.max(angle) - np.min(angle)))


def write_metrics(
    path: Path,
    metrics: NeutralMetrics,
    veer_degrees: float | None,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["metric", "value"])
        for name, value in metrics.__dict__.items():
            writer.writerow([name, value])
        if veer_degrees is not None:
            writer.writerow(["ekman_veer_degrees", veer_degrees])


def plot_profiles(
    path: Path,
    neutral: Profile,
    phi_values: np.ndarray,
    height: float,
    ekman: Profile | None,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    columns = 3 if ekman is not None else 2
    figure, axes = plt.subplots(1, columns, figsize=(5*columns, 5))
    axes[0].plot(neutral.u, neutral.y/height, label="x3d2 mean")
    axes[0].plot(neutral.u_log, neutral.y/height, "k--", label="log law")
    axes[0].set(xlabel=r"$\langle u\rangle$", ylabel=r"$y/H$")
    axes[0].legend()

    axes[1].plot(phi_values, neutral.y/height)
    axes[1].axvline(1.0, color="k", linestyle="--")
    axes[1].axvspan(1.2, 1.5, color="0.8", alpha=0.5)
    axes[1].set(xlabel=r"$\Phi$", ylabel=r"$y/H$")

    if ekman is not None:
        axes[2].plot(ekman.u, ekman.y/height, label=r"$\langle u\rangle$")
        axes[2].plot(ekman.w, ekman.y/height, label=r"$\langle w\rangle$")
        axes[2].set(xlabel="mean velocity", ylabel=r"$y/H$")
        axes[2].legend()
    for axis in axes:
        axis.grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(path, dpi=160)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--ekman-profile", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("validation"))
    parser.add_argument("--height", type=float, default=1000.0)
    parser.add_argument("--kappa", type=float, default=0.4)
    parser.add_argument("--surface-layer", type=float, nargs=2, default=(0.02, 0.2))
    parser.add_argument("--outer-layer", type=float, nargs=2, default=(0.2, 0.5))
    parser.add_argument("--max-log-error", type=float, default=0.15)
    parser.add_argument("--max-ustar-error", type=float, default=0.05)
    parser.add_argument("--min-ekman-veer", type=float, default=1.0)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    neutral = load_profile(args.profile)
    metrics = neutral_metrics(
        neutral, args.height, args.kappa,
        tuple(args.surface_layer), tuple(args.outer_layer),
    )
    phi_values = phi(neutral, args.kappa, neutral.metadata["imposed_u_star"])
    ekman = load_profile(args.ekman_profile) if args.ekman_profile else None
    veer = ekman_veer_degrees(ekman, args.height) if ekman else None

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_metrics(args.output_dir / "metrics.csv", metrics, veer)
    if not args.no_plot:
        try:
            plot_profiles(
                args.output_dir / "profiles.png",
                neutral, phi_values, args.height, ekman,
            )
        except ModuleNotFoundError as error:
            if error.name != "matplotlib":
                raise
            print("Matplotlib is unavailable; metrics were written without a plot.")

    print(f"log-law relative L2: {metrics.log_law_relative_l2:.4e}")
    print(f"outer-layer mean Phi: {metrics.phi_outer_mean:.4f}")
    print(f"surface-layer max Phi: {metrics.phi_surface_max:.4f}")
    print(
        "friction-velocity relative error: "
        f"{metrics.friction_velocity_relative_error:.4e}"
    )
    if veer is not None:
        print(f"Ekman veer: {veer:.3f} degrees")

    if args.check:
        passed = (
            metrics.log_law_relative_l2 <= args.max_log_error
            and metrics.friction_velocity_relative_error <= args.max_ustar_error
            and 0.8 <= metrics.phi_outer_mean <= 1.2
            and 1.2 <= metrics.phi_surface_max <= 1.5
            and (veer is None or veer >= args.min_ekman_veer)
        )
        if not passed:
            raise SystemExit("neutral ABL validation criteria failed")


if __name__ == "__main__":
    main()
