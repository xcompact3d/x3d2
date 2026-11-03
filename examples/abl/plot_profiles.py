#!/usr/bin/env python3
"""Plot the mean profiles written by the ``abl`` examples.

``abl_profile.csv`` (``input.x3d``) gives the streamwise velocity against the
log law and the non-dimensional shear Phi = kappa y / u_* dU/dy, with the
imposed u_*. ``abl_ekman_profile.csv`` (``input_ekman.x3d``) gives the wind
speed and direction against height.

usage: plot_profiles.py [--profile abl_profile.csv]
                        [--ekman-profile abl_ekman_profile.csv]
                        [--kappa 0.4] [--output-dir .]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SERIES = "#2a78d6"
REFERENCE = "#6b6b66"


def load(path: Path) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    """Header values and columns of a profile CSV."""
    meta: dict[str, float] = {}
    with path.open(encoding="utf-8") as f:
        for line in f:
            if not line.startswith("#"):
                break
            key, sep, value = line[1:].partition("=")
            if sep:
                meta[key.strip()] = float(value)
            else:
                names = [n.strip() for n in key.split(",")]
    data = np.loadtxt(path, delimiter=",", comments="#")
    return meta, {n: data[:, i] for i, n in enumerate(names)}


def style(ax: plt.Axes) -> None:
    ax.grid(color="#e4e3dd", linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def plot_neutral(path: Path, kappa: float, out: Path) -> None:
    meta, p = load(path)
    u_star = meta["imposed_u_star"]
    height = p["y"][-1]
    y = p["y"][1:]
    # The central difference at the first vertex spans the no-slip floor
    # (u = 0), so Phi is shown from the second vertex up.
    phi = (kappa*p["y"]/u_star
           * np.gradient(p["u_mean"], p["y"], edge_order=2))[2:]

    fig, (a, b) = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)
    a.semilogx(y/height, p["u_mean"][1:]/u_star, color=SERIES, lw=2,
               marker="o", ms=4, label="x3d2")
    a.semilogx(y/height, p["u_log"][1:]/u_star, color=REFERENCE, lw=1.5,
               ls="--", label=r"log law $(1/\kappa)\ln(y/z_0)$")
    a.set_xlabel("$y/H$")
    a.set_ylabel(r"$\langle u\rangle/u_*$")
    a.legend(frameon=False)
    style(a)

    b.plot(phi, y[1:]/height, color=SERIES, lw=2, marker="o", ms=4)
    b.axvline(1, color=REFERENCE, lw=1.5, ls="--")
    b.set_ylim(0, 0.5)
    b.set_xlabel(r"$\Phi = (\kappa y/u_*)\,\partial\langle u\rangle/\partial y$")
    b.set_ylabel("$y/H$")
    style(b)
    fig.savefig(out, dpi=150)
    print(f"{out}: diagnosed u_* = {meta['diagnosed_u_star']:.4f} m/s "
          f"(imposed {u_star:.4f}), peak Phi below 0.2H = "
          f"{phi[y[1:] <= 0.2*height].max():.2f}")


def plot_ekman(path: Path, out: Path) -> None:
    meta, p = load(path)
    height = p["y"][-1]
    y = p["y"][1:]
    speed = np.hypot(p["u_mean"], p["w_mean"])[1:]
    direction = np.degrees(np.arctan2(p["w_mean"], p["u_mean"]))[1:]

    fig, (a, b) = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)
    a.plot(speed, y/height, color=SERIES, lw=2, marker="o", ms=4)
    a.set_xlabel(r"wind speed $\sqrt{\langle u\rangle^2 + \langle w\rangle^2}$ (m/s)")
    a.set_ylabel("$y/H$")
    style(a)

    b.plot(direction, y/height, color=SERIES, lw=2, marker="o", ms=4)
    b.axvline(0, color=REFERENCE, lw=1.5, ls="--")
    b.set_xlabel(r"wind direction $\arctan(\langle w\rangle/\langle u\rangle)$ (deg)")
    b.set_ylabel("$y/H$")
    style(b)
    fig.savefig(out, dpi=150)
    layer = (p["y"] >= 0.05*height) & (p["y"] <= 0.8*height)
    turn = np.degrees(np.ptp(np.unwrap(np.arctan2(p["w_mean"][layer],
                                                  p["u_mean"][layer]))))
    print(f"{out}: diagnosed u_* = {meta['diagnosed_u_star']:.4f} m/s, "
          f"veer over 0.05-0.8H = {turn:.2f} deg")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--profile", type=Path)
    parser.add_argument("--ekman-profile", type=Path)
    parser.add_argument("--kappa", type=float, default=0.4)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    args = parser.parse_args()
    if args.profile is None and args.ekman_profile is None:
        parser.error("give --profile and/or --ekman-profile")
    if args.profile:
        plot_neutral(args.profile, args.kappa,
                     args.output_dir/"abl_neutral_profiles.png")
    if args.ekman_profile:
        plot_ekman(args.ekman_profile, args.output_dir/"abl_ekman_profiles.png")


if __name__ == "__main__":
    main()
