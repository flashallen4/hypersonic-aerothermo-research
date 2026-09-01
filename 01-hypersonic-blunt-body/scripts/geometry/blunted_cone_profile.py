"""
Project 01 - Hypersonic Blunt Body
Parametric spherically-blunted-cone profile generator.

Given nose radius (R_n), cone half-angle (theta_c), and base radius (R_b),
computes:
  - tangent point (x_t, r_t) where the spherical cap meets the cone tangentially
  - overall body length L
  - the full axisymmetric profile r(x) as a discretized polyline

This module is the single source of truth for the body geometry. It is
imported by the Gmsh .geo generation script (scripts/mesh/) so the profile
math is defined in exactly one place.

Usage (standalone):
    python3 blunted_cone_profile.py
Produces a validation plot of the baseline profile in results/figures/.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from dataclasses import dataclass


@dataclass
class BluntedConeParams:
    R_n: float       # nose radius [m]
    theta_c_deg: float  # cone half-angle [deg]
    R_b: float       # base radius [m]

    @property
    def theta_c(self) -> float:
        return np.radians(self.theta_c_deg)


@dataclass
class BluntedConeGeometry:
    params: BluntedConeParams
    x_t: float   # tangent point, axial [m]
    r_t: float   # tangent point, radial [m]
    L: float     # overall body length [m]
    valid: bool  # False if r_t >= R_b (no valid conical afterbody)


def compute_tangent_and_length(params: BluntedConeParams) -> BluntedConeGeometry:
    """
    Compute tangent point and body length for a spherically-blunted cone.

    Raises no exceptions on invalid geometry; instead flags `valid=False`
    so the caller can decide how to handle it (this project treats invalid
    sweep points as a hard error at the sweep-configuration stage, not here).
    """
    R_n = params.R_n
    theta_c = params.theta_c
    R_b = params.R_b

    r_t = R_n * np.cos(theta_c)
    x_t = R_n * (1.0 - np.sin(theta_c))

    valid = r_t < R_b

    if valid:
        L = x_t + (R_b - r_t) / np.tan(theta_c)
    else:
        # Geometry invalid: sphere alone exceeds base radius before tangency.
        # L is undefined in the blunted-cone sense; report NaN.
        L = float("nan")

    return BluntedConeGeometry(params=params, x_t=x_t, r_t=r_t, L=L, valid=valid)


def profile_r_of_x(geom: BluntedConeGeometry, n_nose: int = 200, n_cone: int = 200) -> tuple[np.ndarray, np.ndarray]:
    """
    Discretize the axisymmetric profile r(x) over [0, L].

    Returns (x, r) arrays, nose cap followed by conical afterbody,
    sharing the tangent point exactly once (no duplicate point).
    """
    if not geom.valid:
        raise ValueError(
            f"Cannot generate profile: geometry invalid "
            f"(r_t={geom.r_t:.5f} m >= R_b={geom.params.R_b:.5f} m). "
            f"No valid tangent conical afterbody exists for these parameters."
        )

    R_n = geom.params.R_n
    theta_c = geom.params.theta_c
    x_t = geom.x_t
    r_t = geom.r_t
    L = geom.L

    # Spherical cap: x in [0, x_t]
    x_nose = np.linspace(0.0, x_t, n_nose)
    r_nose = np.sqrt(np.clip(R_n**2 - (x_nose - R_n) ** 2, 0.0, None))

    # Conical afterbody: x in (x_t, L], skip duplicate tangent point
    x_cone = np.linspace(x_t, L, n_cone)[1:]
    r_cone = r_t + (x_cone - x_t) * np.tan(theta_c)

    x = np.concatenate([x_nose, x_cone])
    r = np.concatenate([r_nose, r_cone])

    return x, r


def print_summary(geom: BluntedConeGeometry) -> None:
    p = geom.params
    print(f"--- Blunted Cone Geometry Summary ---")
    print(f"  R_n (nose radius)     = {p.R_n:.5f} m")
    print(f"  theta_c (half-angle)  = {p.theta_c_deg:.3f} deg")
    print(f"  R_b (base radius)     = {p.R_b:.5f} m")
    print(f"  valid                 = {geom.valid}")
    if geom.valid:
        print(f"  x_t (tangent, axial)  = {geom.x_t:.5f} m")
        print(f"  r_t (tangent, radial) = {geom.r_t:.5f} m")
        print(f"  L (body length)       = {geom.L:.5f} m")
    else:
        print(f"  r_t = {geom.r_t:.5f} m >= R_b = {p.R_b:.5f} m -> INVALID (no conical afterbody)")


def plot_profile(geom: BluntedConeGeometry, out_path: Path) -> None:
    x, r = profile_r_of_x(geom)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, r, color="tab:blue", lw=2, label="Body profile r(x)")
    ax.plot(x, -r, color="tab:blue", lw=2)  # mirror for a full-body visual
    ax.axvline(geom.x_t, color="tab:red", ls="--", lw=1, label=f"Tangent point (x_t={geom.x_t:.4f} m)")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("x [m] (axial, nose at x=0)")
    ax.set_ylabel("r [m] (radial)")
    ax.set_title(
        f"Spherically-Blunted Cone Profile\n"
        f"R_n={geom.params.R_n} m, theta_c={geom.params.theta_c_deg} deg, R_b={geom.params.R_b} m, L={geom.L:.4f} m"
    )
    ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved profile plot: {out_path}")


if __name__ == "__main__":
    # Baseline case
    baseline = BluntedConeParams(R_n=0.05, theta_c_deg=15.0, R_b=0.15)
    geom = compute_tangent_and_length(baseline)
    print_summary(geom)

    if not geom.valid:
        raise SystemExit("Baseline geometry is invalid - check parameters.")

    repo_root = Path(__file__).resolve().parents[3]  # .../01-hypersonic-blunt-body/scripts/geometry -> up to repo root
    out_path = repo_root / "01-hypersonic-blunt-body" / "results" / "figures" / "baseline_profile.png"
    plot_profile(geom, out_path)

    print()
    print("--- Sweep validity check ---")
    for R_n in [0.02, 0.05, 0.10, 0.14]:
        p = BluntedConeParams(R_n=R_n, theta_c_deg=15.0, R_b=0.15)
        g = compute_tangent_and_length(p)
        print_summary(g)
        print()
