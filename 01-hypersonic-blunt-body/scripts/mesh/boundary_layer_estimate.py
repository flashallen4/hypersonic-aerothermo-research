"""
Project 01 - Hypersonic Blunt Body
Boundary-layer thickness estimation for near-wall mesh sizing.

Computes:
  - Freestream dynamic viscosity via Sutherland's law
  - Freestream Reynolds number per unit length
  - Estimated compressible laminar boundary-layer thickness at representative
    body stations (using a compressible flat-plate correlation with a
    reference-temperature method — standard engineering approximation,
    NOT a substitute for resolving the BL in the actual CFD solution)
  - Recommended first-cell wall-normal height and geometric growth ratio

This is a SIZING ESTIMATE to inform mesh generation, not a flow solution.
The actual boundary layer will be resolved (and checked) by the CFD run itself.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Freestream conditions (from CASE_SPECIFICATION.md, 30 km altitude, M=7)
# ---------------------------------------------------------------------------
T_inf = 226.5        # K
p_inf = 1197.0        # Pa
rho_inf = 1.841e-2    # kg/m^3
gamma = 1.4
R_gas = 287.0         # J/(kg K)
M_inf = 7.0

a_inf = np.sqrt(gamma * R_gas * T_inf)
U_inf = M_inf * a_inf

# Sutherland's law for air (dynamic viscosity)
def sutherland_mu(T, mu_ref=1.716e-5, T_ref=273.15, S=110.4):
    """Dynamic viscosity of air [Pa.s] via Sutherland's law."""
    return mu_ref * (T / T_ref) ** 1.5 * (T_ref + S) / (T + S)

mu_inf = sutherland_mu(T_inf)
nu_inf = mu_inf / rho_inf

# Reynolds number per unit length
Re_per_m = rho_inf * U_inf / mu_inf

print("--- Freestream Properties (30 km, M=7) ---")
print(f"  T_inf     = {T_inf:.2f} K")
print(f"  a_inf     = {a_inf:.2f} m/s")
print(f"  U_inf     = {U_inf:.2f} m/s")
print(f"  mu_inf    = {mu_inf:.4e} Pa.s  (Sutherland's law)")
print(f"  nu_inf    = {nu_inf:.4e} m^2/s")
print(f"  Re_per_m  = {Re_per_m:.4e} 1/m")
print()

# ---------------------------------------------------------------------------
# Compressible laminar boundary-layer thickness estimate
#
# Approach: reference-temperature (Eckert) method applied to a laminar
# flat-plate correlation, delta/x ~ 5.0 / sqrt(Re_x), with Re_x evaluated
# using reference-temperature properties to account for compressibility
# (density/viscosity variation across the hot shock layer / BL).
#
# This is a first-order engineering estimate suitable for mesh sizing only.
# It is NOT used anywhere as a validation or physics result.
# ---------------------------------------------------------------------------

def reference_temperature(T_inf, T_wall, M_inf, gamma=1.4, r=0.85):
    """
    Eckert reference temperature for compressible boundary layers.
    r = recovery factor (~0.85 for laminar Pr~0.71 -> Pr^0.5)
    T_aw = adiabatic wall temperature (recovery temperature)
    """
    T_aw = T_inf * (1.0 + r * (gamma - 1.0) / 2.0 * M_inf**2)
    T_ref = T_inf + 0.5 * (T_wall - T_inf) + 0.22 * (T_aw - T_inf)
    return T_ref, T_aw

T_wall = 300.0  # K, fixed isothermal wall (from CASE_SPECIFICATION.md)
T_ref, T_aw = reference_temperature(T_inf, T_wall, M_inf)

print("--- Reference Temperature Method ---")
print(f"  T_wall (fixed)     = {T_wall:.2f} K")
print(f"  T_aw (adiabatic)   = {T_aw:.2f} K   <- shows severity of aerodynamic heating")
print(f"  T_ref (Eckert)     = {T_ref:.2f} K")
print()

mu_ref = sutherland_mu(T_ref)
rho_ref = p_inf / (R_gas * T_ref)   # ideal gas, same p (boundary-layer edge pressure approx)
nu_ref = mu_ref / rho_ref

print(f"  mu_ref  = {mu_ref:.4e} Pa.s")
print(f"  rho_ref = {rho_ref:.4e} kg/m^3")
print(f"  nu_ref  = {nu_ref:.4e} m^2/s")
print()

def bl_thickness_estimate(x, U_e, nu_ref):
    """Laminar compressible BL thickness estimate at station x [m] from
    effective boundary-layer origin, using reference-temperature Re_x."""
    if x <= 0:
        return 0.0
    Re_x_ref = U_e * x / nu_ref
    delta = 5.0 * x / np.sqrt(Re_x_ref)
    return delta

# Representative stations along the baseline body (from blunted_cone_profile.py)
# Using arc length approx by axial distance for this first-order estimate.
stations = {
    "near stagnation (x=0.005 m)": 0.005,
    "tangent point (x=0.0371 m)": 0.0371,
    "mid-cone (x=0.20 m)": 0.20,
    "near base (x=0.40 m)": 0.40,
}

print("--- Estimated Boundary-Layer Thickness (reference-temperature method) ---")
deltas = {}
for label, x in stations.items():
    delta = bl_thickness_estimate(x, U_inf, nu_ref)
    deltas[label] = delta
    print(f"  {label:32s}: delta ~ {delta*1000:.4f} mm")
print()

# ---------------------------------------------------------------------------
# Near-wall mesh sizing recommendation
# ---------------------------------------------------------------------------
delta_min = min(d for d in deltas.values() if d > 0)

# First cell height: aim for delta_min / 100 (conservative, ~30+ cells
# across even the thinnest BL region with a growth ratio ~1.12-1.2)
first_cell_height = delta_min / 100.0
n_cells_across_thinnest_bl = delta_min / first_cell_height  # by construction ~100 at growth=1 baseline

print("--- Near-Wall Mesh Sizing Recommendation ---")
print(f"  Thinnest estimated BL region : {delta_min*1000:.4f} mm (near stagnation)")
print(f"  Recommended first-cell height: {first_cell_height*1e6:.2f} micrometers ({first_cell_height:.3e} m)")
print(f"  Suggested growth ratio        : 1.12 - 1.20 (geometric expansion normal to wall)")
print(f"  With growth ratio 1.12, cells needed to span thinnest BL (~{delta_min*1000:.3f} mm):")

# Estimate number of cells to span delta_min with geometric growth
def n_cells_for_growth(delta, first_cell, growth):
    total = 0.0
    n = 0
    h = first_cell
    while total < delta and n < 200:
        total += h
        h *= growth
        n += 1
    return n

n_needed = n_cells_for_growth(delta_min, first_cell_height, 1.12)
print(f"    -> approx {n_needed} cells (target: >=20-30 for adequate BL resolution)")
