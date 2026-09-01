# Hypersonic & Aerothermodynamic Computational Research Program

## Overview

This repository contains an independent, research-oriented computational aerospace program consisting of 12 standalone projects. It is separate from conventional CFD/FEA software-demonstration portfolios. The purpose here is not to showcase tool usage, but to investigate real aerospace engineering problems through the full research cycle:

Engineering problem → Research question → Hypothesis → Theory →
Computational model → Verification → Validation →
Parametric investigation → Physical interpretation →
Engineering decision → Limitations

## Research Interests

- Hypersonic aerodynamics
- Aerothermodynamics
- Hypersonic boundary-layer heating and transition
- Thermal protection systems (TPS)
- Active/effusion cooling
- Fluid–structure interaction (FSI)
- Aeroelasticity and aerothermoelasticity
- Scramjet inlet/isolator/combustor physics
- Multidisciplinary aerospace design and optimization

## Computational Toolchain

| Category | Tools |
|---|---|
| CFD | OpenFOAM 11, ANSYS Fluent |
| FEA / Structures | CalculiX, ANSYS Mechanical |
| Geometry | OpenSCAD, Fusion 360 |
| Meshing | Gmsh |
| Programming / Analysis | Python (NumPy, SciPy, Pandas, Matplotlib) |
| Visualization | PyVista, Matplotlib, ParaView (secondary) |
| Environment | Windows + WSL2, Ubuntu 24.04 LTS |
| Reproducibility | Git, GitHub, Bash |

## Research Philosophy

Correctness, reproducibility, verification, and validation take priority over visual complexity. Additional physics (multiphysics coupling, reacting flow, real-gas effects, etc.) is introduced only when justified by the research question — never to make a project appear more advanced than its evidence supports. Negative or inconclusive results are documented as such rather than adjusted to force agreement with expectation or literature.

Each project is **independently executable** — no project depends on another's files, geometry, or solver configuration. Methodology and experience carry forward; files do not.

## Standard Project Workflow
Engineering problem → Literature review → Research question → Hypothesis →
Engineering objective → Theory/governing equations → Assumptions →
Computational model → Geometry → Mesh → Solver → Baseline case →
Verification → Validation → Parametric study → Physical interpretation →
Engineering decision → Limitations → Final conclusion

## Projects

| # | Project | Research Area | Status |
|---|---|---|---|
| 01 | Hypersonic Blunt Body | Hypersonic aerodynamics / aerothermodynamics | Not started |
| 02 | Shock–Boundary-Layer Interaction Control | Hypersonic aerodynamics / flow control | Not started |
| 03 | Hypersonic Boundary-Layer Transition and Heating | Hypersonic boundary layers / aerothermodynamics | Not started |
| 04 | Thermal Protection System Mass Optimization | Aerothermodynamics / TPS | Not started |
| 05 | Active/Effusion Cooling of Hypersonic Structures | Aerothermodynamics / thermal management | Not started |
| 06 | Hypersonic Aerothermoelastic Panel | Fluid–structure interaction / aerothermoelasticity | Not started |
| 07 | Flexible Hypersonic Control Surface | FSI / aeroelasticity | Not started |
| 08 | Scramjet Isolator Shock-Train Dynamics | Hypersonic propulsion | Not started |
| 09 | Scramjet Inlet–Isolator–Combustor | Hypersonic propulsion | Not started |
| 10 | Hypersonic Blunt-Body Aerothermal Optimization | Hypersonic vehicle design / optimization | Not started |
| 11 | Aero-Thermo-Structural Optimization | Multiphysics / multidisciplinary design | Not started |
| 12 | Integrated Hypersonic System | Integrated aerospace systems | Not started |

Status will be updated as each project progresses (Not started → In progress → Complete). No project is marked complete until it satisfies the full completion standard: research question answered, verification performed, validation performed or its absence documented, parametric investigation carried out, physical interpretation given, engineering conclusion stated, and limitations documented.

## Reproducibility

Each project directory contains its own `README.md`, `PROJECT_DEFINITION.md`, and complete instructions for reproducing all geometry, mesh, solver, and post-processing steps from a clean environment.

## Author

Aerospace engineering research portfolio — computational aerothermodynamics and multidisciplinary design.
