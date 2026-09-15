# Hypersonic & Aerothermodynamic Computational Research Program

## Overview

An independent computational research program focused on the physics and engineering challenges of high-speed flight.

The program consists of **12 standalone projects**, progressing from fundamental hypersonic aerodynamics and heat transfer toward thermal protection, active cooling, fluid–structure interaction, aeroelasticity, hypersonic propulsion, and multidisciplinary design.

Unlike the separate CFD and FEA portfolios, these projects are organized around **research questions rather than software demonstrations**.

The objective is to investigate problems through a complete computational research cycle:

**Problem → Question → Model → Verification → Validation → Investigation → Physical Interpretation → Engineering Decision**

Each project is treated as an independent study, with its own assumptions, computational setup, evidence, conclusions, and limitations.

---

## Research Scope

The program explores:

* Hypersonic aerodynamics
* Aerothermodynamics and surface heating
* Boundary-layer transition
* Thermal protection systems
* Active and effusion cooling
* Fluid–structure interaction
* Aeroelasticity and aerothermoelasticity
* Scramjet inlet, isolator, and combustor physics
* Multidisciplinary aerospace design and optimization

The projects are deliberately sequenced so that increasingly coupled physics are introduced only after the underlying computational methods are sufficiently established.

---

## Research Standards

A project is not considered successful simply because a simulation converges or produces an expected-looking result.

Each study is expected to establish, where applicable:

* A clearly defined research question
* Relevant physical theory and governing equations
* Explicit modelling assumptions
* A justified computational model
* Numerical verification
* Validation against credible analytical, experimental, or published data where available
* A meaningful parametric investigation
* Physical interpretation of the observed behaviour
* An engineering conclusion
* Explicit limitations and sources of uncertainty

**Verification and validation are treated separately.**

Where credible validation data are unavailable, that limitation will be stated rather than presented as evidence of physical accuracy.

Negative, inconclusive, or contradictory results are retained when they provide useful information about the model or physical problem.

---

## Project Roadmap

| # | Project | Research Area | Status |
|---|---|---|---|
| 01 | [**Hypersonic Blunt Body**](./01-hypersonic-blunt-body/) | Hypersonic aerodynamics / aerothermodynamics | 🟡 In progress |
| 02 | [**Shock–Boundary-Layer Interaction Control**](./02-shock-boundary-layer-interaction/) | Hypersonic aerodynamics / flow control | ⚪ Not started |
| 03 | [**Hypersonic Boundary-Layer Transition and Heating**](./03-hypersonic-boundary-layer-transition-and-heating/) | Boundary layers / aerothermodynamics | ⚪ Not started |
| 04 | [**Thermal Protection System Mass Optimization**](./04-thermal-protection-system-mass-optimization/) | Aerothermodynamics / TPS | ⚪ Not started |
| 05 | [**Active/Effusion Cooling of Hypersonic Structures**](./05-active-effusion-cooling/) | Thermal management | ⚪ Not started |
| 06 | [**Hypersonic Aerothermoelastic Panel**](./06-hypersonic-aerothermoelastic-panel/) | FSI / aerothermoelasticity | ⚪ Not started |
| 07 | [**Flexible Hypersonic Control Surface**](./07-flexible-hypersonic-control-surface/) | FSI / aeroelasticity | ⚪ Not started |
| 08 | [**Scramjet Isolator Shock-Train Dynamics**](./08-scramjet-isolator-shock-train/) | Hypersonic propulsion | ⚪ Not started |
| 09 | [**Scramjet Inlet–Isolator–Combustor**](./09-scramjet-inlet-isolator-combustor/) | Hypersonic propulsion | ⚪ Not started |
| 10 | [**Hypersonic Blunt-Body Aerothermal Optimization**](./10-hypersonic-blunt-body-aerothermal-optimization/) | Vehicle design / optimization | ⚪ Not started |
| 11 | [**Aero-Thermo-Structural Optimization**](./11-aero-thermo-structural-optimization/) | Multiphysics / multidisciplinary design | ⚪ Not started |
| 12 | [**Integrated Hypersonic System**](./12-integrated-hypersonic-system/) | Integrated aerospace systems | ⚪ Not started |

### Program progression

**Hypersonic aerodynamics**
→ **Boundary-layer physics & heating**
→ **Thermal protection & cooling**
→ **Fluid–structure interaction**
→ **Aeroelasticity**
→ **Hypersonic propulsion**
→ **Multidisciplinary optimization**
→ **Integrated hypersonic systems**

The progression is intentional: later projects introduce greater physical coupling and computational complexity rather than simply increasing geometric complexity.

---

## Computational Methods

The program uses a combination of:

`OpenFOAM` · `ANSYS Fluent` · `CalculiX` · `Gmsh` · `Python` · `OpenSCAD` · `Fusion 360` · `PyVista`

Specific tools are selected according to the physics and research question of each project rather than forcing every study into the same software workflow.

---

## Project Structure

Each project is maintained as a self-contained research study.

```text
project/
├── README.md
├── PROJECT_DEFINITION.md
├── geometry/
├── mesh/
├── case/
├── scripts/
├── postprocessing/
├── validation/
├── results/
└── figures/
```

The exact structure may vary with the problem, but a clean environment should be sufficient to reproduce the documented computational workflow.

Projects do **not** depend on files from other projects. Methods and experience carry forward; project-specific files do not.

---

## Project Status

Projects progress through three stages:

* 🟢 **Complete** — the documented completion standard has been satisfied
* 🟡 **In progress** — active development or investigation is underway
* ⚪ **Not started** — work has not yet begun

**In progress does not mean validated or complete.**

A project may contain working simulations, preliminary results, or partial investigations while still remaining in this state.

A project is marked **Complete** only when the available evidence supports the following:

* Research question addressed
* Computational model documented
* Verification performed
* Validation performed, or its absence clearly justified
* Parametric investigation completed where relevant
* Results physically interpreted
* Engineering conclusion stated
* Limitations documented
* Reproduction instructions completed

A visually impressive result is not sufficient for completion.

---

## Long-Term Direction

This program forms the hypersonic and aerothermodynamic core of a broader computational aerospace research portfolio.

The long-term objective is to develop the capability to work on **coupled high-speed flight problems**, where aerodynamics, heat transfer, materials, structures, propulsion, and control cannot be treated independently.

The individual projects are therefore not intended to represent isolated simulations. They are steps toward understanding the coupled physics that ultimately constrain high-speed aerospace systems.

---

## Author

**Harsh Panchal**
Aerospace Engineering · Computational Aerothermodynamics · Hypersonic Research
