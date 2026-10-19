# FOSSEE (OpenFOAM CFD & heat-transfer screening tasks)

Submission for the **FOSSEE Autumn Semester-Long Internship 2026** screening
tasks. Both tasks model laminar air flow through an asymmetric
diverging–converging planar cavity and are solved with **OpenFOAM-13** on
structured `blockMesh` grids, including three-level grid-independence studies.

## Task 1 — Fluid mechanics (`FOSSEE_Task1/`)

Steady, incompressible, laminar flow at **Re = 60, 80, 100**. Quantifies the
reattachment length, primary recirculation vortex centre, distance `Z` and wall
contact length `C = C1 + C2`.

- Report: [`FOSSEE_Task1_Fluid_Mechanics_Report.md`](FOSSEE_Task1_Fluid_Mechanics_Report.md)
- Cases: `Re60/`, `Re80/`, `Re100/`
- Grid study: `grid_independence/` (coarse 6,300 / medium 25,200 / fine 50,400 cells)
- Figures and report sources: `figures/`, `report/`

## Task 2 — Heat transfer (`FOSSEE_Task2/`, `fosseetask2/`)

Laminar convective heat transfer at **Re = 50** with a heated bottom wall at
**77 °C, 127 °C and 177 °C** and 27 °C inlet air. Reports outlet bulk
temperature, wall heat flux, average Nusselt number and vortex/contact lengths.

- Report: [`FOSSEE_Task2_Heat_Transfer_Report.md`](FOSSEE_Task2_Heat_Transfer_Report.md)
- Cases: `Case1_77C/`, `Case2_127C/`, `Case3_177C/`
- Grid study: `grid_independence/`

## Automation scripts

| Script | Purpose |
| --- | --- |
| `generate_cases.py`, `generate_task2_cases.py` | Generate the OpenFOAM case trees |
| `run_master_simulation_and_analysis.py` | Run Task 1 cases + post-processing |
| `run_master_task2_simulation_and_analysis.py` | Run Task 2 cases + post-processing |
| `analyze_cases.py` | Extract vortex/reattachment metrics |
| `measure_vortex_walls.py` | Wall contact-length measurements |
| `verify_geometry.py` | Geometry checks |
| `update_report_vortex.py` | Refresh report numbers |

## Requirements

- OpenFOAM 13 (`blockMesh`, `simpleFoam`, `buoyantSimpleFoam`-style workflow)
- Python 3 with NumPy/SciPy/Matplotlib for the analysis scripts

```bash
# example
python generate_cases.py
python run_master_simulation_and_analysis.py
```

## Repository scope

Raw OpenFOAM time-step fields for every iteration are large and were **not**
committed; the `0/`, `constant/` and `system/` inputs, final iterations,
figures, reports and analysis scripts are. The original `.zip`/`.rar` bundles
are also excluded.
