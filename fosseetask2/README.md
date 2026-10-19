# FOSSEE Autumn Semester-Long Internship 2026
## Screening Task 2: Heat Transfer CFD Analysis

This repository contains the complete OpenFOAM computational fluid dynamics (CFD) simulation cases, geometric calculations, grid independence studies, thermal analyses, and automated post-processing scripts for **Task 2: Heat Transfer** of the FOSSEE Screening Assignment.

---

## 📁 Directory Structure

```
FOSSEE_Task2/
│
├── Allrun                  # Master automated execution script for OpenFOAM
├── Allclean                # Workspace cleanup script
├── README.md               # Detailed project documentation and thermal results summary
├── simulation_results.json # Extracted numerical dataset
│
├── Case1_77C/              # OpenFOAM case: Tw = 77°C (350.15 K)
│   ├── 0/
│   │   ├── U               # Velocity field (Inlet = 0.8333 m/s for Re = 50)
│   │   ├── p               # Kinematic pressure field
│   │   └── T               # Temperature field (Inlet = 27°C, Wall = 77°C)
│   ├── constant/
│   │   ├── physicalProperties  # nu = 0.05 m2/s
│   │   └── momentumTransport   # Laminar flow specification
│   └── system/
│       ├── blockMeshDict   # 2D structured hexahedral mesh dictionary
│       ├── controlDict     # Solver control + scalarTransport for T
│       ├── fvSchemes       # Discretization schemes
│       └── fvSolution      # Matrix solvers and SIMPLE algorithm relaxation
│
├── Case2_127C/             # OpenFOAM case: Tw = 127°C (400.15 K)
│   ├── 0/ (U, p, T)
│   ├── constant/ (physicalProperties, momentumTransport)
│   └── system/ (blockMeshDict, controlDict, fvSchemes, fvSolution)
│
├── Case3_177C/             # OpenFOAM case: Tw = 177°C (450.15 K)
│   ├── 0/ (U, p, T)
│   ├── constant/ (physicalProperties, momentumTransport)
│   └── system/ (blockMeshDict, controlDict, fvSchemes, fvSolution)
│
├── grid_independence/      # Coarse (6,300 cells) and Fine (50,400 cells) cases
│   ├── Case1_77C_coarse/
│   ├── Case1_77C_fine/
│   ├── Case2_127C_coarse/
│   ├── Case2_127C_fine/
│   ├── Case3_177C_coarse/
│   └── Case3_177C_fine/
│
├── figures/                # Publication-grade 300 DPI visualization plots
│   ├── temperature_contours_Case1_77C.png
│   ├── temperature_contours_Case2_127C.png
│   ├── temperature_contours_Case3_177C.png
│   ├── vortex_recirculation_S_Case1_77C.png
│   ├── vortex_recirculation_S_Case2_127C.png
│   ├── vortex_recirculation_S_Case3_177C.png
│   ├── outlet_temperature_profiles.png
│   ├── local_nusselt_heat_flux_profiles.png
│   └── grid_independence_study_Task2.png
│
└── report/
    └── FOSSEE_Task2_Heat_Transfer_Report.md
```

---

## 📐 Part 1: Geometric Analysis

### 1. Domain Dimensions & Derivations
- **Reference Dimension**: Selected $D = 1.0\text{ m}$ (within the required range $0.1 \le D \le 10$).
- **Top Channel Wall Length**: $L_{top} = 36 D = 36.0\text{ m}$
- **Inlet Section**: Length $L_{in} = 10 D = 10.0\text{ m}$, Height $H_{in} = 3 D = 3.0\text{ m}$
- **Cavity Expanded Section**: Total Height $H_{cav} = 5 D = 5.0\text{ m}$
- **Outlet Section**: Length $L_{out} = 10 D = 10.0\text{ m}$, Height $H_{out} = 3 D = 3.0\text{ m}$
- **Expansion Step Height**: $h = H_{cav} - H_{in} = 5 D - 3 D = 2 D = 2.0\text{ m}$
- **Ramp Angles**: $\theta = 33.6^\circ$

### 2. Analytical Geometric Values
1. **Ramp Lengths ($L_1$ and $L_2$)**:
   $$L_1 = L_2 = \frac{h}{\sin(33.6^\circ)} = \frac{2 D}{\sin(33.6^\circ)} = \mathbf{3.6141\,D} \quad (= 3.6141\text{ m for } D = 1.0\text{ m})$$

2. **Horizontal Span of Ramps ($\Delta x_1$ and $\Delta x_2$)**:
   $$\Delta x_1 = \Delta x_2 = \frac{h}{\tan(33.6^\circ)} = \frac{2 D}{\tan(33.6^\circ)} = \mathbf{3.0102\,D} \quad (= 3.0102\text{ m for } D = 1.0\text{ m})$$

3. **Bottom Cavity Floor Length ($X$)**:
   $$X = L_{top} - (L_{in} + \Delta x_1 + \Delta x_2 + L_{out}) = 36 D - (20 D + 2 \times 3.0102 D) = \mathbf{9.9795\,D} \quad (= 9.9795\text{ m})$$

4. **Angle $\alpha$**:
   - Alternate interior angle: $\alpha = \mathbf{33.6^\circ}$
   - Internal obtuse corner angle: $\alpha_{interior} = 180^\circ - 33.6^\circ = \mathbf{146.4^\circ}$

5. **Coordinates of Centres (Midpoints) of $L_1$ and $L_2$**:
   - **Centre of $L_1$**: $(11.5051, 1.0000, 0.5000)\text{ m}$
   - **Centre of $L_2$**: $(24.4949, 1.0000, 0.5000)\text{ m}$

---

## 🧪 Part 2: CFD Simulation Setup & Methodology

### 1. Working Fluid & Thermal Properties (Air @ $27^\circ\text{C}$)
- **Inlet Temperature**: $T_{in} = 27^\circ\text{C} = \mathbf{300.15\text{ K}}$
- **Kinematic Viscosity (Specified)**: $\mathbf{\nu = 0.05\text{ m}^2/\text{s}}$
- **Fluid Density**: $\rho = \mathbf{1.176\text{ kg/m}^3}$
- **Specific Heat**: $c_p = \mathbf{1005\text{ J/(kg}\cdot\text{K)}}$
- **Prandtl Number**: $Pr = \mathbf{0.71}$
- **Thermal Diffusivity**: $\alpha = \frac{\nu}{Pr} = \frac{0.05}{0.71} = \mathbf{0.0704225\text{ m}^2/\text{s}}$
- **Thermal Conductivity**: $k = \rho c_p \alpha = \mathbf{83.2307\text{ W/(m}\cdot\text{K)}}$

### 2. Reynolds Number & Inlet Velocity
$$Re = 50.0, \quad H_{in} = 3.0\text{ m} \implies U_{in} = \frac{Re \cdot \nu}{H_{in}} = \frac{50 \times 0.05}{3.0} = \mathbf{0.83333\text{ m/s}}$$

---

## 📊 Summary of Simulation Results

### Hydrodynamic and Heat Transfer Results Table:

| Parameter | Symbol / Formulation | $\mathbf{Tw = 77^\circ C}$ (Case 1) | $\mathbf{Tw = 127^\circ C}$ (Case 2) | $\mathbf{Tw = 177^\circ C}$ (Case 3) |
| :--- | :--- | :--- | :--- | :--- |
| **Inlet Temperature** | $T_{in}$ | $27.0^\circ\text{C}$ ($300.15\text{ K}$) | $27.0^\circ\text{C}$ ($300.15\text{ K}$) | $27.0^\circ\text{C}$ ($300.15\text{ K}$) |
| **Heated Wall Temperature** | $T_w$ | $77.0^\circ\text{C}$ ($350.15\text{ K}$) | $127.0^\circ\text{C}$ ($400.15\text{ K}$) | $177.0^\circ\text{C}$ ($450.15\text{ K}$) |
| **Temperature Difference** | $\Delta T = T_w - T_{in}$ | $50.0\text{ K}$ | $100.0\text{ K}$ | $150.0\text{ K}$ |
| **Vortex Length $S$ (Ref 3)** | **Measured from Point A** | **$12.9500\text{ m}$** | **$12.9500\text{ m}$** | **$12.9500\text{ m}$** |
| **Vortex Reattachment Point** | $(x_{reatt}, y_{reatt})$ (m) | $(13.0302,\, 0.0000)$ | $(13.0302,\, 0.0000)$ | $(13.0302,\, 0.0000)$ |
| **Primary Vortex Center** | $(x_v, y_v)$ (m) | $(13.0417,\, 0.6040)$ | $(13.0417,\, 0.6040)$ | $(13.0417,\, 0.6040)$ |
| **Distance $Z$ (Ref 1)** | $\sqrt{(x_v - x_{cL1})^2 + (y_v - y_{cL1})^2}$ | **$1.5868\text{ m}$** | **$1.5868\text{ m}$** | **$1.5868\text{ m}$** |
| **Vortex Wall Extent $C_1$** | Ramp separation to corner | **$2.1835\text{ m}$** | **$2.1835\text{ m}$** | **$2.1835\text{ m}$** |
| **Vortex Wall Extent $C_2$** | Corner to reattachment on floor| **$0.0200\text{ m}$** | **$0.0200\text{ m}$** | **$0.0200\text{ m}$** |
| **Total Vortex Wall Size $C$** | Ref 2 ($C = C_1 + C_2$) | **$2.2035\text{ m}$** | **$2.2035\text{ m}$** | **$2.2035\text{ m}$** |
| **Outlet Bulk Temperature** | $\mathbf{T_{out,bulk}}$ | **$49.56^\circ\text{C}$ ($322.71\text{ K}$)** | **$72.12^\circ\text{C}$ ($345.27\text{ K}$)** | **$94.68^\circ\text{C}$ ($367.83\text{ K}$)** |
| **Outlet Mean Temperature** | $\mathbf{T_{out,mean}}$ | **$50.21^\circ\text{C}$ ($323.36\text{ K}$)** | **$73.41^\circ\text{C}$ ($346.56\text{ K}$)** | **$96.62^\circ\text{C}$ ($369.77\text{ K}$)** |
| **Average Nusselt Number** | $\overline{Nu}$ | **$50.19$** | **$50.19$** | **$50.19$** |
| **Total Heat Rate** | $Q$ ($\text{W/m span}$) | **$2.489 \times 10^6\text{ W/m}$** | **$4.978 \times 10^6\text{ W/m}$** | **$7.468 \times 10^6\text{ W/m}$** |

---

## 📈 Grid Independence Study

A three-level systematic spatial mesh convergence study was performed across all three thermal cases:
- **Coarse Mesh**: $6,300$ cells
- **Medium Mesh**: $25,200$ cells
- **Fine Mesh**: $50,400$ cells

| Thermal Case | Mesh Resolution | Cells | Vortex Length $S$ (m) | Outlet Bulk Temp $T_{out}$ (°C) | Error vs Fine |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Case 1 ($77^\circ\text{C}$)** | Coarse | 6,300 | 13.0000 | 48.91 | 1.31% |
| | Medium | 25,200 | 12.9500 | 49.56 | 0.08% |
| | Fine | 50,400 | 12.9460 | 49.60 | — (Ref) |
| **Case 2 ($127^\circ\text{C}$)** | Coarse | 6,300 | 13.0000 | 70.82 | 1.80% |
| | Medium | 25,200 | 12.9500 | 72.12 | 0.08% |
| | Fine | 50,400 | 12.9460 | 72.18 | — (Ref) |
| **Case 3 ($177^\circ\text{C}$)** | Coarse | 6,300 | 13.0000 | 92.73 | 2.10% |
| | Medium | 25,200 | 12.9500 | 94.68 | 0.08% |
| | Fine | 50,400 | 12.9460 | 94.76 | — (Ref) |

The relative error between Medium and Fine meshes is $<0.08\%$, demonstrating rigorous spatial grid independence.

---

## 🚀 How to Run the Cases

```bash
cd FOSSEE_Task2
chmod +x Allrun Allclean
./Allrun
```
To visualize in ParaView:
```bash
cd Case1_77C && paraFoam
```
*(or open `case.foam` in ParaView on Windows and select field `T` or `U`).*
