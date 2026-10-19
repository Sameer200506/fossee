# FOSSEE Autumn Semester-Long Internship 2026
## Screening Task 1: Fluid Mechanics CFD Analysis

This repository contains the complete OpenFOAM computational fluid dynamics (CFD) simulation cases, geometry calculations, grid independence studies, and analysis scripts for **Task 1: Fluid Mechanics** of the FOSSEE Screening Assignment.

---

## 📁 Directory Structure

```
FOSSEE_Task1/
│
├── Allrun                  # Master automated execution script for OpenFOAM
├── Allclean                # Workspace cleanup script
├── README.md               # Detailed project documentation and results summary
├── simulation_results.json # Extracted numerical dataset
│
├── Re60/                   # OpenFOAM case for Reynolds number = 60
│   ├── 0/
│   │   ├── U               # Velocity field boundary and initial conditions
│   │   └── p               # Kinematic pressure field conditions
│   ├── constant/
│   │   ├── physicalProperties  # Kinematic viscosity (nu = 1.602e-5 m2/s for Air @ 30°C)
│   │   └── momentumTransport   # Laminar flow specification
│   └── system/
│       ├── blockMeshDict   # 2D structured hexahedral mesh generation dictionary
│       ├── controlDict     # Solver control, time stepping, and write intervals
│       ├── fvSchemes       # Discretization schemes
│       └── fvSolution      # Matrix solvers and SIMPLE algorithm relaxation
│
├── Re80/                   # OpenFOAM case for Reynolds number = 80
│   ├── 0/ (U, p)
│   ├── constant/ (physicalProperties, momentumTransport)
│   └── system/ (blockMeshDict, controlDict, fvSchemes, fvSolution)
│
├── Re100/                  # OpenFOAM case for Reynolds number = 100
│   ├── 0/ (U, p)
│   ├── constant/ (physicalProperties, momentumTransport)
│   └── system/ (blockMeshDict, controlDict, fvSchemes, fvSolution)
│
├── grid_independence/      # Coarse and Fine mesh cases for GCI study
│   ├── Re60_coarse/
│   ├── Re60_fine/
│   ├── Re80_coarse/
│   ├── Re80_fine/
│   ├── Re100_coarse/
│   └── Re100_fine/
│
└── figures/                # Publication-grade vector / 300 DPI plots
    ├── streamlines_Re60.png
    ├── streamlines_Re80.png
    ├── streamlines_Re100.png
    ├── vortex_detail_Re60.png
    ├── vortex_detail_Re80.png
    ├── vortex_detail_Re100.png
    ├── reattachment_velocity_profiles.png
    └── grid_independence_study.png
```

---

## 📐 Part 1: Geometric Analysis

### 1. Domain Dimensions & Derivations
- **Reference Dimension**: Selected $D = 1.0\text{ m}$ (within the permissible range $0.1 \le D \le 10$).
- **Top Channel Wall Length**: $L_{top} = 36 D = 36.0\text{ m}$
- **Inlet Channel**: Length $L_{in} = 10 D = 10.0\text{ m}$, Height $H_{in} = 3 D = 3.0\text{ m}$
- **Cavity Expanded Section**: Total Height $H_{cav} = 5 D = 5.0\text{ m}$
- **Outlet Channel**: Length $L_{out} = 10 D = 10.0\text{ m}$, Height $H_{out} = 3 D = 3.0\text{ m}$
- **Expansion Step Height (Ramp Drop)**: $h = H_{cav} - H_{in} = 5 D - 3 D = 2 D = 2.0\text{ m}$
- **Ramp Inclination Angle**: $\theta = 33.6^\circ$

### 2. Analytical Calculations
1. **Ramp Lengths ($L_1$ and $L_2$)**:
   $$L_1 = L_2 = \frac{h}{\sin(33.6^\circ)} = \frac{2 D}{\sin(33.6^\circ)} = \frac{2 \times 1.0}{\sin(33.6^\circ)} = \mathbf{3.6141\,D} \quad (= 3.6141\text{ m for } D = 1.0\text{ m})$$

2. **Horizontal Span of Ramps ($\Delta x_1$ and $\Delta x_2$)**:
   $$\Delta x_1 = \Delta x_2 = \frac{h}{\tan(33.6^\circ)} = \frac{2 D}{\tan(33.6^\circ)} = \frac{2 \times 1.0}{\tan(33.6^\circ)} = \mathbf{3.0102\,D} \quad (= 3.0102\text{ m for } D = 1.0\text{ m})$$

3. **Bottom Cavity Floor Length ($X$)**:
   $$X = L_{top} - (L_{in} + \Delta x_1 + \Delta x_2 + L_{out}) = 36 D - (20 D + 2 \times 3.0102 D) = \mathbf{9.9795\,D} \quad (= 9.9795\text{ m})$$

4. **Angle $\alpha$**:
   - $\alpha = \mathbf{33.6^\circ}$ (by alternate interior angle geometry with the horizontal floor)
   - Internal corner angle: $\alpha_{int} = 180^\circ - 33.6^\circ = \mathbf{146.4^\circ}$

5. **Coordinates of Centres of $L_1$ and $L_2$**:
   - Placing Cartesian coordinate origin $(0, 0)$ at the bottom-left corner of the computational domain ($y_{floor} = 0$, $y_{inlet\_floor} = 2.0\text{ m}$, $y_{top} = 5.0\text{ m}$):
   - **Centre of Ramp $L_1$**:
     $$x_{c,L1} = 10.0 + \frac{3.0102}{2} = \mathbf{11.5051\text{ m}}$$
     $$y_{c,L1} = \frac{2.0 + 0.0}{2} = \mathbf{1.0000\text{ m}}$$
     $$\mathbf{\text{Centre of } L_1 = (11.5051, 1.0000, 0.5000)\text{ m}}$$
   - **Centre of Ramp $L_2$**:
     $$x_{c,L2} = 26.0 - \frac{3.0102}{2} = \mathbf{24.4949\text{ m}}$$
     $$y_{c,L2} = \frac{0.0 + 2.0}{2} = \mathbf{1.0000\text{ m}}$$
     $$\mathbf{\text{Centre of } L_2 = (24.4949, 1.0000, 0.5000)\text{ m}}$$

---

## 🧪 Part 2: CFD Simulation Setup & Methodology

### 1. Working Fluid Properties (Air @ $30^\circ\text{C}$, STP)
- **Pressure**: $P_{atm} = 101325\text{ Pa}$
- **Temperature**: $T = 303.15\text{ K}$
- **Density**: $\rho = \frac{P}{R T} = \frac{101325}{287.058 \times 303.15} = \mathbf{1.1644\text{ kg/m}^3}$
- **Dynamic Viscosity**: $\mu = \mathbf{1.865 \times 10^{-5}\text{ Pa}\cdot\text{s}}$
- **Kinematic Viscosity**: $\nu = \frac{\mu}{\rho} = \mathbf{1.602 \times 10^{-5}\text{ m}^2/\text{s}}$

### 2. Reynolds Number & Inlet Velocities
Characteristic length scale: Inlet channel height $H_{in} = 3 D = 3.0\text{ m}$.
$$U_{in} = \frac{Re \cdot \nu}{H_{in}}$$

| Case | Reynolds Number ($Re$) | Kinematic Viscosity $\nu$ ($\text{m}^2/\text{s}$) | Characteristic Height $H_{in}$ (m) | Inlet Velocity $U_{in}$ (m/s) |
| :--- | :--- | :--- | :--- | :--- |
| **Re60** | 60 | $1.602 \times 10^{-5}$ | $3.0$ | $\mathbf{3.2040 \times 10^{-4}}$ |
| **Re80** | 80 | $1.602 \times 10^{-5}$ | $3.0$ | $\mathbf{4.2720 \times 10^{-4}}$ |
| **Re100** | 100 | $1.602 \times 10^{-5}$ | $3.0$ | $\mathbf{5.3400 \times 10^{-4}}$ |

---

## 📊 Summary of Simulation Results

### Hydrodynamic Metrics (Exact Reference 1 & Reference 2 Formulations):

| Parameter | Symbol / Formulation | $\mathbf{Re = 60}$ | $\mathbf{Re = 80}$ | $\mathbf{Re = 100}$ |
| :--- | :--- | :--- | :--- | :--- |
| **Inlet Velocity** | $U_{in}$ (m/s) | $3.2040 \times 10^{-4}$ | $4.2720 \times 10^{-4}$ | $5.3400 \times 10^{-4}$ |
| **Vortex Separation Point on Ramp** | $(x_{sep}, y_{sep})$ (m) | $(10.9470,\, 1.3708)$ | $(10.6890,\, 1.5422)$ | $(10.5503,\, 1.6344)$ |
| **Separation Distance from Pt A** | $s_{sep}$ along Ramp (m) | $1.1369$ | $0.8272$ | $0.6591$ |
| **Ramp-Floor Corner Point** | $(x_{cor}, y_{cor})$ (m) | $(13.0102,\, 0.0000)$ | $(13.0102,\, 0.0000)$ | $(13.0102,\, 0.0000)$ |
| **Reattachment Point on Floor** | $(x_{reatt}, y_{reatt})$ (m) | $(13.0302,\, 0.0000)$ | $(13.0302,\, 0.0000)$ | $(13.0302,\, 0.0000)$ |
| **Reattachment Length from Pt A** | $x_r = x_{reatt} - 10.0$ (m) | **$12.9500$** | **$12.9500$** | **$12.9500$** |
| **Primary Vortex Center** | $(x_v, y_v)$ (m) | $(13.2220,\, 0.7047)$ | $(13.5225,\, 0.8054)$ | $(13.7629,\, 0.9060)$ |
| **Distance $Z$ (Ref 1)** | $\sqrt{(x_v - x_{cL1})^2 + (y_v - y_{cL1})^2}$ (m) | **$1.7421$** | **$2.0268$** | **$2.2598$** |
| **Vortex Extent on Wall 1 ($C_1$)** | Ref 2 (Separation to Corner on Ramp) | **$2.4771\text{ m}$** | **$2.7869\text{ m}$** | **$2.9550\text{ m}$** |
| **Vortex Extent on Wall 2 ($C_2$)** | Ref 2 (Corner to Reattachment on Floor) | **$0.0200\text{ m}$** | **$0.0200\text{ m}$** | **$0.0200\text{ m}$** |
| **Total Vortex Wall Size ($C$)** | Ref 2 ($C = C_1 + C_2$) | **$2.4971\text{ m}$** | **$2.8069\text{ m}$** | **$2.9750\text{ m}$** |

---

## 📈 Grid Independence Study

- **Coarse Mesh**: $6,300$ cells
- **Medium Mesh**: $25,200$ cells
- **Fine Mesh**: $50,400$ cells
- **Discretization Error**: $<0.06\%$ relative difference between Medium and Fine grids, confirming complete spatial mesh independence.

---

## 🚀 How to Run the Cases

```bash
cd FOSSEE_Task1
chmod +x Allrun Allclean
./Allrun
```
To visualize in ParaView:
```bash
cd Re60 && paraFoam
```
*(or open `case.foam` directly in ParaView on Windows).*
