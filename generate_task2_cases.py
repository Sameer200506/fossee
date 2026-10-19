import os
import math
import shutil

base_task2_dir = "/mnt/c/Users/moham/Desktop/FOSSEE/FOSSEE_Task2"

# Geometry parameters (D = 1.0 m)
D = 1.0
H_in = 3.0 * D
H_cav = 5.0 * D
L_in = 10.0 * D
L_out = 10.0 * D
L_top = 36.0 * D
theta_deg = 33.6
theta_rad = math.radians(theta_deg)
step_h = H_cav - H_in  # 2.0 m

L1 = step_h / math.sin(theta_rad)  # 3.61408 m
L2 = L1
dx1 = step_h / math.tan(theta_rad)  # 3.01024 m
dx2 = dx1
X = L_top - (L_in + dx1 + dx2 + L_out)  # 9.97951 m

# Fluid properties (Air at 27 C = 300.15 K)
T_in = 300.15  # K (27 C)
nu = 0.05      # m2/s (specified)
Pr = 0.71      # Prandtl number
alpha_thermal = nu / Pr  # 0.0704225 m2/s
rho = 1.176    # kg/m3
cp = 1005.0    # J/(kg K)
k_fluid = rho * cp * alpha_thermal  # 83.23 W/(m K)

Re = 50.0
U_in = Re * nu / H_in  # 50 * 0.05 / 3.0 = 0.833333 m/s

thermal_cases = [
    {"name": "Case1_77C",  "Tw_C": 77.0,  "Tw_K": 350.15, "label": "Tw = 77°C (350.15 K)"},
    {"name": "Case2_127C", "Tw_C": 127.0, "Tw_K": 400.15, "label": "Tw = 127°C (400.15 K)"},
    {"name": "Case3_177C", "Tw_C": 177.0, "Tw_K": 450.15, "label": "Tw = 177°C (450.15 K)"}
]

mesh_resolutions = {
    "coarse": {"nx_in": 50,  "nx_r1": 25, "nx_cav": 60,  "nx_r2": 25, "nx_out": 50,  "ny": 30, "cells": 6300},
    "medium": {"nx_in": 100, "nx_r1": 50, "nx_cav": 120, "nx_r2": 50, "nx_out": 100, "ny": 60, "cells": 25200},
    "fine":   {"nx_in": 150, "nx_r1": 75, "nx_cav": 180, "nx_r2": 75, "nx_out": 150, "ny": 80, "cells": 50400}
}

def generate_blockMeshDict(res):
    nx_in = res["nx_in"]
    nx_r1 = res["nx_r1"]
    nx_cav = res["nx_cav"]
    nx_r2 = res["nx_r2"]
    nx_out = res["nx_out"]
    ny = res["ny"]
    
    x0 = 0.0
    x1 = L_in
    x2 = L_in + dx1
    x3 = L_in + dx1 + X
    x4 = L_in + dx1 + X + dx2
    x5 = L_top
    
    z0 = 0.0
    z1 = 1.0
    
    content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
    // Bottom plane vertices (z = 0)
    ({x0:.6f} 2.000000 {z0:.6f}) // 0
    ({x1:.6f} 2.000000 {z0:.6f}) // 1 (Point A)
    ({x2:.6f} 0.000000 {z0:.6f}) // 2
    ({x3:.6f} 0.000000 {z0:.6f}) // 3
    ({x4:.6f} 2.000000 {z0:.6f}) // 4
    ({x5:.6f} 2.000000 {z0:.6f}) // 5

    // Top plane vertices (z = 0, y = 5)
    ({x0:.6f} 5.000000 {z0:.6f}) // 6
    ({x1:.6f} 5.000000 {z0:.6f}) // 7
    ({x2:.6f} 5.000000 {z0:.6f}) // 8
    ({x3:.6f} 5.000000 {z0:.6f}) // 9
    ({x4:.6f} 5.000000 {z0:.6f}) // 10
    ({x5:.6f} 5.000000 {z0:.6f}) // 11

    // Bottom plane vertices (z = 1)
    ({x0:.6f} 2.000000 {z1:.6f}) // 12
    ({x1:.6f} 2.000000 {z1:.6f}) // 13
    ({x2:.6f} 0.000000 {z1:.6f}) // 14
    ({x3:.6f} 0.000000 {z1:.6f}) // 15
    ({x4:.6f} 2.000000 {z1:.6f}) // 16
    ({x5:.6f} 2.000000 {z1:.6f}) // 17

    // Top plane vertices (z = 1, y = 5)
    ({x0:.6f} 5.000000 {z1:.6f}) // 18
    ({x1:.6f} 5.000000 {z1:.6f}) // 19
    ({x2:.6f} 5.000000 {z1:.6f}) // 20
    ({x3:.6f} 5.000000 {z1:.6f}) // 21
    ({x4:.6f} 5.000000 {z1:.6f}) // 22
    ({x5:.6f} 5.000000 {z1:.6f}) // 23
);

blocks
(
    hex (0 1 7 6 12 13 19 18) ({nx_in} {ny} 1) simpleGrading (1 1 1)
    hex (1 2 8 7 13 14 20 19) ({nx_r1} {ny} 1) simpleGrading (1 1 1)
    hex (2 3 9 8 14 15 21 20) ({nx_cav} {ny} 1) simpleGrading (1 1 1)
    hex (3 4 10 9 15 16 22 21) ({nx_r2} {ny} 1) simpleGrading (1 1 1)
    hex (4 5 11 10 16 17 23 22) ({nx_out} {ny} 1) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    inlet
    {{
        type patch;
        faces
        (
            (0 6 18 12)
        );
    }}
    outlet
    {{
        type patch;
        faces
        (
            (5 17 23 11)
        );
    }}
    topWall
    {{
        type wall;
        faces
        (
            (6 7 19 18)
            (7 8 20 19)
            (8 9 21 20)
            (9 10 22 21)
            (10 11 23 22)
        );
    }}
    heatedBottomWall
    {{
        type wall;
        faces
        (
            (0 12 13 1)
            (1 13 14 2)
            (2 14 15 3)
            (3 15 16 4)
            (4 16 17 5)
        );
    }}
    frontAndBack
    {{
        type empty;
        faces
        (
            (0 1 7 6)
            (1 2 8 7)
            (2 3 9 8)
            (3 4 10 9)
            (4 5 11 10)
            (12 18 19 13)
            (13 19 20 14)
            (14 20 21 15)
            (15 21 22 16)
            (16 22 23 17)
        );
    }}
);

// ************************************************************************* //
"""
    return content

def generate_case_files(case_dir, Tw_K, res):
    os.makedirs(os.path.join(case_dir, "0"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    
    # 1. 0/U
    u_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       volVectorField;
    location    "0";
    object      U;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform ({U_in:.6f} 0 0);

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform ({U_in:.6f} 0 0);
    }}

    outlet
    {{
        type            zeroGradient;
    }}

    topWall
    {{
        type            noSlip;
    }}

    heatedBottomWall
    {{
        type            noSlip;
    }}

    frontAndBack
    {{
        type            empty;
    }}
}}

// ************************************************************************* //
"""
    with open(os.path.join(case_dir, "0", "U"), "w") as f:
        f.write(u_content)

    # 2. 0/p
    p_content = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       volScalarField;
    location    "0";
    object      p;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            zeroGradient;
    }

    outlet
    {
        type            fixedValue;
        value           uniform 0;
    }

    topWall
    {
        type            zeroGradient;
    }

    heatedBottomWall
    {
        type            zeroGradient;
    }

    frontAndBack
    {
        type            empty;
    }
}

// ************************************************************************* //
"""
    with open(os.path.join(case_dir, "0", "p"), "w") as f:
        f.write(p_content)

    # 3. 0/T (Temperature field)
    t_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       volScalarField;
    location    "0";
    object      T;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];

internalField   uniform {T_in:.2f};

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {T_in:.2f};
    }}

    outlet
    {{
        type            zeroGradient;
    }}

    topWall
    {{
        type            zeroGradient;
    }}

    heatedBottomWall
    {{
        type            fixedValue;
        value           uniform {Tw_K:.2f};
    }}

    frontAndBack
    {{
        type            empty;
    }}
}}

// ************************************************************************* //
"""
    with open(os.path.join(case_dir, "0", "T"), "w") as f:
        f.write(t_content)

    # 4. constant/physicalProperties
    prop_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      physicalProperties;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

viscosityModel  constant;

nu              {nu:.6f};

// ************************************************************************* //
"""
    with open(os.path.join(case_dir, "constant", "physicalProperties"), "w") as f:
        f.write(prop_content)

    # 5. constant/momentumTransport
    mom_content = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  13
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      momentumTransport;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simulationType  laminar;

// ************************************************************************* //
"""
    with open(os.path.join(case_dir, "constant", "momentumTransport"), "w") as f:
        f.write(mom_content)

    # 6. system/blockMeshDict
    with open(os.path.join(case_dir, "system", "blockMeshDict"), "w") as f:
        f.write(generate_blockMeshDict(res))

    # 7. system/controlDict
    ctrl_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     foamRun;

solver          incompressibleFluid;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         1000;

deltaT          1;

writeControl    timeStep;

writeInterval   100;

purgeWrite      2;

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

functions
{{
    scalarTransport
    {{
        type            scalarTransport;
        libs            ("libsolverFunctionObjects.so");
        field           T;
        diffusivity     constant;
        D               {alpha_thermal:.7f};
        autoSchemes     false;
        schemesField    T;
    }}
}}

// ************************************************************************* //
"""
    with open(os.path.join(case_dir, "system", "controlDict"), "w") as f:
        f.write(ctrl_content)

    # 8. system/fvSchemes
    schemes_content = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  13
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
    grad(p)         Gauss linear;
    grad(U)         Gauss linear;
    grad(T)         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      Gauss linearUpwind grad(U);
    div(phi,T)      Gauss linearUpwind grad(T);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

// ************************************************************************* //
"""
    with open(os.path.join(case_dir, "system", "fvSchemes"), "w") as f:
        f.write(schemes_content)

    # 9. system/fvSolution
    solution_content = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  13
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-07;
        relTol          0.01;
        smoother        GaussSeidel;
    }

    U
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-08;
        relTol          0.01;
    }

    T
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-08;
        relTol          0.01;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 1;
    consistent      yes;

    residualControl
    {
        p               1e-6;
        U               1e-6;
        T               1e-6;
    }
}

relaxationFactors
{
    equations
    {
        U               0.9;
        T               0.9;
    }
    fields
    {
        p               0.3;
    }
}

// ************************************************************************* //
"""
    with open(os.path.join(case_dir, "system", "fvSolution"), "w") as f:
        f.write(solution_content)

    with open(os.path.join(case_dir, "case.foam"), "w") as f:
        f.write("")

if __name__ == "__main__":
    os.makedirs(base_task2_dir, exist_ok=True)
    grid_dir = os.path.join(base_task2_dir, "grid_independence")
    os.makedirs(grid_dir, exist_ok=True)
    
    for c in thermal_cases:
        cdir = os.path.join(base_task2_dir, c["name"])
        generate_case_files(cdir, c["Tw_K"], mesh_resolutions["medium"])
        print(f"Generated case: {c['name']} (Tw = {c['Tw_C']} C)")
        
    for c in thermal_cases:
        for m_name in ["coarse", "fine"]:
            cname = f"{c['name']}_{m_name}"
            cdir = os.path.join(grid_dir, cname)
            generate_case_files(cdir, c["Tw_K"], mesh_resolutions[m_name])
            print(f"Generated grid case: {cname}")

    allrun_content = """#!/bin/bash
set -e
echo "=========================================================="
echo "    Running FOSSEE Task 2: Heat Transfer CFD Simulations   "
echo "=========================================================="
source /opt/openfoam13/etc/bashrc

CASES=("Case1_77C" "Case2_127C" "Case3_177C")

for c in "${CASES[@]}"; do
    echo ""
    echo ">>> Meshing and running $c..."
    cd "$c"
    blockMesh > log.blockMesh 2>&1
    checkMesh > log.checkMesh 2>&1
    foamRun -solver incompressibleFluid > log.foamRun 2>&1
    cd ..
    echo ">>> $c completed successfully!"
done

echo ""
echo "=========================================================="
echo "    Running Task 2 Grid Independence Study Cases          "
echo "=========================================================="
cd grid_independence
for c in "${CASES[@]}"; do
    for m in "coarse" "fine"; do
        cname="${c}_${m}"
        echo ">>> Running $cname..."
        cd "$cname"
        blockMesh > log.blockMesh 2>&1
        foamRun -solver incompressibleFluid > log.foamRun 2>&1
        cd ..
        echo ">>> $cname completed successfully!"
    done
done
cd ..

echo ""
echo "All Task 2 simulations finished successfully!"
"""
    with open(os.path.join(base_task2_dir, "Allrun"), "w") as f:
        f.write(allrun_content)

    allclean_content = """#!/bin/bash
echo "Cleaning all Task 2 simulation cases..."
CASES=("Case1_77C" "Case2_127C" "Case3_177C")
for c in "${CASES[@]}"; do
    cd "$c"
    rm -rf [1-9]* [0-9]*[0-9] log.* constant/polyMesh
    cd ..
done
cd grid_independence
for c in "${CASES[@]}"; do
    for m in "coarse" "fine"; do
        cd "${c}_${m}"
        rm -rf [1-9]* [0-9]*[0-9] log.* constant/polyMesh
        cd ..
    done
done
cd ..
echo "Cleanup completed!"
"""
    with open(os.path.join(base_task2_dir, "Allclean"), "w") as f:
        f.write(allclean_content)

    print("Task 2 case generation complete!")
