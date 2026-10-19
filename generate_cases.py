import os
import math

def get_blockmesh_dict(nx_in=100, nx_r1=50, nx_cav=120, nx_r2=50, nx_out=100, ny=60):
    dx1 = 2.0 / math.tan(math.radians(33.6)) # 3.010205...
    dx2 = dx1
    X = 36.0 - (10.0 + dx1 + dx2 + 10.0) # 9.979589...
    
    x0 = 0.0
    x1 = 10.0
    x2 = x1 + dx1
    x3 = x2 + X
    x4 = x3 + dx2
    x5 = 36.0
    
    y0 = 0.0
    y1 = 2.0
    y2 = 5.0
    
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
    location    "system";
    object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
    // z = {z0}
    ({x0:.6f} {y1:.6f} {z0:.6f}) // 0: Inlet bottom-left
    ({x1:.6f} {y1:.6f} {z0:.6f}) // 1: Inlet bottom-right (Point A)
    ({x2:.6f} {y0:.6f} {z0:.6f}) // 2: Ramp 1 bottom
    ({x3:.6f} {y0:.6f} {z0:.6f}) // 3: Ramp 2 bottom
    ({x4:.6f} {y1:.6f} {z0:.6f}) // 4: Ramp 2 top / Outlet bottom-left
    ({x5:.6f} {y1:.6f} {z0:.6f}) // 5: Outlet bottom-right
    ({x0:.6f} {y2:.6f} {z0:.6f}) // 6: Inlet top-left
    ({x1:.6f} {y2:.6f} {z0:.6f}) // 7: Inlet top-right
    ({x2:.6f} {y2:.6f} {z0:.6f}) // 8: Cavity top-left
    ({x3:.6f} {y2:.6f} {z0:.6f}) // 9: Cavity top-right
    ({x4:.6f} {y2:.6f} {z0:.6f}) // 10: Outlet top-left
    ({x5:.6f} {y2:.6f} {z0:.6f}) // 11: Outlet top-right

    // z = {z1}
    ({x0:.6f} {y1:.6f} {z1:.6f}) // 12
    ({x1:.6f} {y1:.6f} {z1:.6f}) // 13
    ({x2:.6f} {y0:.6f} {z1:.6f}) // 14
    ({x3:.6f} {y0:.6f} {z1:.6f}) // 15
    ({x4:.6f} {y1:.6f} {z1:.6f}) // 16
    ({x5:.6f} {y1:.6f} {z1:.6f}) // 17
    ({x0:.6f} {y2:.6f} {z1:.6f}) // 18
    ({x1:.6f} {y2:.6f} {z1:.6f}) // 19
    ({x2:.6f} {y2:.6f} {z1:.6f}) // 20
    ({x3:.6f} {y2:.6f} {z1:.6f}) // 21
    ({x4:.6f} {y2:.6f} {z1:.6f}) // 22
    ({x5:.6f} {y2:.6f} {z1:.6f}) // 23
);

blocks
(
    // Block 0: Inlet
    hex (0 1 7 6 12 13 19 18) ({nx_in} {ny} 1) simpleGrading (1 1 1)
    
    // Block 1: Ramp 1 (Wall 1 expansion)
    hex (1 2 8 7 13 14 20 19) ({nx_r1} {ny} 1) simpleGrading (1 1 1)
    
    // Block 2: Cavity Floor (Wall 2)
    hex (2 3 9 8 14 15 21 20) ({nx_cav} {ny} 1) simpleGrading (1 1 1)
    
    // Block 3: Ramp 2 (Contraction)
    hex (3 4 10 9 15 16 22 21) ({nx_r2} {ny} 1) simpleGrading (1 1 1)
    
    // Block 4: Outlet
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
    bottomWall
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
            // Front (z = 1)
            (12 13 19 18)
            (13 14 20 19)
            (14 15 21 20)
            (15 16 22 21)
            (16 17 23 22)
            // Back (z = 0)
            (0 6 7 1)
            (1 7 8 2)
            (2 8 9 3)
            (3 9 10 4)
            (4 10 11 5)
        );
    }}
);

// ************************************************************************* //
"""
    return content

def get_control_dict(end_time=2000, write_interval=500):
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
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

solver          incompressibleFluid;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         {end_time};

deltaT          1;

writeControl    timeStep;

writeInterval   {write_interval};

purgeWrite      2;

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

// ************************************************************************* //
"""

def get_fv_schemes():
    return """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
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
}

divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
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

def get_fv_solution():
    return """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
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

    "(U|k|epsilon|omega|f|v2)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-07;
        relTol          0.01;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent      yes;

    residualControl
    {
        p               1e-6;
        U               1e-6;
    }
}

relaxationFactors
{
    equations
    {
        U               0.9;
        p               0.3;
    }
}

// ************************************************************************* //
"""

def get_physical_properties(nu=1.602e-05):
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
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

nu              {nu:.6e};

// ************************************************************************* //
"""

def get_momentum_transport():
    return """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
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

def get_u_field(u_in):
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    location    "0";
    object      U;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform ({u_in:.8e} 0 0);

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform ({u_in:.8e} 0 0);
    }}

    outlet
    {{
        type            zeroGradient;
    }}

    topWall
    {{
        type            noSlip;
    }}

    bottomWall
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

def get_p_field():
    return """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
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

    bottomWall
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

def create_case(case_dir, u_in, nx_in=100, nx_r1=50, nx_cav=120, nx_r2=50, nx_out=100, ny=60, end_time=2000):
    os.makedirs(os.path.join(case_dir, "0"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)

    with open(os.path.join(case_dir, "0", "U"), "w") as f:
        f.write(get_u_field(u_in))
    with open(os.path.join(case_dir, "0", "p"), "w") as f:
        f.write(get_p_field())

    with open(os.path.join(case_dir, "constant", "physicalProperties"), "w") as f:
        f.write(get_physical_properties())
    with open(os.path.join(case_dir, "constant", "momentumTransport"), "w") as f:
        f.write(get_momentum_transport())

    with open(os.path.join(case_dir, "system", "blockMeshDict"), "w") as f:
        f.write(get_blockmesh_dict(nx_in, nx_r1, nx_cav, nx_r2, nx_out, ny))
    with open(os.path.join(case_dir, "system", "controlDict"), "w") as f:
        f.write(get_control_dict(end_time=end_time))
    with open(os.path.join(case_dir, "system", "fvSchemes"), "w") as f:
        f.write(get_fv_schemes())
    with open(os.path.join(case_dir, "system", "fvSolution"), "w") as f:
        f.write(get_fv_solution())

if __name__ == "__main__":
    base_dir = "/mnt/c/Users/moham/Desktop/FOSSEE/FOSSEE_Task1"
    os.makedirs(base_dir, exist_ok=True)
    
    nu = 1.602e-05
    H_in = 3.0
    cases = {
        "Re60": 60 * nu / H_in,
        "Re80": 80 * nu / H_in,
        "Re100": 100 * nu / H_in
    }
    
    for case_name, u_val in cases.items():
        cdir = os.path.join(base_dir, case_name)
        print(f"Creating case {case_name} with U_in = {u_val:.6e} m/s in {cdir}")
        create_case(cdir, u_val)
        
    print("Cases created successfully!")
