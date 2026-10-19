import os
import glob
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def parse_openfoam_mesh_and_fields(case_dir):
    """
    Parses OpenFOAM polyMesh (points, faces, owner, neighbour) and latest time step fields (U, p)
    to compute cell center coordinates and values.
    """
    poly_dir = os.path.join(case_dir, "constant", "polyMesh")
    if not os.path.exists(poly_dir):
        poly_dir = os.path.join(case_dir, "0", "polyMesh")
    
    # Read points
    points = []
    with open(os.path.join(poly_dir, "points"), "r") as f:
        lines = f.readlines()
        in_list = False
        for line in lines:
            line = line.strip()
            if line == "(":
                in_list = True
                continue
            if in_list:
                if line == ")":
                    break
                if line.startswith("(") and line.endswith(")"):
                    coords = [float(v) for v in line[1:-1].split()]
                    points.append(coords)
    points = np.array(points)

    # Read faces
    faces = []
    with open(os.path.join(poly_dir, "faces"), "r") as f:
        lines = f.readlines()
        in_list = False
        for line in lines:
            line = line.strip()
            if line == "(":
                in_list = True
                continue
            if in_list:
                if line == ")":
                    break
                if "(" in line and ")" in line:
                    idx_str = line[line.find("(")+1:line.find(")")]
                    f_pts = [int(v) for v in idx_str.split()]
                    faces.append(f_pts)

    # Read owner
    owner = []
    with open(os.path.join(poly_dir, "owner"), "r") as f:
        lines = f.readlines()
        in_list = False
        for line in lines:
            line = line.strip()
            if line == "(":
                in_list = True
                continue
            if in_list:
                if line == ")":
                    break
                if line.isdigit():
                    owner.append(int(line))
    owner = np.array(owner)

    num_cells = int(np.max(owner) + 1)
    
    # Calculate cell centers
    # For each cell, average the coordinates of its vertices
    cell_points = [set() for _ in range(num_cells)]
    for f_idx, c_idx in enumerate(owner):
        for p_idx in faces[f_idx]:
            cell_points[c_idx].add(p_idx)
            
    # Also check neighbour for internal faces
    neigh_file = os.path.join(poly_dir, "neighbour")
    if os.path.exists(neigh_file):
        with open(neigh_file, "r") as f:
            lines = f.readlines()
            in_list = False
            f_idx = 0
            for line in lines:
                line = line.strip()
                if line == "(":
                    in_list = True
                    continue
                if in_list:
                    if line == ")":
                        break
                    if line.isdigit():
                        c_idx = int(line)
                        for p_idx in faces[f_idx]:
                            cell_points[c_idx].add(p_idx)
                        f_idx += 1

    cell_centers = np.zeros((num_cells, 3))
    for i in range(num_cells):
        pts = points[list(cell_points[i])]
        cell_centers[i] = np.mean(pts, axis=0)

    # Find latest time directory
    time_dirs = []
    for d in os.listdir(case_dir):
        if os.path.isdir(os.path.join(case_dir, d)):
            try:
                t_val = float(d)
                if t_val > 0:
                    time_dirs.append((t_val, d))
            except ValueError:
                pass
    time_dirs.sort(key=lambda x: x[0])
    latest_time_dir = os.path.join(case_dir, time_dirs[-1][1])
    print(f"Reading fields from latest time: {time_dirs[-1][1]} in {case_dir}")

    # Read U
    U = []
    with open(os.path.join(latest_time_dir, "U"), "r") as f:
        lines = f.readlines()
        in_internal = False
        in_list = False
        for line in lines:
            line = line.strip()
            if "internalField" in line:
                in_internal = True
                if "nonuniform" not in line and "uniform" in line:
                    u_str = line[line.find("(")+1:line.find(")")]
                    u_val = [float(v) for v in u_str.split()]
                    U = [u_val] * num_cells
                    break
            if in_internal and line == "(":
                in_list = True
                continue
            if in_list:
                if line == ")":
                    break
                if line.startswith("(") and line.endswith(")"):
                    u_val = [float(v) for v in line[1:-1].split()]
                    U.append(u_val)
    U = np.array(U)

    # Read p
    p = []
    with open(os.path.join(latest_time_dir, "p"), "r") as f:
        lines = f.readlines()
        in_internal = False
        in_list = False
        for line in lines:
            line = line.strip()
            if "internalField" in line:
                in_internal = True
                if "nonuniform" not in line and "uniform" in line:
                    p_val = float(line.split()[-1].replace(";", ""))
                    p = [p_val] * num_cells
                    break
            if in_internal and line == "(":
                in_list = True
                continue
            if in_list:
                if line == ")":
                    break
                try:
                    p.append(float(line))
                except ValueError:
                    pass
    p = np.array(p)

    return cell_centers, U, p, time_dirs[-1][1]

def analyze_vortex_and_reattachment(cell_centers, U, p, Re_val, nu=1.602e-05):
    """
    Computes reattachment length (xr), vortex center, distance Z, and wall size C.
    """
    # 2D coordinates
    x = cell_centers[:, 0]
    y = cell_centers[:, 1]
    ux = U[:, 0]
    uy = U[:, 1]

    # Interpolate onto a fine regular 2D grid for streamline & vortex analysis
    grid_x, grid_y = np.meshgrid(
        np.linspace(0.0, 36.0, 1000),
        np.linspace(0.0, 5.0, 250)
    )

    grid_ux = griddata((x, y), ux, (grid_x, grid_y), method='linear', fill_value=0.0)
    grid_uy = griddata((x, y), uy, (grid_x, grid_y), method='linear', fill_value=0.0)

    # Mask out regions outside the domain:
    # Inlet bottom: x in [0, 10], y < 2.0
    # Ramp 1: x in [10, 13.0104], y < 2.0 - (x - 10)/tan(33.6) * 2 / dx1 = 2 - 2*(x-10)/3.0102
    dx1 = 2.0 / math.tan(math.radians(33.6))
    dx2 = dx1
    x_r1_end = 10.0 + dx1
    x_r2_start = 26.0 - dx2
    
    mask = np.zeros_like(grid_x, dtype=bool)
    mask[(grid_x <= 10.0) & (grid_y < 2.0)] = True
    
    mask_r1 = (grid_x > 10.0) & (grid_x < x_r1_end)
    y_r1 = 2.0 - (grid_x[mask_r1] - 10.0) * (2.0 / dx1)
    mask[mask_r1 & (grid_y < y_r1)] = True

    mask_r2 = (grid_x > x_r2_start) & (grid_x < 26.0)
    y_r2 = (grid_x[mask_r2] - x_r2_start) * (2.0 / dx2)
    mask[mask_r2 & (grid_y < y_r2)] = True

    mask[(grid_x >= 26.0) & (grid_y < 2.0)] = True

    grid_ux[mask] = 0.0
    grid_uy[mask] = 0.0

    # 1. Reattachment point along the bottom wall
    # In the cavity region (x in [10, 26]), sample velocity near the bottom wall
    # Bottom wall y_wall(x):
    # For x in [10, 13.0102]: y = 2.0 - (x-10)*(2/dx1)
    # For x in [13.0102, 22.9898]: y = 0.0
    # For x in [22.9898, 26.0]: y = (x - 22.9898)*(2/dx2)
    
    x_wall_samples = np.linspace(10.05, 25.0, 500)
    u_near_wall = []
    
    for xw in x_wall_samples:
        if xw <= x_r1_end:
            yw = 2.0 - (xw - 10.0) * (2.0 / dx1)
            # Sample slightly above wall normal
            y_sample = yw + 0.04
        elif xw <= x_r2_start:
            yw = 0.0
            y_sample = 0.04
        else:
            yw = (xw - x_r2_start) * (2.0 / dx2)
            y_sample = yw + 0.04
            
        u_val = griddata((x, y), ux, (xw, y_sample), method='linear', fill_value=0.0)
        u_near_wall.append(u_val)
        
    u_near_wall = np.array(u_near_wall)

    # Reattachment occurs where near-wall tangential velocity changes from negative (recirculation) to positive (forward flow)
    # Find zero crossing after the recirculation vortex minimum
    min_idx = np.argmin(u_near_wall)
    reattach_idx = None
    for i in range(min_idx, len(u_near_wall) - 1):
        if u_near_wall[i] <= 0 and u_near_wall[i+1] > 0:
            reattach_idx = i
            break
            
    if reattach_idx is not None:
        x0, x1_pt = x_wall_samples[reattach_idx], x_wall_samples[reattach_idx+1]
        u0, u1 = u_near_wall[reattach_idx], u_near_wall[reattach_idx+1]
        x_reattach = x0 - u0 * (x1_pt - x0) / (u1 - u0)
    else:
        # Fallback to closest point to 0
        x_reattach = x_wall_samples[np.argmin(np.abs(u_near_wall[min_idx:])) + min_idx]

    # Reattachment length measured from Point A (x_A = 10.0) in direction of flow
    x_r = x_reattach - 10.0

    # 2. Vortex Center (stream function / minimum ux in recirculation zone)
    # Compute 2D streamfunction psi: d(psi)/dy = ux, d(psi)/dx = -uy
    # psi(x, y) = int_0^y ux dy
    dy = grid_y[1, 0] - grid_y[0, 0]
    psi = np.cumsum(grid_ux, axis=0) * dy

    # Find vortex core in cavity region (x in [10, x_reattach], y in [0, 2.5])
    cav_mask = (grid_x >= 10.0) & (grid_x <= x_reattach + 1.0) & (grid_y >= 0.0) & (grid_y <= 2.5) & (~mask)
    psi_cav = np.where(cav_mask, psi, 0.0)
    
    # In the recirculation bubble, psi is negative (clockwise vortex)
    min_psi_idx = np.unravel_index(np.argmin(psi_cav), psi_cav.shape)
    vortex_center_x = grid_x[min_psi_idx]
    vortex_center_y = grid_y[min_psi_idx]

    # 3. Distance Z from midpoint of L1 to vortex center (Reference 1)
    # Midpoint of L1:
    cL1_x = 10.0 + 0.5 * dx1 # 11.5051
    cL1_y = 1.0
    Z = math.sqrt((vortex_center_x - cL1_x)**2 + (vortex_center_y - cL1_y)**2)

    # 4. Vortex size relative to channel walls C = C1 + C2 (Reference 2)
    # C1: length along Wall 1 (Ramp L1) covered by the vortex
    # Ramp 1 starts at Point A (10, 2) and ends at (10+dx1, 0)
    # Total ramp length L1 = sqrt(dx1^2 + 2^2) = 3.6141 m
    # Vortex covers ramp from separation point (Point A) to bottom of ramp:
    C1 = math.sqrt(dx1**2 + 2.0**2) # Full ramp 1 = 3.6141 m (or distance to detachment)
    
    # C2: distance along Wall 2 (Cavity bottom floor) covered by the primary vortex
    # Wall 2 starts at x = 10 + dx1 (13.0102 m). The vortex extends along Wall 2 up to x_reattach
    if x_reattach >= x_r1_end:
        C2 = x_reattach - x_r1_end
    else:
        C2 = 0.0
        # If reattachment happens on the ramp itself:
        C1 = math.sqrt((x_reattach - 10.0)**2 + (2.0 - (2.0 - (x_reattach - 10.0)*(2.0/dx1)))**2)

    C = C1 + C2

    metrics = {
        "Re": Re_val,
        "x_reattach": x_reattach,
        "x_r": x_r,
        "vortex_center": (vortex_center_x, vortex_center_y),
        "Z": Z,
        "C1": C1,
        "C2": C2,
        "C": C,
        "cL1": (cL1_x, cL1_y)
    }

    return metrics, (grid_x, grid_y, grid_ux, grid_uy, psi, mask, x_wall_samples, u_near_wall)

if __name__ == "__main__":
    base_dir = "/mnt/c/Users/moham/Desktop/FOSSEE/FOSSEE_Task1"
    re_cases = [("Re60", 60), ("Re80", 80), ("Re100", 100)]
    
    results = {}
    for cname, re_val in re_cases:
        cdir = os.path.join(base_dir, cname)
        if os.path.exists(cdir):
            centers, U, p, ltime = parse_openfoam_mesh_and_fields(cdir)
            metrics, grid_data = analyze_vortex_and_reattachment(centers, U, p, re_val)
            results[cname] = metrics
            print(f"\n================ Results for {cname} (Re = {re_val}) ================")
            print(f"  Reattachment Location (x): {metrics['x_reattach']:.4f} m")
            print(f"  Reattachment Length (x_r): {metrics['x_r']:.4f} m (measured from Point A)")
            print(f"  Primary Vortex Center: ({metrics['vortex_center'][0]:.4f}, {metrics['vortex_center'][1]:.4f}) m")
            print(f"  Distance Z (Midpoint L1 to Vortex Center): {metrics['Z']:.4f} m")
            print(f"  Vortex Wall Dimension C1 (along Wall 1): {metrics['C1']:.4f} m")
            print(f"  Vortex Wall Dimension C2 (along Wall 2): {metrics['C2']:.4f} m")
            print(f"  Total Vortex Wall Size C (C1 + C2): {metrics['C']:.4f} m")
