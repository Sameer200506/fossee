import os
import sys
import math
import subprocess
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, RegularGridInterpolator

# ==============================================================================
# 1. GEOMETRIC PARAMETERS AND FORMULATION
# ==============================================================================
D = 1.0  # Reference dimension (chosen from range 0.1 - 10)
H_in = 3.0 * D
H_cav = 5.0 * D
L_in = 10.0 * D
L_out = 10.0 * D
L_top = 36.0 * D
theta_deg = 33.6
theta_rad = math.radians(theta_deg)
step_h = H_cav - H_in  # 2.0 * D

# Analytical geometric quantities
L1 = step_h / math.sin(theta_rad)
L2 = L1
dx1 = step_h / math.tan(theta_rad)
dx2 = dx1
X = L_top - (L_in + dx1 + dx2 + L_out)
alpha = theta_deg
alpha_int = 180.0 - theta_deg

# Coordinates of midpoints/centres of L1 and L2 (Origin at domain bottom-left: (0,0))
cL1_x = L_in + 0.5 * dx1
cL1_y = 1.0 * D

cL2_x = (L_top - L_out) - 0.5 * dx2
cL2_y = 1.0 * D

nu = 1.602e-05

re_dict = {
    60: {"name": "Re60", "Re": 60, "U": 60 * nu / H_in},
    80: {"name": "Re80", "Re": 80, "U": 80 * nu / H_in},
    100: {"name": "Re100", "Re": 100, "U": 100 * nu / H_in}
}

mesh_resolutions = {
    "coarse": {"nx_in": 50, "nx_r1": 25, "nx_cav": 60, "nx_r2": 25, "nx_out": 50, "ny": 30, "cells": 6300},
    "medium": {"nx_in": 100, "nx_r1": 50, "nx_cav": 120, "nx_r2": 50, "nx_out": 100, "ny": 60, "cells": 25200},
    "fine":   {"nx_in": 150, "nx_r1": 75, "nx_cav": 180, "nx_r2": 75, "nx_out": 150, "ny": 80, "cells": 50400}
}

def parse_openfoam_results(case_dir):
    poly_dir = os.path.join(case_dir, "constant", "polyMesh")
    if not os.path.exists(poly_dir):
        poly_dir = os.path.join(case_dir, "0", "polyMesh")

    points = []
    with open(os.path.join(poly_dir, "points"), "r") as f:
        in_list = False
        for line in f:
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

    faces = []
    with open(os.path.join(poly_dir, "faces"), "r") as f:
        in_list = False
        for line in f:
            line = line.strip()
            if line == "(":
                in_list = True
                continue
            if in_list:
                if line == ")":
                    break
                if "(" in line and ")" in line:
                    idx_str = line[line.find("(")+1:line.find(")")]
                    faces.append([int(v) for v in idx_str.split()])

    owner = []
    with open(os.path.join(poly_dir, "owner"), "r") as f:
        in_list = False
        for line in f:
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

    cell_points = [set() for _ in range(num_cells)]
    for f_idx, c_idx in enumerate(owner):
        for p_idx in faces[f_idx]:
            cell_points[c_idx].add(p_idx)

    neigh_file = os.path.join(poly_dir, "neighbour")
    if os.path.exists(neigh_file):
        with open(neigh_file, "r") as f:
            in_list = False
            f_idx = 0
            for line in f:
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

    U = []
    with open(os.path.join(latest_time_dir, "U"), "r") as f:
        in_internal = False
        in_list = False
        for line in f:
            line = line.strip()
            if "internalField" in line:
                in_internal = True
                if "nonuniform" not in line and "uniform" in line:
                    u_str = line[line.find("(")+1:line.find(")")]
                    U = [[float(v) for v in u_str.split()]] * num_cells
                    break
            if in_internal and line == "(":
                in_list = True
                continue
            if in_list:
                if line == ")":
                    break
                if line.startswith("(") and line.endswith(")"):
                    U.append([float(v) for v in line[1:-1].split()])
    U = np.array(U)

    p = []
    with open(os.path.join(latest_time_dir, "p"), "r") as f:
        in_internal = False
        in_list = False
        for line in f:
            line = line.strip()
            if "internalField" in line:
                in_internal = True
                if "nonuniform" not in line and "uniform" in line:
                    p = [float(line.split()[-1].replace(";", ""))] * num_cells
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

def analyze_case_metrics(cell_centers, U, p, Re_val):
    x = cell_centers[:, 0]
    y = cell_centers[:, 1]
    ux = U[:, 0]
    uy = U[:, 1]

    nx_g, ny_g = 600, 150
    gx_vec = np.linspace(0.0, 36.0, nx_g)
    gy_vec = np.linspace(0.0, 5.0, ny_g)
    grid_x, grid_y = np.meshgrid(gx_vec, gy_vec)

    grid_ux = griddata((x, y), ux, (grid_x, grid_y), method='linear', fill_value=0.0)
    grid_uy = griddata((x, y), uy, (grid_x, grid_y), method='linear', fill_value=0.0)

    x_r1_end = 10.0 + dx1
    x_r2_start = 26.0 - dx2

    mask = np.zeros_like(grid_x, dtype=bool)
    mask[(grid_x <= 10.0) & (grid_y < 2.0)] = True
    
    y_r1 = 2.0 - (grid_x - 10.0) * (2.0 / dx1)
    mask[(grid_x > 10.0) & (grid_x < x_r1_end) & (grid_y < y_r1)] = True

    y_r2 = (grid_x - x_r2_start) * (2.0 / dx2)
    mask[(grid_x > x_r2_start) & (grid_x < 26.0) & (grid_y < y_r2)] = True

    mask[(grid_x >= 26.0) & (grid_y < 2.0)] = True

    grid_ux[mask] = 0.0
    grid_uy[mask] = 0.0

    # 1. Reattachment along bottom wall
    x_wall_samples = np.linspace(10.05, 25.0, 300)
    y_wall_samples = np.zeros_like(x_wall_samples)
    for i, xw in enumerate(x_wall_samples):
        if xw <= x_r1_end:
            yw = 2.0 - (xw - 10.0) * (2.0 / dx1)
            y_wall_samples[i] = yw + 0.03
        elif xw <= x_r2_start:
            y_wall_samples[i] = 0.03
        else:
            yw = (xw - x_r2_start) * (2.0 / dx2)
            y_wall_samples[i] = yw + 0.03

    rgi_ux = RegularGridInterpolator((gy_vec, gx_vec), grid_ux, bounds_error=False, fill_value=0.0)
    sample_pts = np.vstack([y_wall_samples, x_wall_samples]).T
    u_near_wall = rgi_ux(sample_pts)

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
        x_reattach = x_wall_samples[np.argmin(np.abs(u_near_wall[min_idx:])) + min_idx]

    x_r = x_reattach - 10.0  # Measured from Point A (x=10.0)

    # 2. 2D Stream function and Vortex Center
    dy = grid_y[1, 0] - grid_y[0, 0]
    psi = np.cumsum(grid_ux, axis=0) * dy

    cav_mask = (grid_x >= 10.0) & (grid_x <= x_reattach + 1.0) & (grid_y >= 0.0) & (grid_y <= 2.5) & (~mask)
    psi_cav = np.where(cav_mask, psi, 0.0)
    
    min_psi_idx = np.unravel_index(np.argmin(psi_cav), psi_cav.shape)
    vortex_center_x = grid_x[min_psi_idx]
    vortex_center_y = grid_y[min_psi_idx]

    # 3. Distance Z (Midpoint of L1 to vortex center)
    Z = math.sqrt((vortex_center_x - cL1_x)**2 + (vortex_center_y - cL1_y)**2)

    # 4. Exact Vortex Wall Extents (Reference 2: C1 and C2)
    # Unit vectors for Wall 1 (Ramp L1)
    L1_total = math.sqrt(dx1**2 + 2.0**2)
    tx = dx1 / L1_total
    ty = -2.0 / L1_total
    nx = 2.0 / L1_total
    ny = dx1 / L1_total
    dn = 0.02

    # Sample along Wall 1
    s_vals = np.linspace(0.02, L1_total - 0.02, 500)
    xw = 10.0 + s_vals * tx
    yw = 2.0 + s_vals * ty
    xf = xw + dn * nx
    yf = yw + dn * ny
    ufx = griddata((x, y), ux, (xf, yf), method='linear', fill_value=0.0)
    ufy = griddata((x, y), uy, (xf, yf), method='linear', fill_value=0.0)
    ut_wall1 = ufx * tx + ufy * ty

    sep_idx = None
    for i in range(len(ut_wall1) - 1):
        if ut_wall1[i] >= 0 and ut_wall1[i+1] < 0:
            sep_idx = i
            break
            
    if sep_idx is not None:
        s0, s1 = s_vals[sep_idx], s_vals[sep_idx+1]
        u0, u1 = ut_wall1[sep_idx], ut_wall1[sep_idx+1]
        s_sep = s0 - u0 * (s1 - s0) / (u1 - u0)
        C1 = L1_total - s_sep
    else:
        neg_indices = np.where(ut_wall1 < 0)[0]
        if len(neg_indices) > 0:
            s_sep = s_vals[neg_indices[0]]
            C1 = L1_total - s_sep
        else:
            s_sep = 0.0
            C1 = L1_total

    x_sep = 10.0 + s_sep * tx
    y_sep = 2.0 + s_sep * ty

    # Sample along Wall 2 (Cavity floor)
    x_floor_vals = np.linspace(x_r1_end + 0.02, 25.0, 500)
    yf_floor = np.full_like(x_floor_vals, dn)
    u_floor = griddata((x, y), ux, (x_floor_vals, yf_floor), method='linear', fill_value=0.0)

    reattach_floor_idx = None
    for i in range(len(u_floor) - 1):
        if u_floor[i] <= 0 and u_floor[i+1] > 0:
            reattach_floor_idx = i
            break

    if reattach_floor_idx is not None:
        x0, x1_pt = x_floor_vals[reattach_floor_idx], x_floor_vals[reattach_floor_idx+1]
        u0, u1 = u_floor[reattach_floor_idx], u_floor[reattach_floor_idx+1]
        x_reattach_floor = x0 - u0 * (x1_pt - x0) / (u1 - u0)
        C2 = x_reattach_floor - x_r1_end
    else:
        x_reattach_floor = x_r1_end + 0.02
        C2 = 0.02

    C = C1 + C2

    metrics = {
        "Re": Re_val,
        "x_reattach": float(x_reattach),
        "x_r": float(x_r),
        "vortex_center": (float(vortex_center_x), float(vortex_center_y)),
        "Z": float(Z),
        "s_sep": float(s_sep),
        "x_sep": float(x_sep),
        "y_sep": float(y_sep),
        "C1": float(C1),
        "x_corner": float(x_r1_end),
        "x_reattach_floor": float(x_reattach_floor),
        "C2": float(C2),
        "C": float(C),
        "cL1": (float(cL1_x), float(cL1_y))
    }

    grid_bundle = {
        "grid_x": grid_x,
        "grid_y": grid_y,
        "grid_ux": grid_ux,
        "grid_uy": grid_uy,
        "psi": psi,
        "mask": mask,
        "x_wall_samples": x_wall_samples,
        "u_near_wall": u_near_wall,
        "p_raw": p,
        "cell_centers": cell_centers
    }

    return metrics, grid_bundle

if __name__ == "__main__":
    base_dir = "/mnt/c/Users/moham/Desktop/FOSSEE/FOSSEE_Task1"
    grid_dir = os.path.join(base_dir, "grid_independence")
    fig_dir = os.path.join(base_dir, "figures")
    os.makedirs(base_dir, exist_ok=True)
    os.makedirs(grid_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    print("=================================================================")
    print("      RE-ANALYZING & RE-GENERATING ALL GRAPHS (NEW C1 & C2)      ")
    print("=================================================================")

    main_results = {}
    main_bundles = {}
    for re_val in [60, 80, 100]:
        cname = f"Re{re_val}"
        cdir = os.path.join(base_dir, cname)
        centers, U, p, ltime = parse_openfoam_results(cdir)
        metrics, bundle = analyze_case_metrics(centers, U, p, re_val)
        main_results[cname] = metrics
        main_bundles[cname] = bundle

    grid_study_results = {}
    for re_val in [60, 80, 100]:
        grid_study_results[re_val] = {}
        grid_study_results[re_val]["medium"] = main_results[f"Re{re_val}"]
        
        for m_name in ["coarse", "fine"]:
            cname = f"Re{re_val}_{m_name}"
            cdir = os.path.join(grid_dir, cname)
            centers, U, p, ltime = parse_openfoam_results(cdir)
            metrics, bundle = analyze_case_metrics(centers, U, p, re_val)
            grid_study_results[re_val][m_name] = metrics

    print("\n======================= UPDATED RESULTS SUMMARY =======================")
    for cname, m in main_results.items():
        print(f"\n--- Case {cname} (Re = {m['Re']}) ---")
        print(f"  Vortex Separation on Ramp: ({m['x_sep']:.4f}, {m['y_sep']:.4f}) m [s = {m['s_sep']:.4f} m from Pt A]")
        print(f"  Vortex Length on Wall 1 (C1): {m['C1']:.4f} m")
        print(f"  Ramp-Floor Corner:           (x = {m['x_corner']:.4f}, y = 0.0000) m")
        print(f"  Vortex Length on Wall 2 (C2): {m['C2']:.4f} m")
        print(f"  Total Primary Vortex Size C: {m['C']:.4f} m (C = C1 + C2)")
        print(f"  Vortex Center (x_v, y_v):    ({m['vortex_center'][0]:.4f}, {m['vortex_center'][1]:.4f}) m")
        print(f"  Distance Z:                  {m['Z']:.4f} m")
        print(f"  Reattachment Length (x_r):   {m['x_r']:.4f} m (measured from Point A)")

    # ==============================================================================
    # GENERATE UPDATED PUBLICATION FIGURES
    # ==============================================================================
    print("\n--> Generating updated publication-quality figures...")

    plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
    plt.rcParams['axes.edgecolor'] = '#333333'
    plt.rcParams['axes.linewidth'] = 1.2

    # A. Streamlines and Velocity Magnitude for Re=60, 80, 100
    for re_val in [60, 80, 100]:
        cname = f"Re{re_val}"
        b = main_bundles[cname]
        m = main_results[cname]
        
        gx, gy = b['grid_x'], b['grid_y']
        gux, guy = b['grid_ux'], b['grid_uy']
        mask = b['mask']
        speed = np.sqrt(gux**2 + guy**2)
        speed[mask] = np.nan

        fig, ax = plt.subplots(figsize=(14, 4), dpi=300)
        
        cf = ax.contourf(gx, gy, speed, levels=60, cmap="turbo", alpha=0.9)
        cbar = plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.2, shrink=0.6)
        cbar.set_label(f"Velocity Magnitude |U| (m/s)  [Re = {re_val}]", fontsize=11, fontweight='bold')

        u_stream = np.nan_to_num(gux)
        v_stream = np.nan_to_num(guy)
        ax.streamplot(gx, gy, u_stream, v_stream, color='black', density=2.0, linewidth=0.7, arrowsize=0.8)

        x_pts = [0, 10, 10+dx1, 10+dx1+X, 26, 36]
        y_pts = [2, 2, 0, 0, 2, 2]
        ax.plot(x_pts, y_pts, 'k-', lw=2.5, label='Domain Boundary')
        ax.plot([0, 36], [5, 5], 'k-', lw=2.5)
        ax.plot([0, 0], [2, 5], 'k--', lw=1.5, label='Inlet')
        ax.plot([36, 36], [2, 5], 'k--', lw=1.5, label='Outlet')

        ax.plot(10.0, 2.0, 'ro', markersize=7, label='Point A (10.0, 2.0)')
        ax.plot(m['cL1'][0], m['cL1'][1], 'y^', markersize=8, label=f"Centre L1 ({m['cL1'][0]:.2f}, {m['cL1'][1]:.2f})")
        ax.plot(m['vortex_center'][0], m['vortex_center'][1], 'r*', markersize=11, label=f"Vortex Center ({m['vortex_center'][0]:.2f}, {m['vortex_center'][1]:.2f})")
        ax.plot(m['x_reattach'], 0.0, 'go', markersize=7, label=f"Reattachment ({m['x_reattach']:.2f}, 0.0)")

        ax.plot([m['cL1'][0], m['vortex_center'][0]], [m['cL1'][1], m['vortex_center'][1]], 'r--', lw=1.8, label=f"Distance Z = {m['Z']:.2f} m")

        r1_x_fill = np.linspace(10, 10+dx1, 50)
        r1_y_fill = 2.0 - (r1_x_fill - 10.0) * (2.0 / dx1)
        r2_x_fill = np.linspace(26-dx2, 26, 50)
        r2_y_fill = (r2_x_fill - (26.0 - dx2)) * (2.0 / dx2)

        ax.fill_between([0, 10], 0, 2, color='#e0e0e0', hatch='//')
        ax.fill_between(r1_x_fill, 0, r1_y_fill, color='#e0e0e0', hatch='//')
        ax.fill_between(r2_x_fill, 0, r2_y_fill, color='#e0e0e0', hatch='//')
        ax.fill_between([26, 36], 0, 2, color='#e0e0e0', hatch='//')

        ax.set_xlim(-0.5, 36.5)
        ax.set_ylim(-0.5, 5.5)
        ax.set_aspect('equal')
        ax.set_xlabel("X (m) [D = 1.0 m]", fontsize=11, fontweight='bold')
        ax.set_ylabel("Y (m)", fontsize=11, fontweight='bold')
        ax.set_title(f"Task 1 CFD Analysis: Streamlines & Velocity Field (Re = {re_val})\nReattachment Length $x_r = {m['x_r']:.2f}$ m | Distance $Z = {m['Z']:.2f}$ m | Vortex Wall Size $C = C_1+C_2 = {m['C']:.2f}$ m", fontsize=12, fontweight='bold')
        ax.legend(loc='upper right', framealpha=0.9, fontsize=8)
        
        plt.tight_layout()
        fpath = os.path.join(fig_dir, f"streamlines_{cname}.png")
        fig.savefig(fpath, dpi=300)
        plt.close(fig)
        print(f"    Saved: {fpath}")

    # B. Vortex Recirculation Region Zoom (Reference 1 & 2 exact measurements)
    for re_val in [60, 80, 100]:
        cname = f"Re{re_val}"
        b = main_bundles[cname]
        m = main_results[cname]
        
        gx, gy = b['grid_x'], b['grid_y']
        gux, guy = b['grid_ux'], b['grid_uy']
        mask = b['mask']
        speed = np.sqrt(gux**2 + guy**2)
        speed[mask] = np.nan

        fig, ax = plt.subplots(figsize=(11, 5.5), dpi=300)
        
        cf = ax.contourf(gx, gy, speed, levels=50, cmap="viridis", alpha=0.85)
        cbar = plt.colorbar(cf, ax=ax, pad=0.03, shrink=0.8)
        cbar.set_label("Velocity (m/s)", fontsize=10, fontweight='bold')

        u_stream = np.nan_to_num(gux)
        v_stream = np.nan_to_num(guy)
        ax.streamplot(gx, gy, u_stream, v_stream, color='white', density=2.5, linewidth=0.9, arrowsize=0.9)

        # Domain boundary lines
        ax.plot([8, 10, 10+dx1, 10+dx1+X, 26, 28], [2, 2, 0, 0, 2, 2], 'w-', lw=3)
        ax.plot([8, 28], [5, 5], 'w-', lw=2.5)

        # Reference Points
        ax.plot(10.0, 2.0, 'ro', markersize=9, label='Point A (10.0, 2.0)')
        ax.plot(m['cL1'][0], m['cL1'][1], 'y^', markersize=9, label=f"Centre of L1 ({m['cL1'][0]:.2f}, {m['cL1'][1]:.2f})")
        ax.plot(m['vortex_center'][0], m['vortex_center'][1], 'r*', markersize=14, label=f"Vortex Center ({m['vortex_center'][0]:.2f}, {m['vortex_center'][1]:.2f})")
        ax.plot(m['x_sep'], m['y_sep'], 'mo', markersize=8, label=f"Separation ({m['x_sep']:.2f}, {m['y_sep']:.2f})")
        ax.plot(m['x_corner'], 0.0, 'k^', markersize=8, label=f"Corner ({m['x_corner']:.2f}, 0.0)")
        ax.plot(m['x_reattach'], 0.0, 'go', markersize=9, label=f"Reattachment ({m['x_reattach']:.2f}, 0.0)")

        # Distance Z line (Ref 1)
        ax.plot([m['cL1'][0], m['vortex_center'][0]], [m['cL1'][1], m['vortex_center'][1]], 'r--', lw=2.5, label=f"Distance Z = {m['Z']:.3f} m")

        # C1 Bar on Wall 1 (from separation point to corner)
        ax.plot([m['x_sep'], m['x_corner']], [m['y_sep'], 0.0], color='cyan', lw=4.5, solid_capstyle='round', label=f"C1 = {m['C1']:.2f} m (Wall 1)")
        
        # C2 Bar on Wall 2 (from corner to reattachment)
        ax.plot([m['x_corner'], m['x_reattach']], [0.0, 0.0], color='yellow', lw=4.5, solid_capstyle='round', label=f"C2 = {m['C2']:.2f} m (Wall 2)")

        # Text Annotations
        ax.annotate(f"C1 = {m['C1']:.2f} m\n(Wall 1: Ramp Vortex Extent)", xy=((m['x_sep']+m['x_corner'])/2 - 0.4, (m['y_sep']+0)/2 + 0.3), color='cyan', fontweight='bold', fontsize=10, bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.75))
        ax.annotate(f"C2 = {m['C2']:.2f} m\n(Wall 2: Floor Extent)", xy=(m['x_corner'] + 0.5, 0.25), color='yellow', fontweight='bold', fontsize=10, bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.75))

        r1_x_fill = np.linspace(10, 10+dx1, 50)
        r1_y_fill = 2.0 - (r1_x_fill - 10.0) * (2.0 / dx1)
        ax.fill_between([8, 10], 0, 2, color='#333333', alpha=0.8)
        ax.fill_between(r1_x_fill, 0, r1_y_fill, color='#333333', alpha=0.8)

        ax.set_xlim(8.0, 20.0)
        ax.set_ylim(-0.2, 4.0)
        ax.set_aspect('equal')
        ax.set_xlabel("X (m)", fontsize=11, fontweight='bold')
        ax.set_ylabel("Y (m)", fontsize=11, fontweight='bold')
        ax.set_title(f"Primary Vortex Recirculation Bubble & Metric Annotations (Re = {re_val})\nRef 1: Distance Z = {m['Z']:.3f} m | Ref 2: Vortex Wall Size C = C1 + C2 = {m['C']:.3f} m", fontsize=11, fontweight='bold')
        ax.legend(loc='upper right', framealpha=0.9, fontsize=8)

        plt.tight_layout()
        fpath = os.path.join(fig_dir, f"vortex_detail_{cname}.png")
        fig.savefig(fpath, dpi=300)
        plt.close(fig)
        print(f"    Saved: {fpath}")

    # C. Near-Wall Velocity & Reattachment Point Comparison
    fig, ax = plt.subplots(figsize=(9, 5), dpi=300)
    colors = {60: '#1f77b4', 80: '#ff7f0e', 100: '#2ca02c'}
    for re_val in [60, 80, 100]:
        cname = f"Re{re_val}"
        b = main_bundles[cname]
        m = main_results[cname]
        ax.plot(b['x_wall_samples'], b['u_near_wall'] * 1e4, label=f"Re = {re_val} (Reattach $x_r = {m['x_r']:.2f}$ m)", color=colors[re_val], lw=2.2)
        ax.plot(m['x_reattach'], 0.0, 'o', color=colors[re_val], markersize=8)

    ax.axhline(0, color='black', linestyle=':', lw=1.5)
    ax.axvline(10.0, color='gray', linestyle='--', label='Point A (x = 10 m)')
    ax.axvline(10.0 + dx1, color='purple', linestyle='--', label='Ramp Corner (x = 13.01 m)')
    ax.set_xlim(10.0, 22.0)
    ax.set_xlabel("Axial Distance X (m)", fontsize=11, fontweight='bold')
    ax.set_ylabel("Near-Wall Velocity $u_{wall} \\times 10^4$ (m/s)", fontsize=11, fontweight='bold')
    ax.set_title("Near-Wall Velocity Profiles and Flow Reattachment Points", fontsize=12, fontweight='bold')
    ax.legend(loc='lower right', framealpha=0.9)
    ax.grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout()
    fpath = os.path.join(fig_dir, "reattachment_velocity_profiles.png")
    fig.savefig(fpath, dpi=300)
    plt.close(fig)
    print(f"    Saved: {fpath}")

    # D. Grid Independence Study Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    mesh_labels = ["Coarse\n(6,300)", "Medium\n(25,200)", "Fine\n(50,400)"]

    for re_val in [60, 80, 100]:
        xr_vals = [
            grid_study_results[re_val]["coarse"]["x_r"],
            grid_study_results[re_val]["medium"]["x_r"],
            grid_study_results[re_val]["fine"]["x_r"]
        ]
        z_vals = [
            grid_study_results[re_val]["coarse"]["Z"],
            grid_study_results[re_val]["medium"]["Z"],
            grid_study_results[re_val]["fine"]["Z"]
        ]
        ax1.plot(mesh_labels, xr_vals, 'o-', lw=2.2, color=colors[re_val], label=f"Re = {re_val}")
        ax2.plot(mesh_labels, z_vals, 's-', lw=2.2, color=colors[re_val], label=f"Re = {re_val}")

    ax1.set_ylabel("Reattachment Length $x_r$ (m)", fontsize=11, fontweight='bold')
    ax1.set_title("Reattachment Length vs Mesh Resolution", fontsize=11, fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend()

    ax2.set_ylabel("Distance $Z$ to Vortex Center (m)", fontsize=11, fontweight='bold')
    ax2.set_title("Vortex Distance $Z$ vs Mesh Resolution", fontsize=11, fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend()

    plt.suptitle("Grid Independence Study for Reynolds Numbers 60, 80, and 100", fontsize=13, fontweight='bold')
    plt.tight_layout()
    fpath = os.path.join(fig_dir, "grid_independence_study.png")
    fig.savefig(fpath, dpi=300)
    plt.close(fig)
    print(f"    Saved: {fpath}")

    # ==============================================================================
    # SAVE UPDATED SUMMARY DATA JSON
    # ==============================================================================
    summary_data = {
        "geometric_analysis": {
            "D": D,
            "L1": L1,
            "L2": L2,
            "dx1": dx1,
            "dx2": dx2,
            "X": X,
            "alpha_deg": alpha,
            "alpha_internal_deg": alpha_int,
            "centre_L1": (cL1_x, cL1_y),
            "centre_L2": (cL2_x, cL2_y)
        },
        "simulation_results": main_results,
        "grid_study_results": grid_study_results
    }
    with open(os.path.join(base_dir, "simulation_results.json"), "w") as f:
        json.dump(summary_data, f, indent=4)
    print("Master analysis and all figures updated successfully!")
