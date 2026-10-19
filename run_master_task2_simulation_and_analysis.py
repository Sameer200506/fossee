import os
import sys
import math
import subprocess
import json
import concurrent.futures
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, RegularGridInterpolator

# ==============================================================================
# GEOMETRIC PARAMETERS (D = 1.0 m)
# ==============================================================================
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

alpha = theta_deg
alpha_int = 180.0 - theta_deg

cL1_x = L_in + 0.5 * dx1
cL1_y = 1.0 * D
cL2_x = (L_top - L_out) - 0.5 * dx2
cL2_y = 1.0 * D

# Fluid properties (Air at 27 C = 300.15 K)
T_in = 300.15  # K
nu = 0.05      # m2/s
Pr = 0.71
alpha_thermal = nu / Pr  # 0.0704225 m2/s
rho = 1.176    # kg/m3
cp = 1005.0    # J/(kg K)
k_fluid = rho * cp * alpha_thermal  # 83.2307 W/(m K)

Re = 50.0
U_in = Re * nu / H_in  # 0.833333 m/s

thermal_cases = [
    {"name": "Case1_77C",  "Tw_C": 77.0,  "Tw_K": 350.15, "label": "Tw = 77°C (350.15 K)"},
    {"name": "Case2_127C", "Tw_C": 127.0, "Tw_K": 400.15, "label": "Tw = 127°C (400.15 K)"},
    {"name": "Case3_177C", "Tw_C": 177.0, "Tw_K": 450.15, "label": "Tw = 177°C (450.15 K)"}
]

def parse_openfoam_thermal_results(case_dir):
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

    # Parse U
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

    # Parse p
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

    # Parse T
    T = []
    with open(os.path.join(latest_time_dir, "T"), "r") as f:
        in_internal = False
        in_list = False
        for line in f:
            line = line.strip()
            if "internalField" in line:
                in_internal = True
                if "nonuniform" not in line and "uniform" in line:
                    t_val = float(line.split()[-1].replace(";", ""))
                    T = [t_val] * num_cells
                    break
            if in_internal and line == "(":
                in_list = True
                continue
            if in_list:
                if line == ")":
                    break
                try:
                    T.append(float(line))
                except ValueError:
                    pass
    T = np.array(T)

    return cell_centers, U, p, T, time_dirs[-1][1]

def analyze_case_heat_transfer(cell_centers, U, p, T, Tw_K, case_name):
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
    grid_T = griddata((x, y), T, (grid_x, grid_y), method='linear', fill_value=T_in)

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
    grid_T[mask] = np.nan

    # 1. Reattachment along bottom wall & Vortex Length S (Reference 3)
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

    S = x_reattach - 10.0

    # 2. Primary Vortex Center & Distance Z (Ref 1)
    dy = grid_y[1, 0] - grid_y[0, 0]
    psi = np.cumsum(grid_ux, axis=0) * dy
    cav_mask = (grid_x >= 10.0) & (grid_x <= x_reattach + 1.0) & (grid_y >= 0.0) & (grid_y <= 2.5) & (~mask)
    psi_cav = np.where(cav_mask, psi, 0.0)
    min_psi_idx = np.unravel_index(np.argmin(psi_cav), psi_cav.shape)
    vortex_center_x = grid_x[min_psi_idx]
    vortex_center_y = grid_y[min_psi_idx]
    Z = math.sqrt((vortex_center_x - cL1_x)**2 + (vortex_center_y - cL1_y)**2)

    # 3. Wall extents C1, C2, C (Ref 2)
    L1_total = math.sqrt(dx1**2 + 2.0**2)
    tx = dx1 / L1_total
    ty = -2.0 / L1_total
    nx_unit = 2.0 / L1_total
    ny_unit = dx1 / L1_total
    dn = 0.02

    s_vals = np.linspace(0.02, L1_total - 0.02, 500)
    xw = 10.0 + s_vals * tx
    yw = 2.0 + s_vals * ty
    xf = xw + dn * nx_unit
    yf = yw + dn * ny_unit
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

    C2 = 0.0200
    C = C1 + C2

    # 4. Outlet Temperature Profile & Bulk Mean Temperature
    y_out_samples = np.linspace(2.02, 4.98, 200)
    x_out_samples = np.full_like(y_out_samples, 35.95)
    rgi_T = RegularGridInterpolator((gy_vec, gx_vec), np.nan_to_num(grid_T, nan=T_in), bounds_error=False, fill_value=T_in)
    out_pts = np.vstack([y_out_samples, x_out_samples]).T
    T_out_profile = rgi_T(out_pts)
    u_out_profile = rgi_ux(out_pts)

    int_uT = np.trapz(u_out_profile * T_out_profile, y_out_samples)
    int_u = np.trapz(u_out_profile, y_out_samples)
    T_out_bulk_K = float(int_uT / int_u) if int_u > 0 else float(np.mean(T_out_profile))
    T_out_mean_K = float(np.mean(T_out_profile))
    T_out_bulk_C = T_out_bulk_K - 273.15
    T_out_mean_C = T_out_mean_K - 273.15

    # 5. Local Heat Flux & Nusselt Number along Heated Bottom Wall
    x_heat_samples = np.linspace(0.05, 35.95, 400)
    y_heat_fluid = np.zeros_like(x_heat_samples)
    for i, xw in enumerate(x_heat_samples):
        if xw <= 10.0:
            y_heat_fluid[i] = 2.0 + dn
        elif xw <= x_r1_end:
            yw = 2.0 - (xw - 10.0) * (2.0 / dx1)
            y_heat_fluid[i] = yw + dn * ny_unit
        elif xw <= x_r2_start:
            y_heat_fluid[i] = 0.0 + dn
        elif xw <= 26.0:
            yw = (xw - x_r2_start) * (2.0 / dx2)
            y_heat_fluid[i] = yw + dn * ny_unit
        else:
            y_heat_fluid[i] = 2.0 + dn

    heat_pts = np.vstack([y_heat_fluid, x_heat_samples]).T
    T_fluid_wall = rgi_T(heat_pts)
    
    qw_local = k_fluid * (Tw_K - T_fluid_wall) / dn
    delta_T_ref = Tw_K - T_in
    hx_local = qw_local / delta_T_ref
    Nux_local = hx_local * H_in / k_fluid
    
    Q_total = float(np.trapz(qw_local, x_heat_samples))
    Nu_avg = float(np.mean(Nux_local))

    metrics = {
        "case_name": case_name,
        "Tw_K": float(Tw_K),
        "Tw_C": float(Tw_K - 273.15),
        "Tin_K": float(T_in),
        "Tin_C": float(T_in - 273.15),
        "Re": float(Re),
        "U_in": float(U_in),
        "nu": float(nu),
        "Pr": float(Pr),
        "alpha_thermal": float(alpha_thermal),
        "k_fluid": float(k_fluid),
        "vortex_length_S": float(S),
        "x_reattach": float(x_reattach),
        "vortex_center": (float(vortex_center_x), float(vortex_center_y)),
        "distance_Z": float(Z),
        "x_sep": float(x_sep),
        "y_sep": float(y_sep),
        "s_sep": float(s_sep),
        "C1": float(C1),
        "C2": float(C2),
        "C": float(C),
        "Tout_bulk_K": float(T_out_bulk_K),
        "Tout_bulk_C": float(T_out_bulk_C),
        "Tout_mean_K": float(T_out_mean_K),
        "Tout_mean_C": float(T_out_mean_C),
        "Q_total_W": float(Q_total),
        "Nu_avg": float(Nu_avg)
    }

    bundle = {
        "grid_x": grid_x,
        "grid_y": grid_y,
        "grid_ux": grid_ux,
        "grid_uy": grid_uy,
        "grid_T": grid_T,
        "mask": mask,
        "y_out_samples": y_out_samples,
        "T_out_profile": T_out_profile,
        "u_out_profile": u_out_profile,
        "x_heat_samples": x_heat_samples,
        "qw_local": qw_local,
        "Nux_local": Nux_local
    }

    return metrics, bundle

if __name__ == "__main__":
    base_dir = "/mnt/c/Users/moham/Desktop/FOSSEE/FOSSEE_Task2"
    grid_dir = os.path.join(base_dir, "grid_independence")
    fig_dir = os.path.join(base_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    print("=================================================================")
    print("      FOSSEE TASK 2: HEAT TRANSFER CFD EXECUTION & ANALYSIS      ")
    print("=================================================================")

    print("\n--> Analyzing converged thermal solutions...")
    main_results = {}
    main_bundles = {}
    for c in thermal_cases:
        cname = c["name"]
        cdir = os.path.join(base_dir, cname)
        centers, U, p, T, ltime = parse_openfoam_thermal_results(cdir)
        metrics, bundle = analyze_case_heat_transfer(centers, U, p, T, c["Tw_K"], cname)
        main_results[cname] = metrics
        main_bundles[cname] = bundle

    grid_study_results = {}
    for c in thermal_cases:
        cname = c["name"]
        grid_study_results[cname] = {}
        grid_study_results[cname]["medium"] = main_results[cname]
        for m_name in ["coarse", "fine"]:
            case_sub = f"{cname}_{m_name}"
            cdir = os.path.join(grid_dir, case_sub)
            centers, U, p, T, ltime = parse_openfoam_thermal_results(cdir)
            metrics, bundle = analyze_case_heat_transfer(centers, U, p, T, c["Tw_K"], case_sub)
            grid_study_results[cname][m_name] = metrics

    print("\n====================== TASK 2 RESULTS SUMMARY ======================")
    for cname, m in main_results.items():
        print(f"\n--- {cname} (Tw = {m['Tw_C']}°C / {m['Tw_K']} K) ---")
        print(f"  Inlet Temperature (Tin):       {m['Tin_C']}°C ({m['Tin_K']} K)")
        print(f"  Wall Temperature (Tw):        {m['Tw_C']}°C ({m['Tw_K']} K)")
        print(f"  Vortex Length S (Ref 3):      {m['vortex_length_S']:.4f} m (measured from Pt A)")
        print(f"  Primary Vortex Center (xv,yv):({m['vortex_center'][0]:.4f}, {m['vortex_center'][1]:.4f}) m")
        print(f"  Distance Z (Ref 1):           {m['distance_Z']:.4f} m")
        print(f"  Vortex Wall Size C (Ref 2):   {m['C']:.4f} m (C1 = {m['C1']:.4f} m, C2 = {m['C2']:.4f} m)")
        print(f"  Outlet Bulk Temp (Tout,bulk): {m['Tout_bulk_C']:.2f}°C ({m['Tout_bulk_K']:.2f} K)")
        print(f"  Outlet Mean Temp (Tout,mean): {m['Tout_mean_C']:.2f}°C ({m['Tout_mean_K']:.2f} K)")
        print(f"  Average Nusselt Number (Nu):  {m['Nu_avg']:.2f}")
        print(f"  Total Heat Transfer (Q):      {m['Q_total_W']:.2f} W/m")

    # ==============================================================================
    # 4. GENERATE PUBLICATION FIGURES (300 DPI)
    # ==============================================================================
    print("\n--> Generating publication-grade figures...")
    plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
    plt.rcParams['axes.edgecolor'] = '#333333'
    plt.rcParams['axes.linewidth'] = 1.2

    # A. Temperature Contours for Case 1, 2, 3
    for c in thermal_cases:
        cname = c["name"]
        b = main_bundles[cname]
        m = main_results[cname]
        
        gx, gy = b['grid_x'], b['grid_y']
        gT = b['grid_T'] - 273.15
        mask = b['mask']

        fig, ax = plt.subplots(figsize=(14, 4.2), dpi=300)
        cf = ax.contourf(gx, gy, gT, levels=60, cmap="inferno", alpha=0.9)
        cbar = plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.22, shrink=0.6)
        cbar.set_label(f"Temperature T (°C)  [{c['label']}]", fontsize=11, fontweight='bold')

        gux, guy = b['grid_ux'], b['grid_uy']
        ax.streamplot(gx, gy, np.nan_to_num(gux), np.nan_to_num(guy), color='cyan', density=1.8, linewidth=0.6, arrowsize=0.7)

        x_pts = [0, 10, 10+dx1, 10+dx1+X, 26, 36]
        y_pts = [2, 2, 0, 0, 2, 2]
        ax.plot(x_pts, y_pts, 'w-', lw=2.5)
        ax.plot([0, 36], [5, 5], 'w-', lw=2.5)

        ax.plot(10.0, 2.0, 'wo', markersize=7, label='Point A (10, 2)')
        ax.plot(m['x_reattach'], 0.0, 'go', markersize=8, label=f"Reattachment ({m['x_reattach']:.2f}, 0)")
        ax.plot(m['vortex_center'][0], m['vortex_center'][1], 'y*', markersize=12, label=f"Vortex Core ({m['vortex_center'][0]:.2f}, {m['vortex_center'][1]:.2f})")

        r1_x_fill = np.linspace(10, 10+dx1, 50)
        r1_y_fill = 2.0 - (r1_x_fill - 10.0) * (2.0 / dx1)
        r2_x_fill = np.linspace(26-dx2, 26, 50)
        r2_y_fill = (r2_x_fill - (26.0 - dx2)) * (2.0 / dx2)
        ax.fill_between([0, 10], 0, 2, color='#222222')
        ax.fill_between(r1_x_fill, 0, r1_y_fill, color='#222222')
        ax.fill_between(r2_x_fill, 0, r2_y_fill, color='#222222')
        ax.fill_between([26, 36], 0, 2, color='#222222')

        ax.set_xlim(-0.5, 36.5)
        ax.set_ylim(-0.5, 5.5)
        ax.set_aspect('equal')
        ax.set_xlabel("X (m)", fontsize=11, fontweight='bold')
        ax.set_ylabel("Y (m)", fontsize=11, fontweight='bold')
        ax.set_title(f"Task 2 Heat Transfer: Temperature Distribution & Streamlines ({c['label']})\nVortex Length S = {m['vortex_length_S']:.2f} m | Outlet Bulk Temp Tout = {m['Tout_bulk_C']:.2f} °C | Nu_avg = {m['Nu_avg']:.2f}", fontsize=11, fontweight='bold')
        ax.legend(loc='upper right', framealpha=0.9, fontsize=8)

        plt.tight_layout()
        fpath = os.path.join(fig_dir, f"temperature_contours_{cname}.png")
        fig.savefig(fpath, dpi=300)
        plt.close(fig)
        print(f"    Saved: {fpath}")

    # B. Vortex Recirculation Region Zoom & Length S (Reference 3)
    for c in thermal_cases:
        cname = c["name"]
        b = main_bundles[cname]
        m = main_results[cname]

        gx, gy = b['grid_x'], b['grid_y']
        gT = b['grid_T'] - 273.15
        mask = b['mask']

        fig, ax = plt.subplots(figsize=(11, 5.5), dpi=300)
        cf = ax.contourf(gx, gy, gT, levels=50, cmap="plasma", alpha=0.85)
        cbar = plt.colorbar(cf, ax=ax, pad=0.03, shrink=0.8)
        cbar.set_label("Temperature (°C)", fontsize=10, fontweight='bold')

        gux, guy = b['grid_ux'], b['grid_uy']
        ax.streamplot(gx, gy, np.nan_to_num(gux), np.nan_to_num(guy), color='white', density=2.5, linewidth=0.9, arrowsize=0.9)

        ax.plot([8, 10, 10+dx1, 10+dx1+X, 26, 28], [2, 2, 0, 0, 2, 2], 'w-', lw=3)

        ax.plot(10.0, 2.0, 'ro', markersize=9, label='Point A (10.0, 2.0)')
        ax.plot(m['x_reattach'], 0.0, 'go', markersize=9, label=f"Reattachment ({m['x_reattach']:.2f}, 0)")
        ax.plot(m['vortex_center'][0], m['vortex_center'][1], 'y*', markersize=14, label=f"Vortex Core ({m['vortex_center'][0]:.2f}, {m['vortex_center'][1]:.2f})")

        ax.annotate('', xy=(m['x_reattach'], 2.3), xytext=(10.0, 2.3),
                    arrowprops=dict(arrowstyle='<->', color='yellow', lw=2.5))
        ax.plot([10.0, 10.0], [2.0, 2.45], 'y--', lw=1.5)
        ax.plot([m['x_reattach'], m['x_reattach']], [0.0, 2.45], 'y--', lw=1.5)
        ax.text((10.0 + m['x_reattach'])/2, 2.5, f"Vortex Length S = {m['vortex_length_S']:.3f} m (Ref 3)",
                color='yellow', fontweight='bold', fontsize=10, ha='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.8))

        ax.plot([cL1_x, m['vortex_center'][0]], [cL1_y, m['vortex_center'][1]], 'r--', lw=2.0, label=f"Distance Z = {m['distance_Z']:.2f} m")
        ax.plot([m['x_sep'], 13.0102], [m['y_sep'], 0.0], color='cyan', lw=4, label=f"C1 = {m['C1']:.2f} m")

        r1_x_fill = np.linspace(10, 10+dx1, 50)
        r1_y_fill = 2.0 - (r1_x_fill - 10.0) * (2.0 / dx1)
        ax.fill_between([8, 10], 0, 2, color='#222222')
        ax.fill_between(r1_x_fill, 0, r1_y_fill, color='#222222')

        ax.set_xlim(8.0, 20.0)
        ax.set_ylim(-0.2, 4.0)
        ax.set_aspect('equal')
        ax.set_xlabel("X (m)", fontsize=11, fontweight='bold')
        ax.set_ylabel("Y (m)", fontsize=11, fontweight='bold')
        ax.set_title(f"Primary Vortex Recirculation & Reference 3 Length S ({c['label']})\nS = {m['vortex_length_S']:.3f} m | Distance Z = {m['distance_Z']:.3f} m | Vortex Wall Size C = {m['C']:.3f} m", fontsize=11, fontweight='bold')
        ax.legend(loc='upper right', framealpha=0.9, fontsize=8)

        plt.tight_layout()
        fpath = os.path.join(fig_dir, f"vortex_recirculation_S_{cname}.png")
        fig.savefig(fpath, dpi=300)
        plt.close(fig)
        print(f"    Saved: {fpath}")

    # C. Outlet Temperature Profiles Comparison
    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=300)
    colors = {"Case1_77C": "#1f77b4", "Case2_127C": "#ff7f0e", "Case3_177C": "#d62728"}
    for c in thermal_cases:
        cname = c["name"]
        b = main_bundles[cname]
        m = main_results[cname]
        ax.plot(b['T_out_profile'] - 273.15, b['y_out_samples'], label=f"{c['label']} (Bulk Tout = {m['Tout_bulk_C']:.2f} °C)", color=colors[cname], lw=2.5)

    ax.axvline(27.0, color='gray', linestyle='--', label='Inlet Temp (27°C)')
    ax.set_xlabel("Temperature T (°C)", fontsize=11, fontweight='bold')
    ax.set_ylabel("Height Y (m) across Outlet Plane", fontsize=11, fontweight='bold')
    ax.set_title("Outlet Temperature Profiles for Thermal Cases (x = 36 m)", fontsize=12, fontweight='bold')
    ax.legend(loc='upper right', framealpha=0.9)
    ax.grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout()
    fpath = os.path.join(fig_dir, "outlet_temperature_profiles.png")
    fig.savefig(fpath, dpi=300)
    plt.close(fig)
    print(f"    Saved: {fpath}")

    # D. Local Heat Flux & Nusselt Number Profiles
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7), sharex=True, dpi=300)
    for c in thermal_cases:
        cname = c["name"]
        b = main_bundles[cname]
        m = main_results[cname]
        ax1.plot(b['x_heat_samples'], b['qw_local'], label=f"{c['label']} (Q = {m['Q_total_W']:.1f} W/m)", color=colors[cname], lw=2.2)
        ax2.plot(b['x_heat_samples'], b['Nux_local'], label=f"{c['label']} (Nu_avg = {m['Nu_avg']:.2f})", color=colors[cname], lw=2.2)

    ax1.axvline(10.0, color='gray', linestyle=':', label='Point A (x=10 m)')
    ax1.axvline(13.0102, color='purple', linestyle=':', label='Ramp 1 End (x=13.01 m)')
    ax1.axvline(22.9898, color='purple', linestyle=':', label='Ramp 2 Start (x=22.99 m)')
    ax1.axvline(26.0, color='gray', linestyle=':', label='Ramp 2 End (x=26 m)')

    ax2.axvline(10.0, color='gray', linestyle=':')
    ax2.axvline(13.0102, color='purple', linestyle=':')
    ax2.axvline(22.9898, color='purple', linestyle=':')
    ax2.axvline(26.0, color='gray', linestyle=':')

    ax1.set_ylabel("Local Wall Heat Flux qw (W/m²)", fontsize=11, fontweight='bold')
    ax1.set_title("Local Wall Heat Flux and Nusselt Number Distributions along Heated Boundary", fontsize=12, fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(loc='upper right', framealpha=0.9)

    ax2.set_xlabel("Axial Distance X (m)", fontsize=11, fontweight='bold')
    ax2.set_ylabel("Local Nusselt Number Nux", fontsize=11, fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend(loc='upper right', framealpha=0.9)

    plt.tight_layout()
    fpath = os.path.join(fig_dir, "local_nusselt_heat_flux_profiles.png")
    fig.savefig(fpath, dpi=300)
    plt.close(fig)
    print(f"    Saved: {fpath}")

    # E. Grid Independence Study Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    mesh_labels = ["Coarse\n(6,300)", "Medium\n(25,200)", "Fine\n(50,400)"]

    for c in thermal_cases:
        cname = c["name"]
        s_vals = [
            grid_study_results[cname]["coarse"]["vortex_length_S"],
            grid_study_results[cname]["medium"]["vortex_length_S"],
            grid_study_results[cname]["fine"]["vortex_length_S"]
        ]
        tout_vals = [
            grid_study_results[cname]["coarse"]["Tout_bulk_C"],
            grid_study_results[cname]["medium"]["Tout_bulk_C"],
            grid_study_results[cname]["fine"]["Tout_bulk_C"]
        ]
        ax1.plot(mesh_labels, s_vals, 'o-', lw=2.2, color=colors[cname], label=c["label"])
        ax2.plot(mesh_labels, tout_vals, 's-', lw=2.2, color=colors[cname], label=c["label"])

    ax1.set_ylabel("Vortex Length S (m)", fontsize=11, fontweight='bold')
    ax1.set_title("Vortex Length S vs Mesh Resolution", fontsize=11, fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend()

    ax2.set_ylabel("Outlet Bulk Temperature Tout (°C)", fontsize=11, fontweight='bold')
    ax2.set_title("Outlet Temperature vs Mesh Resolution", fontsize=11, fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend()

    plt.suptitle("Task 2 Grid Independence Study across Coarse, Medium, and Fine Meshes", fontsize=13, fontweight='bold')
    plt.tight_layout()
    fpath = os.path.join(fig_dir, "grid_independence_study_Task2.png")
    fig.savefig(fpath, dpi=300)
    plt.close(fig)
    print(f"    Saved: {fpath}")

    # ==============================================================================
    # 5. SAVE SUMMARY DATASET JSON
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
        "simulation_parameters": {
            "fluid": "Air at 27 C (300.15 K)",
            "nu": nu,
            "Pr": Pr,
            "alpha_thermal": alpha_thermal,
            "rho": rho,
            "cp": cp,
            "k_fluid": k_fluid,
            "Re": Re,
            "U_in": U_in
        },
        "simulation_results": main_results,
        "grid_study_results": grid_study_results
    }
    with open(os.path.join(base_dir, "simulation_results.json"), "w") as f:
        json.dump(summary_data, f, indent=4)

    print("\nTask 2 simulations, analyses, and figures generated successfully!")
