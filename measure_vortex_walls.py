import os
import math
import numpy as np
from scipy.interpolate import griddata

base_dir = "/mnt/c/Users/moham/Desktop/FOSSEE/FOSSEE_Task1"

def extract_wall_vortex_extents_fast(case_dir, Re_val):
    from run_master_simulation_and_analysis import parse_openfoam_results
    centers, U, p, ltime = parse_openfoam_results(case_dir)
    
    x = centers[:, 0]
    y = centers[:, 1]
    ux = U[:, 0]
    uy = U[:, 1]
    
    # Ramp geometry
    dx1 = 2.0 / math.tan(math.radians(33.6)) # 3.01024 m
    L1_total = math.sqrt(dx1**2 + 2.0**2) # 3.61408 m
    x_corner = 10.0 + dx1 # 13.01024 m
    
    tx = dx1 / L1_total
    ty = -2.0 / L1_total
    nx = 2.0 / L1_total
    ny = dx1 / L1_total
    
    dn = 0.02
    
    # 1. Sample along Wall 1 (Ramp L1)
    s_vals = np.linspace(0.02, L1_total - 0.02, 500)
    xw = 10.0 + s_vals * tx
    yw = 2.0 + s_vals * ty
    xf = xw + dn * nx
    yf = yw + dn * ny
    
    # Single vectorized griddata call
    ufx = griddata((x, y), ux, (xf, yf), method='linear', fill_value=0.0)
    ufy = griddata((x, y), uy, (xf, yf), method='linear', fill_value=0.0)
    ut_wall1 = ufx * tx + ufy * ty
    
    # Separation on Wall 1: where ut changes from >0 to <0
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

    # Separation coordinates in Cartesian (x, y)
    x_sep = 10.0 + s_sep * tx
    y_sep = 2.0 + s_sep * ty
    
    # 2. Sample along Wall 2 (Cavity Floor y = 0)
    x_floor_vals = np.linspace(x_corner + 0.02, 25.0, 500)
    yf_floor = np.full_like(x_floor_vals, dn)
    
    u_floor = griddata((x, y), ux, (x_floor_vals, yf_floor), method='linear', fill_value=0.0)
    
    reattach_idx = None
    for i in range(len(u_floor) - 1):
        if u_floor[i] <= 0 and u_floor[i+1] > 0:
            reattach_idx = i
            break
            
    if reattach_idx is not None:
        x0, x1_pt = x_floor_vals[reattach_idx], x_floor_vals[reattach_idx+1]
        u0, u1 = u_floor[reattach_idx], u_floor[reattach_idx+1]
        x_reattach = x0 - u0 * (x1_pt - x0) / (u1 - u0)
        C2 = x_reattach - x_corner
    else:
        x_reattach = x_floor_vals[np.argmin(np.abs(u_floor))]
        C2 = x_reattach - x_corner
        
    C = C1 + C2
    
    return {
        "Re": Re_val,
        "s_separation": s_sep,
        "x_sep": x_sep,
        "y_sep": y_sep,
        "L1_total": L1_total,
        "C1": C1,
        "x_corner": x_corner,
        "x_reattach": x_reattach,
        "C2": C2,
        "C": C
    }

if __name__ == "__main__":
    for cname, re_val in [("Re60", 60), ("Re80", 80), ("Re100", 100)]:
        cdir = os.path.join(base_dir, cname)
        res = extract_wall_vortex_extents_fast(cdir, re_val)
        print(f"\n================ Case {cname} (Re = {re_val}) ================")
        print(f"  Total Ramp Length (L1):             {res['L1_total']:.4f} m")
        print(f"  Vortex Separation Point on Ramp:    (x = {res['x_sep']:.4f}, y = {res['y_sep']:.4f}) m (s = {res['s_separation']:.4f} m from Point A)")
        print(f"  Vortex Wall Extent C1 (on Wall 1):  {res['C1']:.4f} m")
        print(f"  Ramp-Floor Corner Vertex:           (x = {res['x_corner']:.4f}, y = 0.0000) m")
        print(f"  Vortex Reattachment Point on W2:    (x = {res['x_reattach']:.4f}, y = 0.0000) m")
        print(f"  Vortex Wall Extent C2 (on Wall 2):  {res['C2']:.4f} m")
        print(f"  Total Primary Vortex Size (C=C1+C2):{res['C']:.4f} m")
