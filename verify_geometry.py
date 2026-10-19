import os
import sys
import math
import subprocess

# Geometric and physical parameters
D = 1.0 # reference dimension (between 0.1 and 10)
H_in = 3.0 * D
H_cav = 5.0 * D
L_in = 10.0 * D
L_out = 10.0 * D
L_top = 36.0 * D
theta_deg = 33.6
theta_rad = math.radians(theta_deg)
step_h = H_cav - H_in # 2.0 * D

# Geometric analysis values
L1 = step_h / math.sin(theta_rad)
L2 = L1
dx1 = step_h / math.tan(theta_rad)
dx2 = dx1
X = L_top - (L_in + dx1 + dx2 + L_out)
alpha = theta_deg
alpha_int = 180.0 - theta_deg

# Coordinates of centres (Origin at domain bottom-left: inlet floor y=2.0, cavity floor y=0.0, top wall y=5.0)
# Ramp 1 extends from (L_in, 2.0) to (L_in + dx1, 0.0)
c_L1_x = L_in + 0.5 * dx1
c_L1_y = 0.5 * (2.0 * D + 0.0)

# Ramp 2 extends from (L_in + dx1 + X, 0.0) to (L_in + dx1 + X + dx2, 2.0) = (26.0, 2.0)
c_L2_x = (L_top - L_out) - 0.5 * dx2
c_L2_y = 0.5 * (2.0 * D + 0.0)

print(f"Geometric Analysis:")
print(f"  Reference D: {D} m")
print(f"  L1 = L2 = {L1:.4f} m ({L1/D:.4f} D)")
print(f"  dx1 = dx2 = {dx1:.4f} m ({dx1/D:.4f} D)")
print(f"  X = {X:.4f} m ({X/D:.4f} D)")
print(f"  Angle alpha = {alpha:.1f} deg (Internal angle = {alpha_int:.1f} deg)")
print(f"  Centre of L1: ({c_L1_x:.4f}, {c_L1_y:.4f}) m")
print(f"  Centre of L2: ({c_L2_x:.4f}, {c_L2_y:.4f}) m")

# Fluid properties: Air at 30 C, STP
# T = 303.15 K, P = 101325 Pa
# rho = P / (R * T) = 101325 / (287.058 * 303.15) = 1.1644 kg/m3
# mu = 1.865e-5 Pa.s -> nu = mu / rho = 1.602e-5 m2/s
nu = 1.602e-05

# Reynolds numbers and inlet velocities (based on inlet height H_in = 3.0 m)
# Re = U_in * H_in / nu  =>  U_in = Re * nu / H_in
reynolds_cases = {
    "Re60": {"Re": 60, "U": 60 * nu / H_in},
    "Re80": {"Re": 80, "U": 80 * nu / H_in},
    "Re100": {"Re": 100, "U": 100 * nu / H_in}
}

for name, data in reynolds_cases.items():
    print(f"  {name}: Re = {data['Re']}, U_in = {data['U']:.6e} m/s")
