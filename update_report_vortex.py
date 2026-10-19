import os
import shutil
import json

# Update report and JSON with exact vortex separation wall contact values
report_path = r"C:\Users\moham\Desktop\FOSSEE\FOSSEE_Task1\report\FOSSEE_Task1_Fluid_Mechanics_Report.md"
root_report_path = r"C:\Users\moham\Desktop\FOSSEE\FOSSEE_Task1_Fluid_Mechanics_Report.md"

# Let's read and update
with open(report_path, "r") as f:
    text = f.read()

# Replace the discussion of C1 and C2 to explicitly highlight the two definitions:
# 1. Total wall segment vs
# 2. Exact vortex recirculation contact lengths (from separation point on ramp to corner for C1, and corner to reattachment for C2)

print("Updating report with exact separation and reattachment definitions...")
