import numpy as np
import os
import glob

systems = ["9GHZ", "9GUI", "PUT_bound", "SPD_bound", "SPM_bound"]
results = {}

print(f"{ 'System':<12} | {'PC1 Mean':<10} | {'PC2 Mean':<10} | {'Spread (Area)':<10}")
print("-" * 50)

for sys in systems:
    fpath = f"systems/{sys}/proj2d.xvg"
    if not os.path.exists(fpath):
        continue
        
    # Read XVG (skip @ and # lines)
    data = []
    with open(fpath, 'r') as f:
        for line in f:
            if not line.startswith(("@", "#")):
                parts = line.split()
                if len(parts) >= 2:
                    # XVG usually: Time, PC1, PC2
                    # Assuming cols 1 and 2 (0-indexed) are PC1, PC2
                    data.append([float(parts[0]), float(parts[1])])
                    
    data = np.array(data)
    
    # Calculate Stats
    pc1 = data[:, 0]
    pc2 = data[:, 1]
    
    centroid = (np.mean(pc1), np.mean(pc2))
    std_dev = (np.std(pc1), np.std(pc2))
    
    # "Spread" as Area of ellipse = pi * std_pc1 * std_pc2
    spread_area = np.pi * std_dev[0] * std_dev[1]
    
    results[sys] = {"centroid": centroid, "spread": spread_area}
    
    print(f"{sys:<12} | {centroid[0]:<10.2f} | {centroid[1]:<10.2f} | {spread_area:<10.2f}")

print("-" * 50)

# Separation Analysis
apo_center = results["9GHZ"]["centroid"]
print("\nDistance from APO (9GHZ) Centroid:")
for sys in systems:
    if sys == "9GHZ": continue
    c = results[sys]["centroid"]
    dist = np.sqrt((c[0] - apo_center[0])**2 + (c[1] - apo_center[1])**2)
    print(f"{sys}: {dist:.2f} nm")
