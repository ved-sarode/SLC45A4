import os
import csv
import glob

systems = ['PUT', 'SPD', 'SPM']
lipids = ['CHL1', 'POPC', 'POPE']

results = []

for sys in systems:
    for lipid in lipids:
        path = f"PyLipID/{sys}/{lipid}/Interaction_{lipid}/Dataset_{lipid}/Dataset.csv"
        
        if not os.path.exists(path):
            print(f"Missing: {path}")
            continue
            
        try:
            with open(path, 'r') as f:
                reader = csv.DictReader(f)
                rows = list(reader)
                
                # Sort by Residence Time (need to cast to float)
                # Handle empty strings or invalid numbers just in case
                valid_rows = []
                for row in rows:
                    try:
                        row['Residence Time'] = float(row['Residence Time'])
                        row['Occupancy'] = float(row['Occupancy'])
                        valid_rows.append(row)
                    except ValueError:
                        continue
                        
                # Calculate mean
                if valid_rows:
                    mean_time = sum(r['Residence Time'] for r in valid_rows) / len(valid_rows)
                    print(f"Mean Residence Time for {sys} {lipid}: {mean_time:.3f} us")

                valid_rows.sort(key=lambda x: x['Residence Time'], reverse=True)
                top_hits = valid_rows[:5]
                
                for row in top_hits:
                    results.append({
                        'System': sys,
                        'Lipid': lipid,
                        'Residue': row['Residue'],
                        'Residence_Time_us': row['Residence Time'],
                        'Occupancy_pct': row['Occupancy']
                    })
        except Exception as e:
            print(f"Error reading {path}: {e}")

# Print results
print(f"{ 'System':<10} { 'Lipid':<6} { 'Residue':<10} { 'Time(us)':<10} { 'Occ(%)':<10}")
print("-" * 50)

print("\n--- CHOLESTEROL HOTSPOTS (> 0.5 us) ---")
for r in results:
    if r['Lipid'] == 'CHL1' and r['Residence_Time_us'] > 0.5:
        print(f"{r['System']:<10} {r['Lipid']:<6} {r['Residue']:<10} {r['Residence_Time_us']:<10.3f} {r['Occupancy_pct']:<10.1f}")

print("\n--- POPC TOP HITS ---")
count = 0
for r in results:
    if r['Lipid'] == 'POPC':
        print(f"{r['System']:<10} {r['Lipid']:<6} {r['Residue']:<10} {r['Residence_Time_us']:<10.3f} {r['Occupancy_pct']:<10.1f}")
        count += 1
        if count >= 10: break 

print("\n--- POPE TOP HITS ---")
count = 0
for r in results:
    if r['Lipid'] == 'POPE':
        print(f"{r['System']:<10} {r['Lipid']:<6} {r['Residue']:<10} {r['Residence_Time_us']:<10.3f} {r['Occupancy_pct']:<10.1f}")
        count += 1
        if count >= 10: break
