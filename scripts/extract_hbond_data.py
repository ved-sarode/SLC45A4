import csv
import os
import re

SYSTEMS = {
    'PUT-bound': 'PUT_bound',
    'SPD-bound': 'SPD_bound',
    'SPM-bound': 'SPM_bound'
}

def clean_residue_name(res_str):
    match = re.match(r'([A-Z]+)(\d+)', res_str)
    if match:
        code = match.group(1)
        num = match.group(2)
        return f"{code.title()}{num}"
    return res_str.split('-')[0]

def get_ligand_atom(lig_str):
    return lig_str.split('-')[-1]

def is_ligand(res_string):
    # Check for ligand codes
    return any(x in res_string for x in ["UNL", "SPD", "SPM", "PUT"])

output_file = "hbond_interaction_summary.csv"

print(f"Generating {output_file}...")

with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["System", "Residue ID", "Ligand Atom", "Occupancy (%)"])

    for sys_label, folder_name in SYSTEMS.items():
        dat_path = f"systems/{folder_name}/hbonds-details.dat"
        
        if not os.path.exists(dat_path):
            continue
            
        rows = []
        with open(dat_path, 'r') as f:
            lines = f.readlines()
            # Skip header if it exists (starts with Found or donor)
            data_lines = [l for l in lines if not l.startswith("Found") and not l.startswith("donor")]
            
            for line in data_lines:
                parts = line.strip().split()
                if len(parts) < 3: continue
                
                donor = parts[0]
                acceptor = parts[1]
                occupancy_str = parts[2]
                
                try:
                    occupancy = float(occupancy_str.replace('%', ''))
                except ValueError:
                    continue
                
                # Identify ligand and protein
                if is_ligand(donor):
                    lig_atom = get_ligand_atom(donor)
                    res_raw = acceptor
                elif is_ligand(acceptor):
                    lig_atom = get_ligand_atom(acceptor)
                    res_raw = donor
                else:
                    continue # Both protein?
                
                residue_id = clean_residue_name(res_raw)
                rows.append((residue_id, lig_atom, occupancy))
            
        # Sort by Occupancy desc
        rows.sort(key=lambda x: x[2], reverse=True)
        
        for r in rows:
            writer.writerow([sys_label, r[0], r[1], f"{r[2]:.2f}"])

print(f"Finished. Check {output_file}")