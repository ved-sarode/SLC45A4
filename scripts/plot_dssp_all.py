import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import os
import sys
import MDAnalysis as mda
from matplotlib import gridspec

# Configuration
SYSTEMS = [
    {"label": "9GHZ", "path": "secondary_structure/9ghz/ss.dat"},
    {"label": "9GUI", "path": "secondary_structure/9gui/ss.dat"},
    {"label": "PUT_bound", "path": "secondary_structure/put/ss.dat"},
    {"label": "SPD_bound", "path": "secondary_structure/spd/ss.dat"},
    {"label": "SPM_bound", "path": "secondary_structure/spm/ss.dat"}
]
PDB_PATH = "systems/9GHZ/calpha_traj.pdb"
OUTPUT_DIR = "plots/dssp_analysis"

SS_MAP = {
    'H': 'H', 'G': 'H', 'I': 'H',
    'E': 'E', 'B': 'E',
    'T': 'C', 'S': 'C', ' ': 'C', '~': 'C', 'P': 'C'
}

RES_BREAK_START = 259
RES_BREAK_END = 415

def get_residue_indices(pdb_file):
    u = mda.Universe(pdb_file)
    return u.residues.resids

def parse_dssp_file(filepath):
    if not os.path.exists(filepath):
        print(f"Error: File not found {filepath}")
        return None
    matrix = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            row = []
            for char in line:
                ss_type = SS_MAP.get(char, 'C')
                if ss_type == 'H': val = 2
                elif ss_type == 'E': val = 1
                else: val = 0
                row.append(val)
            matrix.append(row)
    return np.array(matrix)

def analyze_occupancy(matrix):
    n_frames, n_res = matrix.shape
    counts_c = np.sum(matrix == 0, axis=0)
    counts_e = np.sum(matrix == 1, axis=0)
    counts_h = np.sum(matrix == 2, axis=0)
    total = n_frames
    pct_c = (counts_c / total) * 100
    pct_e = (counts_e / total) * 100
    pct_h = (counts_h / total) * 100
    return pct_h, pct_e, pct_c

def split_data_by_break(data_array, residues):
    mask_left = residues <= RES_BREAK_START
    mask_right = residues >= RES_BREAK_END
    return (residues[mask_left], data_array[mask_left]), (residues[mask_right], data_array[mask_right])

def plot_all_systems_stacked(systems_data, residues, output_file):
    n_systems = len(systems_data)
    
    # Determine the minimum length across ALL datasets + residues
    # This handles the mismatch (511 vs 522 etc.)
    min_len = len(residues)
    for sys_data in systems_data:
        h, _, _ = sys_data['data']
        min_len = min(min_len, len(h))
    
    print(f"Trimming all data to {min_len} residues.")
    residues = residues[:min_len]
    
    # Width Ratios based on residues
    len_left = np.sum(residues <= RES_BREAK_START)
    len_right = np.sum(residues >= RES_BREAK_END)
    if len_right == 0: width_ratios = [len_left, 1]
    else: width_ratios = [len_left, len_right]
    
    # Figure Setup: n_systems rows, 2 columns each
    fig = plt.figure(figsize=(14, 3 * n_systems))
    gs = gridspec.GridSpec(n_systems, 2, width_ratios=width_ratios, wspace=0.05, hspace=0.3)
    
    colors = {'H': '#0072B2', 'E': '#D55E00', 'C': '#999999'}
    
    def plot_panel(ax, res, h, e, c, label, hide_y=False):
        if len(res) == 0: return
        
        # Stackplot is much more efficient for vectors than thousands of bars
        ax.stackplot(res, h, e, c, labels=['Helix', 'Sheet', 'Coil'], 
                     colors=[colors['H'], colors['E'], colors['C']], step='mid')
        
        ax.set_ylim(0, 100)
        ax.set_xlim(res[0], res[-1])
        
        if hide_y:
            ax.set_yticks([])
            ax.spines['left'].set_visible(False)
        else:
            ax.set_ylabel(f"{label}\nStructure %", fontweight='bold', fontsize=12)
            ax.spines['right'].set_visible(False)

    for i, sys_data in enumerate(systems_data):
        label = sys_data['label']
        h, e, c = sys_data['data']
        
        # Trim this system's data
        h, e, c = h[:min_len], e[:min_len], c[:min_len]
        
        # Split
        (res_l, h_l), (res_r, h_r) = split_data_by_break(h, residues)
        (_, e_l), (_, e_r)         = split_data_by_break(e, residues)
        (_, c_l), (_, c_r)         = split_data_by_break(c, residues)
        
        # Left Panel
        ax_l = fig.add_subplot(gs[i, 0])
        plot_panel(ax_l, res_l, h_l, e_l, c_l, label, hide_y=False)
        
        # Enforce uniform ticks (step 25) for Left Panel
        if len(res_l) > 0:
            start_l = (res_l[0] // 25) * 25
            ticks_l = np.arange(start_l, res_l[-1] + 1, 25)
            ticks_l = ticks_l[ticks_l >= res_l[0]] # Filter out of bounds
            ax_l.set_xticks(ticks_l)

        # Right Panel
        ax_r = fig.add_subplot(gs[i, 1])
        plot_panel(ax_r, res_r, h_r, e_r, c_r, label, hide_y=True)
        
        # Enforce uniform ticks (step 25) for Right Panel
        if len(res_r) > 0:
            start_r = (res_r[0] // 25) * 25
            if start_r < res_r[0]: start_r += 25 # Ensure start is within or after first residue
            ticks_r = np.arange(start_r, res_r[-1] + 1, 25)
            ax_r.set_xticks(ticks_r)
        
        # Legend only on top plot
        if i == 0:
            ax_l.legend(loc='upper left', bbox_to_anchor=(0, 1.25), ncol=3, frameon=False)
            
        # X Label only on bottom plot
        if i == n_systems - 1:
            ax_l.set_xlabel("Residue Index", fontweight='bold', fontsize=14)
            if len(res_r) > 0:
                ax_r.set_xlabel("Residue Index", fontweight='bold', fontsize=14)
    
    # Save PDF
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    print(f"Saved: {output_file}")
    
    # Save SVG as backup
    svg_file = output_file.replace(".pdf", ".svg")
    plt.savefig(svg_file, bbox_inches='tight')
    print(f"Saved: {svg_file}")
    
    plt.close()

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    
    residues = get_residue_indices(PDB_PATH)
    
    processed_data = []
    
    for sys_info in SYSTEMS:
        print(f"Processing {sys_info['label']}...")
        mat = parse_dssp_file(sys_info['path'])
        if mat is not None:
            h, e, c = analyze_occupancy(mat)
            processed_data.append({
                'label': sys_info['label'],
                'data': (h, e, c)
            })
    
    out_file = os.path.join(OUTPUT_DIR, "dssp_stacked_comparison_all.pdf")
    plot_all_systems_stacked(processed_data, residues, out_file)

if __name__ == "__main__":
    main()