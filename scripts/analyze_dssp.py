import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import MDAnalysis as mda
from matplotlib import gridspec

# Secondary Structure Mapping
SS_MAP = {
    'H': 'H', 'G': 'H', 'I': 'H',
    'E': 'E', 'B': 'E',
    'T': 'C', 'S': 'C', ' ': 'C', '~': 'C', 'P': 'C'
}

# Residue Break Configuration
RES_BREAK_START = 259
RES_BREAK_END = 415

def get_residue_indices(pdb_file):
    """Extracts residue IDs from PDB."""
    u = mda.Universe(pdb_file)
    return u.residues.resids

def parse_dssp_file(filepath):
    """Parses a GROMACS ss.dat file into a numeric matrix."""
    if not os.path.exists(filepath):
        print(f"Error: File not found {filepath}")
        sys.exit(1)
    
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
    """Calculates % H, E, C for each residue."""
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
    """Splits data and residues into Left and Right segments based on break."""
    mask_left = residues <= RES_BREAK_START
    mask_right = residues >= RES_BREAK_END
    
    return (residues[mask_left], data_array[mask_left]), (residues[mask_right], data_array[mask_right])

def plot_stacked_distribution(apo_data, lig_data, residues, apo_label, lig_label, output_file):
    """Generates the Stacked Bar Plot with Axis Break."""
    h_apo, e_apo, c_apo = apo_data
    h_lig, e_lig, c_lig = lig_data
    
    # Ensure data length matches residue length
    min_len = min(len(h_apo), len(h_lig), len(residues))
    
    # Trim all to match shortest
    residues = residues[:min_len]
    h_apo, e_apo, c_apo = h_apo[:min_len], e_apo[:min_len], c_apo[:min_len]
    h_lig, e_lig, c_lig = h_lig[:min_len], e_lig[:min_len], c_lig[:min_len]
    
    # Calculate width ratios
    len_left = np.sum(residues <= RES_BREAK_START)
    len_right = np.sum(residues >= RES_BREAK_END)
    width_ratios = [len_left, len_right]
    
    # Setup Figure
    fig = plt.figure(figsize=(14, 8))
    gs = gridspec.GridSpec(2, 2, width_ratios=width_ratios, wspace=0.05, hspace=0.3)
    
    colors = {'H': '#0072B2', 'E': '#D55E00', 'C': '#999999'}
    
    # Helper to plot one panel
    def plot_panel(ax, res, h, e, c, label, title=None, y_label=False, hide_y=False):
        if len(res) == 0: return # Skip empty panels
        ax.bar(res, h, color=colors['H'], label='Helix', width=1.0)
        ax.bar(res, e, bottom=h, color=colors['E'], label='Sheet', width=1.0)
        ax.bar(res, c, bottom=h+e, color=colors['C'], label='Coil', width=1.0)
        
        ax.set_ylim(0, 100)
        ax.set_xlim(res[0], res[-1]) # Tight X-axis
        
        # Uniform X-tick step
        tick_step = 25
        if len(res) > 0:
            start_tick = (res[0] // tick_step) * tick_step
            end_tick = (res[-1] // tick_step + 1) * tick_step
            ax.set_xticks(np.arange(start_tick, end_tick, tick_step))

        if title: ax.set_title(title, fontweight='bold')
        if y_label: ax.set_ylabel(f"{label}\nStructure %", fontweight='bold', fontsize=12)
        
        if hide_y:
            ax.set_yticks([])
            ax.spines['left'].set_visible(False)
        else:
            ax.spines['right'].set_visible(False)
            
    # --- Row 1: APO ---
    (res_l, h_l), (res_r, h_r) = split_data_by_break(h_apo, residues)
    (_, e_l), (_, e_r)         = split_data_by_break(e_apo, residues)
    (_, c_l), (_, c_r)         = split_data_by_break(c_apo, residues)
    
    ax1_l = fig.add_subplot(gs[0, 0])
    plot_panel(ax1_l, res_l, h_l, e_l, c_l, apo_label, title=f"{apo_label}", y_label=True)
    ax1_l.legend(loc='upper left', bbox_to_anchor=(0, 1.25), ncol=3, frameon=False) 
    
    ax1_r = fig.add_subplot(gs[0, 1])
    plot_panel(ax1_r, res_r, h_r, e_r, c_r, apo_label, title=None, hide_y=True)
    
    # --- Row 2: Ligand ---
    (res_l, h_l), (res_r, h_r) = split_data_by_break(h_lig, residues)
    (_, e_l), (_, e_r)         = split_data_by_break(e_lig, residues)
    (_, c_l), (_, c_r)         = split_data_by_break(c_lig, residues)
    
    ax2_l = fig.add_subplot(gs[1, 0])
    plot_panel(ax2_l, res_l, h_l, e_l, c_l, lig_label, y_label=True)
    ax2_l.set_xlabel("Residue Index", fontweight='bold', fontsize=14)
    
    ax2_r = fig.add_subplot(gs[1, 1])
    plot_panel(ax2_r, res_r, h_r, e_r, c_r, lig_label, hide_y=True)
    if len(res_r) > 0: # Check if right segment exists
        ax2_r.set_xlabel("Residue Index", fontweight='bold', fontsize=14) 
    
    plt.savefig(output_file, dpi=600)
    print(f"Saved: {output_file}")
    plt.close()

def plot_transition_bubbles(apo_data, lig_data, residues, apo_label, lig_label, output_file):
    """Generates the Bubble Plot with Axis Break."""
    h_apo, e_apo, c_apo = apo_data
    h_lig, e_lig, c_lig = lig_data
    
    # Trim to match residues
    min_len = min(len(h_apo), len(h_lig), len(residues))
    residues = residues[:min_len]
    h_apo, e_apo, c_apo = h_apo[:min_len], e_apo[:min_len], c_apo[:min_len]
    h_lig, e_lig, c_lig = h_lig[:min_len], e_lig[:min_len], c_lig[:min_len]
    
    delta_h = h_lig - h_apo
    delta_c = c_lig - c_apo
    delta_e = e_lig - e_apo
    
    # Filter Threshold
    threshold = 15.0
    mask_sig = np.abs(delta_h) > threshold
    
    if np.sum(mask_sig) == 0:
        print("No significant transitions found. Skipping bubble plot.")
        return

    # Colors
    colors = []
    for i in range(len(delta_h)):
        val = delta_h[i]
        if val > 0: colors.append('green')
        else:
            if delta_c[i] > delta_e[i]: colors.append('orange')
            else: colors.append('red')
    colors = np.array(colors)
    sizes = np.abs(delta_h) * 5
    
    # Filter data
    res_sig = residues[mask_sig]
    dh_sig = delta_h[mask_sig]
    sz_sig = sizes[mask_sig]
    col_sig = colors[mask_sig]
    
    # Split significant data
    (res_l, dh_l), (res_r, dh_r) = split_data_by_break(dh_sig, res_sig)
    (_, sz_l), (_, sz_r)         = split_data_by_break(sz_sig, res_sig)
    (_, col_l), (_, col_r)       = split_data_by_break(col_sig, res_sig)
    
    # Width Ratios
    len_left = np.sum(residues <= RES_BREAK_START)
    len_right = np.sum(residues >= RES_BREAK_END)
    if len_right == 0: width_ratios = [len_left, 1]
    else: width_ratios = [len_left, len_right]
    
    fig = plt.figure(figsize=(14, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=width_ratios, wspace=0.05)
    
    # Helper
    def plot_bubble_panel(ax, res, dh, sz, col, show_y=False):
        ax.axhline(0, color='black', linewidth=1, linestyle='-')
        ax.axhline(threshold, color='grey', linewidth=0.5, linestyle='--')
        ax.axhline(-threshold, color='grey', linewidth=0.5, linestyle='--')
        
        if len(res) > 0:
            ax.scatter(res, dh, s=sz, c=col, alpha=0.7, edgecolors='black')
            
        ax.set_xlabel("Residue Index", fontweight='bold', fontsize=14)
        
        # Uniform X-tick step
        tick_step = 25
        if len(res) > 0:
            start_tick = (res[0] // tick_step) * tick_step
            end_tick = (res[-1] // tick_step + 1) * tick_step
            ax.set_xticks(np.arange(start_tick, end_tick, tick_step))
        
        if show_y:
            ax.set_ylabel(f"Δ Helicity (%)\n({lig_label} - {apo_label})", fontweight='bold', fontsize=12)
            # Set Limits for Left Panel: Start of residues to Break Start
            ax.set_xlim(residues[0], RES_BREAK_START)
            ax.spines['right'].set_visible(False)
        else:
            ax.set_yticks([])
            ax.spines['left'].set_visible(False)
            # Set Limits for Right Panel: Break End to End of residues
            if len(residues) > 0 and residues[-1] >= RES_BREAK_END:
                 ax.set_xlim(RES_BREAK_END, residues[-1])
            
        ax.grid(True, alpha=0.3)

    ax_l = fig.add_subplot(gs[0])
    plot_bubble_panel(ax_l, res_l, dh_l, sz_l, col_l, show_y=True)
    
    ax_r = fig.add_subplot(gs[1])
    plot_bubble_panel(ax_r, res_r, dh_r, sz_r, col_r, show_y=False)
    
    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='green', label='Helix Formation', markersize=10),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', label='Unfolding -> Coil', markersize=10),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red', label='Unfolding -> Sheet', markersize=10)
    ]
    fig.legend(handles=legend_elements, loc='upper center', ncol=3, frameon=False)
    
    plt.savefig(output_file, dpi=600)
    print(f"Saved: {output_file}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Analyze DSSP transitions.")
    parser.add_argument("--apo", required=True, help="Path to APO ss.dat")
    parser.add_argument("--ligand", required=True, help="Path to Ligand ss.dat")
    parser.add_argument("--pdb", required=True, help="Path to PDB file for residue mapping")
    parser.add_argument("--label_apo", required=True, help="Label for APO system")
    parser.add_argument("--label_lig", required=True, help="Label for Ligand system")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    
    args = parser.parse_args()
    
    print(f"Analyzing {args.label_apo} vs {args.label_lig}...")
    
    # Get Residues
    residues = get_residue_indices(args.pdb)
    
    mat_apo = parse_dssp_file(args.apo)
    mat_lig = parse_dssp_file(args.ligand)
    
    data_apo = analyze_occupancy(mat_apo)
    data_lig = analyze_occupancy(mat_lig)
    
    out1 = os.path.join(args.output_dir, f"dssp_stacked_{args.label_apo}_vs_{args.label_lig}.png")
    plot_stacked_distribution(data_apo, data_lig, residues, args.label_apo, args.label_lig, out1)
    
    out2 = os.path.join(args.output_dir, f"dssp_bubble_{args.label_apo}_vs_{args.label_lig}.png")
    plot_transition_bubbles(data_apo, data_lig, residues, args.label_apo, args.label_lig, out2)

if __name__ == "__main__":
    main()