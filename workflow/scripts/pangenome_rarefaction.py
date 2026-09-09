import gzip
import itertools
import matplotlib.pyplot as plt
import re
import sys
import csv
import numpy as np
from collections import defaultdict

def build_growth_curve(gfa_path):
    node_presence = defaultdict(set)
    samples = set()
    
    print(f"Reading {gfa_path}...")
    with gzip.open(gfa_path, 'rt') as f:
        for line in f:
            if line.startswith('W\t'):
                parts = line.strip().split('\t')
                sample_name = parts[1]
                samples.add(sample_name)
                
                path_str = parts[6]
                nodes = re.findall(r'[><](\d+)', path_str)
                for node in nodes:
                    node_presence[node].add(sample_name)

    samples = list(samples)
    total_genomes = len(samples)
    print(f"Found {total_genomes} genomes across {len(node_presence)} total nodes.")

    growth_data = [] 
    node_sets = list(node_presence.values())
    
    for k in range(1, total_genomes + 1):
        combinations = list(itertools.combinations(samples, k))
        print(f"Processing k={k} ({len(combinations)} combinations)...")
        
        for combo in combinations:
            combo_set = set(combo)
            node_count = sum(1 for n_set in node_sets if not n_set.isdisjoint(combo_set))
            growth_data.append((k, node_count, ",".join(combo)))
            
    return growth_data

def plot_growth(growth_data, output_file):
    # Explicitly casting to Numpy arrays of floats to bypass the __array__ casting bug
    x_vals = np.array([d[0] for d in growth_data], dtype=float)
    y_vals = np.array([d[1] for d in growth_data], dtype=float)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.scatter(x_vals, y_vals, alpha=0.5, color='#1f77b4', edgecolor='none', s=30)
    
    unique_k = np.unique(x_vals)
    medians = []
    for k in unique_k:
        k_values = y_vals[x_vals == k]
        medians.append(np.median(k_values))
        
    ax.plot(unique_k, medians, color='black', linewidth=2, linestyle='--', label='Median Node Count')
    
    ax.set_xlabel('Number of Genomes')
    ax.set_ylabel('Absolute Node Count')
    ax.set_title('Pangenome Graph Growth')
    ax.set_xticks(unique_k)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rarefaction_curve.py <input.gfa.gz> <output_plot.png>")
        sys.exit(1)
        
    gfa_input = sys.argv[1]
    plot_output = sys.argv[2]
    tsv_output = plot_output.replace(".png", "_data.tsv")
    
    data = build_growth_curve(gfa_input)
    
    # Checkpoint: Save the data so you don't have to re-parse if plotting fails
    print(f"Saving raw data to {tsv_output}...")
    with open(tsv_output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['k_genomes', 'node_count', 'genomes_in_subset'])
        writer.writerows(data)
        
    print("Generating plot...")
    plot_growth(data, plot_output)