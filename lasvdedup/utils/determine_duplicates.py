#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_pairwise_distances(alignment_file, output_csv, outlier_file, threshold=0.1):
    """
    Calculate pairwise distances between sequences in an alignment
    and identify outliers with unusually high average distances
    """
    print(f"Processing alignment: {alignment_file}")

    # Read the alignment
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
    except Exception as e:
        print(f"Error reading alignment: {e}")
        return

    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')

    try:
        distance_matrix = calculator.get_distance(alignment)
    except Exception as e:
        print(f"Error calculating distances: {e}")
        return

    # Convert to pandas DataFrame for easier manipulation
    seq_names = [seq.id for seq in alignment]
    n_seqs = len(seq_names)

    distances_df = pd.DataFrame(np.zeros((n_seqs, n_seqs)),
                               index=seq_names,
                               columns=seq_names)

    # Fill the dataframe with distances
    for i in range(n_seqs):
        for j in range(n_seqs):
            if i != j:  # Skip diagonal (self-comparisons)
                distances_df.iloc[i, j] = distance_matrix[i, j]

    # Save the full distance matrix
    distances_df.to_csv(output_csv)

    # Calculate average distance for each sequence to all others
    avg_distances = distances_df.mean(axis=1)

    # Identify outliers (sequences with unusually high average distances)
    overall_mean = avg_distances.mean()
    overall_std = avg_distances.std()
    threshold_value = overall_mean + (overall_std * threshold)

    outliers = avg_distances[avg_distances > threshold_value]

    # Save outliers
    with open(outlier_file, 'w') as f:
        f.write(f"# Average distance threshold: {threshold_value:.4f}\n")
        f.write(f"# Overall mean: {overall_mean:.4f}, Std: {overall_std:.4f}\n\n")
        f.write("Sequence\tAverage_Distance\n")
        for seq_id, distance in outliers.items():
            f.write(f"{seq_id}\t{distance:.4f}\n")

    # Create a histogram of the average distances
    plt.figure(figsize=(10, 6))
    sns.histplot(avg_distances, kde=True)
    plt.axvline(threshold_value, color='red', linestyle='--',
                label=f'Threshold ({threshold_value:.4f})')
    plt.title(f"Distribution of Average Pairwise Distances\n{os.path.basename(alignment_file)}")
    plt.xlabel("Average Distance")
    plt.ylabel("Count")
    plt.legend()

    # Save the plot
    plot_file = os.path.join(os.path.dirname(outlier_file),
                            f"{os.path.basename(alignment_file)}.distance_histogram.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')

    print(f"Analysis complete. Found {len(outliers)} outliers.")
    return distances_df, outliers

if __name__ == "__main__" or "snakemake" in locals():
    # For Snakemake integration
    if "snakemake" in locals():
        alignment_file = snakemake.input.alignment
        output_csv = snakemake.output.distances
        outlier_file = snakemake.output.outliers
        threshold = snakemake.params.get("threshold", 0.1)
    else:
        # For standalone execution
        if len(sys.argv) < 3:
            print("Usage: python calculate_distances.py <alignment_file> <output_csv> <outlier_file> [threshold]")
            sys.exit(1)

        alignment_file = sys.argv[1]
        output_csv = sys.argv[2]
        outlier_file = sys.argv[3]
        threshold = float(sys.argv[4]) if len(sys.argv) > 4 else 0.1

    calculate_pairwise_distances(alignment_file, output_csv, outlier_file, threshold)
