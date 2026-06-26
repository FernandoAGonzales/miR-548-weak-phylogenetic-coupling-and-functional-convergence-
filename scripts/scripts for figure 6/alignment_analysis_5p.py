#!/usr/bin/env python3

"""
Figure 6B reproducibility script
================================

Manuscript
----------
Functional convergence of the miR-548 family reveals weak coupling
between evolution and regulatory roles in human disease

Purpose
-------
This script performs the sequence variability analysis used to generate
Figure 6B of the manuscript. It analyzes a multiple sequence alignment
of mature human miR-548-5p sequences and calculates positional sequence
variability metrics.

Analyses performed
------------------
• Nucleotide frequencies at each alignment position
• Shannon entropy (H) per nucleotide position
• Conservation score per nucleotide position
• Gap frequency per nucleotide position
• Statistical comparison of entropy values between canonical seed and
  non-seed regions using the Mann–Whitney U test

Seed region
-----------
Canonical miRNA seed positions:
6–12 (alignment-normalized positions)

Input
-----
FASTA multiple sequence alignment of mature miR-548-5p sequences.

Output
------
CSV file containing:

- Position
- Shannon entropy
- Conservation score
- Gap frequency

Software
--------
Python 3.12

Dependencies
------------
NumPy
Pandas
Biopython
SciPy

Reproducibility
---------------
This script forms part of the reproducible computational workflow
used to generate Figure 6 of the manuscript. Together with
alignment_analysis_3p.py and figure6_reproducibility.R, it enables
complete reproduction of the entropy analyses reported in the study.

Author
------
Damaris Asanza-Sanmartin, Enrique A. Solorzano-Toala, Fernando A. Gonzales-Zubiate et al.
"""

import argparse
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
from Bio import AlignIO
from scipy.stats import mannwhitneyu


def nucleotide_frequencies(alignment):
    freq_list = []

    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]

        counts = Counter(column)
        total = sum(counts.values())

        freqs = {
            nt: counts.get(nt, 0) / total
            for nt in ["A", "U", "G", "C", "-"]
        }

        freq_list.append(freqs)

    return pd.DataFrame(freq_list)


def shannon_entropy(column):

    counts = Counter(column)
    total = sum(counts.values())

    entropy = 0

    for nt in counts:
        p = counts[nt] / total

        if p > 0:
            entropy -= p * np.log2(p)

    return entropy


def conservation_score(column):

    counts = Counter(column)
    total = sum(counts.values())

    most_common = counts.most_common(1)[0][1]

    return most_common / total


def analyze_alignment(alignment):

    freq_df = nucleotide_frequencies(alignment)

    entropy_values = []
    conservation_values = []

    for i in range(alignment.get_alignment_length()):

        column = alignment[:, i]

        entropy_values.append(
            shannon_entropy(column)
        )

        conservation_values.append(
            conservation_score(column)
        )

    results = pd.DataFrame({
        "Position": range(
            1,
            len(entropy_values) + 1
        ),
        "Entropy": entropy_values,
        "Conservation": conservation_values
    })

    results["Gap_frequency"] = freq_df["-"]

    return results


def perform_seed_test(results):

    # Canonical miRNA seed region (biological numbering)
    seed_start = 6
    seed_end = 12

    # Convert biological positions (1-based) to Python indices (0-based)
    seed_positions = list(range(seed_start - 1, seed_end))

    seed_entropy = results.iloc[seed_positions]["Entropy"]
    non_seed_entropy = results.drop(seed_positions)["Entropy"]

    stat, p_value = mannwhitneyu(
        seed_entropy,
        non_seed_entropy,
        alternative="two-sided",
        method="auto",
        use_continuity=True
    )

    return stat, p_value

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Calculate nucleotide frequencies, "
            "Shannon entropy, conservation scores "
            "and seed-region statistics from "
            "a FASTA alignment."
        )
    )

    parser.add_argument(
        "--alignment",
        required=True,
        type=Path,
        help="Input FASTA alignment"
    )

    parser.add_argument(
        "--output",
        default=Path(
            "sequence_analysis_results.csv"
        ),
        type=Path,
        help="Output CSV file"
    )

    args = parser.parse_args()

    alignment = AlignIO.read(
        args.alignment,
        "fasta"
    )

    results = analyze_alignment(
        alignment
    )

    stat, p_value = perform_seed_test(
        results
    )

    print("\nMann-Whitney U test")
    print("-------------------")
    print(f"Statistic: {stat}")
    print(f"P-value: {p_value}")

    results.to_csv(
        args.output,
        index=False
    )

    print(
        f"\nAnalysis completed. Results saved to: {args.output}"
    )


if __name__ == "__main__":
    main()