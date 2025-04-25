import random
import json
import numpy as np


def sample_microbial_library(representative_data_path):
    """
    Simulates sampling of a large microbial library from a list of phylogenetically
    representative strains. The sampling probability for each strain is proportional
    to its phylum size.

    Since directly sampling 10^11 cells can cause integer overflows and memory issues,
    this function samples a smaller number (10^6) and scales the results up.

    Parameters:
        representative_data_path (str): Path to JSON file with representative strain data.

    Returns:
        dict: Mapping of strain ID to estimated cell count (total ≈ 10^11).
    """
    # Define realistic and safe sampling numbers
    N_realistic = random.randint(10**9, 10**11)  # Target total number of cells (scaled)
    N_sim = 10**6                                # Safe sample size for multinomial draw

    # Load representative strain data from JSON
    with open(representative_data_path, "r") as f:
        representative_data = json.load(f)

    # Prepare strain list and corresponding weights
    strain_ids = []
    weights = []

    for phylum, info in representative_data.items():
        if info["phylum_size"] > 2:
            strain_ids.append(info["representative_strain"])
            weights.append(info["phylum_size"])

    # Convert phylum sizes into probabilities
    probs = np.array(weights) / np.sum(weights)

    # Sample N_sim cells and scale to N_realistic
    sampled_counts = np.random.multinomial(N_sim, probs)
    scaled_counts = (sampled_counts / N_sim * N_realistic).round().astype(np.int64)

    # Construct the final result as a dictionary
    return {strain+".xml": int(count) for strain, count in zip(strain_ids, scaled_counts) if count > 0}
