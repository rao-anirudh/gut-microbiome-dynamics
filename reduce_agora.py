from cobra.io import read_sbml_model
from ete3 import NCBITaxa
from collections import defaultdict
import os
import numpy as np
import json

# Directory containing SBML model files
model_dir = "AGORA_1_03_sbml"
# List of all model files with XML extension in the directory
model_files = [f for f in os.listdir(model_dir) if f.endswith(".xml")]
# Initialize the NCBITaxa object to fetch taxonomic information
ncbi = NCBITaxa()

# Map strain IDs to their corresponding model file names
strain_to_file = {f.split(".")[0]: f for f in model_files}
# Dictionary to store phylum information for each strain
strain_to_phylum = {}
# Dictionary to store strains grouped by phylum
phylum_to_strains = defaultdict(list)

# Loop over each strain to gather phylum information
for strain_id in strain_to_file:
    try:
        # Extract genus and species from strain_id
        name = strain_id.replace("_", " ").split(" ")[0:2]
        # Get the taxonomic ID from the strain name using NCBITaxa
        taxid = ncbi.get_name_translator([" ".join(name)])
        if not taxid:
            continue  # Skip if taxonomic ID is not found
        # Extract taxonomic ID and lineage
        taxid = list(taxid.values())[0][0]
        lineage = ncbi.get_lineage(taxid)
        # Get taxonomic names and ranks along the lineage
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)
        # Extract the phylum name from the lineage
        phylum = [names[t] for t in lineage if ranks[t] == "phylum"]
        if phylum:
            # Assign phylum to strain and add strain to the corresponding phylum list
            strain_to_phylum[strain_id] = phylum[0]
            phylum_to_strains[phylum[0]].append(strain_id)
    except Exception as e:
        # Log any error that occurs during processing of the strain
        print(f"Skipping {strain_id}: {e}")

# Dictionary to store the representative strain for each phylum
representatives = {}

# Loop over each phylum and calculate representative strain based on Jaccard distance
for phylum, strains in phylum_to_strains.items():

    print("\n", phylum, len(strains))

    # If only one strain exists for the phylum, it is automatically the representative
    if len(strains) == 1:
        representatives[phylum] = (strains[0], 0, 0)
        continue

    reaction_sets = []
    strain_ids = []

    # Loop over strains in the current phylum to gather their reactions
    for sid in strains:
        print(sid)
        model_path = os.path.join(model_dir, strain_to_file[sid])
        try:
            # Read the SBML model for the strain
            model = read_sbml_model(model_path)
            # Collect reactions for this model
            reaction_sets.append(set(r.id for r in model.reactions))
            strain_ids.append(sid)
        except:
            continue  # Skip strains that cannot be read

    # If fewer than two reaction sets were found, assign the first strain as representative
    if len(reaction_sets) < 2:
        representatives[phylum] = (strain_ids[0], 0, 0)
        continue

    n = len(reaction_sets)
    # Initialize a matrix to store pairwise Jaccard distances between reaction sets
    dist_matrix = np.zeros((n, n))

    # Calculate pairwise Jaccard distances between reaction sets
    for i in range(n):
        for j in range(i+1, n):
            # Jaccard distance = 1 - (intersection size / union size)
            d = 1 - len(reaction_sets[i] & reaction_sets[j]) / len(reaction_sets[i] | reaction_sets[j])
            dist_matrix[i, j] = dist_matrix[j, i] = d

    # Calculate average distance for each strain
    avg_dist = dist_matrix.mean(axis=1)
    # Find the index of the strain with the smallest average distance
    rep_idx = np.argmin(avg_dist)
    # Store the representative strain and distance metrics for the phylum
    representatives[phylum] = (strain_ids[rep_idx], np.min(avg_dist), np.mean(avg_dist))

# Print the selected representative strain for each phylum
print("\nRepresentative strains per phylum:")
for phylum, (strain, min_avg_dist, mean_avg_dist) in representatives.items():
    print(f"{phylum}: {strain}")

# Define output file path to save the representative strains data
output_path = "representative_strains.json"

# Prepare data in a format suitable for JSON serialization
serialisable_dict = {
    phylum: {
        "phylum_size": len(phylum_to_strains[phylum]),
        "representative_strain": strain,
        "min_jaccard_distance": float(min_dist),
        "avg_jaccard_distance": float(avg_dist)
    }
    for phylum, (strain, min_dist, avg_dist) in representatives.items()
}

# Write the serializable dictionary to a JSON file
with open(output_path, "w") as f:
    json.dump(serialisable_dict, f, indent=4)

print(f"\nSaved representative strains to {output_path}")
