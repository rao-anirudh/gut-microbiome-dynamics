import os
from cobra.io import read_sbml_model
import csv

# Path to the folder containing all AGORA SBML models
path_to_agora = "AGORA_1_03_sbml"
agora_files = os.listdir(path_to_agora)

# Dictionaries to map metabolite names to IDs and vice versa
metabolites_names_id = dict()
metabolites_id_names = dict()

# Loop over each AGORA model file
for file in agora_files:
    filepath = os.path.join(path_to_agora, file)

    # Load the metabolic model from SBML
    model = read_sbml_model(filepath)
    print(file, len(model.exchanges))  # Log the number of exchange reactions

    # For each exchange reaction, extract the metabolite being exchanged
    for exchange in model.exchanges:
        metabolite = list(exchange.metabolites.keys())[0]  # Only one metabolite per exchange

        # Map metabolite name → ID(s)
        if metabolite.name not in metabolites_names_id:
            metabolites_names_id[metabolite.name] = {metabolite.id}
        else:
            metabolites_names_id[metabolite.name].add(metabolite.id)

        # Map metabolite ID → name(s)
        if metabolite.id not in metabolites_id_names:
            metabolites_id_names[metabolite.id] = {metabolite.name}
        else:
            metabolites_id_names[metabolite.id].add(metabolite.name)

# Load and process the small intestine model
small_intestine_model = "MODEL1310110020_url_small.xml"
model = read_sbml_model(small_intestine_model)
for exchange in model.exchanges:
    metabolite = list(exchange.metabolites.keys())[0]

    if metabolite.name not in metabolites_names_id:
        metabolites_names_id[metabolite.name] = {metabolite.id}
    else:
        metabolites_names_id[metabolite.name].add(metabolite.id)

    if metabolite.id not in metabolites_id_names:
        metabolites_id_names[metabolite.id] = {metabolite.name}
    else:
        metabolites_id_names[metabolite.id].add(metabolite.name)

# Load and process the large intestine model
large_intestine_model = "MODEL1310110043_url_large_cleaned.xml"
model = read_sbml_model(large_intestine_model)
for exchange in model.exchanges:
    metabolite = list(exchange.metabolites.keys())[0]

    if metabolite.name not in metabolites_names_id:
        metabolites_names_id[metabolite.name] = {metabolite.id}
    else:
        metabolites_names_id[metabolite.name].add(metabolite.id)

    if metabolite.id not in metabolites_id_names:
        metabolites_id_names[metabolite.id] = {metabolite.name}
    else:
        metabolites_id_names[metabolite.id].add(metabolite.name)

# Write the name → ID mapping to a CSV file
with open("metabolites_names_to_ids.csv", mode="w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Metabolite Name", "Metabolites IDs"])
    for met_name, met_ids in metabolites_names_id.items():
        writer.writerow([met_name, "; ".join(sorted(met_ids))])

# Write the ID → name mapping to another CSV file
with open("metabolites_ids_to_names.csv", mode="w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Metabolite ID", "Metabolites Names"])
    for met_id, met_names in metabolites_id_names.items():
        writer.writerow([met_id, "; ".join(sorted(met_names))])
