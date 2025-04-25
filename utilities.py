import numpy as np
import os
from cobra.io import read_sbml_model
import random
import concurrent.futures
import multiprocessing
import copy


def read_sbml_with_timeout(filepath, timeout=5):
    import concurrent.futures
    import cobra.io

    def load_model():
        return cobra.io.read_sbml_model(filepath)

    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        future = executor.submit(load_model)
        try:
            model = future.result(timeout=timeout)
            return model
        except concurrent.futures.TimeoutError:
            return None


class SmallIntestine:

    def __init__(self):
        self.metabolome = dict()  # in mmol
        self.microbiome = dict()  # in cell counts
        self.model = read_sbml_model("MODEL1310110020_url_small.xml")
        self.growth_rate = float
        self.input_frequency = 24  # in hours
        self.output_frequency = 4  # in hours
        self.biomass = 640  # in gDCW

    def add_to_metabolome(self, metabolites):
        self.metabolome = metabolites

    def add_to_microbiome(self, microbes):
        for microbe in microbes.keys():
            if microbe in self.microbiome.keys():
                self.microbiome[microbe] = self.microbiome[microbe] + microbes[microbe]
            else:
                self.microbiome[microbe] = microbes[microbe]

    def process_species(self, metabolome, species, total_biomass):
        path_to_agora = "AGORA_1_03_sbml"
        filepath = os.path.join(path_to_agora, species)

        try:
            model = read_sbml_with_timeout(filepath)
            bacterial_cell_volume = 1e-12  # in cm^3
            dry_weight_per_unit_volume = 0.33  # in gDCW/cm^3
            biomass = self.microbiome[species] * bacterial_cell_volume * dry_weight_per_unit_volume  # in gDCW

            for exchange in model.exchanges:
                metabolite = list(exchange.metabolites.keys())[0].id
                if metabolite in metabolome.keys():
                    availability = metabolome[metabolite]
                    species_share = availability * (biomass / total_biomass)
                    exchange.lower_bound = min(-1e-6, round(-species_share / (biomass * self.output_frequency), 3))
                else:
                    exchange.lower_bound = -1e-6

            solution = model.optimize()
            growth_rate = solution.objective_value
            new_cell_count = int(self.microbiome[species] * np.exp(growth_rate * self.output_frequency))
            self.microbiome[species] = new_cell_count

            exchanges = dict()
            for exchange in model.exchanges:
                metabolite = list(exchange.metabolites.keys())[0].id
                exchange_flux = solution[exchange.id]
                exchange_amount = exchange_flux * biomass * self.output_frequency
                exchanges[metabolite] = exchange_amount

            return species, growth_rate, exchanges

        except:
            return species, 0, dict()

    def metabolise(self):

        total_biomass = 0
        for species in self.microbiome.keys():
            bacterial_cell_volume = 1e-12  # in cm^3
            dry_weight_per_unit_volume = 0.33  # in gDCW/cm^3
            biomass = self.microbiome[species] * bacterial_cell_volume * dry_weight_per_unit_volume  # in gDCW
            total_biomass += biomass

        growth_rates = {species: 0 for species in self.microbiome.keys()}
        combined_exchanges = {}

        num_cpus = os.cpu_count()
        current_metabolome = copy.deepcopy(self.metabolome)
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpus) as executor:
            futures = {executor.submit(self.process_species, current_metabolome, species, total_biomass): species for
                       species in self.microbiome.keys()}
            for future in concurrent.futures.as_completed(futures):
                try:
                    species, growth_rate, exchanges = future.result()
                    if growth_rate is not None:
                        growth_rates[species] = growth_rate
                    for metabolite, amount in exchanges.items():
                        combined_exchanges[metabolite] = combined_exchanges.get(metabolite, 0) + amount
                except:
                    continue

        for metabolite, amount in combined_exchanges.items():
            if metabolite in self.metabolome and self.metabolome[metabolite] + amount != 0:
                self.metabolome[metabolite] += amount
            elif metabolite not in self.metabolome and amount != 0:
                self.metabolome[metabolite] = amount

        model = self.model
        for exchange in model.exchanges:
            metabolite = list(exchange.metabolites.keys())[0].id
            if metabolite in self.metabolome:
                availability = self.metabolome[metabolite]
                exchange.lower_bound = min(-1e-6, round(-availability / (self.biomass * self.output_frequency), 3))
            else:
                exchange.lower_bound = -1e-6
        solution = model.optimize()
        self.growth_rate = solution.objective_value
        for exchange in model.exchanges:
            metabolite = list(exchange.metabolites.keys())[0].id
            exchange_flux = solution[exchange.id]
            exchange_amount = exchange_flux * self.biomass * self.output_frequency
            if metabolite in self.metabolome:
                self.metabolome[metabolite] += exchange_amount
            elif metabolite not in self.metabolome and exchange_amount != 0:
                self.metabolome[metabolite] = exchange_amount

        return growth_rates

    def transfer(self, LargeIntestine, growth_rates):

        LargeIntestine.add_to_metabolome(self.metabolome)
        self.metabolome = dict()

        total_microbiome_growth = sum(growth_rates.values())
        transfer_probability = {species: 1 - (gr / total_microbiome_growth) for species, gr in growth_rates.items()}
        num_cells_to_transfer = sum(self.microbiome.values()) - random.randint(10 ** 3, 10 ** 8)
        cells_transferred = {species: 0 for species in self.microbiome.keys()}
        num_cells_transferred = 0
        while num_cells_transferred < num_cells_to_transfer:
            chosen_species = \
                random.choices(list(transfer_probability.keys()), weights=list(transfer_probability.values()), k=1)[0]
            num_species_cells_transferred = random.randint(0, self.microbiome[chosen_species])
            self.microbiome[chosen_species] -= num_species_cells_transferred
            num_cells_transferred += num_species_cells_transferred
            cells_transferred[chosen_species] += num_species_cells_transferred
        cells_transferred = {species: count for species, count in cells_transferred.items() if count != 0}
        LargeIntestine.add_to_microbiome(cells_transferred)
        self.microbiome = {species: count for species, count in self.microbiome.items() if count != 0}


class LargeIntestine:

    def __init__(self):
        self.metabolome = dict()  # in mmol
        self.microbiome = dict()  # in cell counts
        self.model = read_sbml_model("MODEL1310110043_url_large_cleaned.xml")
        self.growth_rate = float
        self.input_frequency = 4  # in hours
        self.output_frequency = 24  # in hours
        self.biomass = 370  # in gDCW

    def add_to_metabolome(self, metabolites):
        for metabolite in metabolites.keys():
            if metabolite in self.metabolome.keys():
                self.metabolome[metabolite] = self.metabolome[metabolite] + metabolites[metabolite]
            else:
                self.metabolome[metabolite] = metabolites[metabolite]

    def add_to_microbiome(self, microbes):
        for microbe in microbes.keys():
            if microbe in self.microbiome.keys():
                self.microbiome[microbe] = self.microbiome[microbe] + microbes[microbe]
            else:
                self.microbiome[microbe] = microbes[microbe]

    def process_species(self, metabolome, species, total_biomass):
        path_to_agora = "AGORA_1_03_sbml"
        filepath = os.path.join(path_to_agora, species)

        try:
            model = read_sbml_with_timeout(filepath)
            bacterial_cell_volume = 1e-12  # in cm^3
            dry_weight_per_unit_volume = 0.33  # in gDCW/cm^3
            biomass = self.microbiome[species] * bacterial_cell_volume * dry_weight_per_unit_volume  # in gDCW

            exchanges = dict()
            for exchange in model.exchanges:
                metabolite = list(exchange.metabolites.keys())[0].id
                if metabolite in metabolome.keys():
                    availability = metabolome[metabolite]
                    species_share = availability * (biomass / total_biomass)
                    exchange.lower_bound = min(-1e-6, round(
                        -species_share / (biomass * (self.output_frequency - self.input_frequency)), 3))
                else:
                    exchange.lower_bound = -1e-6

            solution = model.optimize()
            growth_rate = solution.objective_value
            new_cell_count = int(self.microbiome[species] * np.exp(
                growth_rate * (self.output_frequency - self.input_frequency)))
            self.microbiome[species] = new_cell_count

            for exchange in model.exchanges:
                metabolite = list(exchange.metabolites.keys())[0].id
                exchange_flux = solution[exchange.id]
                exchange_amount = exchange_flux * biomass * (self.output_frequency - self.input_frequency)
                exchanges[metabolite] = exchanges.get(metabolite, 0) + exchange_amount

            return species, growth_rate, exchanges

        except concurrent.futures.TimeoutError:
            return species, 0, dict()

    def metabolise(self):
        total_biomass = 0
        for species in self.microbiome.keys():
            bacterial_cell_volume = 1e-12  # in cm^3
            dry_weight_per_unit_volume = 0.33  # in gDCW/cm^3
            biomass = self.microbiome[species] * bacterial_cell_volume * dry_weight_per_unit_volume  # in gDCW
            total_biomass += biomass

        growth_rates = {species: 0 for species in self.microbiome.keys()}
        combined_exchanges = {}

        num_cpus = os.cpu_count()
        current_metabolome = copy.deepcopy(self.metabolome)
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpus) as executor:
            futures = {executor.submit(self.process_species, current_metabolome, species, total_biomass): species for
                       species in self.microbiome.keys()}
            for future in concurrent.futures.as_completed(futures):
                try:
                    species, growth_rate, exchanges = future.result()
                    if growth_rate is not None:
                        growth_rates[species] = growth_rate
                    for metabolite, amount in exchanges.items():
                        combined_exchanges[metabolite] = combined_exchanges.get(metabolite, 0) + amount
                except:
                    continue

        for metabolite, amount in combined_exchanges.items():
            if metabolite in self.metabolome and self.metabolome[metabolite] + amount != 0:
                self.metabolome[metabolite] += amount
            elif metabolite not in self.metabolome and amount != 0:
                self.metabolome[metabolite] = amount

        model = self.model
        for exchange in model.exchanges:
            metabolite = list(exchange.metabolites.keys())[0].id
            if metabolite in self.metabolome:
                availability = self.metabolome[metabolite]
                exchange.lower_bound = min(-1e-6, round(
                    -availability / (self.biomass * (self.output_frequency - self.input_frequency)), 3))
            else:
                exchange.lower_bound = -1e-6
        solution = model.optimize()
        self.growth_rate = solution.objective_value
        for exchange in model.exchanges:
            metabolite = list(exchange.metabolites.keys())[0].id
            exchange_flux = solution[exchange.id]
            exchange_amount = exchange_flux * self.biomass * (self.output_frequency - self.input_frequency)
            if metabolite in self.metabolome:
                self.metabolome[metabolite] += exchange_amount
            elif metabolite not in self.metabolome and exchange_amount != 0:
                self.metabolome[metabolite] = exchange_amount

        return growth_rates

    def transfer(self, growth_rates):
        self.metabolome = dict()

        total_microbiome_growth = sum(growth_rates.values())
        transfer_probability = {species: 1 - (gr / total_microbiome_growth) for species, gr in growth_rates.items()}
        num_cells_to_transfer = sum(self.microbiome.values()) - random.randint(10 ** 8, 10 ** 10)
        num_cells_transferred = 0
        while num_cells_transferred < num_cells_to_transfer:
            chosen_species = \
                random.choices(list(transfer_probability.keys()), weights=list(transfer_probability.values()), k=1)[0]
            num_species_cells_transferred = random.randint(0, self.microbiome[chosen_species])
            self.microbiome[chosen_species] -= num_species_cells_transferred
            num_cells_transferred += num_species_cells_transferred
        self.microbiome = {species: count for species, count in self.microbiome.items() if count != 0}
