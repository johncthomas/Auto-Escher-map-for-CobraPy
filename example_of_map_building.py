from __future__ import print_function
__author__ = 'JT'

"""
Call generate_escher_map.gen_map() passing the Cobra model, a list of
the reactions you wish to include and a list of common intermediates you do
not wish to be used to create links between reactions. This returns an
EscherMap object. Use map_obj.dump_json() to get a JSON that can be passed
to escher.Builder()

The layout of the metabolites is produced by the igraph module using the
Fruchterman-Reingold force-directed algorithm. It may be possible to alter
the settings in the generate_escher_map.assign_positions() function;
see the IGraph documentation for more info. http://igraph.org/python/

Requires cobra, escher, igraph, json and prerequisites for those packages.
Information on installing Escher:
https://github.com/zakandrewking/escher

Known issues:
- Reactions that share the same primary metabolites are layed over the top of
each other and can end up with secondary metabolites further away.
- Unlinked reactions are doubled for some reason
- Label placing could be better
- unnecesary midnodes

"""

import cobra, cobra.test
import generate_escher_map
import escher.validate
import os

# Load a model.
#ecoli_model = selected_model = cobra.io.load_json_model(cobra.test.ecoli_json)
ecoli_model = cobra.io.read_sbml_model(cobra.test.ecoli_sbml)
# Escher requires cobra files in JSON
cobra.io.save_json_model(ecoli_model, "test.json")
# Set the objective reaction and optimise
objective = 'ALCD2x'
ecoli_model.change_objective(objective)
ecoli_model.reactions.get_by_id(objective).upper_bound = 1000
ecoli_model.optimize()
print('Flux of',objective, ecoli_model.solution)

# The flux for each reaction is stored as a dictionary (x_dict) in the model
flux = ecoli_model.solution.x_dict

# # The below code can be used to save the model
# file_name = 'temp_cobra_model.json'
# cobra.io.save_json_model(ecoli_model, file_name)

# To draw a map with generate_escher_map you need:
# * a list of reactions to be included
# * a dictionary of metabolite frequency (produced with
# generate_escher_map.metabolite_occurance()
# * a list of common metabolites that will not be used to link reactions

# Cobra id of reactions to be excluded from the map
#excluded_reactions = ('biomass', 'PROTEIN', 'DNA', 'PEPTIDO')
excluded_reactions = ()

common_mets = set()
rxn_with_flux = set()
rxn_with_greater_flux = set()

threshold_flux = 0.0001

# Get a list of reactions with flux and those where flux exceeds the threshold
# The rxn_with_flux will be used to determine what's a common metabolite
# rxn_with_greater flux used to draw the map
for rxn_name in flux:
    if flux[rxn_name] != 0 and rxn_name not in excluded_reactions:
        rxn = ecoli_model.reactions.get_by_id(rxn_name)
        rxn_with_flux.add(rxn)
        if 0-threshold_flux > flux[rxn_name] or threshold_flux < flux[rxn_name]:
            rxn_with_greater_flux.add(rxn)

# Get a dictionary of metabolite frequency
# (keys are the Cobra metabolite objects)
met_count = generate_escher_map.metabolite_occurence(rxn_with_flux)

# Get common metabolites.
# included_common_m will be used be treated as uncommon whatever their occurence
included_common_m = ('accoa_c', 'pry_c', 'g3p_c', 'acald_c', 'f6p_c', 'pep_c')
for met in met_count:
    if met_count[met] > 4 and met.id not in included_common_m:
        common_mets.add(met)

# Convert the sets to lists so that they can be iterated over
common_mets, rxn_with_flux, rxn_with_greater_flux = \
    list(common_mets), list(rxn_with_flux), list(rxn_with_greater_flux)


# Get the Escher JSON
escher_map = generate_escher_map.gen_map(
     ecoli_model, rxn_with_greater_flux,  common_mets, 350, met_count
)
escher_json = escher_map.dump_json()
# Currently Escher doesn't accept json strings, needs a file, i assume that's a bug
with open('escher.json', 'w') as f:
    f.write(escher_json)
# See escher documentation for more info on escher.
import escher
ecoli_metabolite_map = escher.Builder(
    #map_json=escher_json,
    map_json='escher.json',
    reaction_data=flux,
    model_json='test.json'
)
import json
parsed_json = json.loads(escher_json)
escher.validate.check_map(parsed_json)
# Use display in browser to get an interactive version of the map.
ecoli_metabolite_map.display_in_browser()
