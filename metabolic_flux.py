#############################################################################################################
# This file contains a Python implementation and adaptation of the Matlab algorithm used by Granata et al.
# (Granata et al. BMC Bioinformatics 2019, 20(Suppl 4):162 https://doi.org/10.1186/s12859-019-2685-9)
# The paper itself uses the Lee-12 algorithm (Matlab) for data integration and model simulation
# (Lee et al. BMC Systems Biology 2012, 6:73 https://doi.org/10.1186/1752-0509-6-73)
#
# COBRApy, the openCOBRA project (https://opencobra.github.io/cobrapy/)
# COBRA Toolbox (https://opencobra.github.io/cobratoolbox/stable/; https://doi.org/10.1038/s41596-018-0098-2)
# Gurobipy obtained form Gurobi and used with an academic license (https://www.gurobi.com/)
#
# Daniel Ribeiro 2020
#############################################################################################################


# System / Python
import os
import sys
import re
import tempfile
from pprint import pprint
from dataclasses import dataclass
import argparse
# Pandas and numpy
import pandas as pd
import numpy as np
# COBRA - cobrapy
import cobra
# Import extensions to Cobrapy
import cobrapy_extensions
# import 'SMOP' converted and reviewed 'Lee12' script
import lee12py
from libsmop import *


######
#Functions
######

def genesToCobraGenes(gene_names: list, cobra_model: cobra.core.Model, attribute: str) -> tuple:
    '''
    Returns tuple(dict_list, dict_names)
    dict_list: Cobrapy DictList of Cobra Genes mapped from our data into the model)
    dict_names: dict of key(gene names): value(Cobra Gene)

    The same name can be mappted to several Cobra Genes (eg. 31_AT1, 31_AT2, 31_AT3, 31_AT4, 31_AT5).
    The number in the ID string, corresponds to the NCBI Entrez ID. Therefore, they signify the same gene (possibly different isoforms).
    Therefore, if a gene name is found, collect only the 1st instance of the Cobra Gene.
    '''
    result = []
    dict_names = {}
    for gene in gene_names:
        regex = re.compile('^{match}$'.format(match=gene))
        query = cobra_model.genes.query(regex, attribute=attribute)
        if len(query) > 0:
            result.append(query[0])
            dict_names[gene] = query[0]
    dict_list = cobra.core.DictList(result)
    return dict_list, dict_names


######
#CODE
######
pd.options.mode.use_inf_as_na = True
idx=pd.IndexSlice


# Define arguments
parser = argparse.ArgumentParser(description="Perform flux balance analysis (FBA) on RNA-seq data.\n\
                                            ")
parser.add_argument("-f", "--file", required=True, help="path/to/rnaseq/file. Requires a Gene per line, a Mean expression (norm counts) column and a SD column.")
parser.add_argument("-g", "--gene_col", required=True, help="Gene column name.")
parser.add_argument("-e", "--expression_col", required=True, help="Expression column name.")
parser.add_argument("-s", "--sd_col", required=True, help="SD column name.")
parser.add_argument("-F", "--measured_fluxes", required=True, help="path/to/measured/fluxes/file. Format:\n\
                    First column are reaction names in the chosen metabolic model.\n\
                    Second column is the measured flux. \n\
                    By convention, exchange reactions are written as export reactions (e.g. ‘glc[e] <==>’),\n\
                    so import of a metabolite is a negative flux and so export of a metabolite is a positive flux\n")
parser.add_argument("-m", "--model", required=True, help="path/to/metabolic/model/file, in XML format.")
parser.add_argument("-o", "--output_path", help="path/to/output/folder.", default="./")
parser.add_argument("-n", "--name", help="Optional name for the analysis", default="flux")
args = parser.parse_args()

# Import tables
print("\nImport tables")
df_rna = pd.read_csv(args.file, sep='\t', na_filter=False, index_col=0, header=0)
df_fluxes = pd.read_csv(args.measured_fluxes, sep='\t', na_filter=False, index_col=0, header=0)
print("Imported RNA-seq:")
print(df_rna)

# Build experimental fluxes dict
experimental_fluxes = {}
for i in range(len(df_fluxes.index)):
    experimental_fluxes[df_fluxes.index[i]] = (df_fluxes.index[i], df_fluxes.iloc[i, 0])
print("Imported fluxes:")
print(experimental_fluxes)
# Build condition dict
condition = {}
condition = {'condition': args.name,
            'df_rna': df_rna.copy(),
            'mean_col': args.expression_col,
            'sd_col': args.sd_col,
            'gene_name': df_rna.index.values.copy(),
            'model': args.model,
            'output': args.output_path,
            'experimental_fluxes': experimental_fluxes
            }


# Create Cobra configuration
print("\n***Load COBRApy***")
cobra_config = cobra.Configuration()

# Display relevant Cobra configuration
print("InstalledPackage versions:")
cobra.show_versions()
cobra_config.solver = 'gurobi'	#Set solver
print("\n***COBRApy configuration***")
print(cobra_config)

# Load model in Cobra
print(f'''\nLoading model '{condition["model"]}' into COBRApy...''')
model = cobra.core.Model()
model = cobra.io.read_sbml_model(condition['model'])

#Model Properties
print("\n***Model description***")
print("Model summary: ", model)
print("ID: ", model.id)
print("Reactions no.: ", len(model.reactions))
print("Genes no.: ", len(model.genes))
print("Subcellular compartment no.: ", len(model.compartments))
#print("Model solver: ", model.solver) #WARNING: will print a lot of info
print("Model objective expression: ", model.objective.expression)
print("Model objective direction: ", model.objective.direction)
objective_reaction = cobra.util.solver.linear_reaction_coefficients(model)
print("Model objective reaction: ")
for key, value in objective_reaction.items():
    print(f"\t{key}", value)
    print("\tName: ", key.name)
    print("\tID: ", key.id)

print("\n***Experimental fluxes to constrain ***")
for key, val in condition['experimental_fluxes'].items():
    measured_flux = val[0] # eg. 'D-glucose exchange' (extracellular consumpion  of glucose)
    measured_flux = model.reactions.get_by_id(measured_flux)
    print("Reaction name: ", measured_flux.name)
    print("Bounds (default): ", measured_flux.bounds)
    print("Reversibility (default): ", measured_flux.reversibility)
#Prepare our data for Lee12 algorithm
print("\nPrepare data for Lee12 algorithm...")
print(f"Using condition: {condition['condition']}")
print(f"Using model: {str(model)}")
print("Clean gene names of data...")
clean_name = [name for name in condition['gene_name'] if name]
print("Lenght of gene names: ", len(clean_name))

print("Convert gene names into model Gene objects...")
gene_dictlist, gene_name_dict = genesToCobraGenes(clean_name, model, attribute='name')
print("Select Mean and SD of genes existing in the model")
boolean = condition['df_rna'].index.isin(gene_dictlist.list_attr('name'))
df_filter = condition['df_rna'].loc[idx[boolean], :].copy()
print("Select the corresponding Gene objects")
gene_name, *_ = genesToCobraGenes(df_filter.index.values, model, attribute='name')
condition['gene_name'] = gene_name
gene_mean = df_filter[condition['mean_col']].values
gene_sd = df_filter[condition['sd_col']].values
print("Filtered genes: ", len(df_filter.index))
##Uncomment to verify gene order
#for i in range(len(df_filter.index)):
#    print("Index: ", df_filter.iloc[i, [2, 3]],
#          "Cobra gene: ", gene_name.list_attr('name')[i])
print(f"\n***Execute Lee12 algorithm with: {condition['condition']} - WIP***")
#Run the gene expression constraint FBA, and get the resulting flux distribution
print("Run the gene expression constraint FBA, and get the resulting flux distribution")
fluxes = lee12py.call_Lee12(model=model,
                            gene_names=condition['gene_name'],
                            gene_exp=gene_mean,
                            gene_exp_sd=gene_sd,
                            #gene_names_cobra=gene_name_dict,
                            #gene_scale_rxn=gene_to_scale,
                            #flux_scale_rxn=flux_to_scale,
                            flux=condition['experimental_fluxes'],
                            model_id=condition['condition'])
#Match Reaction to Flux
#DEBUG
print(f"len(fluxes) ==  len(model.reactions): {fluxes.size == len(model.reactions)}")

print("Save predicted reaction flux table...")
flux_table = {'Reaction.Name': [],
              'Reaction.Id': [],
              'Flux': []
              }
for k in range(len(model.reactions)):
    flux_table['Reaction.Name'].append(model.reactions[k].name)
    flux_table['Reaction.Id'].append(model.reactions[k].id)
    flux_table['Flux'].append(fluxes[k])
flux_table = pd.DataFrame.from_dict(flux_table)
flux_table.to_csv(f"./{condition['output']}/flux_table_{condition['condition']}.txt", sep='\t', index=False)
print("Save the model as COBRApy model in SMBL and JSON formats...")
cobra.io.write_sbml_model(model, f"./{condition['output']}/{condition['condition']}.xml")
cobra.io.save_json_model(model, f"./{condition['output']}/{condition['condition']}.json", pretty=True)

print("******", "\nDONE", "\n******")
