import json
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
import textwrap
import matplotlib
import pandas as pd
matplotlib.use('TkAgg')

from helper_functions import *



"""
Based on Kemp et al., 2012, I want to check the enrichment of each target codon (the codon corresponding to each of the 
tRNAs that were KOd in the library) in the beginning (5') area (30-50 first codons) of the CDS for each of the proteins
in S. cerevisiae. Because seemingly it's expected that those proteins which have a high concentration of the target 
codon at the beginning are subject to ribosomal queueing, which will prevent more ribosomes from joining and will 
strongly decrease the expression of the protein. In order for ribosomal queueing to happen we also need the mRNA to have
high initiation rates, so I should also keep that in mind. 
"""




###############################################################
# 0. Set parameters
###############################################################
initial_codons_num = 30
working_from = "charite"
if working_from == "charite":
    base_dir = "C:\\MyStuff\\"
elif working_from == "home":
    base_dir = "\\home\\alvaro\\MyStuff\\"




###############################################################
# 1. Load data
###############################################################
# 1.1. Unique target codons
path = os.path.join(base_dir, "tRNA_KOs\\Data\\check_codon_enrichment_in_protein_sequences\\check_in_5_prime_end\\unique_target_codons.json")
with open(path, encoding='utf-8') as f:
    unique_target_codons = json.load(f)


# 1.2. SGD reference coding regions
path = os.path.join(base_dir, "tRNA_KOs\\Data\\databases\\sgdref_cds.fa")
with open(path) as handle:
    ref_dict = {}
    all_prot_ids = []
    repeated = {}

    # In each file, iterate over the proteins
    for seq_id, seq in SimpleFastaParser(handle):
        # For each protein, want to get a list of its codons as the entry in the dictionary
        strain_list = []

        # Get what is going to be the protein ID. Also append it to the "repeated" list if we´ve seen that ID before in this file
        limit = seq_id.rfind("|")
        last_chunk = seq_id[limit + 1:len(seq_id)]
        first_chunk = seq_id[0:(len(seq_id) - limit - 1)]

        if last_chunk != first_chunk:
            prot_id = last_chunk
            if prot_id in all_prot_ids:
                if prot_id not in list(repeated.keys()):
                    repeated[prot_id] = [full_protein_seqs_strain[prot_id], seq]
                else:
                    repeated[prot_id].append(seq)
            else:
                # Get the codons for this protein
                strain_list = textwrap.wrap(seq, 3)

                # Add this entry to the output dictionary
                ref_dict[prot_id] = strain_list


# 1.3. For proteins with multiple possible names separated by /, keep only the first one to make it cleaner
for protein_name in list(ref_dict.keys()):
    if "/" in protein_name:
        new_protein_name = protein_name[:protein_name.find("/")]
        ref_dict[new_protein_name] = ref_dict.pop(protein_name)



###############################################################
# 2. Iterate over unique target codons and come up with the proportion of each of them in each protein
###############################################################
"""
Iterate first over target codons, then for each of them over all proteins and get the proportion of that codon among the
first initial_codons_num codons - put this into a dataframe where target codons are columns and proteins are rows
"""

data_for_df = {"Protein": list(ref_dict.keys())}

for target_codon in unique_target_codons:                                                                   # Iterate over target codons
    proportions_for_this_codon = []                                                                         # Create list where I'll save the proportion of target codon in first codons for each protein (for this codon)
    for protein in list(ref_dict.keys()):                                                                   # Iterate over proteins
        first_codons = ref_dict[protein][:initial_codons_num]                                               # Grab the first few codons in this protein (exact number defined at the beginning of the script)
        proportions_for_this_codon.append(first_codons.count(target_codon)/initial_codons_num)              # Get the ratio (appearances of the target codon / codons I'm looking at) for this protein and save it
    data_for_df[target_codon] = proportions_for_this_codon

df = pd.DataFrame(data_for_df)


# Save resulting dataframe as .csv to read it in R
filename = f"codon_proportions_first_{initial_codons_num}_codons.csv"
output_path = os.path.join(base_dir, "tRNA_KOs", "Data", "check_codon_enrichment_in_protein_sequences", "check_in_5_prime_end", filename)

df.to_csv(output_path, index = False)




