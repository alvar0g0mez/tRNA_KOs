import json
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
import textwrap
import matplotlib
matplotlib.use('TkAgg')

from helper_functions import *




"""
I have an R list containing vectors, each of them with the proteins that were found to be DE in each of the KO strains. 
The goal here is to, for each KO strain:
    - Go through the DE proteins, and in the FASTA file with all protein sequences, check how many times the codon in the KOd tRNA was present in the protein sequence
    - Do the same for the proteins that were not DE
    - Turn these into ratios! This is, the % of codons in the protein that are the codon of interest - this corrects for protein length
    - Then we can test these vectors against each other, to see if the mean of the first one is significantly higher than the second one - do this in R
    - Hence proving that the KO of a tRNA causes differential expression of, particularly, proteins containing more of the corresponding codon
"""


###############################################################
# 0. Set parameters
###############################################################
alpha = 0.01
lfc_threshold = 1.5
working_from = "charite"
if working_from == "charite":
    base_dir = "C:/MyStuff/tRNA_KOs/"
elif working_from == "home":
    base_dir = "/home/alvaro/MyStuff/tRNA_KOs/"




###############################################################
# 1. Load data
###############################################################
# 1.1. DE proteins at alpha = 0.01
json_file = os.path.join(base_dir, "Data\\Other\\enrichment_analysis\\de_proteins_list_001.json")
with open(json_file, 'r') as f:
    de_proteins_dict_001 = json.load(f)
del de_proteins_dict_001["WT"]

# 1.2. SGD reference coding regions
path = os.path.join(base_dir, "Data\\sgdref_cds.fa")
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




###############################################################
# 2. Iterate over KO strains and come up with vectors with counts of codon appearances for each protein
###############################################################
'''I am getting 2 versions of the dictionary because one contains the names of the proteins, which will make it easier in case I want to check which is the protein which has a certain count value,
but the other one just has count vectors, which will be much easier to plot from'''

final_count_dict_with_protein_names = {}
final_count_dict_just_counts = {}
for strain in list(de_proteins_dict_001.keys()):
    strain_dict_DE = {}
    strain_dict_non_DE = {}
    strain_dict_just_counts = {"DE": [],
                              "non_DE": []}
    # Come up with the codon
    anticodon =  strain[strain.rfind("(")+1:strain.rfind(")")]
    codon = anticodon_to_codon(anticodon)

    # Iterate over all proteins in reference file and get the count of the appearances of the codon in it
    for protein in list(ref_dict.keys()):
        ratio = ref_dict[protein].count(codon)/len(ref_dict[protein])

        # Figure out where to save this count: is this protein DE or not? - a bit annoying bc of "/" in some of the reference protein names
        if "/" in protein:
            proteins = protein.split("/")
            if any(element in de_proteins_dict_001[strain] for element in proteins):
                protein_new = next((element for element in de_proteins_dict_001[strain] if element in proteins), None)
                strain_dict_DE[protein_new] = ratio
                strain_dict_just_counts["DE"].append(ratio)
            else:
                strain_dict_non_DE[protein] = ratio
                strain_dict_just_counts["non_DE"].append(ratio)

        else:
            if protein in de_proteins_dict_001[strain]:
                strain_dict_DE[protein] = ratio
                strain_dict_just_counts["DE"].append(ratio)
            else:
                strain_dict_non_DE[protein] = ratio
                strain_dict_just_counts["non_DE"].append(ratio)

    # Add the strain dictionaries to the full dictionary
    final_count_dict_with_protein_names[strain] = {"DE": strain_dict_DE,
                                                   "non_DE": strain_dict_non_DE}
    final_count_dict_just_counts[strain] = strain_dict_just_counts




###############################################################
# 3. Save these as JSON files to analyze in R
###############################################################
## 3.1. With protein names
json_file = os.path.join(base_dir, "Data\\Other\\check_codon_enrichment_in_protein_sequences\\codon_counts_with_protein_names.json")
with open(json_file, 'w') as fp:
    json.dump(final_count_dict_with_protein_names, fp)

## 3.2. Without protein names
json_file = os.path.join(base_dir, "Data\\Other\\check_codon_enrichment_in_protein_sequences\\codon_counts.json")
with open(json_file, 'w') as fp:
    json.dump(final_count_dict_just_counts, fp)



































