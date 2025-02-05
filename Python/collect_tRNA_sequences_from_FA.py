import re

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser



"""
The plan here is to, for each of the tRNA genes, grab its sequence and, from the header, come up with an ID of the same
shape as those I have in the R files, so that I can match them and add the sequence to those files. I will do this for 
both the DNA and mature sequences. 
"""



#-----------------------------------------------------------------------------------------------------------------------

# 1. For the DNA (non-mature) sequences
path = "C:\\MyStuff\\tRNAs\\Data\\GtRNAdb\\sacCer3-tRNAs.fa"
with open(path) as handle:
    # Create 2 empty lists to which to add IDs and sequences and then create a dataframe from them
    IDs_gene_seqs = []
    seqs_gene = []

    # In each file, iterate over the different tRNA genes
    for seq_id, seq in SimpleFastaParser(handle):
        # Extract the different parts that we are going to need to come up with the ID and put them together
        pattern = r'tRNAscan-SE ID:\s*([^)]*)'
        first_part = re.findall(pattern, seq_id)[0]

        pattern = r'trna\s*\d*\) ([^ (]*)'
        second_part = re.findall(pattern, seq_id)[0]

        pattern = rf'{re.escape(second_part)} \((\w+)\)'
        third_part = re.findall(pattern, seq_id)[0]

        ID = first_part+"-"+second_part+third_part
        IDs_gene_seqs.append(ID)

        # Add the sequence to the corresponding list
        seqs_gene.append(seq)

# Create a dictionary and a dataframe from it
my_dict_gene_seq = {"GtRNADB_name": IDs_gene_seqs, "DNA_sequence": seqs_gene}
df_gene_seq = pd.DataFrame(my_dict_gene_seq)




#-----------------------------------------------------------------------------------------------------------------------

# 2. For the mature sequences
path = "C:\\MyStuff\\tRNAs\\Data\\GtRNAdb\\sacCer3-mature-tRNAs.fa"
with open(path) as handle:
    # Create 2 empty lists to which to add IDs and sequences and then create a dataframe from them
    IDs_mature_seqs = []
    seqs_mature = []

    # In each file, iterate over the different tRNA genes
    for seq_id, seq in SimpleFastaParser(handle):
        # Extract the different parts that we are going to need to come up with the ID and put them together
        pattern = r'tRNAscan-SE ID:\s*([^)]*)'
        first_part = re.findall(pattern, seq_id)[0]

        pattern = r'trna\s*\d*\) ([^ (]*)'
        second_part = re.findall(pattern, seq_id)[0]

        pattern = rf'{re.escape(second_part)} \((\w+)\)'
        third_part = re.findall(pattern, seq_id)[0]

        ID = first_part+"-"+second_part+third_part
        IDs_mature_seqs.append(ID)

        # Add the sequence to the corresponding list
        seqs_mature.append(seq)

# Create a dictionary and a dataframe from it
my_dict_mature_seq = {"GtRNADB_name": IDs_mature_seqs, "mature_sequence": seqs_mature}
df_mature_seq = pd.DataFrame(my_dict_mature_seq)



#-----------------------------------------------------------------------------------------------------------------------


# 3. Final checks and export the dataframe for use in R
# Compare the gene IDs from the 2 files and see if we have exactly the same ones
IDs_gene_seqs == IDs_mature_seqs

# They are exactly the same, so let's join the 2 datasets based on the gene IDs
df_merged = df_gene_seq.merge(df_mature_seq, on="GtRNADB_name")

# Check if the mature sequences are just the DNA sequences changing T for U - THEY ARE NOT! AT LEAST THEY ARE NOT EXACTLY EQUAL
matured_sequences = [sequence.replace("T", "U") for sequence in seqs_gene]
matured_sequences[0] == seqs_mature[0]

# Save to .csv
df_merged.to_csv("C:\\MyStuff\\tRNAs\\Data\\GtRNAdb\\gene_and_mature_tRNA_seqs.csv", sep=",", index=False, header=True)




































