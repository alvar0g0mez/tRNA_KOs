# This script creates a master dataset with all information about all tRNA genes
# The 3 main sources of information (for now) are: 
#   - GtRNAdb_gene_list.xlsx --> directly from the website, I created it by copy-
#   pasting the table to a .xlsx document
#   - phenotypic_results.tsv --> generated before somewhere else (maybe by hand?),
#   slight modification of the 1st dataset in the SI from the original article
#   - amino_acids.csv --> super quick dataset I created to just have the amino 
#   acid name, 3-letter code and 1-letter code mapped to each other



# Some of the stuff that happens in this file - THIS DOES NOT COVER EVERYTHING THAT HAPPENS:
#   - Separate "Features" column into 2 columns: Intron and Mismatch, where I 
#   keep just the number of each that corresponds to each row. Then delete
#   Features column
#   - For the columns with 3-letter amino acid codes, create a corresponding 
#   column with the 1-letter code.
#   - For some reason, in GtRNAdb they use "TGC" as the anticodon name,
#   which doesn't make any sense, that is neither the codon (which would be "GCA")
#   nor the anticodon (which should be "UGC"), so I had to write a function to 
#   go through these and turn all the Ts to Us for it to make sense
#   - Output it as a .csv and delete the .xlsx original file


# Load libraries
library(data.table)
library(dplyr)
library(xlsx)
library(stringr)



#-----------------------------------------------------------------------------------------

# Set directories to be used
working_from = "charite"

if (working_from == "home") {
  base_dir = "/home/alvaro/MyStuff/tRNA_KOs/"
  source("/home/alvaro/MyStuff/tRNA_KOs/Code/R/Mine/0.general_use_functions.R")
} else
  if (working_from == "charite") {
    base_dir = "C:/MyStuff/tRNA_KOs/"
    source("C:/MyStuff/tRNA_KOs/Code/R/Mine/0.general_use_functions.R")
  }

#-----------------------------------------------------------------------------------------




# 1. Dealing with the database from GtRNAdb

# Load data
db <- read.xlsx(paste(base_dir, "Data/Other/GtRNAdb/GtRNAdb_gene_list.xlsx", sep=""), 1)
source("C:/MyStuff/tRNA_KOs/Code/R/Mine/0.general_use_functions.R")


# 1.1. Come up with Intron and Mismatch columns from the Features column
## Create empty columns
db$Intron <- NA
db$Mismatch <- NA

for (i in 1:nrow(db)) {
  # If we have mismatch information from the row above
  if (grepl("mismatch", db$Features[i]) & is.na(db$Anticodon_and_isotype_model_agreement[i])) {
    db$Mismatch[i-1] <- as.numeric(regmatches(db$Features[i], regexpr("(?<=: ).*", db$Features[i], perl = TRUE)))
    #db <- db[-c(i),]
  } 
  
  # If we have mismatch information from this row
  else if (grepl("mismatch", db$Features[i])) {
    db$Mismatch[i] <- as.numeric(regmatches(db$Features[i], regexpr("(?<=: ).*", db$Features[i], perl = TRUE)))
  }
  
  # If we have intron information (only possible for this row)
  else if (grepl("intron", db$Features[i])) {
    db$Intron[i] <- as.numeric(regmatches(db$Features[i], regexpr("(?<=: ).*", db$Features[i], perl = TRUE)))
  }
}

## Remove rows that only contained mismatch information for the previous row
db <- subset(db, !is.na(Anticodon))

## Remove Features column
db <- db %>% dplyr::select(-Features)

## Turn NAs in these 2 columns into 0s
db <- db %>%
  mutate(Intron = case_when(is.na(Intron) ~ 0,
                            TRUE ~ Intron),
         Mismatch = case_when(is.na(Mismatch) ~ 0,
                              TRUE ~ Mismatch))


# 1.2. Match 3-letter codes to 1-letter codes
## Load dataframe with amino acid information
aas <- fread(paste(base_dir, "Data/Other/GtRNAdb/amino_acids.csv",sep=""))
aas <- as.data.frame(aas)

## Add new columns
aas_temp <- aas %>% 
  dplyr::select(Amino_acid_3_letter, Amino_acid_1_letter) %>%
  dplyr::rename(Isotype_from_anticodon = Amino_acid_3_letter)
db <- left_join(db, aas_temp, by = "Isotype_from_anticodon") %>%
  dplyr::rename(Amino_acid_3_letter = Isotype_from_anticodon)

colnames(aas_temp) <- c("Best_isotype_model", "Best_isotype_model_1_letter")
db <- left_join(db, aas_temp, by = "Best_isotype_model")

## Remove unnecessary variables
rm(aas, aas_temp, i)


# 1.3. Come up with an ID like the one they use in the original article, which they 
# call "GtRNADB_name", to be able to join both datasets based on it

## Create empty column
db$GtRNADB_name <- NA

## Iterate over every row and add the corresponding name to the new column
for (i in 1:nrow(db)) {
  og_name <- db$tRNAscan_SE_ID[i]
  
  # Special case - these two tRNAs whose codons are not known for some reason?
  if (grepl("Und", db$GtRNAdb_gene_symbol[i])) {
    new_name <- paste(og_name, "-Undet???", sep = "")
    db$GtRNADB_name[i] <- new_name 
  }
  
  # All other cases
  else {
    new_name <- paste(og_name, "-", db$Amino_acid_3_letter[i], db$Anticodon[i], sep = "")
    db$GtRNADB_name[i] <- new_name 
  }
}


# 1.4. Turn those "iMet" IDs into just "Met" so they can match with the phenotypic dataset, and create an 
# extra column where I keep that info of which Methionines are initiators
db <- db %>%
  mutate(iMet = case_when(grepl("iMet", GtRNADB_name) ~ "Yes",
                          TRUE ~ "No")) %>%
  mutate(GtRNADB_name = case_when(grepl("iMet", GtRNADB_name) ~ gsub("iMetCAT", "MetCAT", GtRNADB_name),
                                  TRUE ~ GtRNADB_name))



## Remove unnecessary variables
rm(og_name, new_name, i)



#-----------------------------------------------------------------------------------------



# 2. Add a KOd column which represents whether that gene could be KOd in the original study or not (so represents 
# if it is included in our data or not)
# Also, add the column from the phenotypic results dataframe "gene_name", it's the only one that will allow me to 
# match the sample_layout to the master dataframe.


# 2.0. Load phenotypic results dataset
phenotypic_results <- fread(paste(base_dir, "Data/Other/Articles/bloom_ackermann_2014/phenotypic_results_2014.tsv", sep=""))

# 2.1. Add columns
master_dataset <- db %>%
  mutate(KOd = case_when(GtRNADB_name %in% phenotypic_results$GtRNAdb_name ~ "Yes",
                         TRUE ~ "No")) %>%
  dplyr::rename(GtRNAdb_name = GtRNADB_name) %>%
  left_join(phenotypic_results[, c("GtRNAdb_name", "Strain.Name")], by="GtRNAdb_name")


## Remove unnecessary variables
rm(phenotypic_results, db)



#-----------------------------------------------------------------------------------------



# 3. Add genetic and mature tRNA sequences from the FASTA files
trna_seqs <- read.csv(paste(base_dir, "Data/Other/GtRNAdb/gene_and_mature_tRNA_seqs.csv", sep="")) %>%
  dplyr::rename(GtRNAdb_name = GtRNADB_name)
master_dataset <- left_join(master_dataset, trna_seqs, by = "GtRNAdb_name")

## Remove unnecessary variables
rm(trna_seqs)



#-----------------------------------------------------------------------------------------



# 4. Add things I forgot to add and realized later

## What we have in the Anticodon column is not exactly an Anticodon: it's the Anticodon sequence
## but where there should be a U there's a T for some reason. So first I am going to fix that,
## then create a codon column from the anticodon one
master_dataset <- master_dataset %>%
  dplyr::mutate(Anticodon = sapply(Anticodon, turn_t_to_u_in_wrong_GtRNAdb_anticodons), 
                Codon = sapply(Anticodon, anticodon_to_codon))

## Make the order of the columns a bit nicer
master_dataset <- master_dataset %>%
  relocate(Strain.Name, .after = Locus) %>%
  relocate(Anticodon, .after = Strain.Name) %>%
  relocate(Codon, .after = Anticodon) %>%
  relocate(Amino_acid_3_letter, .after = Codon) %>%
  relocate(Amino_acid_1_letter, .after = Amino_acid_3_letter)

## Create a chromosome column based on tRNAscan_SE_ID 
master_dataset <- master_dataset %>%
  dplyr::mutate(chromosome = LETTERS[as.numeric(as.roman(str_extract(tRNAscan_SE_ID, "(?<=chr)[^\\.]+(?=\\.)")))]) %>%
  relocate(chromosome, .after = Codon)

## Add columns with the length of each of the sequences
master_dataset <- master_dataset %>%
  mutate(length_DNA_seq = nchar(DNA_sequence),
         length_mature_seq = nchar(mature_sequence))

## Add column with family size
master_dataset <- master_dataset %>%
  group_by(Anticodon) %>%
  mutate(Family_size = n())

## Add a column which contains the number of tRNAs loading each amino acid (like family_size, but for amino acid loaded, not anticodon)
master_dataset <- master_dataset %>%
  group_by(Amino_acid_3_letter) %>%
  mutate(Number_of_tRNAs_loading_this_aa = n())

## Add a column which contains the number of tRNA genes that could and couldn't be KOd which load a certain amino acid
master_dataset <- master_dataset %>%
  group_by(Amino_acid_3_letter) %>%
  mutate(Number_of_not_KOd_genes_per_aa = sum(KOd == "No"),
         Number_of_KOd_genes_per_aa = sum(KOd == "Yes"))

## Add a column which contains the number of tRNA genes that could and couldn't be KOd with a certain anticodon
master_dataset <- master_dataset %>%
  group_by(Anticodon) %>%
  mutate(Number_of_not_KOd_genes_per_anticodon = sum(KOd == "No"),
         Number_of_KOd_genes_per_anticodon = sum(KOd == "Yes"))

## Add a column with percentage of KOd genes per amino acid and per anticodon
master_dataset <- master_dataset %>%
  mutate(Perc_KOd_genes_per_anticodon = Number_of_KOd_genes_per_anticodon/Family_size,
         Perc_KOd_genes_per_amino_acid = Number_of_KOd_genes_per_aa/Number_of_tRNAs_loading_this_aa)

## Column saying if there is a U in position 34 (the one that binds the last nucleotide in the codon), and a couple things more
master_dataset <- master_dataset %>%
  mutate(U_34 = case_when(grepl("\\(U", Strain.Name) ~ T,
                          TRUE ~ F),
         A_34 = case_when(grepl("\\(A", Strain.Name) ~ T,
                          TRUE ~ F),
         Nt_at_1 = substr(mature_sequence, 1, 1),
         Nt_at_2 = substr(mature_sequence, 2, 2),
         Nt_at_3 = substr(mature_sequence, 3, 3),
         Nt_at_4 = substr(mature_sequence, 4, 4),
         Nt_at_5 = substr(mature_sequence, 5, 5),
         Nt_at_6 = substr(mature_sequence, 6, 6),
         Nt_at_7 = substr(mature_sequence, 7, 7),
         Nt_at_8 = substr(mature_sequence, 8, 8),
         Nt_at_9 = substr(mature_sequence, 9, 9),
         Nt_at_10 = substr(mature_sequence, 10, 10),
         Nt_at_11 = substr(mature_sequence, 11, 11),
         Nt_at_12 = substr(mature_sequence, 12, 12),
         Nt_at_13 = substr(mature_sequence, 13, 13),
         Nt_at_14 = substr(mature_sequence, 14, 14),
         Nt_at_15 = substr(mature_sequence, 15, 15),
         Nt_at_16 = substr(mature_sequence, 16, 16),
         Nt_at_17 = substr(mature_sequence, 17, 17),
         Nt_at_18 = substr(mature_sequence, 18, 18),
         Nt_at_19 = substr(mature_sequence, 19, 19),
         Nt_at_20 = substr(mature_sequence, 20, 20),
         Nt_at_21 = substr(mature_sequence, 21, 21),
         Nt_at_22 = substr(mature_sequence, 22, 22),
         Nt_at_23 = substr(mature_sequence, 23, 23),
         Nt_at_24 = substr(mature_sequence, 24, 24),
         Nt_at_25 = substr(mature_sequence, 25, 25),
         Nt_at_26 = substr(mature_sequence, 26, 26),
         Nt_at_27 = substr(mature_sequence, 27, 27),
         Nt_at_28 = substr(mature_sequence, 28, 28),
         Nt_at_29 = substr(mature_sequence, 29, 29),
         Nt_at_30 = substr(mature_sequence, 30, 30),
         Nt_at_31 = substr(mature_sequence, 31, 31),
         Nt_at_32 = substr(mature_sequence, 32, 32),
         Nt_at_33 = substr(mature_sequence, 33, 33),
         Nt_at_34 = substr(mature_sequence, 34, 34),
         Nt_at_35 = substr(mature_sequence, 35, 35),
         Nt_at_36 = substr(mature_sequence, 36, 36),
         Nt_at_37 = substr(mature_sequence, 37, 37),
         Nt_at_38 = substr(mature_sequence, 38, 38),
         Nt_at_39 = substr(mature_sequence, 39, 39),
         Nt_at_40 = substr(mature_sequence, 40, 40),
         Nt_at_41 = substr(mature_sequence, 41, 41),
         Nt_at_42 = substr(mature_sequence, 42, 42),
         Nt_at_43 = substr(mature_sequence, 43, 43),
         Nt_at_44 = substr(mature_sequence, 44, 44),
         Nt_at_45 = substr(mature_sequence, 45, 45),
         Nt_at_46 = substr(mature_sequence, 46, 46),
         Nt_at_47 = substr(mature_sequence, 47, 47),
         Nt_at_48 = substr(mature_sequence, 48, 48),
         Nt_at_49 = substr(mature_sequence, 49, 49),
         Nt_at_50 = substr(mature_sequence, 50, 50),
         Nt_at_51 = substr(mature_sequence, 51, 51),
         Nt_at_52 = substr(mature_sequence, 52, 52),
         Nt_at_53 = substr(mature_sequence, 53, 53),
         Nt_at_54 = substr(mature_sequence, 54, 54),
         Nt_at_55 = substr(mature_sequence, 55, 55),
         Nt_at_56 = substr(mature_sequence, 56, 56),
         Nt_at_57 = substr(mature_sequence, 57, 57),
         Nt_at_58 = substr(mature_sequence, 58, 58),
         Nt_at_59 = substr(mature_sequence, 59, 59),
         Nt_at_60 = substr(mature_sequence, 60, 60),
         Nt_at_61 = substr(mature_sequence, 61, 61),
         Nt_at_62 = substr(mature_sequence, 62, 62),
         Nt_at_63 = substr(mature_sequence, 63, 63),
         Nt_at_64 = substr(mature_sequence, 64, 64),
         Nt_at_65 = substr(mature_sequence, 65, 65),
         Nt_at_66 = substr(mature_sequence, 66, 66),
         Nt_at_67 = substr(mature_sequence, 67, 67),
         Nt_at_68 = substr(mature_sequence, 68, 68),
         Nt_at_69 = substr(mature_sequence, 69, 69),
         Nt_at_70 = substr(mature_sequence, 70, 70))





#-----------------------------------------------------------------------------------------



# 5. That's it for now, just save this resulting master dataframe
fwrite(master_dataset, paste(base_dir, "Data/Other/GtRNAdb/master_tRNA_dataset.csv", sep=""))





# ----> This continues in "C:/MyStuff/tRNAs/Scripts/Python/complete_master_trna_dataset.py"














