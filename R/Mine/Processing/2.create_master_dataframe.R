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
library(tidyr)
library(xlsx)
library(stringr)



#-----------------------------------------------------------------------------------------

# Set directories to be used
working_from = "charite"

if (working_from == "home") {
  base_dir = "/home/alvaro/MyStuff/"
} else
  if (working_from == "charite") {
    base_dir = "C:/MyStuff/"
  }

#-----------------------------------------------------------------------------------------




# 1. Dealing with the database from GtRNAdb

# Load data
db <- read.xlsx(paste(base_dir, "tRNA_KOs/Data/databases/GtRNAdb/GtRNAdb_gene_list.xlsx", sep=""), 1)
phenotypic_data <- as.data.frame(fread(paste(base_dir, "tRNA_KOs/Data/Articles/bloom_ackermann_2014/phenotypic_results_full.tsv", sep="")))
source(paste(base_dir, "tRNA_KOs/Code/R/Mine/0.general_use_functions.R", sep=""))


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
aas <- fread(paste(base_dir, "tRNA_KOs/Data/databases/GtRNAdb/amino_acids.csv",sep=""))
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
phenotypic_results <- fread(paste(base_dir, "tRNA_KOs/Data/Articles/bloom_ackermann_2014/phenotypic_results_2014.tsv", sep=""))

# 2.1. Add columns
master_dataset <- db %>%
  mutate(KOd = case_when(GtRNADB_name %in% phenotypic_results$GtRNAdb_name ~ "Yes",
                         TRUE ~ "No")) %>%
  dplyr::rename(GtRNAdb_name = GtRNADB_name) %>%
  left_join(phenotypic_results[, c("GtRNAdb_name", "Strain.Name")], by="GtRNAdb_name")


# 2.2. Add a column with the information on whether each strain grew better or worse than the WT in DTT
# I am using 2 and -2 as a threshold here because it's what they effectively use in the article, although they say they use the 
# WT mean +- the WT SD. The thing is that I still haven't been able to figure out how they obtained these values, how the 
# calculated them or normalized them or anything like that, so for now I'm going with this. 
## Need to get rid of the lethal strains before I select here
phenotypic_data <- phenotypic_data %>%
  filter(lethal != "Yes",
         !(is.na(GR_DTT_2014)),
         !is.na(GY_DTT_2014))

## Actually perform the selection
strains_better_in_DTT_GR <- phenotypic_data$Strain.Name[phenotypic_data$GR_DTT_2014 > 2]
strains_better_in_DTT_GY <- phenotypic_data$Strain.Name[phenotypic_data$GY_DTT_2014 > 2]
strains_worse_in_DTT_GR <- phenotypic_data$Strain.Name[phenotypic_data$GR_DTT_2014 < -2]
strains_worse_in_DTT_GY <- phenotypic_data$Strain.Name[phenotypic_data$GY_DTT_2014 < -2]

master_dataset <- master_dataset %>%
  dplyr::mutate(GR_in_DTT_compared_to_WT = case_when(master_dataset$Strain.Name %in% strains_better_in_DTT_GR ~ "Better",
                                                     master_dataset$Strain.Name %in% strains_worse_in_DTT_GR ~ "Worse",
                                                     TRUE ~ "Un-affected"),
                GY_in_DTT_compared_to_WT = case_when(master_dataset$Strain.Name %in% strains_better_in_DTT_GY ~ "Better",
                                                     master_dataset$Strain.Name %in% strains_worse_in_DTT_GY ~ "Worse",
                                                     TRUE ~ "Un-affected"))

master_dataset <- master_dataset %>%
  dplyr::relocate(GR_in_DTT_compared_to_WT, .after = Amino_acid_1_letter) %>%
  dplyr::relocate(GY_in_DTT_compared_to_WT, .after = GR_in_DTT_compared_to_WT)


## Remove unnecessary variables
rm(phenotypic_results, db)



#-----------------------------------------------------------------------------------------



# 3. Add genetic and mature tRNA sequences from the FASTA files
trna_seqs <- read.csv(paste(base_dir, "tRNA_KOs/Data/databases/GtRNAdb/gene_and_mature_tRNA_seqs.csv", sep="")) %>%
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

## Column saying if there is a U or A in position 34 (the one that binds the last nucleotide in the codon)
master_dataset <- master_dataset %>%
  mutate(U_34 = case_when(grepl("\\(U", Strain.Name) ~ T,
                          TRUE ~ F),
         A_34 = case_when(grepl("\\(A", Strain.Name) ~ T,
                          TRUE ~ F))

## Create a column for each Nt in the tRNA structure, making the value the letter for the Nt in that position in each tRNA
## HAVEN'T TESTED THIS, I THINK IT SHOULD BE WORKING FINE BUT IDK
master_dataset <- master_dataset %>%
  dplyr::mutate(chars = strsplit(mature_sequence, "")) %>%
  unnest_wider(chars, names_sep = "_") %>%
  rename_with(~ paste0("Nt_at_", seq_along(.)), starts_with("chars_"))


## Change the strain names from using parenthesis to using dots
master_dataset <- master_dataset %>%
  dplyr::mutate(Strain.Name = str_replace(Strain.Name, "\\(", ".")) %>%
  dplyr::mutate(Strain.Name = str_replace(Strain.Name, "\\)", "."))


#-----------------------------------------------------------------------------------------



# 5. That's it for now, just save this resulting master dataframe
fwrite(master_dataset, paste(base_dir, "tRNA_KOs/Data/basic/master_tRNA_dataset.csv", sep=""))





# ----> This continues in "C:/MyStuff/tRNA_KOs/Code/Python/complete_master_trna_dataset.py"














