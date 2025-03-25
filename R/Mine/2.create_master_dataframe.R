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
#   - Output it as a .csv and delete the .xlsx original file


# Load libraries
library(data.table)
library(dplyr)
library(xlsx)
library(stringr)



#-----------------------------------------------------------------------------------------

# Set directories to be used
working_from = "home"

if (working_from == "home") {
  base_dir = "/home/alvaro/MyStuff/tRNA_KOs/"
} else
  if (working_from == "charite") {
    base_dir = "C:/MyStuff/tRNA_KOs/"
  }

#-----------------------------------------------------------------------------------------




# 1. Dealing with the database from GtRNAdb

# Load data
db <- read.xlsx(paste(base_dir, "Data/Other/GtRNAdb/GtRNAdb_gene_list.xlsx", sep=""), 1)


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
aas_temp <- aas %>% dplyr::select(Amino_acid_3_letter, Amino_acid_1_letter)
colnames(aas_temp) <- c("Isotype_from_anticodon", "Isotype_from_anticodon_1_letter")
db <- left_join(db, aas_temp, by = "Isotype_from_anticodon")

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
    new_name <- paste(og_name, "-", db$Isotype_from_anticodon[i], db$Anticodon[i], sep = "")
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
phenotypic_results <- fread(paste(base_dir, "Data/Other/Articles/bloom_ackermann_2014/phenotypic_results.tsv", sep=""))

# 2.1. Add columns
master_dataset <- db %>%
  mutate(KOd = case_when(GtRNADB_name %in% phenotypic_results$GtRNADB_name ~ "Yes",
                         TRUE ~ "No")) %>%
  left_join(phenotypic_results[, c("GtRNADB_name", "gene_name")], by="GtRNADB_name")


## Remove unnecessary variables
rm(phenotypic_results, db)



#-----------------------------------------------------------------------------------------



# 3. Add genetic and mature tRNA sequences from the FASTA files
trna_seqs <- read.csv(paste(base_dir, "Data/Other/GtRNAdb/gene_and_mature_tRNA_seqs.csv", sep=""))
master_dataset <- left_join(master_dataset, trna_seqs, by = "GtRNADB_name")

## Remove unnecessary variables
rm(trna_seqs)



#-----------------------------------------------------------------------------------------



# 4. Add things I forgot to add and realized later

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
  group_by(Isotype_from_anticodon) %>%
  mutate(Number_of_tRNAs_loading_this_aa = n())

## Add a column which contains the number of tRNA genes that could and couldn't be KOd which load a certain amino acid
master_dataset <- master_dataset %>%
  group_by(Isotype_from_anticodon) %>%
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
         Nt_at_32 = substr(mature_sequence, 32, 32),
         Nt_at_38 = substr(mature_sequence, 38, 38))



#-----------------------------------------------------------------------------------------



# 5. That's it for now, just save this resulting master dataframe
fwrite(master_dataset, paste(base_dir, "Data/Other/GtRNAdb/master_tRNA_dataset.csv", sep=""))





# ----> This continues in "C:/MyStuff/tRNAs/Scripts/Python/complete_master_trna_dataset.py"












































