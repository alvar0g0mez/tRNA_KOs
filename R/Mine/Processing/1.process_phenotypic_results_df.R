# PART 1
# The goal of this file is to process the Excel files obtained from the original
# article's SI, one of them originally named: "List of all tRNA deletion strains
# in the library and their respective phenotypes across conditions", so as to 
# make it easier to work with in R, and turn it into a .tsv. 
# The other one, originally named "4.Microarray Fold change measurements for 
# selected tRNA deletion strains.", just to change colnames and turn it into a 
#.tsv. 


# PART 2
# Added on 30.04.2025 - I am also processing here the GRs and GYs in SDC and YPD
# that came with the library when we received it in 2020 - I am actually going to
# save both these and the original phenotypic resutls to the same file, I think 
# that makes sense






################################################################################
#################################### PART 1 ####################################
################################################################################


# FIRST FILE
# Things to do here: 
#   - Fix colnames
#   - Create a "lethal", a "MC/SC", and a "UCU family" columns, from the 
#     "comments" one, and get rid of the latter. 
#   - Change the arabic numerals in the GtRNADB_name to roman ones to match how
#     they are in the databases online
#   - Write out as "phenotypic_results.tsv"




# Packages
library(dplyr)
library(data.table)
library(xlsx)

# Set directories to be used
working_from = "charite"

if (working_from == "home") {
  base_dir = "/home/alvaro/MyStuff/"
} else
  if (working_from == "charite") {
    base_dir = "C:/MyStuff/"
  }


# Load first file
og_results <- read.xlsx(paste(base_dir, "tRNA_KOs/Data/Articles/bloom_ackermann_2014/phenotypic_results.xlsx", sep=""), sheetIndex = 1)


# Fix colnames
colnames(og_results) <- c("Strain.Name", "GtRNAdb_name", "anticodon", "Amino_acid_1_letter", "codon", "chromosome", "start", "end", "strand", 
                          "family_size", "GR_YPD_2014", "GR_SDC_2014", "GR_low_glucose_2014", "GR_galactose_2014", "GR_NaCl_2014", "GR_DTT_2014", 
                          "GY_YPD_2014", "GY_SDC_2014", "GY_low_glucose_2014", "GY_galactose_2014", "GY_NaCl_2014", "GY_DTT_2014", "comments")


## Create "lethal" and "MC/SC" columns and remove "comments" column, and change 
### () to . in gene_name column
og_results <- og_results %>% 
  mutate(lethal = case_when(str_detect(comments, "lethal") ~ "Yes",
                            TRUE ~ "No"),
         MC_or_SC = case_when(str_detect(comments, "part of the MC") ~ "MC",
                              str_detect(comments, "part of the SC") ~ "SC",
                              TRUE ~ "None"),
         UCU_family = case_when(str_detect(comments, "major of the UCU family") ~ "Major",
                                str_detect(comments, "minor of the UCU family") ~ "Minor",
                                TRUE ~ "No")) %>%
  dplyr::select(-comments)


## Change the arabic numerals in the GtRNADB_name from the phenotypic results to roman numerals, so they match those from the database
### Create empty vector that will be the new column
GtRNADB_names <- c()

### Create the new version of each name and add them to the vector above
for (i in 1:nrow(og_results)) {
  og_name <- og_results$GtRNAdb_name[i]
  roman_num <- str_extract(og_name, "(?<=chr)([^\\.]+)")
  part_to_be_replaced <- sub("\\..*", "", og_name)
  replacement <- paste("chr", as.character(as.roman(as.integer(roman_num))), sep="")
  new_name <- gsub(part_to_be_replaced, replacement, og_name)
  GtRNADB_names <- c(GtRNADB_names, new_name)
}

### Substitute the column in the dataframe
og_results$GtRNAdb_name <- GtRNADB_names

# Change the format tA(AGC)J to the format tA.AGC.J
og_results <- og_results %>%
  dplyr::mutate(Strain.Name = str_replace(Strain.Name, "\\(", ".")) %>%
  dplyr::mutate(Strain.Name = str_replace(Strain.Name, "\\)", "."))


## Write final version
fwrite(og_results, paste(base_dir, "tRNA_KOs/Data/Articles/bloom_ackermann_2014/phenotypic_results_2014.tsv", sep=""))





# SECOND FILE
# Load file
microarray_data <- read.xlsx("C:/MyStuff/tRNA_KOs/Articles/Bloom-Ackermann et al., 2014/Supplementary Information/4.Microarray Fold change measurements for selected tRNA deletion strains.xls", sheetIndex = 1)

# Get new colnames
colnames(microarray_data)[colnames(microarray_data) == "gene.name"] <- "Gene.secondaryIdentifier"

# One of the modifications to the colnames is that I turn "tL.GAG.G1" to "tL.GAG.G", since this is the only one we have in our proteomics data
colnames(microarray_data)[colnames(microarray_data) == "tL.GAG.G1"] <- "tL.GAG.G"

# Write as a .tsv to the appropriate folder
fwrite(microarray_data, paste(base_dir, "tRNA_KOs/Data/Articles/bloom_ackermann_2014/microarray_fold_change_data.tsv", sep=""))










################################################################################
#################################### PART 2 ####################################
################################################################################

# Packages
library(dplyr)
library(data.table)
library(xlsx)

# Load data
growth_2020_1 <- read.xlsx("S:/AG/AG-CF-HTMS/AG-Ralser-Share/30-0092_AndreaLehmann-AlternativeAAUsage-tRNA/01_ProjectManagement/05_Samples/01_SampleMetadata/20200302_tRNA deletion library_96plates_Orna_rearranged.xlsx", sheetName = "Plate 1")
growth_2020_2 <- read.xlsx("S:/AG/AG-CF-HTMS/AG-Ralser-Share/30-0092_AndreaLehmann-AlternativeAAUsage-tRNA/01_ProjectManagement/05_Samples/01_SampleMetadata/20200302_tRNA deletion library_96plates_Orna_rearranged.xlsx", sheetName = "Plate 2")
growth_2020_3 <- read.xlsx("S:/AG/AG-CF-HTMS/AG-Ralser-Share/30-0092_AndreaLehmann-AlternativeAAUsage-tRNA/01_ProjectManagement/05_Samples/01_SampleMetadata/20200302_tRNA deletion library_96plates_Orna_rearranged.xlsx", sheetName = "Plate 3")
og_results <- as.data.frame(fread(paste(base_dir, "tRNA_KOs/Data/Articles/bloom_ackermann_2014/phenotypic_results_2014.tsv", sep="")))


# Fix colnames
colnames(growth_2020_1) <- c("Analysis.Plate.96", "Position.Bad.Format", "Strain.ID", "Strain.Name", "anticodon", "Amino_acid_1_letter", 
                                 "GR_YPD_2020", "GY_YPD_2020", "GR_SDC_2020", "GY_SDC_2020")
colnames(growth_2020_2) <- c("Analysis.Plate.96", "Position.Bad.Format", "Strain.ID", "Strain.Name", "anticodon", "Amino_acid_1_letter", 
                                 "GR_YPD_2020", "GY_YPD_2020", "GR_SDC_2020", "GY_SDC_2020")
colnames(growth_2020_3) <- c("Analysis.Plate.96", "Position.Bad.Format", "Strain.ID", "Strain.Name", "anticodon", "Amino_acid_1_letter", 
                                 "GR_YPD_2020", "GY_YPD_2020", "GR_SDC_2020", "GY_SDC_2020")

# Join them into a single dataframe
growth_2020 <- rbind(rbind(growth_2020_1, growth_2020_2), growth_2020_3)

# Change the format tA(AGC)J to the format tA.AGC.J
growth_2020 <- growth_2020 %>%
  dplyr::mutate(Strain.Name = str_replace(Strain.Name, "\\(", ".")) %>%
  dplyr::mutate(Strain.Name = str_replace(Strain.Name, "\\)", "."))

# Save this dataframe as is just in case
fwrite(growth_2020, file = paste(base_dir, "tRNA_KOs/Data/Articles/bloom_ackermann_2014/phenotypic_results_2020.tsv", sep=""))

# Join this to the original phenotypic results dataframe
temp_growth_2020 <- growth_2020 %>%
  dplyr::select(Strain.Name, starts_with("GR"), starts_with("GY")) %>%
  na.omit()
og_results <- left_join(og_results, temp_growth_2020, by = "Strain.Name")

# Save the dataframe with both phenotypic datas, 2014 and 2020
fwrite(og_results, file = paste(base_dir, "tRNA_KOs/Data/Articles/bloom_ackermann_2014/phenotypic_results_full.tsv", sep=""))




