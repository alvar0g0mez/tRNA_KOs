
# The goal of this file is to process the Excel files obtained from the original
# article's SI, one of them originally named: "List of all tRNA deletion strains
# in the library and their respective phenotypes across conditions", so as to 
# make it easier to work with in R, and turn it into a .tsv. 
# The other one, originally named "4.Microarray Fold change measurements for 
# selected tRNA deletion strains.", just to change colnames and turn it into a 
#.tsv. 



# FIRST FILE
# Note that I already modified some things by hand in the first Excel file prior
# to this!
#   - Add _ instead of space to the column names
#   - Start colnames with small letters except "NaCl", "DTT"...
#   - Leave growth rate metrics with their original name, and change growth 
#     yield ones from "(Y)" to "_GY"


# Next things to do here: 
#   - Create a "lethal", a "MC/SC", and a "UCU family" columns, from the 
#     "comments" one, and get rid of the latter. 
#   - Change the arabic numerals in the GtRNADB_name to roman ones to match how
#     they are in the databases online
#   - Write out as "phenotypic_results.tsv"




# Packages
library(dplyr)
library(data.table)




# Load first file
og_results <- read_xlsx("C:/MyStuff/tRNAs/Articles/Bloom-Ackermann et al., 2014/Supplementary Information/phenotypic_results.xlsx", sheet = 1)


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
  select(-comments)


## Change the arabic numerals in the GtRNADB_name from the phenotypic results to roman numerals, so they match those from the database
### Create empty vector that will be the new column
GtRNADB_names <- c()

### Create the new version of each name and add them to the vector above
for (i in 1:nrow(og_results)) {
  og_name <- og_results$GtRNADB_name[i]
  roman_num <- str_extract(og_name, "(?<=chr)([^\\.]+)")
  part_to_be_replaced <- sub("\\..*", "", og_name)
  replacement <- paste("chr", as.character(as.roman(as.integer(roman_num))), sep="")
  new_name <- gsub(part_to_be_replaced, replacement, og_name)
  GtRNADB_names <- c(GtRNADB_names, new_name)
}

### Substitute the column in the dataframe
og_results$GtRNADB_name <- GtRNADB_names



## Write final version
fwrite(og_results, "C:/MyStuff/tRNAs/data/bloom_ackermann_2014/phenotypic_results.tsv")







# SECOND FILE
# Load file
microarray_data <- read_xls("C:/MyStuff/tRNAs/Original article/Supplementary Information/4.Microarray Fold change measurements for selected tRNA deletion strains.xls")

# Get new colnames
colnames(microarray_data)[colnames(microarray_data) == "gene name"] <- "gene_names"

# One of the modifications to the colnames is that I turn "tL.GAG.G1" to "tL.GAG.G", since this is the only one we have in our proteomics data
colnames(microarray_data)[colnames(microarray_data) == "tL.GAG.G1"] <- "tL.GAG.G"

# Write as a .tsv to the appropriate folder
fwrite(microarray_data, "C:/MyStuff/tRNAs/Data/microarray_fold_change_data.tsv")























