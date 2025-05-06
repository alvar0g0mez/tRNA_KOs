# I did the basic processing in the Excel sheets that I got from Andrea already, since due to the original format that was way easier 
# I am only adding some more stuff here:
#   - Which wells contained samples and which didn't (based on which ones appear in my sample layout and which don't)



# Libraries
library(dplyr)
library(tidyr)
library(xlsx)




# Set directories to be used
working_from = "home"

if (working_from == "home") {
  base_dir = "/home/alvaro/MyStuff/tRNA_KOs/"
  
  # Load the sample layout
  sample_layout <- as.data.frame(fread("/home/alvaro/MyStuff/tRNA_KOs/Data/Other/proteomics_data/sample_layout_alvaro.tsv"))
  
} else
  if (working_from == "charite") {
    base_dir = "C:/MyStuff/tRNA_KOs/"
    
    # Load the sample layout
    sample_layout <- as.data.frame(fread("S:/AG/AG-CF-HTMS/AG-Ralser-Share/30-0092_AndreaLehmann-AlternativeAAUsage-tRNA/05_DataAnalysis/12_Analysis_Alvaro/sample_layout_alvaro.tsv"))
  }




# Load data
#processed_od <- read.xlsx(paste(base_dir, "Data/OD_processed.xlsx", sep=""), sheetIndex = 1)
processed_od <- as.data.frame(fread(paste(base_dir, "Data/OD_processed.csv", sep="")))
master_dataset <- as.data.frame(fread(paste(base_dir, "Data/Other/GtRNAdb/master_tRNA_dataset.csv", sep="")))




# Process OD data - check in the sample layout which wells are there and which ones are not, add a column with this information
replicates <- c(1,2,3)
out <- data.frame(matrix(nrow = 0, ncol = ncol(processed_od)+3))

for (i in replicates) {
  temp_sample_layout <- sample_layout %>%
    filter(Replicate == i & Analysis.Plate.96 == i) %>%
    mutate(Position = paste(Analysis.Row.96, Analysis.Column.96, sep="")) 
  temp_od <- processed_od %>%
    filter(Plate_96 == i)
  
  temp_od <- temp_od %>%
    mutate(Empty_well = case_when(Position %in% temp_sample_layout$Position ~ "No",
                                  TRUE ~ "Yes"))
  
  # Also add 3 columns, with the tRNA KO name and anticodon, based on sample layout
  temp_sample_layout <- temp_sample_layout %>%
    dplyr::select(Position, Strain.Name, Anticodon)
  temp_od <- left_join(temp_od, temp_sample_layout, by = "Position")
  
  out <- rbind(out, temp_od)
}


# Add a column also with the amino acid - need to do this from the master dataframe
temp_master <- master_dataset %>%
  dplyr::select(Strain.Name, Anticodon, Isotype_from_anticodon_1_letter) %>%
  distinct(Strain.Name, Anticodon, .keep_all = T)
out <- left_join(out, temp_master, by = c("Anticodon", "Strain.Name")) %>%
  dplyr::rename(Amino_acid_1_letter = Isotype_from_anticodon_1_letter)


# The 1-letter amino acid names come out empty for iMets, I'm going to just fill them with "M" (I have it saved in the master_dataframe which ones are iMet and which ones are just Met)
out <- out %>%
  mutate(iMet = case_when(Amino_acid_1_letter == "" ~ "Yes",
                                         TRUE ~ "No"),
         Amino_acid_1_letter = case_when(Amino_acid_1_letter == "" ~ "M",
                                         TRUE ~ Amino_acid_1_letter))


# Save the final version with the new column
fwrite(out, paste(base_dir, "Data/OD_final.csv", sep=""))


















