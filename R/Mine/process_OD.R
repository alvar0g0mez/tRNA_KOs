# I did the basic processing in the Excel sheets that I got from Andrea already, since due to the original format that was way easier 
# I am only adding some more stuff here:
#   - Which wells contained samples and which didn't (based on which ones appear in my sample layout and which don't)



# Libraries
library(dplyr)
library(tidyr)
library(xlsx)




# Set directories to be used
working_from = "charite"

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
processed_od <- read.xlsx(paste(base_dir, "Data/OD_processed.xlsx", sep=""), sheetIndex = 1)




# Process OD data - check in the sample layout which wells are there and which ones are not, add a column with this information
replicates <- c(1,2,3)
out <- data.frame(matrix(nrow = 0, ncol = ncol(processed_od)))

for (i in replicates) {
  temp_sample_layout <- sample_layout %>%
    filter(Replicate == i & Analysis.Plate.96 == i) %>%
    mutate(Position.Temp = paste(Analysis.Row.96, Analysis.Column.96, sep="")) 
  temp_od <- processed_od %>%
    filter(Plate_96 == i)
  
  temp_od <- temp_od %>%
    mutate(Empty_well = case_when(Position %in% temp_sample_layout$Position.Temp ~ "No",
                                  TRUE ~ "Yes"))
  out <- rbind(out, temp_od)
}


# Save the final version with the new column
fwrite(out, paste(base_dir, "Data/OD_final.csv", sep=""))


















