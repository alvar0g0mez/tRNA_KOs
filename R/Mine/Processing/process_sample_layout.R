

# This only works when run from the Charité computer, since the version of the sample_layout I am starting from is in the S drive
# This sample_layout was already prepared by Boris from previous ones, I'm just making some small modifications based on how I want to use it

# This file reads the sample layout by Boris from the S drive, modifies slightly and writes it to my directory within 30-0092




# Packages
library(dplyr)
library(data.table)
library(stringr)


# Set up
working_from = "home"

if (working_from == "home") {
  base_dir = "/home/alvaro/MyStuff/"
} else
  if (working_from == "charite") {
    base_dir = "C:/MyStuff/"
  }


# Load data
if (working_from == "home") {
  sample_layout <- as.data.frame(fread("/home/alvaro/MyStuff/tRNA_KOs/Data/Other/proteomics_data/sample_layout_Boris.tsv"))
} else {
  sample_layout <- as.data.frame(fread("S:/AG/AG-CF-HTMS/AG-Ralser-Share/30-0092_AndreaLehmann-AlternativeAAUsage-tRNA/05_DataAnalysis/11_Preprocessing_Boris/AlternativeAAUsage-tRNA/AlternativeAAUsage-tRNA_peptidecentric_PrecursorQuantity_filename_annotations.tsv"))
}

master_dataset <- as.data.frame(fread(paste(base_dir, "tRNA_KOs/Data/Other/GtRNAdb/master_tRNA_dataset.csv", sep="")))





# Turn all dashes to underscores
sample_layout <- sample_layout %>%
  mutate(Sample.ID = str_replace_all(Sample.ID, "-", "_"),
         Sample.ID.unique = str_replace_all(Sample.ID.unique, "-", "_"))

# Grab the position within plate from the Plate.Position column, so that I can use it to obtain plate images with platetools
sample_layout <- sample_layout %>%
  mutate(Position.Within.Plate.384 = case_when(nchar(as.character(Analysis.Column.384)) == 2 ~ paste(Analysis.Row.384, Analysis.Column.384, sep=""),
                                               nchar(as.character(Analysis.Column.384)) == 1 ~ paste(Analysis.Row.384, 0, Analysis.Column.384, sep="")),
         Position.Within.Plate.96 = case_when(nchar(as.character(Analysis.Column.96)) == 2 ~ paste(Analysis.Row.96, Analysis.Column.96, sep=""),
                                              nchar(as.character(Analysis.Column.96)) == 1 ~ paste(Analysis.Row.96, 0, Analysis.Column.96, sep="")))

# Create a column that can match the colnames of the proteomics data as I get it from Boris
sample_layout <- sample_layout %>%
  mutate(raw_proteomics_colnames = case_when(Strain.Name == "WT" ~ paste(str_replace_all(Sample.ID, "_", "."), ".0", Replicate, sep=""),
                                             Strain.Name == "QC" ~ str_replace_all(Sample.ID.unique, "_", "."),
                                             TRUE ~ paste("X", Strain.ID, ".0", Replicate, sep="")))


# Create a column with the column IDs of the shape I want to be working with for my proteomics data
sample_layout <- sample_layout %>%
  mutate(final_proteomics_colnames = case_when(Strain.Name == "WT" ~ Sample.ID.unique,
                                              Strain.Name == "QC" ~ Sample.ID.unique,
                                              TRUE ~ paste(Strain.Name, "_0", Replicate, sep="")))


# Add some columns - I used to do this in my main analysis file
## Fix the QC rows to say "QC" instead of "NA" in the columns that don't apply
sample_layout <- sample_layout %>%
  mutate(Analysis.Plate.384 = as.factor(case_when(Plate.Position == "QC" ~ "QC",
                                                  TRUE ~ as.character(Analysis.Plate.384))),
         Analysis.Row.384 = as.factor(case_when(Plate.Position == "QC" ~ "QC",
                                                TRUE ~ as.character(Analysis.Row.384))),
         Analysis.Column.384 = as.factor(case_when(Plate.Position == "QC" ~ "QC",
                                                   TRUE ~ as.character(Analysis.Column.384))),
         Analysis.Plate.96 = as.factor(case_when(Plate.Position == "QC" ~ "QC",
                                                 TRUE ~ as.character(Analysis.Plate.96))),
         Analysis.Row.96 = as.factor(case_when(Plate.Position == "QC" ~ "QC",
                                               TRUE ~ as.character(Analysis.Row.96))),
         Analysis.Column.96 = as.factor(case_when(Plate.Position == "QC" ~ "QC",
                                                  TRUE ~ as.character(Analysis.Column.96))))

## Extract the date for when each sample was run, as well as Date_Injection
sample_layout <- sample_layout %>% 
  mutate(date = str_extract(File.Name, "(?<=/30-0092/).*?(?=_Z2_KTT_)")) %>%
  mutate(Injection_Order = str_extract(File.Name, "(?<=KTT_).*?(?=_30-0092_tRNA)"), 
         Date_Injection_Order = paste(str_extract(File.Name, "(?<=/30-0092/).*?(?=_Z2_KTT_)"), str_extract(File.Name, "(?<=KTT_).*?(?=_30-0092_tRNA)"),
                                      sep="_"))

## Create a column that simply has the info of if each sample is a KO strain or a WT replicate
sample_layout <- sample_layout %>% 
  mutate(Strain.Type = case_when(Strain.Name == "WT" ~ "WT",
                                 TRUE ~ "KO"))

## Add a column with the following format: Analysis.Plate.96_Replicate
sample_layout <- sample_layout %>%
  mutate(Analysis.Plate.96_Replicate = paste(Analysis.Plate.96, Replicate, sep="_"))

## Add a column which identifies samples in Analysis.Plate.96 = 3, Replicate = 2
sample_layout <- sample_layout %>%
  mutate(Wrong_batch = case_when(Analysis.Plate.96_Replicate == "3_2" ~ "Yes",
                                 TRUE ~ "No"))

# Add a column with Strain.Name but using "." instead of "_", for makeContrasts()
sample_layout <- sample_layout %>%
  mutate(Strain.Name.Dots.Delete = str_replace_all(Strain.Name, "\\(", "."),
         Strain.Name.Dots = str_replace_all(Strain.Name.Dots.Delete, "\\)", ".")) %>%
  dplyr::select(-Strain.Name.Dots.Delete)




# Merge with master dataframe - I used to do this in the main tRNA script
sample_layout <- sample_layout %>% 
  dplyr::select(-Anticodon) %>%
  left_join(master_dataset, by = c("Strain.Name"))








# Save final version
if (working_from == "home") {
  fwrite(sample_layout, "/home/alvaro/MyStuff/tRNA_KOs/Data/Other/proteomics_data/sample_layout_alvaro.tsv")
} else {
  fwrite(sample_layout, "S:/AG/AG-CF-HTMS/AG-Ralser-Share/30-0092_AndreaLehmann-AlternativeAAUsage-tRNA/05_DataAnalysis/12_Analysis_Alvaro/sample_layout_alvaro.tsv")
}





































