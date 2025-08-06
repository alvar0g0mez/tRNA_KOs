
# Libraries
library(helperfunctions)
library(dplyr)
library(data.table)
library(stringr)
library(limpa)
library(jsonlite)





# 0. Set up 
## Significance level to be used for all tests and plots in this file
alpha <- 0.05

## Significance level as a plain string to use when loading or writing files
alpha_plain <- str_replace(as.character(alpha), "\\.", "")

## Minimum threshold for the average log2-expression across all samples in order to keep a protein in the data before DE analysis
mean_log2_across_all_samples_threshold <- 2

## Minimum threshold for the variance of the log2-expression across all samples in order to keep a protein in the data before DE analysis
var_across_log2_all_samples_threshold <- 0.6

## Separate significance level, the one used for the enrichment analysis
alpha_enrichment <- 0.05

## Significance level above as a plain string to use when loading or writing files
alpha_enrichment_plain <- str_replace(as.character(alpha_enrichment), "\\.", "")

## Log-fold change limit to be considered "biologically significant"
lfc_threshold <- 0.5

## Set directories to be used
working_from = "charite"

if (working_from == "home") {
  base_dir = "/home/alvaro/MyStuff/"
} else
  if (working_from == "charite") {
    base_dir = "C:/MyStuff/"
  }




# 0. Load data
## Sample layout
sample_layout <- as.data.frame(fread(paste(base_dir, "tRNA_KOs/Data/basic/sample_layout_alvaro.tsv", sep="")))

## UniProt dataset
uniprot_db <- as.data.frame(fread(paste(base_dir, "tRNA_KOs/Data/databases/uniprotkb_proteome_UP000002311_2025_05_24.tsv", sep="")))

## Messner et al., 2023
protein_abundances <- as.data.frame(fread(paste(base_dir, "tRNA_KOs/Data/5k/yeast5k_impute_wide.csv", sep=""))) %>%
  dplyr::rename(gene_names = Protein.Group)
protein_abundances <- match_systematic_and_standard_protein_names(protein_abundances,
                                                                  uniprot = uniprot_db,
                                                                  input = "uniprot",
                                                                  output = "systematic") %>%
  dplyr::select(-uniprot_IDs)
protein_abundances <- set_column_as_rownames(protein_abundances, 
                                             "systematic_SGD_IDs")
protein_abundances <- as.matrix(protein_abundances)
protein_abundances <- ifelse(protein_abundances == 0, 0, log2(protein_abundances))
my_metadata <- as.data.frame(fread(paste(base_dir, "tRNA_KOs/Data/5k/yeast5k_metadata.csv", sep="")))




# 2. Prepare to perform the DEA
## Need to change protein IDs from using - to _ or the flipping makeContrasts complains
my_metadata <- my_metadata %>%
  dplyr::mutate(ORF_contrast_friendly = case_when(grepl("-", ORF) ~ str_replace(ORF, "-", "_"),
                                                  TRUE ~ ORF))

## Make sure the samples are in exactly the same order in the metadata than in the protein abundance df
sum(colnames(protein_abundances) == my_metadata$Filename) == nrow(my_metadata)

## Get the names of all KO strains we have, and set His3 as the reference
strain_levels <- as.factor(my_metadata$ORF_contrast_friendly)
strain_levels <- relevel(strain_levels, ref = "YOR202W")

## Generate design matrix
mm <- model.matrix(~ 0 + strain_levels)
colnames(mm) <- levels(strain_levels)






# 3. Don't think I can do DE with limpa if I didn't also do protein summarization with limpa - so DE analysis with limma - from ChatGPT
## Ensure your data is a numeric matrix with proteins as rows and samples as columns
y_protein <- as.matrix(protein_abundances)
mode(y_protein) <- "numeric"

## Double-check dimensions match
stopifnot(ncol(y_protein) == nrow(my_metadata))

## Create contrast matrix - compare each strain vs. His3 (reference, YOR202W)
contrast_matrix <- makeContrasts(
  levels = colnames(mm),
  contrasts = paste(setdiff(colnames(mm), "YOR202W"), "-YOR202W", sep = "")
)

## Fit linear model
fit <- lmFit(y_protein, mm)

## Apply contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)

## Empirical Bayes moderation
fit2 <- eBayes(fit2)

## Save topTables to a da list
da <- list()

for (contrast_name in colnames(contrast_matrix)) {
  new_name <- substr(contrast_name, 1, str_locate(contrast_name, "-")-1)
  if (grepl("_", new_name)) {
    new_name <- str_replace(new_name, "_", "-")
  }
  da[[new_name]] <- topTable(fit2, coef = contrast_name, adjust.method = 'BH', number = Inf, sort.by = 'none')
  da[[new_name]]$Strain.Name <- rep(new_name, nrow(da[[new_name]]))
}

## Generate a list of DE proteins per strain, like for the tRNA_KOs
de_proteins_list <- list()
for (i in 1:length(da)) {
  temp <- da[[i]]
  temp <- na.omit(temp)
  
  # Collect names of DE proteins
  general_protein_names <- rownames(temp)[(temp$adj.P.Val < alpha) &
                                          (temp$logFC > lfc_threshold | temp$logFC < -lfc_threshold)]
  
  # Save protein names
  de_proteins_list[[names(da)[i]]] <- general_protein_names
}
output_file <- paste(base_dir, "tRNA_KOs/Data/5k/de_proteins_list_", alpha_plain, "_logFC_", as.character(lfc_threshold), ".json", sep="")
write_json(de_proteins_list, path=output_file, pretty = TRUE)







# 4. Finish processing the DE analysis results
## Format data
da <- bind_rows(da) %>%
  dplyr::mutate(protein = rep(rownames(da[[1]]), times = length(da))) %>%
  dplyr::distinct(Strain.Name, protein, logFC, .keep_all = T)
temp <- sample_layout %>%
  dplyr::distinct(Strain.Name, .keep_all = T)
da <- left_join(da, temp, by = "Strain.Name") %>%
  dplyr::relocate(Strain.Name, .before = logFC) %>%
  dplyr::relocate(protein, .after = Strain.Name) %>%
  filter(Strain.Name != "WT")

da <- da %>%
  dplyr::select(protein, Strain.Name, logFC, P.Value, adj.P.Val) %>%                              # From here on in this function it's added by me
  dplyr::mutate(diffexpressed_adjusted = case_when((logFC > lfc_threshold) & (adj.P.Val < alpha) ~ "Up_regulated",
                                                   (logFC < -lfc_threshold) & (adj.P.Val < alpha) ~ "Down_regulated",
                                                   TRUE ~ "Not_significant"),
                diffexpressed_non_adjusted = case_when((logFC > lfc_threshold) & (P.Value < alpha) ~ "Up_regulated",
                                                       (logFC < -lfc_threshold) & (P.Value < alpha) ~ "Down_regulated",
                                                       TRUE ~ "Not_significant"))

da <- da %>%
  dplyr::distinct(Strain.Name, protein, logFC, .keep_all = T)

## Add a column to da with the number of replicates per KO - by me
unique_KOs <- unique(da$Strain.Name)
replicates <- c()
for (i in 1:length(unique_KOs)) {
  KO <- unique_KOs[i]
  replicates <- c(replicates, sum(grepl(KO, colnames(y_protein), fixed = T)))
}
KOs_replicates <- data.frame(unique_KOs, replicates)
colnames(KOs_replicates) <- c("Strain.Name", "Replicate_num")
da <- merge(da, KOs_replicates, by = "Strain.Name")

## Responsiveness
responsiveness <- da %>%
  group_by(Strain.Name) %>%
  summarise(nDEP = sum(abs(logFC) >= lfc_threshold & adj.P.Val <= alpha, na.rm = T),
            Up_regulated_adjusted = sum(abs(logFC) >= lfc_threshold & adj.P.Val <= alpha & diffexpressed_adjusted == "Up_regulated", na.rm = T),
            Down_regulated_adjusted = sum(abs(logFC) >= lfc_threshold & adj.P.Val <= alpha & diffexpressed_adjusted == "Down_regulated", na.rm = T),
            Up_regulated_non_adjusted = sum(abs(logFC) >= lfc_threshold & P.Value <= alpha & diffexpressed_non_adjusted == "Up_regulated", na.rm = T),
            Down_regulated_non_adjusted = sum(abs(logFC) >= lfc_threshold & P.Value <= alpha & diffexpressed_non_adjusted == "Down_regulated", na.rm = T),
            Replicate_num = mean(Replicate_num)) %>%
  mutate(Amino_acid_1_letter = str_sub(Strain.Name, 2, 2),
         chromosome_letter = substr(Strain.Name, 8, 8),
         anticodon = str_extract(Strain.Name, "(?<=\\()[[:alpha:]]{3}(?=\\))"),
         up_down_regulated_ratio_adjusted = Up_regulated_adjusted/Down_regulated_adjusted,
         up_down_regulated_ratio_non_adjusted = Up_regulated_non_adjusted/Down_regulated_non_adjusted) %>%
  arrange(Amino_acid_1_letter) %>%
  mutate(anticodon = factor(anticodon, levels = unique(anticodon)))

## Add all amino acid names
amino_acids <- as.data.frame(fread(paste(base_dir, "tRNA_KOs/Data/databases/GtRNAdb/amino_acids.csv", sep="")))
responsiveness <- left_join(responsiveness, amino_acids, by = "Amino_acid_1_letter")

## Save results from the DE analysis in the same way as in the file where I do DE separately for each batch
DE <- list(fit = fit,
           fit2 = fit2,
           da = da, 
           responsiveness= responsiveness)


## Remove unnecesary variables
rm(amino_acids, DE, de_proteins_list, fit, fit2, KOs_replicates, mm, temp, trna_levels, contrast.matrix, design, fit3,
   yeastmine, final_protein_names, general_protein_names, i, KO, output_file, replicates, sample_name, standard_protein_names, strain_name,
   systematic_protein_names, unique_KOs, trna_ko, dpcfit)





## Save DE proteins (up- and down-regulated separately)
de_proteins_up_down <- list()
strains <- unique(da$Strain.Name)

for (i in 1:length(strains)) {
  strain <- strains[i]
  temp <- da %>%
    dplyr::filter(Strain.Name == strain)
  
  up_regulated <- temp$protein[temp$diffexpressed_adjusted == "Up_regulated"]
  down_regulated <- temp$protein[temp$diffexpressed_adjusted == "Down_regulated"]
  
  de_proteins_up_down[[strain]] <- list("Up" = up_regulated,
                                        "Down" = down_regulated)
}

output_file <- paste(base_dir, "tRNA_KOs/Data/5k/de_proteins_list_up_down_", alpha_plain, "_logFC_", as.character(lfc_threshold), ".json", sep="")
write_json(de_proteins_up_down, path=output_file)







