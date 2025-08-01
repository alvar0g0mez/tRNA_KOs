
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




















