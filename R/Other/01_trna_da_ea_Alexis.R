# Description ####
# Project: tRNAs, Alternative usage of AAs
# 01_trna_da_ea, Script for Differential and Enrichment Analysis
# author: Alexis García Avilés, date: 27.06.24


#Load packages
library(dplyr)
library(limma)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tibble)
library(stringr)
library(gprofiler2)



# Import data
proteomics_raw <- read.delim2('S:/AG/AG-CF-HTMS/AG-Ralser-Share/30-0092_AndreaLehmann-AlternativeAAUsage-tRNA/05_DataAnalysis/11_Preprocessing_Boris/AlternativeAAUsage-tRNA/AlternativeAAUsage-tRNA_peptidecentric_ProteinMaxLFQ_PCAoutlier_removed_batchcorrected.tsv', header = T)




# Format data ####
# Generate a named character vector to format colnames
sample_names <- colnames(proteomics_raw) # generate list of old colnames
names(sample_names) <- paste(colnames(proteomics_raw), as.character(proteomics_raw[1,]), sep = '_') # generate list of new colnames


# Format data to get a biological protein abundance dataframe with samples names as colnanames and gene names as rownames
trna_ko <- proteomics_raw %>%
  rename(all_of(sample_names)) %>% # rename colnames
  rename(genes = 'X_') %>% # rename genes colnames
  filter(!(genes %in% c('', 'Genes'))) %>% # filter rows with metadata
  select(-sample_group & !contains('QC')) %>% # remove UNIPORT ids and QCs columns
  column_to_rownames(var = 'genes') # convert gene name column to rownames


# convert to a numeric dataframe
trna_ko <- as.data.frame(lapply(trna_ko, as.numeric), row.names =  rownames(trna_ko))
rm(proteomics_raw, sample_names)




# Data Analysis ####
## Differential Analysis ####
sample_trna <- factor(str_extract(colnames(trna_ko), '(?<=\\_)[:graph:]+'))
sample_trna <- relevel(sample_trna, ref = 'WT')


# Generate design matrix
mm <- model.matrix(~sample_trna)
colnames(mm) <- levels(sample_trna)
voom(trna_ko, mm, plot = TRUE)


# Fit linear model
fit <- lmFit(log2(trna_ko), mm)
fit2 <- eBayes(fit, trend = TRUE)


# Identify differentially expressed proteins
da <- list()
for (i in colnames(mm)) {
  da[[i]] <- topTable(fit2, coef = i, adjust = 'BH', n = Inf, sort.by = 'none')
}

# Format data
da <- bind_rows(da) %>%
  mutate(protein = rep(rownames(da[[1]]), times = length(da)),
         ko = rep(names(da), each = nrow(da[[1]])), .before = 1) %>%
  `rownames<-`(NULL) %>%
  select(protein, ko, logFC, adj.P.Val) %>%
  filter(ko != 'WT')

# Remove intermediate objects
rm(mm, fit, fit2, i, sample_trna)


# Plot data
# Responsiveness
responsiveness <- da %>%
  group_by(ko) %>%
  summarise(nDEP = sum(abs(logFC) >= log2(1.5) & adj.P.Val <= 0.05)) %>%
  mutate(AA = str_sub(ko, 2, 2),
         chromosome = str_extract(str_remove(ko, '[:digit:]$'), '[:alpha:]$'),
         anticodon = str_extract(ko, '(?<=\\.)[:alpha:]{3}')) %>%
  arrange(AA) %>%
  mutate(anticodon = factor(anticodon, levels = unique(anticodon)))

p1_1 <- ggplot(responsiveness, aes(ko, nDEP)) +
  geom_bar(aes(fill = AA), stat = 'identity') +
  geom_hline(yintercept = 20, linetype = 'dashed') +
  geom_text(aes(label = nDEP), vjust = -1, size = 3) +
  labs(title = 'tRNA KOs Responsiveness', subtitle = 'Number of Differentially Expressed Proteins (DEPs)') +
  theme(panel.background = element_rect(fill = 'white', colour = 'gray'), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())

p1_2 <- ggplot(responsiveness, aes(ko, 1, groug = AA)) +
  geom_tile(aes(fill = AA)) +
  guides(fill = 'none') +
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_blank(), plot.subtitle = element_blank())
  
p1 <- p1_1 + p1_2 + plot_layout(design = 'A\nB', heights = c(0.9, 0.05), guides = 'collect')
rm(p1_1, p1_2)

p2 <- ggplot(responsiveness, aes(nDEP)) +
  geom_density(fill = '#8AC8CC', color = '#2C5985', alpha = 0.5, linewidth = 1.25) +
  geom_vline(xintercept = 20, linetype = 'dashed') +
  annotate('text', label = paste('>= 20 DEPs:', round(sum(responsiveness$nDEP >= 20)/nrow(responsiveness), 2), 'KOs', sep = ' ') , y = 0.05, x = 21, hjust = 0) +
  labs(title = 'Distribution of number of DEP per KO') +
  theme(panel.background = element_rect(fill = 'white', colour = 'gray'))

p3 <- ggplot(responsiveness, aes(AA, nDEP)) +
  geom_boxplot(aes(fill = AA)) +
  geom_jitter(data = mutate(responsiveness, nDEP = case_when(nDEP > 50 ~ 50, TRUE ~ nDEP))) +
  ylim(0, 50) +
  labs(title = 'tRNA KOs Responsiveness per AA', subtitle = 'Number of Differentially Expressed Proteins (DEPs) - truncated > 50 nDEP') +
  theme(panel.background = element_rect(fill = 'white', colour = 'gray'))

p4 <- ggplot(responsiveness, aes(chromosome, nDEP)) +
  geom_boxplot(aes(fill = chromosome)) +
  geom_jitter(data = mutate(responsiveness, nDEP = case_when(nDEP > 50 ~ 50, TRUE ~ nDEP))) +
  ylim(0, 50) +
  labs(title = 'tRNA KOs Responsiveness per Chromosome', subtitle = 'Number of Differentially Expressed Proteins (DEPs) - truncated > 50 nDEP') +
  theme(panel.background = element_rect(fill = 'white', colour = 'gray'))

p5 <- ggplot(responsiveness, aes(anticodon, nDEP, groupping = AA)) +
  geom_boxplot(aes(fill = AA)) +
  geom_jitter(data = mutate(responsiveness, nDEP = case_when(nDEP > 50 ~ 50, TRUE ~ nDEP))) +
  ylim = c(0, 50) +
  labs(title = 'tRNA KOs Responsiveness per Anticodon', subtitle = 'Number of Differentially Expressed Proteins (DEPs) - truncated > 50 nDEP') +
  theme(panel.background = element_rect(fill = 'white', colour = 'gray'), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



# Enrichment Analysis
ea <- list()
for (i in unique(da$ko)) {
  ea[[i]] <- pull(filter(da, ko %in% i & abs(logFC) >= log2(1.5) & adj.P.Val <= 0.05), protein)
}
rm(i)

ea <- gost(ea, organism = 'scerevisiae', correction_method = 'fdr', domain_scope = 'custom', custom_bg = unique(da$protein), sources = c('GO', 'KEGG', 'TF'))

ea_terms <- ea$result %>%
  group_by(query) %>%
  summarise(n_terms = n()) %>%
  mutate(AA = str_sub(query, 2, 2),
         chromosome = str_extract(str_remove(query, '[:digit:]$'), '[:alpha:]$'),
         anticodon = str_extract(query, '(?<=\\.)[:alpha:]{3}'))

p6 <- ggplot(ea_terms, aes(query, n_terms)) +
  geom_bar(aes(fill = AA), stat = 'identity') +
  labs(title = 'tRNA KOs Enriched Terms', subtitle = 'Number of Enriched Terms') +
  geom_text(aes(label = n_terms), vjust = -1, size = 3) +
  theme(panel.background = element_rect(fill = 'white', colour = 'gray'), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))

p7 <- ggplot(ea_terms, aes(n_terms)) +
  geom_density(fill = '#C2E393', color = '#3C8200', alpha = 0.5, linewidth = 1.25) +
  geom_vline(xintercept = 10, linetype = 'dashed') +
  annotate('text', label = paste('>= 10 terms:', round(sum(ea_terms$n_terms >= 10)/nrow(ea_terms), 2), 'KOs', sep = ' ') , y = 0.01, x = 50, hjust = 0) +
  labs(title = 'Distribution of number of Enriched KEGG/GO terms per KO') +
  theme(panel.background = element_rect(fill = 'white', colour = 'gray'))


ea_terms_ko <- ea$result %>%
  group_by(term_name) %>%
  summarise(n_kos = n(), source = unique(source)) %>%
  ungroup() %>%
  filter(n_kos > 25)

p8 <- ggplot(ea_terms_ko, aes(term_name, n_kos)) +
  geom_bar(aes(fill = source), stat = 'identity') +
  labs(title = 'Number of KOs per enriched Terms') +
  scale_x_discrete(label = function(x) abbreviate(x, minlength = 25)) +
  theme(panel.background = element_rect(fill = 'white', color = 'gray'), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

pdf(file = 'C:/MyStuff/tRNAs/20240702_tRNA_KO.pdf', onefile = T, width = 12, height = 5, paper = 'USr')
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
dev.off()

rm(p1, p2, p3, p4, p5, p6, p7, p8)





p1
p2
p3
p4
p5
p6
p7
p8








############################# My Stuff #############################

# Check how many replicates there are of each tRNA, QCs and WTs
trnas <- levels(sample_trna)
counts <- c()
  
for (i in 1:length(trnas)) {
  trna <- trnas[i]
  counts <- c(counts, sum(grepl(trna, colnames(trna_ko))))
}


# What's going on with the design matrix
apply(mm, 2, sum)
























