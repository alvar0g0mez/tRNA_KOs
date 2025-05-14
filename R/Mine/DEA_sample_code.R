library(limma)

# Example: create a factor for your treatments
treatments <- factor(colnames(fit$coefficients))

# Create a design matrix (without intercept, to get one column per treatment)
design <- model.matrix(~ 0 + treatments)
colnames(design) <- levels(treatments)

# Generate contrast strings for each treatment vs WT
contrast_list <- paste0(levels(treatments)[-1], "-WT")

# Use makeContrasts with these contrasts
contrast.matrix <- makeContrasts(contrasts = contrast_list, levels = design)

















# Example: create a factor for your treatments
treatments <- factor(c("WT", paste0("T", 1:196)))

# Create a design matrix (without intercept, to get one column per treatment)
design <- model.matrix(~ 0 + treatments)
colnames(design) <- levels(treatments)















