library(dplyr)
library(tidyr)
library(seqinr)
library(data.table)
library(Biostrings)
library(stringr)



# Set directories to be used
working_from = "charite"

if (working_from == "home") {
  base_dir = "/home/alvaro/MyStuff/"
} else
  if (working_from == "charite") {
    base_dir = "C:/MyStuff/"
  }


# Load the master tRNA dataset
master_dataset <- as.data.frame(fread(paste(base_dir, "tRNA_KOs/Data/basic/master_tRNA_dataset.csv", sep="")))



# Modify FASTA files (obtained from GtRNAdb) and re-write them
## Raw DNA sequence - keep all sequences, labelled with GtRNAdb IDs
seqs <- readDNAStringSet(paste(base_dir, "tRNA_KOs/Data/sequence_alignments/sequence_fastas/raw/sacCer3-tRNAs.fa", sep=""))

new_names <- c()
for (i in 1:length(names(seqs))) {
  old_name <- names(seqs)[i]
  start <- str_locate(old_name, "_cerevisiae_")[2]+1
  end <- str_locate(old_name, " \\(tRNAscan-SE ID:")[1]-1
  new_name <- substr(old_name, start, end)
  new_names <- c(new_names, new_name)
}

names(seqs) <- new_names
writeXStringSet(seqs, paste(base_dir, "tRNA_KOs/Data/sequence_alignments/sequence_fastas/processed/DNA_seqs_all.fa", sep=""))



## Mature tRNA sequence - keep all sequences, labelled with GtRNAdb IDs
seqs <- readDNAStringSet(paste(base_dir, "tRNA_KOs/Data/sequence_alignments/sequence_fastas/raw/sacCer3-mature-tRNAs.fa", sep=""))

new_names <- c()
for (i in 1:length(names(seqs))) {
  old_name <- names(seqs)[i]
  start <- str_locate(old_name, "_cerevisiae_")[2]+1
  end <- str_locate(old_name, " \\(tRNAscan-SE ID:")[1]-1
  new_name <- substr(old_name, start, end)
  new_names <- c(new_names, new_name)
}

names(seqs) <- new_names
writeXStringSet(seqs, paste(base_dir, "tRNA_KOs/Data/sequence_alignments/sequence_fastas/processed/mature_seqs_all.fa", sep=""))




## Raw DNA sequence - keep only tRNAs that were KOd in our library, labelled with our Strain.Name
seqs <- readDNAStringSet(paste(base_dir, "tRNA_KOs/Data/sequence_alignments/sequence_fastas/processed/DNA_seqs_all.fa", sep=""))

### Keep only the sequences which were included as KOs in our library
tRNAs_to_select <- master_dataset %>%
  dplyr::filter(KOd == "Yes") %>%
  dplyr::pull(GtRNAdb_gene_symbol)
seqs_filtered <- seqs[(names(seqs) %in% tRNAs_to_select)]
temp_names <- data.frame(GtRNAdb_gene_symbol = names(seqs_filtered))
master_temp <- master_dataset %>%
  dplyr::select(GtRNAdb_gene_symbol, Strain.Name, KOd)
new_names <- left_join(temp_names, master_temp, by = "GtRNAdb_gene_symbol") %>%
  dplyr::pull(Strain.Name)

### Save with new names
names(seqs) <- new_names
writeXStringSet(seqs, paste(base_dir, "tRNA_KOs/Data/sequence_alignments/sequence_fastas/processed/DNA_seqs_only_KO.fa", sep=""))



## Mature tRNA sequence - keep only tRNAs that were KOd in our library, labelled with our Strain.Name
seqs <- readDNAStringSet(paste(base_dir, "tRNA_KOs/Data/sequence_alignments/sequence_fastas/processed/mature_seqs_all.fa", sep=""))

### Keep only the sequences which were included as KOs in our library
tRNAs_to_select <- master_dataset %>%
  dplyr::filter(KOd == "Yes") %>%
  dplyr::pull(GtRNAdb_gene_symbol)
seqs_filtered <- seqs[(names(seqs) %in% tRNAs_to_select)]
temp_names <- data.frame(GtRNAdb_gene_symbol = names(seqs_filtered))
master_temp <- master_dataset %>%
  dplyr::select(GtRNAdb_gene_symbol, Strain.Name, KOd)
new_names <- left_join(temp_names, master_temp, by = "GtRNAdb_gene_symbol") %>%
  dplyr::pull(Strain.Name)

### Save with new names
names(seqs) <- new_names
writeXStringSet(seqs, paste(base_dir, "tRNA_KOs/Data/sequence_alignments/sequence_fastas/processed/mature_seqs_only_KO.fa", sep=""))






# All tRNAs - using GtRNAdb symbol, since there is no Strain.Name for the tRNAs that were not KOd
## DNA sequence
sequences <- as.list(master_dataset$DNA_sequence)
seq_names <- as.list(master_dataset$GtRNAdb_gene_symbol)
write.fasta(sequences = sequences,
            names = seq_names,
            file.out = paste(base_dir, "tRNA_KOs/Data/sequence_alignments/sequence_fastas/all_DNA_sequences.fasta", sep=""))

## Mature sequence
sequences <- as.list(master_dataset$mature_sequence)
seq_names <- as.list(master_dataset$GtRNAdb_gene_symbol)
write.fasta(sequences = sequences,
            names = seq_names,
            file.out = paste(base_dir, "tRNA_KOs/Data/sequence_alignments/sequence_fastas/all_mature_sequences.fasta", sep=""))




# Only tRNAs that were KOd in our experiment - using Strain.Name 
master_dataset <- master_dataset %>%
  dplyr::filter(Strain.Name != "")
  
## DNA sequence
sequences <- as.list(master_dataset$DNA_sequence)
seq_names <- as.list(master_dataset$Strain.Name)
write.fasta(sequences = sequences,
            names = seq_names,
            file.out = paste(base_dir, "tRNA_KOs/Data/sequence_alignments/sequence_fastas/only_KOd_DNA_sequences.fasta", sep=""))

## Mature sequence
sequences <- as.list(master_dataset$mature_sequence)
seq_names <- as.list(master_dataset$Strain.Name)
write.fasta(sequences = sequences,
            names = seq_names,
            file.out = paste(base_dir, "tRNA_KOs/Data/sequence_alignments/sequence_fastas/only_KOd_mature_sequences.fasta", sep=""))













# CHECKS
## which sequences are exactly the same between DNA sequence and mature sequence, and which ones are not
library(helperfunctions)

ts_to_us <- function(letter) {
  if (letter == "T") {
    return("U")
  }
  else {
    return(letter)
  }
}

DNA_to_mature <- function(sequence) {
  nts <- nchar(sequence)
  out <- ""
  for (i in 1:nts) {
    nt <- substr(sequence, i, i)
    new_nt <- ts_to_us(nt)
    out <- paste(out, new_nt, sep="")
  }
  return(out)
}


transcribed_sequences <- c()
for (i in 1:nrow(master_dataset)) {
  input <- master_dataset$DNA_sequence[i]
  transcribed_sequences <- c(transcribed_sequences, DNA_to_mature(input))
}


sum(transcribed_sequences == master_dataset$mature_sequence)
which(transcribed_sequences != master_dataset$mature_sequence)

print(master_dataset$DNA_sequence[96])
print(transcribed_sequences[96])
print(master_dataset$mature_sequence[96])



## How many tRNAs do we not have sequences for? 7, doesn't sound terrible
sum(master_dataset$DNA_sequence == "" & master_dataset$Strain.Name != "")
sum(master_dataset$mature_sequence == "")

check <- master_dataset %>%
  dplyr::filter(Strain.Name != "")




























