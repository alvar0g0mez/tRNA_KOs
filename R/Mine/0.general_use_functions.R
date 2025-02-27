# Packages
library(data.table)
library(dplyr)
library(readODS)
library(stringr)
library(roxygen2)



################### THIS FILE ONLY CONTAINS FUNCTIONS TO BE USED IN OTHER FILES ###################



#' Borrowed this one from the ScRAP
#' Match systematic to standard protein names and vice versa
#' 
#' We provide a dataframe or a vector with either systematic or standard protein
#' names, and the output is either a vector of (or a dataframe where one of the 
#' columns is) the corresponding standard/systematic protein names. It is 
#' important to notice that when standard names are required as output and there
#' is no standard name in the database for a certain protein, the systematic 
#' name is returned instead.
#' The default parameters take in a dataframe with a column containing
#' systematic names and return the same dataframe with an extra column with the 
#' standard protein names.
#' 
#' @param data This can be a vector with the systematic protein names, or a 
#' dataframe where one column has the systematic protein names. If it is a 
#' dataframe, the name of this column must be "gene_names"
#' @param yeastmine A dataframe with the database information for protein names 
#' in S. cerevisiae, as downloaded from AllianceMine. 
#' @param input A string indicating the type of names that are provided as 
#' input, and hence which type of names that should be returned (the other one).
#' Either "systematic" or "standard". 
#' @param simplify A boolean value indicating if we want the output to be simply
#' a vector with the standard protein names (TRUE), or the input dataframe where
#' the standard protein names are added as a new column (FALSE).
#' @param add_extra_columns A boolean value indicating, if simplify == FALSE, 
#' whether we only want to add to the dataframe the column with the standard 
#' protein names (FALSE) or also all other columns in the provided yeastmine 
#' dataframe.
#' 

match_systematic_and_standard_protein_names <- function(data, 
                                                       yeastmine,
                                                       input = "systematic",
                                                       simplify = FALSE,
                                                       add_extra_columns = FALSE) {
  # FOR SYSTEMATIC TO STANDARD
  if (input == "systematic") {
    
    # First of all, if we have received a vector as input, turn it into a dataframe and work from there
    if (class(data) == "character") {
      data <- data.frame(data)
      colnames(data) <- c("gene_names")
    }
    
    # Set the gene names column to Gene.secondaryIdentifier for merging with database
    colnames(data)[colnames(data) == "gene_names"] <- "Gene.secondaryIdentifier"
    
    # Match the names to the YeastMine ones
    df <- left_join(data, yeastmine, by = join_by(Gene.secondaryIdentifier))
    
    # Create the new column we'll keep as output, where we take standard gene names, but if this is not present, we fill it in with the systematic one
    df <- df %>% 
      mutate(Final.Ids = case_when(Gene.symbol == "" ~ Gene.secondaryIdentifier,
                                   is.na(Gene.symbol) ~ Gene.secondaryIdentifier,
                                   TRUE ~ Gene.symbol))
    
    # Prepare the output according to the specifications provided when calling the function
    if (simplify == TRUE) {
      out <- as.character(df$Final.Ids)
    }
    else {
      if (add_extra_columns == TRUE) {
        out <- df
      }
      else if (class(data) == "data.frame") {
        colnames_to_remove <- colnames(yeastmine)
        colnames_to_remove <- colnames_to_remove[!colnames_to_remove %in% c("Gene.secondaryIdentifier")]
        out <- df %>%
          dplyr::select(-any_of(colnames_to_remove))}
      else {
        out <- df %>%
          dplyr::select(Gene.secondaryIdentifier, Final.Ids)
      }
      
      # Change the final column name to Gene.symbol
      colnames(out)[colnames(out) == "Final.Ids"] <- "Gene.symbol"
    }
    
    # Return output
    return(out)
  }
  
  # FOR STANDARD TO SYSTEMATIC
  else if (input == "standard") {
    
    # First of all, if we have received a vector as input, turn it into a dataframe and work from there
    if (class(data) == "character") {
      data <- data.frame(data)
      colnames(data) <- c("gene_names")
    }
    
    # Set the gene names column to Gene.symbol for merging with database
    colnames(data)[colnames(data) == "gene_names"] <- "Gene.symbol"
    
    # Match the names to the YeastMine ones
    df <- left_join(data, yeastmine, by = join_by(Gene.symbol))
    
    # Create the new column we'll keep as output, where we take systematic gene names
    df <- df %>% 
      mutate(Final.Ids = case_when(Gene.secondaryIdentifier == "" ~ Gene.symbol,
                                   is.na(Gene.secondaryIdentifier) ~ Gene.symbol,
                                   TRUE ~ Gene.secondaryIdentifier))
    
    # Prepare the output according to the specifications provided when calling the function
    if (simplify == TRUE) {
      out <- as.character(df$Final.Ids)
    }
    else {
      if (add_extra_columns == TRUE) {
        out <- df
      }
      else if (class(data) == "data.frame") {
        colnames_to_remove <- colnames(yeastmine)
        colnames_to_remove <- colnames_to_remove[!colnames_to_remove %in% c("Gene.symbol")]
        out <- df %>%
          dplyr::select(-any_of(colnames_to_remove))}
      else {
        out <- df %>%
          dplyr::select(Gene.symbol, Final.Ids)
      }
      
      # Change the final column name to Gene.secondaryIdentifier
      colnames(out)[colnames(out) == "Final.Ids"] <- "Gene.secondaryIdentifier"
    }
    
    # Return output
    return(out)
  }
}






################################################################################
# Functions to turn a codon to its anticodon and viceversa
################################################################################


#' Change a nucleotide for another one following the DNA transcription pattern
#' 
#' Literally just provide one of the nucleotides that are observed in DNA sequences
#' (A, C, G, T) and the corresponding one in an RNA sequence will be returned 
#' (U, G, C, A). If you provide a different letter, an error message will pop up.
#' 
#' @param nt The nucleotide you want to transcribe
#' @return The corresponding nucleotide (a 1-character string)
transcribe_nucleotide <- function(nt) {
  if (nt == "A") {return("U")}
  else if (nt == "C") {return("G")}
  else if (nt == "G") {return("C")}
  else if (nt == "T") {return("A")}
  else if (nt == "U") {return("A")}                     # NEED TO ADD THIS TO THE PYTHON VERSION - FOR SOME REASON IN SOME DATA THE CODON SEQUENCES HAVE U IN THEM
  else {return("Provided letter is not a nucleotide")}
}


#' Change a nucleotide for another one following the opposite pattern to that of 
#' DNA transcription - so to get from an RNA sequence to the DNA sequence that 
#' produced it
#' 
#' Literally just provide one of the nucleotides that are observed in RNA sequences
#' (A, C, G, U) and the corresponding one in a DNA sequence will be returned 
#' (T, G, C, A). If you provide a different letter, an error message will pop up.
#' 
#' @param nt The nucleotide you want to reverse transcribe
#' @return The corresponding nucleotide (a 1-character string)
reverse_transcribe_nucleotide <- function(nt) {
  if (nt == "A") {return("T")}
  else if (nt == "C") {return("G")}
  else if (nt == "G") {return("C")}
  else if (nt == "U") {return("A")}
  else {return("Provided letter is not a nucleotide")}
}


#' Turn a codon into its corresponding anticodon
#' 
#' Provide a codon sequence as a 3-character string, and another 3-character
#' string will be returned with the corresponding anticodon. Makes use of the 
#' above defined function "transcribe_nucleotide".
#' 
#' @param codon A 3-character string composed only of valid nucleotide letters 
#' (as defined by "transcribe_nucleotide()")
#' @return The corresponding anticodon (a 3-character string composed of only 
#' nucleotide-valid letters)
codon_to_anticodon <- function(codon) {
  anticodon <- ""
  for (i in 1:nchar(codon)) {
    nt <- substr(codon, i, i)
    new_nt <- transcribe_nucleotide(nt)
    anticodon <- paste(anticodon, new_nt, sep="")
  }
  return(anticodon)
}


#' Turn an anticodon into the codon that it came from
#' 
#' Provide an anticodon sequence as a 3-character string, and another 3-character
#' string will be returned with the corresponding codon. Makes use of the 
#' above defined function "reverse_transcribe_nucleotide".
#' 
#' @param codon A 3-character string composed only of valid nucleotide letters 
#' (as defined by "reverse_transcribe_nucleotide()")
#' @return The corresponding codon (a 3-character string composed of only 
#' nucleotide-valid letters)
anticodon_to_codon <- function(anticodon) {
  codon <- ""
  for (i in 1:nchar(anticodon)) {
    nt <- substr(anticodon, i, i)
    new_nt <- reverse_transcribe_nucleotide(nt)
    codon <- paste(codon, new_nt, sep="")
  }
  return(codon)
}






################################################################################
# Functions to take log2 or log10 leaving 0 values as 0
################################################################################

#' log2 ignoring 0s - this one is not working or used anywhere, maybe delete it?
#' 
#' Literally just want to be able to take log2 of my data without generating -Inf 
#' values, so when the value is 0 this function returns 0, when it is positive it
#' returns its log2. 
#' 
#' @param value A number 
log2_ignoring_zeros <- function(value) {
  if (value == 0) {return(0)}
  else {return(log2(value))}
}




################################################################################
# Jaccard Similarity function
################################################################################
#' Jaccard similarity index
#' This is a statistic used for evaluating the similarity or diversity of sample
#' sets: it's just the intersection of 2 sets divided by their union. It takes in
#' 2 vectors and returns a numeric value
#' 
#' @param a A vector containing the first set of values
#' @param b A vector containing another set of values
#' 
#' @return A numeric value, the Jaccard index between the 2 provided vectors
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection/union)
}








################################################################################
# Generate Jaccard index matrix
################################################################################
#' Take a list of vectors, and a vector indicating the names of those vectors in
#' the list to be used. For those, calculating the Jaccard index pairwise between
#' all of them, and return the corresponding matrix.
#' 
#' @param list_of_vectors A list of vectors
#' @param vectors_to_be_used A vector containing all or some of the names of the
#' vectors in the list, those to be used
#' 
#' @return A numeric matrix containing Jaccard indexes
get_jaccad_index_matrix <- function(list_of_vectors, vectors_to_be_used) {
  # Create emtpy matrix to be filled and returned
  output_matrix <- matrix(nrow = length(vectors_to_be_used), ncol = length(vectors_to_be_used))
  colnames(output_matrix) <- rownames(output_matrix) <- vectors_to_be_used
  
  # Iterate over the vectors
  for (i in 1:length(vectors_to_be_used)) {
    name_1 <- vectors_to_be_used[i]
    vector_1 <- unlist(list_of_vectors[names(list_of_vectors) == name_1])
    
    # For each of the vectors, iterate over all vectors again, for each pair get
    # the Jaccard index and add it to a vector, which will be a row in the matrix
    new_row <- c()
    for(j in 1:length(vectors_to_be_used)) {
      name_2 <- vectors_to_be_used[j]
      vector_2 <- unlist(list_of_vectors[names(list_of_vectors) == name_2])
      
      new_row <- c(new_row, jaccard(vector_1, vector_2))
    }
    output_matrix[i,] <- new_row
  }
  return(output_matrix)
}











