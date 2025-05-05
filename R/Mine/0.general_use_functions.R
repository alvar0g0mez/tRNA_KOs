# Packages
library(data.table)
library(dplyr)
library(readODS)
library(stringr)
library(roxygen2)
library(stringi)
library(stringr)



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
  
  # First of all, make sure all letters are in uppercase. Then turn data to a dataframe and proceed
  if (class(data) == "data.frame") {
    data$gene_names <- toupper(data$gene_names)
  }
  else {
    data <- toupper(data)
  }
  
  
  
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
#' above defined function "transcribe_nucleotide". This function takes
#' into account that codons and anticodons are both usually read in a 5'-to-3'
#' fashion, although for their pairing, the anticodon is oriented in the opposite
#' direction - so this function reverses the codon and then produces an anticodon
#' that matches this one. 
#' 
#' @param codon A 3-character string composed only of valid nucleotide letters 
#' (as defined by "transcribe_nucleotide()")
#' @return The corresponding anticodon (a 3-character string composed of only 
#' nucleotide-valid letters)
codon_to_anticodon <- function(codon) {
  anticodon <- ""
  codon <- stringi::stri_reverse(codon)
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
#' above defined function "reverse_transcribe_nucleotide". This function takes
#' into account that codons and anticodons are both usually read in a 5'-to-3'
#' fashion, although for their pairing, the anticodon is oriented in the opposite
#' direction - so this function reverses the anticodon and then produces a 
#' codon that matches this one. 
#' 
#' @param codon A 3-character string composed only of valid nucleotide letters 
#' (as defined by "reverse_transcribe_nucleotide()")
#' @return The corresponding codon (a 3-character string composed of only 
#' nucleotide-valid letters)
anticodon_to_codon <- function(anticodon) {
  codon <- ""
  anticodon <- stringi::stri_reverse(anticodon)
  for (i in 1:nchar(anticodon)) {
    nt <- substr(anticodon, i, i)
    new_nt <- reverse_transcribe_nucleotide(nt)
    codon <- paste(codon, new_nt, sep="")
  }
  return(codon)
}






#' Go through wrong anticodons from GtRNAdb and turn any T into U
#' 
#' Provide a wrong anticodon sequence from GtRNAdb as a 3-character string, 
#' and another 3-character string will be returned with the fixed anticodon. This
#' is only the case if the anticodon contains a T, which will then be turned into
#' a U. Otherwise, the inputed anticodon is returned without any modification.
#' 
#' @param codon A 3-character string 
#' @return Another 3-character string
turn_t_to_u_in_wrong_GtRNAdb_anticodons <- function(anticodon) {
  fixed_anticodon <- ""
  for (i in 1:nchar(anticodon)) {
    nt <- substr(anticodon, i, i)
    if (nt == "T") {
      new_nt <- "U"
    } else {
      new_nt <- nt
    }
    fixed_anticodon <- paste(fixed_anticodon, new_nt, sep="")
  }
  return(fixed_anticodon)
}













################################################################################
# Jaccard Similarity function
################################################################################
#' Jaccard similarity index
#' This is a statistic used for evaluating the similarity or diversity of sample
#' sets: it's just the intersection of 2 sets divided by their union. It takes in
#' 2 vectors and if they both contain at least one observation, returns a numeric
#' value. Otherwise, if at least one of them is empty, it returns an NA
#' 
#' @param a A vector containing the first set of values
#' @param b A vector containing another set of values
#' 
#' @return A numeric value, the Jaccard index between the 2 provided vectors, or
#' NA
jaccard <- function(a, b) {
  if (length(a) > 0 & length(b) > 0) {
    intersection <- length(intersect(a, b))
    union <- length(a) + length(b) - intersection
    return(intersection/union)
  }
  else {return(NA)}
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







# From https://stackoverflow.com/questions/57128090/remove-baseline-color-for-geom-histogram
# modified version of StatBin2 inherits from StatBin, except for an
# additional 2nd last line in compute_group() function
StatBin2 <- ggproto(
  "StatBin2", 
  StatBin,
  compute_group = function (data, scales, binwidth = NULL, bins = NULL, 
                            center = NULL, boundary = NULL, 
                            closed = c("right", "left"), pad = FALSE, 
                            breaks = NULL, origin = NULL, right = NULL, 
                            drop = NULL, width = NULL) {
    if (!is.null(breaks)) {
      if (!scales$x$is_discrete()) {
        breaks <- scales$x$transform(breaks)
      }
      bins <- ggplot2:::bin_breaks(breaks, closed)
    }
    else if (!is.null(binwidth)) {
      if (is.function(binwidth)) {
        binwidth <- binwidth(data$x)
      }
      bins <- ggplot2:::bin_breaks_width(scales$x$dimension(), binwidth, 
                                         center = center, boundary = boundary, 
                                         closed = closed)
    }
    else {
      bins <- ggplot2:::bin_breaks_bins(scales$x$dimension(), bins, 
                                        center = center, boundary = boundary, 
                                        closed = closed)
    }
    res <- ggplot2:::bin_vector(data$x, bins, weight = data$weight, pad = pad)
    
    # drop 0-count bins completely before returning the dataframe
    res <- res[res$count > 0, ] 
    
    res
})








################################################################################
# My version of reshape2::melt (not working anymore)
################################################################################
#' Take a correlation matrix and return its long format, where the first 2 columns
#' are "Var1" and "Var2", with all the possible combinations that we have in the 
#' correlation matrix, and the third column is "Value", with the corresponding 
#' correlation, between the entries in the first 2 columns in that row. 
#' 
#' This is necessary because reshape2 doesn't exist anymore for some reason. 
#' This is very easily done with dplyr::pivot_longer, which is what I am doing
#' here, I'm just adding a couple things to make the usage exactly equal to the 
#' old melt. 
#' 
#' @param df A correlation matrix or dataframe (i.e. colnames and rownames are 
#' the same, and the cells only contain correlations).
#' 
#' @return A dataframe where the first 2 columns are "Var1" and "Var2", with all
#' the possible combinations that we have in the correlation matrix, and the 
#' third column is "Value", with the corresponding correlation, between the 
#' entries in the first 2 columns in that row
melt <- function(df) {
  if (is.matrix(df))  {
    df <- as.data.frame(df)
  }
  
  df$Var1 <- as.factor(rownames(df))
  out <- pivot_longer(df, cols = -Var1, names_to = "Var2", values_to = "value")
  out$Var2 <- as.factor(out$Var2)
  return(out)
}








################################################################################
# dplyr::summarize without collapsing
################################################################################
#' Take a dataframe, the name of one of its columns and a name for the new output
#' column that will be produced. Return the same dataframe with an extra column,
#' containing the number of times the term in that row and in the afore 
#' mentioned column appears throughout the whole column. Basically what 
#' dplyr::summarize does, but without actually summarizing or collapsing the 
#' dataset in any way.
#' 
#' @param df Input dataframe
#' @param column_name A string, with the name of the column containing the
#' terms we want to count 
#' @param output_column_name A string, with the desired column name for the new
#' column in the dataframe
#' 
#' @return The input dataframe with a new column, counting the appearances of 
#' the terms in the specified column
summarize_mine <- function(df, column_name, output_column_name) {
  # Perform the count
  temp <- df %>%
    group_by(get(column_name)) %>%
    dplyr::count()
  
  # Match to full dataset
  colnames(temp) <- c(column_name, output_column_name)
  df <- left_join(df, temp, by = column_name)
  
  # Return output
  return(df)
}








################################################################################
# Extract anticodon from tRNA name
################################################################################
#' Just grab the 3 letters in between parenthesis, very simple but by defining 
#' it as a function it's easier to apply it to lists/vectors. 
#' 
#' @param trna A string containing the name of one of the tRNA KO strains, in 
#' the format "tA(AGC)M"
#' 
#' @return The anticodon in this tRNA, this is, the 3 letters in between 
#' parenthesis in the name of the tRNA KO strain - "AGC" in this case
extract_anticodon_from_trna_name <- function(trna) {
  start <- str_locate(trna, "\\(")[1]+1
  end <- str_locate(trna, "\\)")[1]-1
  out <- substr(trna, start, end)
  return(out)
}















