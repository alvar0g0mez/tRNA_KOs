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
          select(-any_of(colnames_to_remove))}
      else {
        out <- df %>%
          select(Gene.secondaryIdentifier, Final.Ids)
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
          select(-any_of(colnames_to_remove))}
      else {
        out <- df %>%
          select(Gene.symbol, Final.Ids)
      }
      
      # Change the final column name to Gene.secondaryIdentifier
      colnames(out)[colnames(out) == "Final.Ids"] <- "Gene.secondaryIdentifier"
    }
    
    # Return output
    return(out)
  }
}
















# Old version just in case
#match_systematic_to_standard_protein_names <- function(data, 
#                                                       yeastmine,
#                                                       simplify = FALSE,
#                                                       add_extra_columns = FALSE) {
#  
#  # First of all, if we have received a vector as input, turn it into a dataframe and work from there
#  if (class(data) == "character") {
#    data <- data.frame(data)
#    colnames(data) <- c("Gene.secondaryIdentifier")
#  }
#  
#  # Match the names to the YeastMine ones
#  df <- left_join(data, yeastmine, by = join_by(Gene.secondaryIdentifier))
#  
#  # Create the new column we'll keep as output, where we take standard gene names, but if this is not present, we fill it in with the systematic one
#  df <- df %>% 
#    mutate(Final.Ids = case_when(Gene.symbol == "" ~ Gene.secondaryIdentifier,
#                                 is.na(Gene.symbol) ~ Gene.secondaryIdentifier,
#                                 TRUE ~ Gene.symbol))
#  
#  # Prepare the output according to the specifications provided when calling the function
#  if (simplify == TRUE) {
#    out <- as.character(df$Final.Ids)
#  }
#  else {
#    if (add_extra_columns == TRUE) {
#      out <- df
#    }
#    else if (class(data) == "data.frame") {
#      colnames_to_remove <- colnames(yeastmine)
#      colnames_to_remove <- colnames_to_remove[!colnames_to_remove %in% c("Gene.secondaryIdentifier")]
#      out <- df %>%
#        select(-c(colnames_to_remove))}
#    else {
#      out <- df %>%
#        select(Gene.secondaryIdentifier, Final.Ids)
#    }
#  }
#  
#  # Return output
#  return(out)
#}




