#' Adaptive column renaming
#' 
#' This function renames columns in case the input data.frame includes any
#' colnames required by internal functions (e.g., \code{"y"}).
#' 
#' @param df Input data.frame.
#' @param old_name Name of column to be renamed.
#' 

col_rename <- function(df, old_name) {
  k <- 1L
  converged <- FALSE
  while (!isTRUE(converged)) {
    new_name <- paste0(old_name, k)
    if (!new_name %in% colnames(df)) {
      coverged <- TRUE
    } else {
      k <- k + 1L
    }
  }
  return(new_name)
}