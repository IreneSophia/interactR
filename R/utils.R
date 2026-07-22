#' Checks column names of dataframe
#'
#' Assesses whether each column named in `colnames` is a column of the dataframe `df`.
#'
#' @param df Dataframe. The dataset to be checked.
#' @param colnames Character vector. The exact name or names of the column(s) in \code{df} 
#' which must be in `df`.
#'
#' @return Throws an error if any of the `colnames` are missing in `df`. 
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' 
checkDF = function(df, colnames) {
  colDiff = paste(setdiff(colnames, colnames(df)), collapse = ", ")
  if (nchar(colDiff) > 0) stop("Missing columns in df: ", colDiff)
}