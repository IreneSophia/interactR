#' Aggregate a dataframe containing the results of WLCC to extract relevant features 
#'
#' Takes the dataframe created by [extractIPS()] with the method WLCC and aggregates the results.
#' Identifier values contain leading by this Identifier, while Dyad values contain the overall
#' WLCC across all Lags. 
#'
#' @param df Dataframe. Dataframe created by [extractIPS()] with method WLCC.
#' @param FUN Function. Function used for the aggregation across Lags for each window.
#' @param absolute Logical. Switch whether to convert all the values to their 
#'   absolute before aggregation. Default is `FALSE`.
#' @param rs.path Character. Path to destination directory for saved files. 
#'   If empty (is.null(rs.path) == TRUE), then nothing is saved. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns dataframe. If provided, the dataframe is saved as rds to `rs.path`.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @import dplyr
#' @export
#' 

featWLCC = function(df, FUN, absolute = F, r2z = F, rs.path = c(), suffix = "", 
                   verbose = T, recompute = F, return = T) {
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename 
    flnm  = file.path(rs.path, sprintf("featWLCC%s.csv", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading WLCC feature dataframe\n")
      df.out = read_csv(flnm)
    }
  } else {
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Extracting WLCC features from dataframe\n")
    
    # potentially use Fisher's z-transformation
    if (r2z) {
      df = df |>
        mutate(
          # winsorise to avoid infinity
          WLCC = pmin(pmax(WLCC, -0.9999), 0.9999),
          # perform the transformation
          WLCC = atanh(WLCC)
        )
    }
    # potentially remove sign from all WLCC values
    if (absolute) df = df |> mutate(WLCC = abs(WLCC))
  
    # get leading by Identifier
    df.indi = df |>
      mutate(
        Identifier = case_when(
          Leading == "Left" ~ gsub("-.*", "", Dyad),
          Leading == "Right" ~ gsub(".*-", "", Dyad)
        )
      ) |> tidyr::drop_na() |>
      # first, aggregate using FUN for each window
      group_by(Dyad, Time, Method, Identifier, Feature, Window) |>
      summarise(
        WLCC = FUN(WLCC),
        .groups = "drop"
      ) |> 
      # then, aggregate for each Dyad over all windows
      group_by(Dyad, Time, Method, Identifier, Feature) |>
      summarise(
        STD   = sd(WLCC),
        AVG = mean(WLCC),
        .groups = "drop"
      ) |>
      tidyr::pivot_wider(names_from = Feature, values_from = c(AVG, STD),
                         names_glue = "WLCC_{Feature}_{.value}")
    
    # get overall WLCC
    df.dyad = df |>
      # first, aggregate using FUN for each window
      group_by(Dyad, Time, Method, Feature, Window) |>
      summarise(
        WLCC = FUN(WLCC),
        .groups = "drop"
      ) |> 
      # then, aggregate for each Dyad over all windows
      group_by(Dyad, Time, Method, Feature) |>
      summarise(
        STD   = sd(WLCC),
        AVG = mean(WLCC),
        .groups = "drop"
      ) |>
      tidyr::pivot_wider(names_from = Feature, values_from = c(AVG, STD),
                         names_glue = "DyadWLCC_{Feature}_{.value}")
    
    df.out = merge(df.indi, df.dyad) |> 
      relocate(Dyad, Time, Identifier)
    
    # save the data
    if (!is.null(rs.path)) saveRDS(df.out, file = flnm)
    
  }
  
  if (return) return(df.out)
  
}