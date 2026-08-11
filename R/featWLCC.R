#' Aggregate Dataframe of WLCC Results
#'
#' Aggregates values extracted by [extractIPS()] with the method WLCC using the
#' standard deviation, mean, median and mean of window peaks. Identifier values
#' contain leading by this Identifier, while Dyad values contain the overall
#' WLCC across all Lags. 
#'
#' @param df Dataframe. Dataframe created by [extractIPS()] with method WLCC.
#' @param absolute Logical. Switch whether to convert all the values to their 
#'   absolute before aggregation. Default is `FALSE`.
#' @param rs.path Character. Path to destination directory for saved files. 
#'   If empty (is.null(rs.path) == TRUE), then nothing is saved. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns dataframe. If provided, the dataframe is saved to `rs.path`.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @import dplyr
#' @export
#' 

featWLCC = function(df, absolute = F, r2z = F, rs.path = c(), suffix = "", 
                   verbose = T, recompute = F, return = T) {
  
  if (verbose) cat("-------------- Extracting IPS features based on WLCC  --------------\n")
  
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
      if (verbose) cat(format(Sys.time(), "%X"), ": Loading WLCC feature dataframe\n")
      df.out = readr::read_csv(flnm)
    }
  } else {
    
    if (verbose) cat(format(Sys.time(), "%X"), ": Extracting WLCC features from dataframe\n")
    
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
    
    # prepare the data frame with Identifier information
    df = df |> 
      tidyr::drop_na() |>
      mutate(
        Identifier = case_when(
          Leading == "Left" ~ gsub("-.*", "", Dyad),
          Leading == "Right" ~ gsub(".*-", "", Dyad),
          T ~ "Simultaneous"
        )
      )
    
    # get the peaks for the Identifier leading and for the Dyad
    df.peak = df |>
      # first, aggregate using peak detection for each window
      group_by(Dyad, Time, Method, Identifier, Feature, Window) |>
      summarise(
        indi = max(WLCC, na.rm = T),
        .groups = "drop"
      ) |>
      group_by(Dyad, Time, Method, Feature, Window) |>
      mutate(
        dyad = max(indi, na.rm = T)
      ) |> filter(Identifier != "Simultaneous") |>
      # then, aggregate for over all windows
      group_by(Dyad, Time, Method, Identifier, Feature) |>
      summarise(
        DyadWLCC = mean(dyad, na.rm = T),
        WLCC     = mean(indi, na.rm = T),
        .groups = "drop"
      ) |>
      tidyr::pivot_wider(values_from = ends_with("WLCC"), 
                         names_from = Feature, 
                         names_glue = "{.value}_{Feature}_PEAK")
  
    # get leading by Identifier
    df.indi = df |> filter(Identifier != "Simultaneous") |>
      # aggregate for each Dyad and Identifier over all windows
      group_by(Dyad, Time, Method, Identifier, Feature) |>
      summarise(
        STD  = sd(WLCC, na.rm = T),
        AVG  = mean(WLCC, na.rm = T),
        MED  = median(WLCC, na.rm = T),
        .groups = "drop"
      ) |>
      tidyr::pivot_wider(names_from = Feature, 
                         values_from = c(AVG, MED, STD),
                         names_glue = "WLCC_{Feature}_{.value}")
    
    # get overall WLCC
    df.dyad = df |>
      # aggregate for each Dyad over all windows
      group_by(Dyad, Time, Method, Feature) |>
      summarise(
        STD  = sd(WLCC, na.rm = T),
        AVG  = mean(WLCC, na.rm = T),
        MED  = median(WLCC, na.rm = T),
        .groups = "drop"
      ) |>
      tidyr::pivot_wider(names_from = Feature, 
                         values_from = c(AVG, MED, STD),
                         names_glue = "DyadWLCC_{Feature}_{.value}")
    
    df.out = merge(df.indi, df.dyad) |> 
      merge(df.peak) |>
      relocate(Dyad, Time, Identifier)
    
    # save the data
    if (!is.null(rs.path)) {
      if (verbose) cat(format(Sys.time(), "%X"), ": Saving the data\n")
      readr::write_csv(df.out, file = flnm)
    }
    
  }
  
  if (verbose) cat(format(Sys.time(), "%X"), ": Done\n")
  
  if (return) return(df.out)
  
}