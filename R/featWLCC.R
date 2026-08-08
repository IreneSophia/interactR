#' Compares pseudo and observed WLCC for each Lag
#'
#' Takes the dataframe created by [extractIPS()] with method WLCC and compares 
#' the pseudo and the observed WLCC values based on the mean. 
#'
#' @param df.observed Dataframe. Dataframe created by [extractIPS()] extracting observed values.
#' @param df.pseudo Dataframe. Dataframe created by [extractIPS()] extracting pseudo values. 
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

compareWTC = function(df.observed, df.pseudo, rs.path = c(), suffix = "", 
                      verbose = T, recompute = F, return = T) {
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename 
    flnm  = file.path(rs.path, 
                      sprintf("featWLCC_pseudo-comp%s.csv", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading WLCC comparison\n")
      df.rsq = read_csv(flnm)
    }
  } else {
    
    # get the shared Features
    features = union(df.pseudo$Feature, df.observed$Feature)
    
    # check that df.pseudo only contains one type of pseudo WLCC
    if (length(unique(df.pseudo$Method)))
    
    # combine the dataframes
    df = rbind(df.pseudo, df.observed) |>
      filter(Feature %in% features)
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Comparing pseudo and observed WLCC from", features, "\n")
    
    # aggregate the values
    df = df |>
      group_by(Lag, Dyad, Method, Feature) |>
      summarise(
        WLCC = mean(WLCC, na.rm = T),
        .groups = "drop"
      )
    
    # use non-parametric tests to assess differences
    df.stat = df |> 
      group_by(Feature, Lag) |> 
      rstatix::wilcox_test(WLCC ~ Method, detailed = T) |> 
      rstatix::adjust_pvalue(method = "BH") |>
      mutate(
        Sig  = if_else(p.adj < 0.05, "*", ""),
        Direction = case_when(estimate > 0 & p.adj < 0.05 ~ "correlation",
                         estimate < 0 & p.adj < 0.05 ~ "hypo-correlation",
                         T ~ "")
      )
    
    # aggregate for plotting
    df.agg = df |>
      group_by(Lag, Method, Feature) |>
      summarise(
        STD = sd(WLCC),
        AVG = mean(WLCC),
        .groups = "drop"
      ) |>
      left_join(df.stat |> select(Feature, Lag, statistic, p, p.adj, estimate, conf.low, conf.high, Sig, Direction),
                by = c("Lag", "Feature"))
    
    # save the data
    if (!is.null(rs.path)) write_csv(df.rsq, file = flnm)
    
  }
  
  if (return) return(df.rsq)
  
}

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

featWTC = function(df, FUN, absolute = F, r2z = F, rs.path = c(), suffix = "", 
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
      group_by(Dyad, Time, Type, Identifier, Feature, Window) |>
      summarise(
        WLCC = FUN(WLCC),
        .groups = "drop"
      ) |> 
      # then, aggregate for each Dyad over all windows
      group_by(Dyad, Time, Type, Identifier, Feature) |>
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
      group_by(Dyad, Time, Type, Feature, Window) |>
      summarise(
        WLCC = FUN(WLCC),
        .groups = "drop"
      ) |> 
      # then, aggregate for each Dyad over all windows
      group_by(Dyad, Time, Type, Feature) |>
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