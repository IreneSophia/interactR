#' Compares pseudo and observed WLCC for each Lag
#'
#' Takes the dataframe created by [extractIPS()] with method WLCC and compares 
#' the pseudo and the observed WLCC values based on the mean. Comparison 
#' can be either made by aggregating the values across Dyads and then comparing
#' the pseudo and the observed values for each lag with a t-test to assess whether 
#' observed is larger than pseudo. Frequentist or Bayesian t-tests can be chosen
#' by using the `Bayesian` logical switch. Alternatively, if the sample of dyads is small, 
#' observed WLCC can be aggregated across all dyads and compared with a larger sample of 
#' pseudo WLCC in a permutation test by setting `perm = TRUE`. 
#' When frequentist t-tests are computed, multiple comparison correction is applied based on 
#' Benjamini & Hochberg's (1995) FDR method for the number of lags. 
#'
#' @param df.observed Dataframe. Dataframe created by [extractIPS()] extracting observed values.
#' @param df.pseudo Dataframe. Dataframe created by [extractIPS()] extracting pseudo values. 
#' @param Bayesian Logical. Switch to perform Bayesian or frequentist testing. Default is `TRUE`.
#' @param perm Logical. Switch to use permutation testing instead of comparison across Dyads.
#'   If the dataset is smaller, then one can create a larger set of pseudo values 
#'   against which the mean observed value per Lag can be compared. Default is `FALSE`.
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

compareWTC = function(df.observed, df.pseudo, Bayesian = T, perm = F, rs.path = c(), suffix = "", 
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
    features = intersect(df.pseudo$Feature, df.observed$Feature)
    
    # check that df.pseudo only contains one type of pseudo WLCC
    if (length(unique(df.pseudo$Method)) != 1) stop("Function performs comparison with one pseudo method.")
    
    # combine the dataframes
    df = rbind(df.pseudo, df.observed) |>
      filter(Feature %in% features)
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Comparing pseudo and observed WLCC from", features, "\n")
    
    # aggregate the values across the windows
    df.lag = df |>
      group_by(Lag, Dyad, Method, Feature) |>
      summarise(
        WLCC = mean(WLCC, na.rm = T),
        .groups = "drop"
      )
    
    # aggregate values for plotting
    df.agg = df.lag |>
      group_by(Lag, Method, Feature) |>
      summarise(
        STD = sd(WLCC),
        AVG = mean(WLCC),
        .groups = "drop"
      )
    
    if (!perm) {
      
      if (!Bayesian) {
      
        # use non-parametric tests to assess differences
        df.stat = df.lag |> 
          group_by(Feature, Lag) |> 
          rstatix::t_test(WLCC ~ Method, detailed = T, alternative = "g") |> 
          group_by(Feature) |>
          rstatix::adjust_pvalue(method = "BH") |>
          mutate(
            Sig  = if_else(p.adj < 0.05, "*", "")
          ) |> rename(p_BH = p.adj)
        
        # add to the aggregated dataframe
        df.agg = df.agg |>
          left_join(df.stat |> select(Feature, Lag, p, p.adj, Sig),
                    by = c("Lag", "Feature"))
        
      } else {
        
        # compute the Bayesian t-tests
        df.stat = df.lag |>
          group_by(Feature, Lag) |>
          group_modify(~ {
            observed = .x |>
              filter(Method == "observed") |>
              summarise(mean = mean(WLCC, na.rm = TRUE)) |> pull(mean)
            pseudo = .x |>
              filter(Method != "observed") |>
              summarise(mean = mean(WLCC, na.rm = TRUE)) |> pull(mean)
            ttest = BayesFactor::ttestBF(formula = WLCC ~ Method, data = as.data.frame(.x))
            data.frame(logBF = ttest@bayesFactor$bf, Direction = if_else(observed > pseudo, "greater", "lesser"))
          }) |>
          ungroup() |>
          mutate(
            Cred = if_else(logBF > log(3) & Direction == "greater", "*", "")
          )
        
        # add to the aggregated dataframe
        df.agg = df.agg |>
          left_join(df.stat |> select(Feature, Lag, logBF, Cred),
                    by = c("Lag", "Feature"))
        
      }
      
    } else {
      
      # aggregate the observed values to get grand average per Lag
      df.perm = df.observed |>
        group_by(Lag, Dyad, Feature) |>
        summarise(
          WLCC = mean(WLCC, na.rm = T),
          .groups = "drop"
        ) |> 
        group_by(Lag, Feature) |>
        summarise(
          observed = mean(WLCC, na.rm = T),
          .groups = "drop"
        ) |> 
        # merge with the pseudo dataframe
        right_join(df.pseudo, by = c("Lag", "Feature"))
      
      # calculate permutation values
      df.stat = df.perm |>
        group_by(Feature, Lag) |>
        summarise(
          Probability = mean(observed > WLCC),
          p = 1 - Probability,
          .groups = "drop"
        ) |> group_by(Feature) |> 
        # adjust the permutation value for the number of lags
        rstatix::adjust_pvalue(method = "BH") |>
        mutate(
          Probability_BH = 1 - p.adj,
          Sig  = if_else(p.adj < 0.05, "*", "")
        ) |> select(-p, -p.adj) |> ungroup()
      
      # add the result to the aggregated dataframe
      df.agg = df.agg |>
        left_join(df.stat, by = c("Lag", "Feature"))
        
    }
    
    # save the data
    if (!is.null(rs.path)) write_csv(df.agg, file = flnm)
    
  }
  
  if (return) return(df.agg)
  
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