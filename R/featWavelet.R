
#' Aggregate a dataframe containing the results of WTC to extract relevant features 
#'
#' Takes the dataframe created by [extractIPS()] and aggregates the results in prespecified Frequency
#' Bins and for Phase also in percentages of Angle Bins. Limits for Bins and aggregation function 
#' may differ between Rsq and Phase. Values outside COI are excluded. 
#'
#' @details Relative phase difference was calculated as Phi_left - Phi_right. Thus, if Phase values
#'   are positive, this indicates the left Identifier was leading, while negative values indicate that
#'   the right Identifier was leading (based on the Dyad ID). Values closer to 0 indicate in-phase 
#'   and values closer to +- pi indicate anti-phase.
#'
#' @param df Dataframe. Dataframe created by [extractIPS()].
#' @param Limits List. List containing one or two dataframes with numeric of frequency limits in Hz.
#'   First entry is used for the Coherence Bins, optionally a second for Phase Bins. If only
#'   one dataframe with columns `upper` and `lower` as limits is provided, it is used for both.
#'   Lower limits are included in the Bins but upper limits are not. 
#' @param Funs List. List of function to be used for the aggregation. First is used for Coherence. 
#'   Optionally, a second is used for Phase. If Funs contains only one function,
#'   then this function is used for both Phase and Coherence. 
#' @param Bins Numeric. Number of Angle Bins for Phases. Default is `8`.
#' @param rs.path Character. Path to destination directory for saved files. 
#'   If empty (is.null(rs.path) == TRUE), then nothing is saved. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns dataframe. If provided, the dataframe is saved as a csv to `rs.path`.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @import dplyr
#' @export
#' 

featWTC = function(df, Limits, Funs, Bins = 8, 
                   rs.path = c(), suffix = "", 
                   verbose = T, recompute = F, return = T) {
  
  if (verbose) cat("--------------- Extracting IPS features based on WTC ---------------\n")
  
  # check if the variables are correct
  if ((class(Limits) != "list")  | (class(Funs) != "list")) stop("Both Limits and Funs need to be lists.")
  if (!(length(Limits) %in% 1:2) | !(length(Limits) %in% 1:2)) stop("Both Limits and Funs need to be of length 1 or 2 - first for coherence and, potentially, second for phase.")
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename 
    flnm  = file.path(rs.path, sprintf("featWTC%s.csv", suffix))
  }
  
  # if no recompute and the CSV file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (verbose) cat(format(Sys.time(), "%X"), ": Loading WTC features\n")
    df.agg = readr::read_csv(flnm, show_col_types = F)
    
  } else {
    # no recompute and the file exists, it is simply loaded
    if (verbose) cat(format(Sys.time(), "%X"), ": Extracting WTC features\n")
  
    # if two Limits but only one Fun, then use the same Fun twice
    if ((length(Limits) == 2) & (length(Funs) == 1)) Funs = c(Funs, Funs)
    
    # if only one limits, use the same for both
    if (length(Limits) == 1) Limits = c(Limits, Limits)
    
    # ensure only within COI
    df = df |> filter(WithinCOI)
    
    # aggregate for Rsq
    df.rsq = Limits[[1]] |>
      left_join(
        df, 
        join_by(lower <= Frequency, upper > Frequency), 
        relationship = "many-to-many"
      ) |>
      group_by(lower, upper, Dyad, Method) |>
      summarise(
        value = Funs[[1]](Rsq, na.rm = TRUE), 
        .groups = "drop"
      ) |>
      mutate(
        name = sprintf("Rsq_%.2fd-%.2fHz", lower, upper)
      ) |> select(Dyad, Method, name, value) |>
      tidyr::pivot_wider()
    
    # aggregate for Phase
    df.phase = Limits[[2]] |>
      left_join(
        df, 
        join_by(lower <= Frequency, upper > Frequency), 
        relationship = "many-to-many"
      ) |>
      group_by(lower, upper, Dyad, Method) |>
      summarise(
        data = list({
          tibble(Phase = Phase) %>%
            mutate(
              Bin = cut(
                Phase, 
                breaks = seq(-pi, pi, length.out = Bins + 1), 
                include.lowest = TRUE
              )
            ) %>%
            count(Bin, name = "count", .drop = FALSE) |>
            mutate(Bin = as.numeric(as.factor(Bin)),
                   count = (count / sum(count)))
          }),
        Phase = Funs[[2]](Phase, na.rm = TRUE), 
        .groups = "drop"
      ) |>
      mutate(
        name = sprintf("%.2f-%.2fHz", lower, upper)
      ) |> tidyr::unnest(data) |> 
      select(Dyad, Method, name, Phase, count, Bin) |>
      tidyr::pivot_wider(values_from = count, names_from = Bin, names_prefix = "PhaseBin") |>
      tidyr::pivot_wider(values_from = starts_with("Phase"))
    
    # combine both
    df.agg = merge(df.rsq, df.phase)
    
    # save the data
    if (!is.null(rs.path)) {
      if (verbose) cat(format(Sys.time(), "%X"), ": Saving the data\n")
      readr::write_csv(df.agg, flnm)
    }
  }
  if (verbose) cat(format(Sys.time(), "%X"), ": Done\n")
  
  if (return) return(df.agg)
  
}