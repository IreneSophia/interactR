#' Extract Zero-Crossing Frequency from Time Series Data
#'
#' Evaluates localised oscillatory characteristics within specified time series 
#' channels by quantifying zero-crossing counts over moving windows. This allows
#' for the extraction of nodding and head shaking from rotational head movement.
#' The function is loosely based on the algorithm presented in Hale et al. (2020) and
#' shifts data dynamically against a running local mean, calculates zero-crossing 
#' frequencies, filters outcomes within a target frequency band and applies a 
#' smoothing threshold filter across the results. 
#' 
#' @details 
#' The calculation executes sequentially across one or more specified columns:
#' 1. **Baseline Centering:** Removes slow-moving baseline drift using a running mean ([aggSlide()]).
#' 2. **Frequency Mapping:** Counts sign inversions using [findZCrossing()] over an explicit window interval.
#' 3. **Bandpass Filtering:** Identifies sequences falling strictly within the bounds defined by `minFreq` and `maxFreq`.
#' 4. **Smoothing:** Enforcing a secondary duration validation window with only segments of at least half its length considered relevant.
#'
#' Newly generated vectors are appended back to the source data frame utilizing a 
#' systematic naming template (e.g., `[column]_centred`, `[column]_sum`, 
#' `[column]_rel`, and `[column]_smooth`).
#'
#' @param df Dataframe. The dataset containing the variables to be processed. Must explicitly feature columns `Dyad`, 
#'   `Identifier`, `Frame`, `Timestamp`, `Time` as well as all columns contained in `colnames`. 
#' @param colnames Character vector. The exact name or names of the column(s) in \code{df} from which
#'   to extract zero-crossing features. 
#' @param fps Numeric. Frame processing rate frequency profile (frames per second) of the dataset.
#' @param minDegree Numeric. How many degree of rotational difference are needed for the movement to be considered relevant. 
#' Depends on the fps and the specific movement. Setting to negative number leads to no thresholding based on degrees. 
#' @param rs.path Path to the directory where the output file will be saved, if empty (is.null(rs.path) == TRUE), 
#'   then nothing is saved. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param win Numeric. Window duration scale evaluated in seconds for the moving frequency summary. Default is \code{2}.
#' @param minFreq Numeric. The lower cutoff boundary of the targeted frequency band in Hz. Default is \code{1.5}.
#' @param maxFreq Numeric. The upper cutoff boundary of the targeted frequency band in Hz. Default is \code{6.5}.
#' @param winCentre Numeric. Seconds for detrending before zero crossings are extracted. 
#' Default is `NULL` translating to same size as `win`. Setting it to \code{0} translates to no centring.
#' @param winSmooth Numeric. Seconds for smoothing to remove noise based on majority presence - only if relevant ZCrossings are 
#' present in more than half of this time window, this translates to true. Default is \code{0} translating to no smoothing.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns the dataframe and saves file to `rs.path` if provided.
#' 
#' @references Hale et al. (2020). Journal of Nonverbal Behavior.
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @seealso \code{\link{aggSlide}}, \code{\link{findZCrossing}}
#' @import dplyr
#' @export

extractZCrossing = function(df, colnames, fps, minDegree, rs.path = c(), suffix = "", 
                            win = 2, minFreq = 1.5, maxFreq = 6.5, 
                            winCentre = NULL, winSmooth = 0, 
                            verbose = T, recompute = F, return = T) {
  
  if (verbose) cat("------------------ Extracting Z Crossing Features ------------------\n")
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename
    flnm  = file.path(rs.path, sprintf("dataZC%s.arrow", suffix))
  }
  
  # adjust winCentre if necessary
  if (is.null(winCentre)) winCentre = win
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%X"), ": Loading Zero Crossing features\n")
      df = arrow::read_feather(flnm)
    }
  } else {
    
    if (verbose) cat(format(Sys.time(), "%X"), ": Exctracting Zero Crossings from ", paste(colnames, collapse = ", "), "\n")
    
    checkDF(df, c("Dyad", "Identifier", "Time", "Frame", "Timestamp", colnames))
    
    # focus on relevant columns 
    df = df |>
      select(Dyad, Identifier, Frame, Time, Timestamp, 
             any_of(c("Speaking", "Listening", "Communication", colnames)))
    
    # check whether detrending
    if (winCentre > 0) {
      df = df |>
        group_by(Dyad, Identifier, Time) |> arrange(Dyad, Identifier, Frame) |>
        mutate(
          # detrend data with local mean (based on winCentre)
          across(.cols = all_of(colnames), .fns = list(centred   = ~ .x - aggSlide(.x, mean, fps * winCentre)),
                 .names = "{.col}_{.fn}")
        )
      colnames = paste0(colnames, "_centred")
    }
    
    # further preprocess the colnames
    df = df |>
      group_by(Dyad, Identifier, Time) |> arrange(Dyad, Identifier, Frame) |>
      mutate(
        # compute zero-crossings and downstream filtering on data
        across(
          .cols = all_of(colnames),
          .fns = list(
            # extract the zero crossings
            zc     = ~ findZCrossing(.x),
            # find zero crossings with difference above minDegree
            zcT    = ~ findZCrossing(.x) & (abs(.x - lag(.x, default = .x[1])) >= minDegree)
            # get the difference in rotation with 0 for first entry
            diff   = ~ abs(.x - lag(.x, default = .x[1]))
            # sum up the above threshold ZCrossings across the window and divide by window for frequency
            sum    = ~ aggSlide(findZCrossing(.x) & (abs(.x - lag(.x, default = .x[1])) >= minDegree), 
                                sum, fps * win) / win,
            # compare sum / window to the minimum and maximum frequencies
            rel    = ~ 
              # larger than the minimum frequency
              (aggSlide(findZCrossing(.x) & (abs(.x - lag(.x, default = .x[1])) >= minDegree), 
                        sum, fps * win) / win > minFreq ) & 
              # smaller than the maximum frequency
              (aggSlide(findZCrossing(.x) & (abs(.x - lag(.x, default = .x[1])) >= minDegree), 
                        sum, fps * win) / win < maxFreq )
          ),
          .names = "{.col}_{.fn}"
        )
      ) |> ungroup()
    
    # check if smoothing
    if (winSmooth > 0) {
      df = df |>
        mutate(
          across(
            .cols = ends_with("_rel"),
            .fns = list(smooth = ~ aggSlide(.x, sum, fps * winSmooth) > (fps * winSmooth/2)), 
            .names = "{.col}_{.fn}"
          )
        )
    }
    
    # save the data frame
    if (!is.null(rs.path)) {
      if (verbose) cat(format(Sys.time(), "%X"), ": Saving the data\n")
      arrow::write_feather(df, flnm, compression = "zstd")
    }
    
  }
  if (verbose) cat(format(Sys.time(), "%X"), ": Done\n")
  
  if (return) return(df)
  
}

#' Agreggate Zero-Crossing Features Extracted from Time Series Data
#'
#' Aggregates the results from \code{\link{extractZCrossing}} or \code{\link{featHeadGestures}} to provide one 
#' absolute and one relative value per Identifier per time series. 
#'
#' @param df Dataframe. The dataset containing the variables to be processed, created by \code{\link{extractZCrossing}}. 
#'   Must explicitly feature columns `Dyad`, `Identifier`, `Frame`, `Time`, all columns contained in `colnames`. 
#'   If `Communication` is a column, zero crossings are also aggregated based on its classification. 
#' @param colnames Character vector. The exact name or names of the column(s) in \code{df} from which
#'   to extract zero-crossing features. 
#' @param rs.path Path to the directory where the output file will be saved, if empty (is.null(rs.path) == TRUE), 
#'   then nothing is saved. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns the dataframe or saves consolidated summary CSV to `rs.path` if provided.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @seealso \code{\link{extractZCrossing}}
#' @import dplyr
#' @export

featZCrossing = function(df, colnames, rs.path = c(), suffix = "",
                        verbose = T, recompute = F, return = T) {
  
  if (verbose) cat("----------------- Aggregating Z Crossing Features  -----------------\n")
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename
    flnm  = file.path(rs.path, sprintf("featZC%s.csv", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%X"), ": Loading Zero Crossing features\n")
      df.out = readr::read_csv(flnm, show_col_types = F)
    }
  } else {
    
    # aggregate the ZC information
    if (verbose) cat(format(Sys.time(), "%X"), ": Aggregating information\n")
    checkDF(df, c("Dyad", "Identifier", "Time", "Frame", colnames))
    
    # overall ZC information
    df.out = df |> 
      group_by(Dyad, Identifier, Time) |>
      mutate(
        # get the total number of frames
        Frames.total = n()
      ) |> 
      # get rid of NAs
      drop_na(all_of(colnames)) |>
      # aggregate the information
      group_by(Dyad, Identifier, Time, Frames.total) |>
      summarise(
        across(matches(colnames), sum),
        .groups = "drop"
      ) |> ungroup() |>
      mutate(
        across(matches(colnames), ~ .x * 100/Frames.total, .names = "RelativeZC_{.col}")
      ) |> select(Dyad, Identifier, Time, Frames.total, matches(colnames))
    
    # potentially add the values depending on Communication
    if ("Communication" %in% colnames(df)) {
      df.out = merge(
        df.out, 
        df |> 
          group_by(Dyad, Identifier, Time) |>
          mutate(
            # get the total number of frames
            Frames.total = n()
          ) |>
          group_by(Dyad, Identifier, Time, Communication, Frames.total) |>
          summarise(
            across(matches(colnames), sum),
            .groups = "drop"
          ) |> ungroup() |>
          mutate(
            across(matches(colnames), ~ .x * 100/Frames.total, .names = "RelativeZC_{.col}")
          ) |> select(Dyad, Identifier, Time, Frames.total, Communication, matches(colnames)) |>
          tidyr::pivot_wider(names_from = Communication, values_from = matches(colnames),
                             names_glue = "{.value}_{Communication}")
      )
    }
    
    # save speech dwell dataframe
    # save the data frame
    if (!is.null(rs.path)) {
      if (verbose) cat(format(Sys.time(), "%X"), ": Saving the data\n")
      readr::write_csv(df.out, flnm)
    }
    
  }
  
  if (verbose) cat(format(Sys.time(), "%X"), ": Done\n")
  
  if (return) return(df.out)
  
}

#' Perform Sliding Window Aggregations
#'
#' A wrapper around `slider::slide_dbl` that matches MATLAB-style centre-left biased 
#' windowing rules to calculate rolling averages, sums or other statistical functions.
#'
#' @param x Numeric vector. The time series data to pass to the sliding window.
#' @param FUN Function. The function object to apply over the window (e.g., `mean`, `sum`, `sd`).
#' @param samples Numeric. The total number of frame samples inside the tracking window.
#'
#' @return A numeric vector of the same length as `x` containing the rolling calculation.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
aggSlide = function(x, FUN, samples) {
  
  # calculate exact offsets to match MATLAB centre-left bias
  lo = ceiling((samples + 1) / 2) - 1
  up = floor((samples + 1) / 2) - 1
  
  # pass FUN into slider pipeline and return the resulting vector
  return(slider::slide_dbl(x, FUN, .before = lo, .after = up, .complete = FALSE))
  
}

#' Find Zero-Crossing Instances in a Time Series
#'
#' This function identifies indices where a vector changes sign (zero-crossings). 
#' Exact zeros are treated as positive values to stay consistent with specific signal 
#' processing conventions.
#'
#' @param data Numeric vector. The time series data to analyse, arranged by time.
#'
#' @return A numeric vector of the same length as `data`, where a `1` indicates 
#'   a zero-crossing occurred *after* that frame index, and `0` indicates no crossing.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
findZCrossing = function(x) {
  
  # get sign of each data point (-1, 0, or 1)
  signum = sign(x)
  
  # convert exact zeros to 1 (positive sign)
  signum[signum == 0] = 1
  
  # use diff() to find where sign changes
  return(as.integer(c(diff(signum) != 0, 0)))
}