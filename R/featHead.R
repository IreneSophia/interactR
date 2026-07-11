#' Extract Zero-Crossing Frequency from Time Series Data
#'
#' Evaluates localised oscillatory characteristics within specified time series 
#' channels by quantifying zero-crossing counts over moving windows. This allows
#' for the extraction of nodding and head shaking from rotational head movement.
#' The function is based on the algorithm presented in Hale et al. (2020) and
#' shifts data dynamically against a running local mean, calculates zero-crossing 
#' frequencies, filters outcomes within a target frequency band and applies a 
#' smoothing threshold filter across the results. 
#' 
#' @details 
#' The calculation executes sequentially across one or more specified columns:
#' 1. **Baseline Centering:** Removes slow-moving baseline drift using a running mean ([aggSlide()]).
#' 2. **Frequency Mapping:** Counts sign inversions using [findZCrossing()] over an explicit window interval.
#' 3. **Bandpass Filtering:** Identifies sequences falling strictly within the bounds defined by `minFreq` and `maxFreq`.
#' 4. **Smoothing:** Stabilizes transient state transitions by enforcing a secondary duration validation window.
#'
#' Newly generated vectors are appended back to the source data frame utilizing a 
#' systematic naming template (e.g., `[column]_centred`, `[column]_sum`, 
#' `[column]_rel`, and `[column]_smooth`).
#'
#' @param df Dataframe. The dataset containing the variables to be processed. [!MISSING]
#' @param rs.path Character. Path to destination directory for saved files. If empty (is_empty(rs.path) == TRUE), then nothing is saved.
#' @param colnames Character vector. The exact name or names of the column(s) in \code{df} from which
#'   to extract zero-crossing features. 
#' @param fps Numeric. Frame processing rate frequency profile (frames per second) of the dataset.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param win Numeric. Window duration scale evaluated in seconds for the moving frequency summary. Default is \code{2}.
#' @param minFreq Numeric. The lower cutoff boundary of the targeted frequency band in Hz. Default is \code{1.5}.
#' @param maxFreq Numeric. The upper cutoff boundary of the targeted frequency band in Hz. Default is \code{6.5}.
#' @param winSmooth Numeric. Secondary window scale evaluated in seconds for state smoothing. Default is \code{5}.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns the dataframe or saves consolidated summary CSV and full RDS file to `rs.path` if provided.
#' 
#' @references Hale et al. (2020). Journal of Nonverbal Behavior.
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @seealso \code{\link{aggSlide}}, \code{\link{findZCrossing}}
#' @import tidyverse
#' @export

featZCrossing = function(df, rs.path, colnames, fps, suffix = "",
                         win = 2, minFreq = 1.5, maxFreq = 6.5, winSmooth = 5, 
                         verbose = T, recompute = F, return = T) {
  
  # check rs.path
  if (is_empty(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename
    flnm  = file.path(rs.path, sprintf("featZC%s.csv", suffix))
  }

  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading Zero Crossing features\n")
      df.out = read_csv(flnm, show_col_types = F)
    }
  } else {
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Exctracting Zero Crossings from ", paste(colnames, collapse = ", "), "\n")
  
    # preprocess the colnames
    df = df %>%
      select(Dyad, Identifier, Frame, any_of(c("Communication", colnames))) %>%
      group_by(Dyad, Identifier) %>% arrange(Frame) %>%
      mutate(
        # detrend data with local mean (1 second window)
        across(.cols = all_of(colnames), .fns = list(centred   = ~ .x - aggSlide(.x, mean, fps)),
               .names = "{.col}_{.fn}")
      ) %>%
      mutate(
        # compute zero-crossings and downstream filtering on centred data
        across(
          .cols = ends_with("_centred"),
          .fns = list(
            sum    = ~ aggSlide(findZCrossing(.x), sum, fps * win) / win,
            rel    = ~ (aggSlide(findZCrossing(.x), sum, fps * win) / win > minFreq * win) & 
              (aggSlide(findZCrossing(.x), sum, fps * win) / win <= maxFreq * win),
            smooth = ~ aggSlide(
              (aggSlide(findZCrossing(.x), sum, fps * win) / win > minFreq * win) & 
                (aggSlide(findZCrossing(.x), sum, fps * win) / win <= maxFreq * win), 
              sum, winSmooth * fps
            ) > (winSmooth * fps / win)
          ),
          .names = "{.col}_{.fn}"
        )
      )
    
    # save the data for plotting
    if (!is_empty(rs.path)) saveRDS(df, file = file.path(rs.path, sprintf("dataZC%s.rds", suffix)))
  
    # aggregate the ZC information
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Aggregating information\n")
    
    # overall ZC information
    df.agg = df %>% 
      group_by(Dyad, Identifier) %>%
      mutate(
        # get the total number of frames
        Frames.total = n()
      ) %>% 
      group_by(Dyad, Identifier, Frames.total) %>%
      summarise(
        across(matches(paste0(colnames, "_centred_smooth")), sum)
      ) %>% ungroup() %>%
      mutate(
        across(matches(colnames), ~ .x * 100/Frames.total, .names = "RelativeZC_{.col}")
      ) %>% select(Dyad, Identifier, matches("Relative.*centred_smooth")) %>%
      rename_with(~ gsub("_centred_smooth", "_Total", .x), .cols = ends_with("centred_smooth"))
    
    # potentially add the values depending on Communication
    if ("Communication" %in% colnames(df)) {
      df.agg = merge(
        df.agg, 
        df %>% 
          group_by(Dyad, Identifier) %>%
          mutate(
            # get the total number of frames
            Frames.total = n()
          ) %>%
          group_by(Dyad, Identifier, Communication, Frames.total) %>%
          summarise(
            across(matches(paste0(colnames, "_centred_smooth")), sum)
          ) %>% ungroup() %>%
          mutate(
            across(matches(colnames), ~ .x * 100/Frames.total, .names = "RelativeZC_{.col}")
          ) %>% select(Dyad, Identifier, Communication, matches("Relative.*centred_smooth")) %>%
          rename_with(~ gsub("_centred_smooth", "", .x), .cols = ends_with("centred_smooth")) %>%
          pivot_wider(names_from = Communication, values_from = matches(colnames),
                      names_glue = "{.value}_{Communication}")
      )
    }
    
    # save speech dwell dataframe
    if (!is_empty(rs.path)) write_csv(df.out, flnm)
    
  }
  
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

#' Correct for Geometric Circularity in Rotational Time Series
#'
#' A tracking data helper function that identifies and removes artificial phase-wrapping 
#' jumps (e.g., flipping between -180 and +180 or 0 and 360 degrees) across consecutive frames.
#' Accumulates correction factors forward via cumulative sums. Inspired by the methods in Hale et al. (2020).
#'
#' @note **Caution:** The input vector `x` must be sorted chronologically (i.e., in the correct order of Frame / Timecourse).
#'
#' @param x Numeric vector. Uncorrected rotational data column entries ordered sequentially by time or frame.
#' @param th Numeric. Boundary threshold angle used to detect an artificial phase jump. Default is `270`.
#'
#' @return A numeric vector of the same length as `x` corrected for wrapping discontinuities.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export

fixCirc = function(x, th = 270) {
  
  # do nothing, if x only has one entry
  if (length(x) <= 1) return(x)
  
  dz = diff(x)
  
  # initialise a vector of zeros for the shifts
  shifts = rep(0, length(x))
  
  # identify jumps and assign the +360 or -360 correction factor
  up   = which(dz > th)
  down = which(dz < -th)
  
  # *subsequent* indices get shifted
  shifts[up   + 1] = -360
  shifts[down + 1] = +360
  
  # accumulate the shifts forward using cumsum
  return(x + cumsum(shifts))
}

#' Calculate Circular Distance Between Two Rotational Time Series
#'
#' A helper function that computes the shortest angular difference between two rotational 
#' positions (e.g., calculating independent head rotation relative to a baseline neck position) 
#' while correcting for modulo 360-degree wrapping.
#'
#' @param x Numeric vector. Target rotational position timecourse (e.g., head).
#' @param y Numeric vector. Baseline or relative comparison rotational position timecourse (e.g., neck).
#'
#' @return A numeric vector representing the wrapped angular difference, bounded strictly between -180 and +180 degrees.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export

rotDiff = function(x, y) {
  return((((x - y) + 180) %% 360) - 180)
}

#' Preprocess Head Motion Coordinates
#'
#' Preprocesses spatial head trajectory streams. Rotational paths are scrubbed of wrap-around geometric circularity 
#' or baseline-adjusted using corresponding body segments. Circularity corrected and translational trajectories are detrended.
#'
#' @param df Dataframe containing head coordinate streams. Requires tracking coordinates alongside core columns `Dyad`, `Identifier` and `Frame`.
#' @param rs.path Character. Path to the directory where the output files will be saved.
#'   If empty (`is_empty(rs.path) == TRUE`), nothing is saved to disk.
#' @param rotnames Character vector. String labels identifying target rotational columns.
#' @param tranames Character vector. String labels identifying target translational columns.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param performFixCirc Logical. Directs execution to correct boundaries for angular circular values. Default is `TRUE`.
#' @param cornames Character vector. Optional mapping labels matched in length to `rotnames` for baseline adjustments 
#'   (e.g., using spine tracking vectors to separate head posture from whole-body sway). Default is `c()`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `FALSE`.
#'
#' @return If `return = TRUE`, returns an adjusted movement dataframe. Saves `dataHead[suffix].rds` to `rs.path` if provided.
#' 
#' @references Hale et al. (2020). Journal of Nonverbal Behavior.
#' @import tidyverse
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 

preproHead = function(df, rs.path, rotnames, tranames, suffix = '',
                      performFixCirc = T, cornames = c(),
                      verbose = T, recompute = F, return = F) {
  
  # check rs.path
  if (is_empty(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = flrds = ''
  } else {
    # create filenames
    flnm = file.path(rs.path, sprintf("dataHead%s.rds", suffix))
  }
  
  # give some info
  if (verbose) cat("----------- Preprocess head motion data -----------\n")
  if (file.exists(flnm) & !recompute) {
    df = readRDS(flnm)
  } else {
    # focus on relevant columns
    df = df %>% 
      select(Dyad, Identifier, Frame, 
             any_of(c("Speaking", "Listening", "Communication")),
             any_of(c(rotnames, tranames, cornames))) %>%
      arrange(Dyad, Identifier, Frame)
    
    if (performFixCirc) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Fixing circularity\n")
      df = df %>%
        group_by(Dyad, Identifier) %>%
        # fix circularity based on the algorithm of Hale et al. (2020),
        # default threshold is 270
        mutate(across(all_of(rotnames), fixCirc, .names = "{.col}_fixCirc")) %>%
        ungroup()
    }
    if (length(cornames) == length(rotnames)) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Rotational difference with to other body part\n")
      # compute rotational difference to a different body part
      for (i in 1:length(cornames)) {
        new_name = paste0(rotnames[i], "_rotDiff")
        df[[new_name]] = rotDiff(df[[rotnames[i]]], df[[cornames[i]]])
      }
    }
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Detrend translational columns\n")
    fixnames = paste0(rotnames, "_fixCirc")
    df = df %>%
      # de-trended translational and fixCirc values by subtracting mean value
      group_by(Dyad, Identifier) %>%
      mutate(across(all_of(c(tranames, fixnames)), ~ .x - mean(.x), .names = "{.col}_detrended")) %>%
      ungroup() %>% arrange(Dyad, Identifier, Frame)
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Save data\n")
    if (!is_empty(rs.path)) saveRDS(df, file = flnm)
    
  }
  
  if (return) return(df)

}
