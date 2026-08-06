#' Extract Head Gestures Extracted from Time Series Data
#'
#' Extracts nodding and head shaking by detecting Zero Crossing using \code{\link{featZCrossing}}. 
#' Gestures are classified exclusively as nodding, head shaking or nothing for each Frame.
#' Input can be data preprocessed using \code{\link{preproHead}}.
#'
#' @param df Dataframe. The dataset containing the variables to be processed, potentially created by \code{\link{preproHead}}. 
#'   Must explicitly feature columns `Dyad`, `Identifier`, `Time`, `Frame` and the columns `colNodding` and `colShaking`. This dataframe
#'   will be processed using \code{\link{featHeadGestures}} to extract nodding and head shaking, assuming that only one can happen at a time. 
#' @param colNodding Character. The exact name of the column in \code{df} from which
#'   to extract, and then plot, nodding. 
#' @param colShaking Character. The exact name of the column in \code{df} from which
#'   to extract, and then plot, head shaking. 
#' @param fps Numeric. Frame processing rate frequency profile (frames per second) of the dataset.
#' @param minDegree Numeric. How many degree of rotational difference are needed for the movement to be considered relevant. 
#'   Depends on the fps and the specific movement. Setting to negative number leads to no thresholding based on degrees. 
#' @param rs.path Character. Path to destination directory for saved files. 
#'   If empty (is.null(rs.path) == TRUE), then nothing is saved. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param win Numeric. Window duration scale evaluated in seconds for the moving frequency summary. Default is \code{2}.
#' @param minFreq Numeric. The lower cutoff boundary of the targeted frequency band in Hz. Default is \code{1.5}.
#' @param maxFreq Numeric. The upper cutoff boundary of the targeted frequency band in Hz. Default is \code{6.5}.
#' @param winCentre Numeric. Seconds for detrending before zero crossings are extracted. 
#' Default is `NULL` translating to same size as `win`. Setting it to \code{0} translates to no centring.
#' @param winSmooth Numeric. Seconds for smoothing to remove noise based on majority presence - only if head gesture is 
#' present in more than half of this time window, this translates to true. Default is \code{0} translating to no smoothing.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns the dataframe or saves RDS file to `rs.path` if provided.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @seealso \code{\link{featZCrossing}} \code{\link{preproHead}}
#' @import dplyr
#' @export
#' 
featHeadGestures = function(df, colNodding, colShaking, fps, minDegree, 
                            rs.path = c(), suffix = "", 
                            win = 2, minFreq = 1.5, maxFreq = 6.5, 
                            winCentre = NULL, winSmooth = 0, 
                            verbose = T, recompute = F, return = T) {
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename
    flnm  = file.path(rs.path, sprintf("dataHeadGestures%s.rds", suffix))
  }
  
  # adjust winCentre if necessary
  if (is.null(winCentre)) winCentre = win
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading Head Gesture features\n")
      df = readRDS(flnm)
    }
  } else {
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Exctracting Head Gestures from ", paste(c(colNodding, colShaking), collapse = ", "), "\n")
    
    checkDF(df, c("Dyad", "Identifier", "Frame", "Time", colNodding, colShaking))
    
    # process the dataframe to extract Zero Crossings
    df = featZCrossing(df, c(), c(colNodding, colShaking), fps, minDegree,
                       win = win, minFreq = minFreq, maxFreq = maxFreq, 
                       winCentre = winCentre, winSmooth = winSmooth, verbose = F)
    
    # if centring was used, adjust the colnames
    if (winCentre > 0) {
      colNodding = paste0(colNodding, "_centred")
      colShaking = paste0(colShaking, "_centred")
    }
    
    # preprocess the extracted z crossings
    df = df |>
      # rename the columns to nodding and shaking
      rename_with(~ gsub(colNodding, "nodding", .x), .cols = matches(colNodding)) |>
      rename_with(~ gsub(colShaking, "shaking", .x), .cols = matches(colShaking)) |>
      # if both are relevant, then use the one with the larger frame-wise difference
      mutate(
        shaking_rel = case_when(
          shaking_rel & nodding_rel & shaking_diff >  nodding_diff ~ TRUE,
          shaking_rel & nodding_rel & shaking_diff <= nodding_diff ~ FALSE,
          T ~ shaking_rel
        ),
        nodding_rel = case_when(
          shaking_rel & nodding_rel & shaking_diff >  nodding_diff ~ FALSE,
          shaking_rel & nodding_rel & shaking_diff <= nodding_diff ~ TRUE,
          T ~ nodding_rel
        )
      )  
    
    # save the data for plotting
    if (!is.null(rs.path)) saveRDS(df, file = flnm)
    
  }
  
  if (return) return(df)
  
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
#' @param df Dataframe containing head coordinate streams. Requires tracking coordinates alongside 
#' core columns `Dyad`, `Identifier`, `Time`, `Timestamp` and `Frame` as well as all columns
#' specified to be processed.
#' @param rotnames Character vector. String labels identifying target rotational columns.
#' @param tranames Character vector. String labels identifying target translational columns.
#' @param rs.path Character. Path to the directory where the output files will be saved.
#'   If empty (`is.null(rs.path) == TRUE`), nothing is saved to disk. Default is `c()`.
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
#' @import dplyr
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 

preproHead = function(df, rotnames, tranames, rs.path = c(), suffix = '',
                      performFixCirc = T, cornames = c(),
                      verbose = T, recompute = F, return = F) {
  
  checkDF(df, c("Dyad", "Identifier", "Frame", "Time", "Timestamp", 
                rotnames, tranames, cornames))
  
  # check rs.path
  if (is.null(rs.path)) {
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
    df = df |> 
      select(Dyad, Identifier, Frame, Timestamp, Time,
             any_of(c("Speaking", "Listening", "Communication")),
             any_of(c(rotnames, tranames, cornames))) |>
      arrange(Dyad, Identifier, Frame)
    
    if (performFixCirc) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Fixing circularity\n")
      df = df |>
        group_by(Dyad, Identifier, Time) |>
        # fix circularity based on the algorithm of Hale et al. (2020),
        # default threshold is 270
        mutate(across(all_of(rotnames), fixCirc, .names = "{.col}_fixCirc")) |>
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
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Detrend non-difference columns\n")
    fixnames = paste0(rotnames, "_fixCirc")
    df = df |>
      # de-trended translational and fixCirc values by subtracting mean value
      group_by(Dyad, Identifier, Time) |>
      mutate(across(all_of(c(tranames, fixnames)), ~ .x - mean(.x), .names = "{.col}_detrended")) |>
      ungroup() |> arrange(Dyad, Identifier, Frame)
    
    if (!is.null(rs.path)) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Save data\n")
      saveRDS(df, file = flnm)
      }
    
  }
  
  if (return) return(df)

}
