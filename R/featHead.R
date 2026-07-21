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
#'   If empty (`is.null(rs.path) == TRUE`), nothing is saved to disk.
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
#' @import dplyr
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 

preproHead = function(df, rs.path, rotnames, tranames, suffix = '',
                      performFixCirc = T, cornames = c(),
                      verbose = T, recompute = F, return = F) {
  
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
        group_by(Dyad, Identifier) |>
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
      group_by(Dyad, Identifier) |>
      mutate(across(all_of(c(tranames, fixnames)), ~ .x - mean(.x), .names = "{.col}_detrended")) |>
      ungroup() |> arrange(Dyad, Identifier, Frame)
    
    if (!is.null(rs.path)) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Save data\n")
      saveRDS(df, file = flnm)
      }
    
  }
  
  if (return) return(df)

}
