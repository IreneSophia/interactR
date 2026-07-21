#' Extract Head Gestures Extracted from Time Series Data
#'
#' Extracts nodding and head shaking by detecting Zero Crossing using \code{\link{featZCrossing}}. 
#' Gestures are classified exclusively as nodding, head shaking or nothing for each Frame.
#'
#' @param df Dataframe. The dataset containing the variables to be processed, potentially created by \code{\link{preproHead}}. 
#'   Must explicitly feature columns `Identifier`, `Frame` and the columns `colNodding` and `colShaking`. This dataframe
#'   will be processed using \code{\link{featHeadGestures}} to extract nodding and head shaking, assuming that only one can happen at a time. 
#' @param rs.path Character. Path to destination directory for saved files. If empty (is.null(rs.path) == TRUE), then nothing is saved.
#' @param colNodding Character. The exact name of the column in \code{df} from which
#'   to extract, and then plot, nodding. 
#' @param colShaking Character. The exact name of the column in \code{df} from which
#'   to extract, and then plot, head shaking. 
#' @param fps Numeric. Frame processing rate frequency profile (frames per second) of the dataset.
#' @param minDegree Numeric. How many degree of rotational difference are needed for the movement to be considered relevant. 
#' Depends on the fps and the specific movement. Setting to negative number leads to no thresholding based on degrees. 
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param win Numeric. Window duration scale evaluated in seconds for the moving frequency summary. Default is \code{2}.
#' @param minFreq Numeric. The lower cutoff boundary of the targeted frequency band in Hz. Default is \code{1.5}.
#' @param maxFreq Numeric. The upper cutoff boundary of the targeted frequency band in Hz. Default is \code{6.5}.
#' @param winCentre Numeric. Seconds for detrending before zero crossings are extracted. Default is \code{0} translating to no detrending
#' @param winSmooth Numeric. Seconds for state smoothing. Default is \code{0} translating to no smoothing.
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
featHeadGestures = function(df, rs.path, colNodding, colShaking, fps, minDegree, 
                            suffix = "", 
                            win = 2, minFreq = 1.5, maxFreq = 6.5, 
                            winCentre = 0, winSmooth = 0, 
                            verbose = T, recompute = F, return = T) {
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename
    flnm  = file.path(rs.path, sprintf("dataHeadGestures%s.rds", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading Head Gesture features\n")
      df = readRDS(flnm)
    }
  } else {
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Exctracting Head Gestures from ", paste(c(colNodding, colShaking), collapse = ", "), "\n")
    
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
