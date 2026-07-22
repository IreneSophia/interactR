#' Compute Wavelet Transformation for Time-course Data
#'
#' Uses the \code{\link{biwavelet::wt}} function to compute the wavelet transformation
#' for a given column. 
#'
#' @param df Dataframe. The dataset containing the variables to be processed. 
#'   Must explicitly feature columns `Dyad`, `Identifier`, `Time`, either `Frame` or `Timestamp` and the column `colname`. 
#' @param rs.path Character. Path to destination directory for saved files. If empty (is.null(rs.path) == TRUE), then nothing is saved.
#' @param colname Character The exact name in \code{df} from which to extract the data. 
#' @param fps Numeric. Frame processing rate frequency profile (frames per second) of the dataset.
#' @param order Numeric. Order for the wavelet transformation. Default is `6`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns a list with each Identifiers wavelet result. If provided, list is saved as rds to `rs.path`.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @import dplyr
#' @export

compWT = function(df, rs.path, colname, fps, order = 6, suffix = "", 
                  verbose = T, recompute = F, return = T) {
  
  # extract the individual Identifiers
  ls.IDs = unique(df$Identifier)
  nos    = length(ls.IDs)
  
  # output structure
  ls.out = list()
  
  # extract the timeseries from either Timestamp or Frame + fps
  if ("Timestamp" %in% colnames(df)) {
    df = df |> group_by(Dyad, Identifier, Time) |> 
      arrange(Timestamp) |>
      mutate(start = min(Timestamp), 
             Timecourse = Timestamp - start,
             Duration = Timecourse - lag(Timecourse)) |>
      select(-start)
  } else {
    df = df |> group_by(Dyad, Identifier, Time) |> 
      arrange(Frame) |>
      mutate(Timecourse = Frame/fps,
             Duration   = 1/fps)
  }
  
  # loop through the Identifiers
  for (i in 1:nos) {
    
    # focus on this Identifier
    df.sel = df |> 
      filter(Identifier == ls.IDs[i])
    
    # get the timestep
    dt = mean(df.sel$Duration, na.rm = T)
    
    # extract data into a vector
    data = df.sel |> pull(colname)
    
    # extract the time into a vector
    time = df.sel |> pull(Timecourse)
    
    # configuration for wavelet transformation
    S0 = 2 * dt # smallest scale, set here to Nyquist limit
    Dj = 1/12
    J1 = round(log2(((length(time) * 0.17) * 2 * dt) / S0) / Dj) # configure number of frequency scales (minus 1)
    
    # compute the wavelet transformation
    # comparable to MATLAB's wavelet() by Torrence & Compo
    wt = biwavelet::wt(cbind(time, data), 
                       pad = T, # padding with zeros to length of power of 2, faster and less cone of influence problems
                       dt  = dt,
                       dj  = Dj, # spacing / resolution between scales, 12 sub-octaves for fine frequency resolution
                       s0  = S0, 
                       J1  = J1, 
                       mother = "morlet", 
                       param = order)
    
    # store results
    ls.out[i] = list(wt)
    names(ls.out)[i] = paste0(df.sel$Dyad[1], "_", ls.IDs[i])
    
  }
  
  return(ls.out)
  
}