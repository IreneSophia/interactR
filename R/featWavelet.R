#' Compute Wavelet Coherence for Time-course Data
#'
#' Uses the \code{\link{biwavelet::wtc}} function to compute the wavelet coherence
#' for a given column for a dyad. Optionally creates pseudo coherence based 
#' on Dyad shuffling. Observed and pseudo coherence can be compared using `[!MISSING]`
#'
#' @param df Dataframe. The dataset containing the variables to be processed. 
#'   Must explicitly feature columns `Dyad`, `Identifier`, `Time`, either `Frame` or `Timestamp`, and the column `colname`. 
#'   For each Dyad, there must be exactly two Identifiers in the data. 
#' @param rs.path Character. Path to destination directory for saved files. If empty (is.null(rs.path) == TRUE), then nothing is saved.
#' @param colname Character The exact name in \code{df} from which to extract the data. 
#' @param fps Numeric. Frame processing rate frequency profile (frames per second) of the dataset.
#' @param order Numeric. Order for the wavelet transformation. Default is `8` based on Issartel et al. (2006).
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param nsim. Numeric. Number of Monte Carlo randomisations for computing WTC. Default is `400`.
#' @param pseudoDyad Logical. Flags whether to generate pseudo-WTC benchmarks using dyad shuffling 
#'   instead of observed WTC Default is `FALSE`.
#' @param nDyad Numeric. Total number of synthetic dyad simulations to execute when `pseudoDyad = TRUE`. 
#'   Default is `NULL` which is converted into the number of real dyads.
#' @param seed Character or Numeric. Seed supporting reproducibility. Takes an integer seed or `"random"`. Default is `"random"`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns dataframe with the results. If provided, the dataframe and the full list are saved as rds to `rs.path`.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @references Issartel et al. (2006): A Practical Guide to Time—Frequency Analysis 
#'    in the Study of Human Motor Behavior: The Contribution of Wavelet Transform, 
#'    Journal of Motor Behavior, 38(2), 139-159.
#'    Veleda, D., Montagne, R., & Araujo, M. (2012). Cross-wavelet bias corrected 
#'    by normalizing scales. Journal of Atmospheric and Oceanic Technology, 29(9), 1401-1408.
#' @import dplyr
#' @export

compWTC = function(df, rs.path, colname, fps, order = 8, suffix = "", 
                   nsim = 400, pseudoDyad = F, nDyads = NULL, seed = "random",
                   verbose = T, recompute = F, return = T) {
  
  # get a random seed
  if (!is.numeric(seed)) {
    seed = sample(1000:9999, 1)
  }
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename depending on whether this is pseudo or not
    if (pseudoDyad) {
      flnm  = file.path(rs.path, sprintf("dataWTC_%s_seed-%d_pseudo%s.rds", 
                                         colname, seed, suffix))
    } else {
      flnm  = file.path(rs.path, sprintf("dataWTC_%s_seed-%d%s.rds", 
                                         colname, seed, suffix))
    }
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading WTC features\n")
      df.out = readRDS(flnm)
    }
  } else {
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Exctracting WTC features from ", colname, "\n")
    
    # check whether all columns in dataframe
    cols = c("Dyad", "Identifier", "Time", colname)
    checkDF(df, cols)
    
    # check if either Frame or Timestamp
    if (!("Frame" %in% colnames(df)) & !("Timestamp" %in% colnames(df))) {
      stop("Dataframe needs to contain either Timestamp or Frame")
    }
    
    # if pseudoDyad, then create a random list of combinations
    if (pseudoDyad) {
      # check if nDyad needs to be set
      if (is.null(nDyad)) nDyad = length(unique(df$Dyad))
      # get a list of shuffled dyads: either out of all options
      if (nDyads > length(unique(df$Dyad))) {
        # randomly draw from all possible options if large nDyads
        df.dyad = shuffleDyads(df |> select(Dyad, Time, Identifier) |> distinct(),
                     seed = seed, nsim = nDyad)
      } else {
        # shuffle the right Identifier
        df.dyad = shuffleIdentifier(df |> select(Dyad, Time, Identifier) |> distinct(), 
                                    seed = seed, side = "right", nsim = nDyad)
      }
    } else {
      # get a list of real dyads
      df.dyad = df |>
        select(Dyad, Time, Identifier) |> 
        distinct() |>
        mutate(
          side  = if_else(gsub("(.+)-.*", "\\1", Dyad) == Identifier, "left", "right")
        )
    }
    
    # set the seed
    set.seed(seed)
    
    # output structure
    ls.out = list()
    df.out = data.frame()
    
    # focus on the relevant columns
    df = df |>
      select(all_of(cols), any_of(c("Frame", "Timestamp")))
    
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
    
    # loop through the Dyads
    ls.dyads = unique(df.dyad$Dyad)
    for (i in 1:length(ls.dyads)) {
      
      # get the info for this dyad
      ID.left = df.dyad |> 
        filter(Dyad == ls.dyads[i] & side == "left") |> pull(Identifier)
      t.left  = df.dyad |> 
        filter(Dyad == ls.dyads[i] & side == "left") |> pull(Time)
      ID.right = df.dyad |> 
        filter(Dyad == ls.dyads[i] & side == "right") |> pull(Identifier)
      t.right  = df.dyad |> 
        filter(Dyad == ls.dyads[i] & side == "right") |> pull(Time)
      
      # check if two valid IDs
      if ((nchar(ID.right) == 0) | (nchar(ID.left) == 0)) stop("This function requires exactly two Identifiers per Dyad.")
      # check if two valid times
      if ((nchar(t.right) == 0) | (nchar(t.left) == 0)) stop("One Identifier is missing a valid Time.")
      
      # check if time is the same or if it is pseudo - both cases take left time
      if ((t.left == t.right) | pseudoDyad) {
        t = t.left
      } else {
        stop("Different Time without dyad shuffling!")
      }
      
      # get the timestep
      dt = mean(df |> filter(Dyad == ls.dyads[i]) |> pull(Duration), na.rm = T)
      
      # extract data into a vector
      data1 = df |> filter(Identifier == ID.left  & Time == t.left)  |> arrange(Timecourse) |> pull(colname)
      data2 = df |> filter(Identifier == ID.right & Time == t.right) |> arrange(Timecourse) |> pull(colname)
      
      # extract the time into a vector
      time1 = df |> filter(Identifier == ID.left  & Time == t.left)  |> arrange(Timecourse) |> pull(Timecourse)
      time2 = df |> filter(Identifier == ID.right & Time == t.right) |> arrange(Timecourse) |> pull(Timecourse)
      
      # configuration for wavelet transformation
      S0 = 2 * dt # smallest scale, set here to Nyquist limit
      Dj = 1/12
      J1 = round(log2(((length(time1) * 0.17) * 2 * dt) / S0) / Dj) # configure number of frequency scales (minus 1)
      
      # compute the wavelet coherence
      wtc = biwavelet::wtc(cbind(time1, data1), cbind(time2, data2), 
                           pad = T, # padding with zeros to length of power of 2, faster and less cone of influence problems
                           dj  = Dj, # spacing / resolution between scales, 12 sub-octaves for fine frequency resolution
                           s0  = S0, 
                           J1  = J1, 
                           nrands = nsim, 
                           mother = "morlet", 
                           param = order)
      
      # convert wavelet coherence into dataframe
      df.out = data.frame(wtc[["rsq"]]) |>
        mutate(Period = wtc[["period"]]) |>
        tidyr::pivot_longer(cols = starts_with("X"), values_to = "Rsq") |>
        merge(data.frame(wtc[["phase"]]) |>
                mutate(Period = wtc[["period"]]) |>
                tidyr::pivot_longer(cols = starts_with("X"), values_to = "Phase")) |>
        arrange(name) |>
        mutate(
          Timecourse = rep(wtc[["t"]], each = length(wtc[["period"]])),
          COI   = rep(wtc[["coi"]], each = length(wtc[["period"]])),
          Frame = as.numeric(gsub("X", "", name)),
          Frequency = 1 / Period,
          WithinCOI = Period < COI 
        ) |> select(-name) |>
        mutate(
          Dyad = ls.dyads[i],
          Time = t,
          pseudoDyad = pseudoDyad
        )
      
    }
    
    # save the data
    if (!is.null(rs.path)) saveRDS(df.out, file = flnm)
    
  }
  
  if (return) return(df.out)
  
}