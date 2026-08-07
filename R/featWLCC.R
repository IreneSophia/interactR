#' Compute Windowed-Lagged Cross-Correlation (WLCC) of two Timecourses
#'
#' Computes windowed-lagged cross-correlations based on the standard `cor` function. 
#' In a Dyad named \[Identifier\]-\[Identifier\], left and right refer to the the Identifier
#' before and after the hyphen, respectively. 
#' 
#' @details
#' This function selects windows so that each lag of each window has the exact same
#' amount of samples to compute the correlation. Thus, the first window starts one 
#' lag from the first index, to allow for one window to be shifted in this direction. 
#' This approach results in some datapoints, at the beginning and the end, not being
#' used as often as other datapoints.
#'
#' @param df Dataframe. Containing the columns `Dyad`, `Frame`, `left` and `right`.
#' @param winSample Numeric. Window size in samples.
#' @param incSample Numeric. Window increment step in samples.
#' @param lagSample Numeric. Evaluated cross-correlation lag in samples. Lag is applied
#'   in both directions, i.e., a lag of 10 samples results in +- 10 samples as lags.
#'
#' @return WLCC dataframe.
#' 
#' @import dplyr
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 
computeWLCC = function(df, winSample, incSample, lagSample) {
  
  # check the dataframe
  checkDF(df, c("Dyad", "Frame", "left", "right"))
  
  # check that it contains the data of one Dyad
  if (length(unique(df$Dyad)) != 1) stop("Dataframe must contain data of one Dyad in columns left and right.")
  
  # ensure it is properly sorted
  df = df |> arrange(Frame)
  
  # extract the data 
  minFrame = min(df$Frame)
  left  = df$left
  right = df$right
  
  # create cross correlation dataframe
  df.wlcc = tidyr::expand_grid(
    # get all starting indices, considering the increment given
    Start = seq(1+lagSample, length(left) - winSample - lagSample + 1, by = incSample), 
    # get all the defined lags
    Lag   = (-lagSample):(lagSample)) |>
    mutate(
      # map two values of this row to a function: Start and Lag
      WLCC = purrr::map2_dbl(Start, Lag, ~ {
        s = .x                  # start index
        e = .x + winSample - 1  # end index
        l = .y                  # lag
        # compute the cross correlation
        cor(left[s:e], right[(s+l):(e+l)], use = "complete.obs")
      })
    ) |>
    # add the window number
    group_by(Lag) |>
    mutate(
      Window = if_else(Lag == min(Lag), paste0("w", row_number()), NA)
    ) |> ungroup() |> tidyr::fill(Window, .direction = "down") |>
    mutate(
      # add the start and end index of each window
      winStart = Start - lagSample,
      winEnd   = Start + winSample + lagSample,
      # interpret the lag as to who was acting before whom
      Direction = case_when(
        Lag < 0 ~ "rightLeading",
        Lag > 0 ~ "leftLeading",
        T ~ "Simultaneous"
      ),
      # add the Dyad
      Dyad = df$Dyad[1]
    ) |> select(-Start)
  
  return(df.wlcc)
  
}

#' Extract Observed or Pseudo Windowed-Lagged Cross-Correlation (WLCC)
#'
#' !ADD DETAILS!
#' 
#' @details
#' If pseudo WLCC is computed, this should be either through setting pseudoSegment to TRUE
#' or by using a df.pseudo to provide shuffled dyads, but not both. 
#' 
#'
#' @param df Dataframe. Must contain the columns `Dyad`, `Time`, `Identifier`, 
#'   `Frame` as well as the column described by `colname`.
#' @param winSample Numeric. Window size in samples.
#' @param incSample Numeric. Window increment step in samples.
#' @param lagSample Numeric. Evaluated cross-correlation lag step in samples.
#' @param colname Character. Name of the column from which WLCC is extracted. 
#' @param rs.path Character. Path to destination directory for saved files. 
#'   If empty (is.null(rs.path) == TRUE), then nothing is saved. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param pseudoSegment Boolean. Switch to create pseudo WLCC instead by shuffling
#'   the windows of the right Identifier before computing WLCC. Default is `FALSE`.
#' @param nSegment Numeric. Number of shuffling to be performed if pseudoSegment is TRUE.
#'   Default is `100`.
#' @param df.pseudo Dataframe. If it contains rows, then instead of using the observed 
#'    dyads listed in `df`, these pairings are tested for pseudo-WTC. Dataframe can be 
#'    created using [shuffleIdentifier()] or [shuffleDyads()]. Default is `data.frame()`.
#' @param seed Numeric. Seed supporting reproducibility. Takes an integer seed or `NA` for random. Default is `NA`.
#'    Only used, if pseudoSegment is TRUE, otherwise there is no random process executed in this function.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return WLCC dataframe.
#' 
#' @import dplyr
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 

extractWLCC = function(df, winSample, incSample, lagSample, colname, 
                       rs.path = c(), suffix = "", 
                       pseudoSegment = F, nSegment = 100, 
                       df.pseudo = data.frame(), seed = NA,
                       verbose = T, recompute = F, return = T) {
  
  # get a random seed
  if (is.na(seed)) {
    seed = sample(1000:9999, 1)
  }
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
    # get whether this is pseudo or not
    if (nrow(df.pseudo) > 0) pseudoDyad = T else pseudoDyad = F
  } else {
    # create filename depending on whether this is pseudo or not
    if (nrow(df.pseudo) > 0) {
      flnm  = file.path(rs.path, sprintf("dataWLCC_%s_seed-%d-pseudo%s.rds", 
                                         colname, seed, suffix))
      pseudoDyad = T
    } else {
      flnm  = file.path(rs.path, sprintf("dataWLCC_%s_seed-%d%s.rds", 
                                         colname, seed, suffix))
      pseudoDyad = F
    }
  }
  
  # check which type is computed
  if (pseudoDyad & pseudoSegment) {
    stop("Only shuffle either by Dyad or by Segment.")
  } else if (pseudoDyad) {
    Type = "pseudoDyad"
  } else if (pseudoSegment) {
    Type = "pseudoSegment"
  } else {
    Type = "observed"
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading WLCC features\n")
      df.wlcc = readRDS(flnm)
    }
  } else {
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Extracting WLCC features from ", colname, "\n")
    
    # check whether all columns in dataframe
    cols = c("Dyad", "Identifier", "Time", "Frame", colname)
    checkDF(df, cols)
    
    # if pseudoDyad, then use the provided list of combinations
    if (pseudoDyad) {
      df.dyad = df.pseudo
    } else {
      # get a list of real dyads - one row per dyad same as with pseudo
      df.dyad = df |> ungroup() |>
        select(Dyad, Time) |> 
        distinct() |>
        mutate(
          left_Identifier  = gsub("(.+)-.*", "\\1", Dyad),
          right_Identifier  = gsub(".*-(.+)", "\\1", Dyad),
          left_Time = Time
        ) |> rename(right_Time = Time)
    }
    
    # set the seed
    set.seed(seed)
    
    # output structure
    df.wlcc = data.frame()
    
    # focus on the relevant columns
    df = df |>
      select(all_of(cols))
    
    # loop through the Dyads
    for (i in 1:nrow(df.dyad)) {
      
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Starting with WLCC from dyad ", i, " of ", nrow(df.dyad), "\n")
      
      # check whether all information here
      if (sum(is.na(df.dyad[i,])) > 0) stop("Dyad ", df.dyad$Dyad[i], " is missing crucial information (Dyad, Time or Identifier)")
      
      # extract data into a vector
      df.sel = rbind(
        df |> filter(Identifier == df.dyad$left_Identifier[i]   & Time == df.dyad$left_Time[i])   |> arrange(Frame) |> mutate(Side = "left"),
        df |> filter(Identifier == df.dyad$right_Identifier[i]  & Time == df.dyad$right_Time[i])  |> arrange(Frame) |> mutate(Side = "right")) |>
        select(-Identifier) |>
        tidyr::pivot_wider(names_from = Side, values_from = all_of(colname))
      
      if (pseudoSegment) {
        # loop through the number of shuffles
        for (j in 1:nSegment) {
          # shuffle the data of the right person
          df.sel = df.sel |>
            mutate(
              right = unlist(sample(split(right, ceiling(seq_along(right) / winSample))), 
                             use.names = FALSE)
              )
          # compute the WLCC
          df.wlcc = rbind(df.wlcc, 
                          computeWLCC(df.sel, winSample, incSample, lagSample) |>
                            mutate(Type = Type, seed = seed, iteration = j))
        }
        
      } else {
        # compute the WLCC
        df.wlcc = rbind(df.wlcc, 
                        computeWLCC(df.sel, winSample, incSample, lagSample) |>
                          mutate(Type = Type, seed = NA, iteration = NA))
      }
      
    }
    
    # save the data
    if (!is.null(rs.path)) saveRDS(df.wlcc, file = flnm)
    
  }
  
  if (return) return(df.wlcc)
  
}