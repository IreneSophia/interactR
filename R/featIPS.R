#' Compares pseudo and observed IPS for each Lag (WLCC) or each Frequency (WTC)
#'
#' Takes the dataframe created by [extractIPS()] and compares the pseudo and the 
#' observed IPS values based on the mean. Comparison can be either made by 
#' aggregating the values across Dyads and then comparing the pseudo and the 
#' observed z-transformed values for each lag or Frequency with a t-test to assess. 
#' Frequentist or Bayesian t-tests can be chosen with the `Bayesian` logical switch. 
#' Alternatively, if the sample of dyads is small, observed IPS can be aggregated 
#' across all dyads and compared with a larger sample of pseudo IPS in a 
#' permutation test by setting `perm = TRUE`. When frequentist t-tests are computed, 
#' multiple comparison correction is applied based on Benjamini & Hochberg's (1995) FDR. 
#'
#' @param df.observed Dataframe. Dataframe created by [extractIPS()] extracting observed values.
#' @param df.pseudo Dataframe. Dataframe created by [extractIPS()] extracting pseudo values. 
#' @param Bayesian Logical. Switch to perform Bayesian or frequentist testing. Default is `TRUE`.
#' @param perm Logical. Switch to use permutation testing instead of comparison across Dyads.
#'   If the dataset is smaller, then one can create a larger set of pseudo values 
#'   against which the mean observed value per Lag can be compared. Default is `FALSE`.
#' @param minBF Numeric. Threshold above which log Bayes Factor is considered credible evidence. Default is `log(3)`.
#' @param alpha Numeric. Threshold above which permutation and Frequentist is considered significant. Default is `0.05`.
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

compareIPS = function(df, Bayesian = T, perm = F, 
                      minBF = log(3), alpha = 0.05, freqLimits = c(0.2, 8),
                      rs.path = c(), suffix = "", 
                      verbose = T, recompute = F, return = T) {
  
  if (verbose) cat("---------------- Comparing pseudo and observed IPS  ----------------\n")
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename 
    flnm  = file.path(rs.path, 
                      sprintf("featIPS_pseudo-comp%s.csv", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading IPS comparison\n")
      df.agg = readr::read_csv(flnm, show_col_types = F)
    }
  } else {
    
    # check that df only contains one type of pseudo WLCC and must contain observed values
    if ((length(unique(df$Method)) != 2) | !("observed" %in% unique(df$Method))) stop("Function performs comparison of observed with one pseudo method.")
    
    # rename the main column depending on whether WTC or WLCC
    if ("Rsq" %in% colnames(df)) {
      df = df |>
        # if WTC, then filter out the frequencies outside of the limits
        filter(Frequency >= freqLimits[1] & Frequency <= freqLimits[2]) |>
        rename(value = Rsq) |>
        select(value, Frame, Frequency, Dyad, Method, Feature, iteration)
      colnm = "Frequency"
    } else {
      df = df |>
        rename(value = WLCC) |>
        select(value, Lag, Dyad, Method, Feature, iteration)
      colnm = "Lag"
    }
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Comparing pseudo and observed IPS from", unique(df$Feature), "\n")
    
    # aggregate the values across the windows
    df.tmp = df |>
      group_by(across(all_of(colnm)), Dyad, Method, Feature) |>
      summarise(
        value = mean(value, na.rm = T),
        .groups = "drop"
      )
    
    # aggregate values for plotting
    df.agg = df.tmp |>
      group_by(across(all_of(colnm)), Method, Feature) |>
      summarise(
        STD = sd(value),
        AVG = mean(value),
        .groups = "drop"
      )
    
    if (!perm) {
      
      if (!Bayesian) {
        
        if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Computing frequentist stats\n")

                # use t-tests to assess differences
        df.stat = df.tmp |> 
          # z-transform to achieve normal distribution
          mutate(zvalue = atanh(pmin(pmax(value, -0.9999), 0.9999))) |>
          group_by(Feature, across(all_of(colnm))) |> 
          rstatix::t_test(zvalue ~ Method, detailed = T, p.adjust.method = "BH") |> 
          mutate(
            Evaluation = if_else(p < alpha, "*", ""), 
            Direction  = if_else(estimate1 > estimate2, "greater", "lesser")
          )
        
        # add to the aggregated dataframe
        df.agg = df.agg |>
          left_join(df.stat |> select(Feature, all_of(colnm), p, Evaluation, Direction),
                    by = c(colnm, "Feature"))
        
      } else {
        
        if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Computing Bayesian stats\n")
        
        # compute the Bayesian t-tests
        df.stat = df.tmp |>
          # z-transform to achieve normal distribution
          mutate(zvalue = atanh(pmin(pmax(value, -0.9999), 0.9999))) |>
          group_by(Feature, across(all_of(colnm))) |>
          group_modify(~ {
            observed = .x |>
              filter(Method == "observed") |>
              summarise(mean = mean(value, na.rm = TRUE)) |> pull(mean)
            pseudo = .x |>
              filter(Method != "observed") |>
              summarise(mean = mean(value, na.rm = TRUE)) |> pull(mean)
            ttest = BayesFactor::ttestBF(formula = zvalue ~ Method, data = as.data.frame(.x))
            data.frame(logBF = ttest@bayesFactor$bf, Direction = if_else(observed > pseudo, "greater", "lesser"))
          }) |>
          ungroup() |>
          mutate(
            Evaluation = if_else(logBF > minBF, "*", "")
          )
        
        # add to the aggregated dataframe
        df.agg = df.agg |>
          left_join(df.stat |> select(Feature, all_of(colnm), logBF, Evaluation, Direction),
                    by = c(colnm, "Feature"))
        
      }
      
    } else {
      
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Computing stats with permutation\n")
      
      # aggregate the observed values to get grand average per Lag or Frequency
      df.perm = df |> filter(Method == "observed") |>
        group_by(across(all_of(colnm)), Dyad, Feature) |>
        summarise(
          value = mean(value, na.rm = T),
          .groups = "drop"
        ) |> 
        group_by(across(all_of(colnm)), Feature) |>
        summarise(
          observed = mean(value, na.rm = T),
          .groups = "drop"
        ) |> 
        # merge with the aggregated pseudo dataframe (across windows)
        right_join(df |> filter(Method != "observed") |>
                     group_by(across(all_of(colnm)), Dyad, Feature, iteration) |>
                     summarise(
                       pseudo = mean(value, na.rm = T),
                       .groups = "drop"
                     ), by = c(colnm, "Feature"))
      
      # calculate permutation values
      df.stat = df.perm |>
        group_by(Feature, across(all_of(colnm))) |>
        summarise(
          Probability = max(c(mean(observed > pseudo), mean(pseudo > observed))),
          Direction   = if_else(mean(observed > pseudo) > mean(pseudo > observed),
                                "greater", "lesser"),
          .groups = "drop"
        ) |> group_by(Feature) |> 
        mutate(
          Evaluation = if_else(Probability > (1 - alpha), "*", "")
        )
      
      # add the result to the aggregated dataframe
      df.agg = df.agg |>
        left_join(df.stat, by = c(colnm, "Feature"))
      
    }
    
    # save the data
    if (!is.null(rs.path)) {
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Saving the data\n")
      readr::write_csv(df.agg, file = flnm)
    }
    
  }
  
  if (return) return(df.agg)
  
  if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Done\n")
  
}

#' Compute Wavelet Coherence for Time-course Data
#'
#' Uses the \code{\link{biwavelet::wtc}} function to compute the Wavelet Coherence
#' and Phase for one dyad. 
#'
#' @details Relative phase differences are calculated as Phi_left - Phi_right. Thus, if Phase values
#'   are positive, this indicates the left Identifier was leading, while negative values indicate that
#'   the right Identifier was leading (based on the Dyad ID). Values closer to 0 indicate in-phase 
#'   and values closer to +- pi indicate anti-phase. No Monte Carlo randomisations are computed. 
#'
#' @param df Dataframe. The dataset containing the variables to be processed. 
#'   Must explicitly feature columns `Dyad`, `Frame`, `left` and `right`. 
#'   For each Dyad, there must be exactly two Identifiers in the data. 
#' @param fps Numeric. Frame processing rate frequency profile (frames per second) of the dataset.
#' @param featname Character. Name of the Feature from which WTC is computed. 
#' @param order Numeric. Order for the wavelet transformation. Default is `8` based on Issartel et al. (2006).
#' @param rs.path Character. Path to destination directory for saved files. 
#'   If empty (is.null(rs.path) == TRUE), then nothing is saved. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param seed Numeric. Seed supporting reproducibility. Takes an integer seed or `NA` for random. Default is `NA`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns dataframe with the results. If provided, the wtc object is saved as rds to `rs.path`.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @references Issartel et al. (2006): A Practical Guide to Time—Frequency Analysis 
#'    in the Study of Human Motor Behavior: The Contribution of Wavelet Transform, 
#'    Journal of Motor Behavior, 38(2), 139-159.
#' @import dplyr
#' @export
#' 
computeWTC = function(df, fps, featname,
                      order = 8, rs.path = c(), suffix = "", seed = NA, 
                      verbose = T, recompute = F, return = T) {
  
  # get a random seed
  if (is.na(seed)) {
    seed = sample(1000:9999, 1)
  }
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename including the dyad
    flnm  = file.path(rs.path, sprintf("dataWTC_%s_%s_seed-%d%s.rds", 
                                       featname, df$Dyad[1], seed, suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading WTC of ", df$Dyad[1], "\n")
      wtc = readRDS(flnm)
    }
  } else {
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Starting WTC for dyad ", df$Dyad[1], "\n")
  
    # check if the correct columns
    checkDF(df, c("Dyad", "Frame", "left", "right"))
    
    # check that this data only contains one Dyad with the same Frames
    if (length(unique(df$Dyad)) != 1) stop("computeWTC() should have data of one Dyad as input")
    
    # check for NAs
    if (any(is.na(df))) stop("Input to computeWTC() cannot contain NAs")
    
    # extract the timeseries from Frame + fps
    df = df |> 
      arrange(Frame) |>
      mutate(Timecourse = Frame/fps,
             Duration   = 1/fps) |>
      ungroup()
    
    # extract data into a vector
    data1 = df |> arrange(Timecourse) |> pull(left)
    data2 = df |> arrange(Timecourse) |> pull(right)
    
    # extract the time into a vector
    time = df |> pull(Timecourse)
    
    # configuration for wavelet transformation
    dt = 1/fps
    S0 = 2 * dt # smallest scale, set here to Nyquist limit
    Dj = 1/12
    J1 = round(log2(((length(time) * 0.17) * 2 * dt) / S0) / Dj) # configure number of frequency scales (minus 1)
    
    # try to compute the wavelet coherence 
    wtc = tryCatch(
      {
        biwavelet::wtc(cbind(time, data1), cbind(time, data2), 
                       pad = T,  # padding with zeros to length of power of 2, faster and less cone of influence problems
                       dj  = Dj, # spacing / resolution between scales, 12 sub-octaves for fine frequency resolution
                       s0  = S0, 
                       J1  = J1, 
                       nrands = 0, 
                       mother = "morlet", 
                       param = order, 
                       quiet = TRUE)
      },
      error = function(e) {
        message(paste("Skipping Dyad", df$Dyad[1], "due to error:", e$message))
        return(e$message)
      }
    )
    
    # save the wtc object
    if (!is.null(rs.path)) saveRDS(wtc, flnm)
  }
  
  # convert to a dataframe and return it
  df.out = convertWTC(wtc, df$Dyad[1])
  
  return(df.out)
  
}

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
      startFrame = Start - lagSample,
      endFrame   = Start + winSample + lagSample,
      # interpret the lag as to who was acting before whom
      Leading = case_when(
        Lag < 0 ~ "Right",
        Lag > 0 ~ "Left",
        T ~ "Simultaneous"
      ),
      # add the Dyad
      Dyad = df$Dyad[1]
    ) |> select(-Start)
  
  return(df.wlcc)
  
}

#' Extract Observed or Pseudo Interpersonal Synchrony (IPS)
#'
#' This function computes IPS either using Wavelect Coherence (WTC, [computeWTC()]) or using
#' Windowed-Lagged Cross-Correlation (WLCC, [computeWLCC()]). It allows for a `df.pseudo`
#' to be provided to compute pseudo IPS based on dyad shuffling as well. Furthermore, 
#' one can compute pseudo IPS based on segment shuffling by setting nSegment > 0. 
#' 
#' @details
#' If pseudo IPS is computed, this should be either through setting pseudoSegment to TRUE
#' or by using a df.pseudo to provide shuffled dyads, but not both. 
#' 
#'
#' @param df Dataframe. Must contain the columns `Dyad`, `Time`, `Identifier`, 
#'   `Frame` as well as the column described by `colname`.
#' @param colname Character. Name of the column from which IPS is extracted. 
#' @param type Character. Name of the method to compute IPS: either "WTC" for 
#'   wavelet coherence or "WLCC" for windowed lagged cross correlation. 
#' @param fps Numeric. Frame processing rate frequency profile (frames per second) of the dataset.
#' @param featname Character. Name of the Feature. If `NA`, then the colname is used. Default is `NA`. 
#' @param rs.path Character. Path to destination directory for saved files. 
#'   If empty (is.null(rs.path) == TRUE), then nothing is saved. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param cores Numeric. How many cores to use for parallelisation. 
#'   Default is `NA`, translating to as many as possible.
#' @param winSample Numeric. Window size in samples for WLCC. Default is `NA`.
#' @param incSample Numeric. Window increment step in samples for WLCC. Default is `NA`.
#' @param lagSample Numeric. Evaluated cross-correlation lag step in samples for WLCC. Default is `NA`.
#' @param order Numeric. Order for the wavelet transformation. Default is `8` based on Issartel et al. (2006).
#' @param nSegment Numeric. If larger than 0, then instead of extracting observed IPS, 
#'   windows are shuffled to created pseudo IPS. Number of shuffling to be performed.
#'   Default is `100`.
#' @param df.pseudo Dataframe. If it contains rows, then instead of using the observed 
#'    dyads listed in `df`, these pairings are tested for pseudo IPS. Dataframe can be 
#'    created using [shuffleIdentifier()] or [shuffleDyads()]. Default is `data.frame()`.
#' @param seed Numeric. Seed supporting reproducibility. Takes an integer seed or `NA` for random. Default is `NA`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If return is TRUE, then IPS dataframe. If rs.path is provided, then the dataframe is also saved there.
#' 
#' @import dplyr
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 

extractIPS = function(df, colname, type, fps, featname = NA,
                      rs.path = c(), suffix = "", cores = NA,
                      winSample = NA, incSample = NA, lagSample = NA, # settings for WLCC 
                      order = 8,                                      # settings for WTC
                      nSegment = 0, df.pseudo = data.frame(), seed = NA,
                      verbose = T, recompute = F, return = T) {

  if (verbose) cat(sprintf("---------------------- Extract IPS using %-4s ----------------------\n", type))
  
  # get a random seed
  if (is.na(seed)) {
    seed = sample(1000:9999, 1)
  }
  
  # potentially get feature name
  if (is.na(featname)) featname = colname
  
  # get whether this is pseudo or not
  if (nrow(df.pseudo) > 0) pseudoDyad = T else pseudoDyad = F
  if (nSegment > 0) pseudoSegment = T else pseudoSegment = F
  
  if (pseudoDyad & pseudoSegment) {
    stop("Only shuffle either by Dyad or by Segment.")
  }
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename depending on whether this is pseudo or not
    if (pseudoDyad) {
      flnm  = file.path(rs.path, sprintf("data%s_%s_seed-%d-pseudoDyad%s.arrow", 
                                         type, featname, seed, suffix))
    } else if (pseudoSegment) {
      flnm  = file.path(rs.path, sprintf("data%s_%s_seed-%d-pseudoSegment%s.arrow", 
                                         type, featname, seed, suffix))
    } else {
      flnm  = file.path(rs.path, sprintf("data%s_%s_seed-%d%s.arrow", 
                                         type, featname, seed, suffix))
    }
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading ",type ," features\n")
      df.out = arrow::read_feather(flnm)
    }
  } else {
  
    # check which Method is computed
    if (pseudoDyad) {
      Method = "pseudoDyad"
    } else if (pseudoSegment) {
      Method = "pseudoSegment"
    } else {
      Method = "observed"
    }
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Extracting ", type, Method, " from ", colname, "\n")
    
    # check whether settings for WLCC are complete, if needed
    if (type == "WLCC") {
      if (any(is.na(c(incSample, lagSample, winSample)))) {
        stop("To use WLCC, incSample, lagSample and winSample need to be determined.")
      }
    } else if (type == "WTC") {
      # for segment shuffling in WTC, one second is used
      winSample = fps
    } else {
          stop("type must be either WTC or WLCC.")
    }
    
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
    
    # focus on the relevant columns
    df = df |>
      select(all_of(cols))
    
    # set up the parallelisation
    if (is.na(cores)) {
      future::plan(future::multisession, workers = future::availableCores() - 1)
    } else {
      future::plan(future::multisession, workers = cores)
    } 
    
    # if pseudoSegment, then increase the rows of df.dyad to reflect the iterations
    if (pseudoSegment) {
      df.dyad = df.dyad |>
        slice(rep(1:n(), each = nSegment)) |>
        mutate(k = rep(1:nSegment, length.out = n()))
    } else {
      df.dyad = df.dyad |>
        mutate(k = 1)
    }
    
    # loop through the df.dyad - parallelised
    ls.out = furrr::future_map(
      1:nrow(df.dyad), 
      \(dyadRow) {
        if (verbose & (dyadRow%%10 == 1)) cat(format(Sys.time(), "%x %X %Z"), ": Starting with dyad ", dyadRow, " of ", nrow(df.dyad), "\n")
        
        # check whether all information here
        if (sum(is.na(df.dyad[dyadRow,])) > 0) stop("Dyad ", df.dyad$Dyad[dyadRow], " is missing crucial information (Dyad, Time or Identifier)")
        
        # extract data into a vector
        df.sel = rbind(
          df |> filter(Identifier == df.dyad$left_Identifier[dyadRow]   & Time == df.dyad$left_Time[dyadRow])   |> arrange(Frame) |> mutate(Side = "left"),
          df |> filter(Identifier == df.dyad$right_Identifier[dyadRow]  & Time == df.dyad$right_Time[dyadRow])  |> arrange(Frame) |> mutate(Side = "right")) |>
          select(-Identifier, -Time) |>
          mutate(Dyad = df.dyad$Dyad[dyadRow]) |>
          tidyr::pivot_wider(names_from = Side, values_from = all_of(colname)) |>
          # drop NA to ensure that we have the same Frames in both interaction partners
          tidyr::drop_na()
        
        # check if a Frame is missing in-between
        if (length(setdiff(min(df.sel$Frame):max(df.sel$Frame), df.sel$Frame)) != 0) {
          warning("Missing Frames in dyad ", df.dyad$Dyad[dyadRow])
          # return an empty dataframe for this one
          return(data.frame())
        }
        
        if (pseudoSegment) {
          # shuffle the data of the right person
          df.sel = df.sel |>
            mutate(
              right = unlist(sample(split(right, ceiling(seq_along(right) / winSample))), 
                             use.names = FALSE)
            )
        } 
        # compute IPS
        if (type == "WLCC") {
          minFrame = min(df.sel$Frame)
          computeWLCC(df.sel, winSample, incSample, lagSample) |>
            mutate(Method = Method, Feature = featname, Time = df.dyad$left_Time[dyadRow],
                   seed = seed, iteration = df.dyad$k[dyadRow],
                   # adjust the WLCC indices with the Frame
                   startFrame = startFrame + minFrame - 1, endFrame = endFrame + minFrame - 1)
        } else if (type == "WTC") {
          if (pseudoSegment) {
            # do not save the iterations
            computeWTC(df.sel, fps, featname, order = order, rs.path = c(), suffix = suffix, 
                       seed = seed, verbose = F, recompute = recompute, return = T) |>
              mutate(Method = Method, Feature = featname, Time = df.dyad$left_Time[dyadRow],
                     seed = seed, iteration = df.dyad$k[dyadRow])
          } else {
            # save the wtc output
            computeWTC(df.sel, fps, featname, order = order, rs.path = rs.path, suffix = suffix, 
                       seed = seed, verbose = F, recompute = recompute, return = T) |>
              mutate(Method = Method, Feature = featname, Time = df.dyad$left_Time[dyadRow],
                     seed = seed, iteration = df.dyad$k[dyadRow])
          }
        }
      },
      .options = furrr::furrr_options(seed = seed)
    )
    
    # clean up parallelisation
    future::plan(future::sequential)
    invisible(gc())
    
    # unpack the data into a dataframe
    df.out = purrr::list_rbind(ls.out)
    
    # if WLCC, then divide by fps
    if (type == "WLCC") {
      df.out = df.out |>
        mutate(
          Lag = Lag / fps
        )
    }
    
    # save the data
    if (!is.null(rs.path)) {
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Saving the data\n")
      arrow::write_feather(df.out, flnm, compression = "zstd")
    }
    
  }
  
  if (return) return(df.out)
  
  if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Done\n")
  
}

#' Convert a wavelet coherence object into a dataframe
#'
#' Take an object created by [computeWTC()] and transforms it into a dataframe, 
#' allowing easy aggregation and adding plotting options.
#'
#' @details Relative phase difference was calculated as Phi_left - Phi_right. Thus, if Phase values
#'   are positive, this indicates the left Identifier was leading, while negative values indicate that
#'   the right Identifier was leading (based on the Dyad ID). Values closer to 0 indicate in-phase 
#'   and values closer to +- pi indicate anti-phase.
#'
#' @param wtc Varying class.If wtc is an object created by \code{\link{biwavelet::wtc}}, 
#'   its data is converted into a dataframe. If not, potentially because wtc failed, 
#'   then a dataframe with the same columns is created with NAs in all but the `Dyad` column.
#' @param Dyad Character. Name of the Dyad of which WTC had been computed or attempted to be computed. 
#' @param withinCOI Logical. Whether only values inside the COI should be included. Default is `TRUE`.
#'
#' @return Returns dataframe.
#' 
#' @seealso \code{\link{featIPS}} \code{\link{computeWTC}}
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @import dplyr
#' @export
#' 

convertWTC = function(wtc, Dyad, withinCOI = T) {
  
  # check if WTC was successful
  if (class(wtc) != "biwavelet") {
    warning("Skipping ", names(ls)[i], ": not biwavelet.") 
    # return a dataframe with NAs
    Period = Rsq = Phase = Timecourse = COI = Frame = Frequency = WithinCOI = NA
    pseudoDyad = grepl("|", Dyad, fixed = T)
    return(data.frame(Period, Rsq, Phase, Timecourse, COI, Frame, Frequency, WithinCOI,
                      Dyad, pseudoDyad))
  }
  # convert wavelet coherence list into dataframe
  df.out = 
    data.frame(wtc[["rsq"]]) |>
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
      Dyad = Dyad,
      pseudoDyad = grepl("|", Dyad, fixed = T)
    )
  
  # potentially filter within COI
  if (withinCOI) df.out = df.out |> filter(WithinCOI)
  
  return(df.out)
  
}