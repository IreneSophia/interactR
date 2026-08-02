
#' Compares pseudo and observed WTC within given Frequencies
#'
#' Takes the dataframe created by [aggWTC()] and compares the pseudo and the
#' observed Rsq values. 
#'
#' @param df Dataframe. Dataframe created by [aggWTC()].
#' @param rs.path Character. Path to destination directory for saved files. If empty (is.null(rs.path) == TRUE), then nothing is saved.
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

compareWTC = function(df.rsq, rs.path, suffix = "", 
                      verbose = T, recompute = F, return = T) {
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename 
    flnm  = file.path(rs.path, 
                      sprintf("featWTC_pseudo-comp%s.csv", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading WTC comparison\n")
      df.rsq = read_csv(flnm)
    }
  } else {
  
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Comparing pseudo and observed WTC from ", unique(df.rsq$Feature), "\n")
    
    # use non-parametric tests to assess differences
    df.stat = df.rsq |> 
      group_by(Bin) |> 
      rstatix::wilcox_test(Rsq_avg ~ pseudoDyad, detailed = T) |> 
      rstatix::adjust_pvalue(method = "BH") |>
      mutate(
        Sig  = if_else(p.adj < 0.05, "*", ""),
        Type = case_when(estimate > 0 & p.adj < 0.05 ~ "coherence",
                         estimate < 0 & p.adj < 0.05 ~ "hypo-coherence",
                         T ~ "")
      )
    
    # merge with the df.rsq
    df.rsq = merge(df.stat |> select(Bin, statistic, p, p.adj, estimate, conf.low, conf.high, Sig, Type),
                   df.rsq)
    
    # save the data
    if (!is.null(rs.path)) write_csv(df.rsq, file = flnm)
    
  }
  
  if (return) return(df.rsq)
  
}

#' Aggregate a dataframe containing the results of WTC to compare pseudo and observed WTC
#'
#' Takes the dataframe created by [convertWTC()] and aggregates the results per Frequency, 
#' within given limits of a minimum and maximum Frequency. Values outside COI can be excluded. 
#'
#' @param df Dataframe. Dataframe created by [convertWTC()].
#' @param rs.path Character. Path to destination directory for saved files. If empty (is.null(rs.path) == TRUE), then nothing is saved.
#' @param minFreq Numeric. Minimum frequency which is included. 
#' @param maxFreq Numeric. Maximum frequency which is included. 
#' @param withinCOI Logical. Whether only values inside the COI should be included. Default is `TRUE`.
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

aggWTC = function(df, rs.path, minFreq, maxFreq, withinCOI = T, suffix = "", 
                  verbose = T, recompute = F, return = T) {
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename 
    flnm  = file.path(rs.path, 
                      sprintf("dataWTC_pseudo-comp-agg%s.rds", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading aggregated WTC for comparison\n")
      df.rsq = readRDS(flnm)
    }
  } else {
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Aggregating pseudo and observed WTC from ", unique(df.rsq$Feature), "\n")
    
    # filter the data if necessary
    if (withinCOI) df = df |> filter(WithinCOI)
    
    # filter based on the Limits
    df = df |>
      filter(Frequency >= minFreq & Frequency <= maxFreq)
    
    # aggregate the data
    df.rsq = df |>
      group_by(Period, Frequency, Dyad, Time, pseudoDyad, Feature) |>
      summarise(
        Rsq_avg = mean(Rsq), 
        Rsq_std = sd(Rsq),
        .groups = "drop"
      ) |>
      mutate(
        Bin = as.numeric(as.factor(Frequency))
      )
  
    # save the data
    if (!is.null(rs.path)) saveRDS(df.rsq, file = flnm)
    
  }
  
  if (return) return(df.rsq)
  
}

#' Compute Wavelet Coherence for Time-course Data
#'
#' Uses the \code{\link{biwavelet::wtc}} function to compute the Wavelet Coherence
#' and Phase for a given column for a dyad. Optionally creates pseudo coherence based 
#' on Dyad shuffling. Resulting list can be transformed to dataframe using [convertWTC()].
#' Observed and pseudo coherence can be aggregated using [featWTC()].
#'
#' @details Relative phase differences are calculated as Phi_left - Phi_right. Thus, if Phase values
#'   are positive, this indicates the left Identifier was leading, while negative values indicate that
#'   the right Identifier was leading (based on the Dyad ID). Values closer to 0 indicate in-phase 
#'   and values closer to +- pi indicate anti-phase.
#'
#' @param df Dataframe. The dataset containing the variables to be processed. 
#'   Must explicitly feature columns `Dyad`, `Identifier`, `Time`, either `Frame` or `Timestamp`, and the column `colname`. 
#'   For each Dyad, there must be exactly two Identifiers in the data. 
#' @param rs.path Character. Path to destination directory for saved files. If empty (is.null(rs.path) == TRUE), then nothing is saved.
#' @param colname Character. The exact name in \code{df} from which to extract the data. 
#' @param fps Numeric. Frame processing rate frequency profile (frames per second) of the dataset.
#' @param order Numeric. Order for the wavelet transformation. Default is `8` based on Issartel et al. (2006).
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param nsim. Numeric. Number of Monte Carlo randomisations for computing WTC. Default is `300`.
#' @param df.pseudo Dataframe. If it contains row, then instead of using the observed 
#'    dyads listed in `df`, these pairings are tested for pseudo-WTC. Dataframe can be 
#'    created using [shuffleIdentifier()] or [shuffleDyads()]. Default is an empty dataframe, 
#'    leading to observed and not pseudo WTC being extracted.
#' @param seed Character or Numeric. Seed supporting reproducibility. Takes an integer seed or `"random"`. Default is `"random"`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns list with the results. If provided, the full list is saved as rds to `rs.path`.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @seealso \code{\link{convertWTC}} \code{\link{featWTC}}
#' @references Issartel et al. (2006): A Practical Guide to Time—Frequency Analysis 
#'    in the Study of Human Motor Behavior: The Contribution of Wavelet Transform, 
#'    Journal of Motor Behavior, 38(2), 139-159.
#' @import dplyr
#' @export
#' 
extractWTC = function(df, rs.path, colname, fps, order = 8, suffix = "", 
                      nsim = 300, df.pseudo = data.frame(), seed = "random",
                      verbose = T, recompute = F, return = T) {
  
  # get a random seed
  if (!is.numeric(seed)) {
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
      flnm  = file.path(rs.path, sprintf("dataWTC_%s_seed-%d-pseudo%s.rds", 
                                         colname, seed, suffix))
      pseudoDyad = T
    } else {
      flnm  = file.path(rs.path, sprintf("dataWTC_%s_seed-%d%s.rds", 
                                         colname, seed, suffix))
      pseudoDyad = F
    }
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading WTC features\n")
      ls.out = readRDS(flnm)
    }
  } else {
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Extracting WTC features from ", colname, "\n")
    
    # check whether all columns in dataframe
    cols = c("Dyad", "Identifier", "Time", "Frame", colname)
    checkDF(df, cols)
    
    # if pseudoDyad, then create a random list of combinations
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
    ls.out = vector(mode = "list", length = nrow(df.dyad))
    names(ls.out) = as.character(1:nrow(df.dyad))
    
    # focus on the relevant columns
    df = df |>
      select(all_of(cols), any_of(c("Frame", "Timestamp")))
    
    # extract the timeseries from Frame + fps
    df = df |> group_by(Dyad, Identifier, Time) |> 
      arrange(Frame) |>
      mutate(Timecourse = Frame/fps,
             Duration   = 1/fps) |>
      ungroup()
    
    # loop through the Dyads
    for (i in 1:nrow(df.dyad)) {
      
      # check whether all information here
      if (sum(is.na(df.dyad[i,])) > 0) stop("Dyad ", df.dyad$Dyad[i], " is missing crucial information (Dyad, Time or Identifier)")
      
      # check if time is the same or if it is pseudo - both cases take left time
      if ((df.dyad$left_Time[i] == df.dyad$right_Time[i]) | pseudoDyad) {
        t = df.dyad$left_Time[i]
      } else {
        stop("Different Time without dyad shuffling!")
      }
      
      # extract data into a vector
      data1 = df |> filter(Identifier == df.dyad$left_Identifier[i]  & Time == df.dyad$left_Time[i])  |> arrange(Timecourse) |> pull(colname)
      data2 = df |> filter(Identifier == df.dyad$right_Identifier[i] & Time == df.dyad$right_Time[i]) |> arrange(Timecourse) |> pull(colname)
      
      # extract the time into a vector
      time1 = df |> filter(Identifier == df.dyad$left_Identifier[i]  & Time == df.dyad$left_Time[i])  |> arrange(Timecourse) |> pull(Timecourse)
      time2 = df |> filter(Identifier == df.dyad$right_Identifier[i] & Time == df.dyad$right_Time[i]) |> arrange(Timecourse) |> pull(Timecourse)
      
      # if pseudoDyad, there might be a difference in timepoints - cut the longer one
      if (pseudoDyad & length(time1) != length(time2)) {
        ls.data = list(data1 = data1, data2 = data2, 
                       time1 = time1, time2 = time2)
        ls.trim = lapply(ls.data, function(v) v[1:min(sapply(ls.data, length))])
        data1 = ls.trim$data1
        data2 = ls.trim$data2
        time1 = ls.trim$time1
        time2 = ls.trim$time2
      }
      
      # configuration for wavelet transformation
      dt = 1/fps
      S0 = 2 * dt # smallest scale, set here to Nyquist limit
      Dj = 1/12
      J1 = round(log2(((length(time1) * 0.17) * 2 * dt) / S0) / Dj) # configure number of frequency scales (minus 1)
      
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Starting WTC for dyad ", df.dyad$Dyad[i], "\n")
      
      # try to compute the wavelet coherence 
      wtc = tryCatch(
        {
          biwavelet::wtc(cbind(time1, data1), cbind(time2, data2), 
                         pad = T, # padding with zeros to length of power of 2, faster and less cone of influence problems
                         dj  = Dj, # spacing / resolution between scales, 12 sub-octaves for fine frequency resolution
                         s0  = S0, 
                         J1  = J1, 
                         nrands = nsim, 
                         mother = "morlet", 
                         param = order, 
                         quiet = !verbose)
        },
        error = function(e) {
          message(paste("Skipping Dyad", df.dyad$Dyad[i], "due to error:", e$message))
          return(e$message)
        }
      )
      
      # add to the output list
      ls.out[[i]] = wtc
      names(ls.out)[i] = paste0(df.dyad$Dyad[i], "_", as.character(t))
      
    }
    
    # save the data
    if (!is.null(rs.path)) saveRDS(ls.out, file = flnm)
    
  }
  
  if (return) return(ls.out)
  
}

#' Convert list of wavelet coherence objects into a dataframe
#'
#' Takes the list created by [extractWTC()] and transforms it into a dataframe, 
#' allowing easy aggregation and adding plotting options.
#'
#' @details Relative phase difference was calculated as Phi_left - Phi_right. Thus, if Phase values
#'   are positive, this indicates the left Identifier was leading, while negative values indicate that
#'   the right Identifier was leading (based on the Dyad ID). Values closer to 0 indicate in-phase 
#'   and values closer to +- pi indicate anti-phase.
#'
#' @param ls List. This list must contain objects created by \code{\link{biwavelet::wtc}}.
#'   Names of the list entries should follow this structure: `[Dyad]_[Time]`.
#' @param rs.path Character. Path to destination directory for saved files. If empty (is.null(rs.path) == TRUE), then nothing is saved.
#' @param featname Character. The name of the feature of which WTC was computed.
#' @param withinCOI Logical. Whether only values inside the COI should be included. Default is `TRUE`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns dataframe. If provided, the full dataframe is saved as rds to `rs.path`.
#' 
#' @seealso \code{\link{featWTC}} \code{\link{extractWTC}} \code{\link{plotPhaseRose}} \code{\link{plotWTCLine}}
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @import dplyr
#' @export
#' 
convertWTC = function(ls, rs.path, featname, withinCOI = T,
                      suffix = "", verbose = T, recompute = F, return = T) {

  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename 
    flnm  = file.path(rs.path, sprintf("dataWTC_%s-df%s.rds", 
                                       featname, suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading WTC feature dataframe\n")
      df.out = readRDS(flnm)
    }
  } else {
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Converting WTC features from list to dataframe\n")
    
    # initialise dataframe
    df.out = data.frame()
    
    for (i in 1:length(ls)) {
      
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Processing", i, " of ", length(ls), "\n")
      # check whether the object in the list is the target - if not skip
      if (class(ls[[i]]) != "biwavelet") {
        warning("Skipping ", names(ls)[i], ": not biwavelet.") 
        next
      }
      # get Dyad and Time
      Time = as.POSIXct(gsub(".*_", "", names(ls)[i]))
      Dyad = gsub("_.*", "", names(ls)[i])
      # convert wavelet coherence list into dataframe
      df.out = rbind(
        df.out,
        data.frame(ls[[i]][["rsq"]]) |>
          mutate(Period = ls[[i]][["period"]]) |>
          tidyr::pivot_longer(cols = starts_with("X"), values_to = "Rsq") |>
          merge(data.frame(ls[[i]][["phase"]]) |>
                  mutate(Period = ls[[i]][["period"]]) |>
                  tidyr::pivot_longer(cols = starts_with("X"), values_to = "Phase")) |>
          arrange(name) |>
          mutate(
            Timecourse = rep(ls[[i]][["t"]], each = length(ls[[i]][["period"]])),
            COI   = rep(ls[[i]][["coi"]], each = length(ls[[i]][["period"]])),
            Frame = as.numeric(gsub("X", "", name)),
            Frequency = 1 / Period,
            WithinCOI = Period < COI
          ) |> select(-name) |>
          mutate(
            Dyad = Dyad,
            Time = Time,
            pseudoDyad = grepl("|", Dyad, fixed = T)
          )
      )
    }
    
    # add the featname
    df.out = df.out |> mutate(Feature = featname)
    
    # potentially filter within COI
    if (withinCOI) df.out = df.out |> filter(WithinCOI)
    
    # save the data
    if (!is.null(rs.path)) saveRDS(df.out, file = flnm)
    
  }
  
  if (return) return(df.out)
  
}

#' Aggregate a dataframe containing the results of WTC to extract relevant features 
#'
#' Takes the dataframe created by [convertWTC()] and aggregates the results in prespecified Bins
#' and across all Bins within the specified limits. Limits may differ between Rsq and Phase. 
#' Values outside COI can be excluded. 
#'
#' @details Relative phase difference was calculated as Phi_left - Phi_right. Thus, if Phase values
#'   are positive, this indicates the left Identifier was leading, while negative values indicate that
#'   the right Identifier was leading (based on the Dyad ID). Values closer to 0 indicate in-phase 
#'   and values closer to +- pi indicate anti-phase.
#'
#' @param df Dataframe. Dataframe created by [convertWTC()].
#' @param rs.path Character. Path to destination directory for saved files. If empty (is.null(rs.path) == TRUE), then nothing is saved.
#' @param phaseLimits Numeric. Frequency limits in Hz for the Phase Bins. All but the last value specify the included lower limit. 
#'   The last value is the excluded maximum period. Must include at least two values.
#' @param Limits List. List containing numeric of frequency limits in Hz.
#'   First entry is used for the Coherence Bins, optionally a second for Phase Bins. 
#'   All but the last value of each numeric specify the included lower limit. 
#'   The last value is the excluded maximum period. List must include at least one value, 
#'   each numeric must include at least two. 
#' @param Funs List. List of function to be used for the aggregation. First is used for Coherence. 
#'   Optionally, a second is used for Phase. If Limits contains two sets of Limits, but Funs only
#'   one function, then this function is used for both Phase and Coherence. 
#' @param labels. Logical. Whether to use labels showing the Bin limits. If false, then Bins are numbered. Default is `TRUE`.
#' @param withinCOI Logical. Whether only values inside the COI should be included. Default is `TRUE`.
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

featWTC = function(df, rs.path, Limits, Funs, labels = T,
                   withinCOI = T, suffix = "", 
                   verbose = T, recompute = F, return = T) {
  
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
  
  # if two Limits but only one Fun, then use the same Fun twice
  if ((length(Limits) == 2) & (length(Funs) == 1)) {
    Funs = c(Funs, Funs)
  }
  
  # create labels
  if (labels) rsqLabels = NULL else rsqLabels = FALSE
  rsqFun      = Funs[[1]]
  rsqLimits   = Limits[[1]]
  if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Extract coherence aggregates\n")
  
  # potentially exclude outside of COI
  if (withinCOI) df = df |> filter(WithinCOI)

  # classify into rsq bins
  df = df |>
    mutate(rsqBin =   cut(Frequency, 
                          breaks = rsqLimits, 
                          labels = rsqLabels, 
                          right = FALSE,
                          include.lowest = TRUE),
           # convert all outside of limits to NA
           Rsq   = if_else(is.na(rsqBin), NA, Rsq)
    )
  
  # aggregate
  df.agg = df |>
    group_by(Dyad, Time, pseudoDyad, Feature) |>
    summarise(
      DyadWTC_Rsq   = rsqFun(Rsq,   na.rm = T),
      .groups = "drop"
    )
  
  # if second Limits, then also for phases
  if (length(Limits) == 2) {
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Extract Phase aggregates\n")
    if (labels) phaseLabels = NULL else phaseLabels = FALSE
    phaseFun    = Funs[[2]]
    phaseLimits = Limits[[2]]
    df = df |>
      mutate(phaseBin = cut(Frequency, 
                            breaks = phaseLimits, 
                            labels = phaseLabels, 
                            right = FALSE,
                            include.lowest = TRUE),
             #convert outside to NA
             Phase = if_else(is.na(phaseBin), NA, Phase)
      )
    df.agg = merge(df.agg, 
                   df |>
                     group_by(Dyad, Time, pseudoDyad, Feature) |>
                     summarise(
                       DyadWTC_Phase = phaseFun(Phase, na.rm = T),
                       .groups = "drop"
                     ))
  } else {
    phaseLimits = c()
  }
  
  if (length(phaseLimits) > 2) {
    # aggregate Phase
    df.phase = df |> filter(!is.na(phaseBin)) |>
      group_by(Dyad, Time, pseudoDyad, Feature, phaseBin) |>
      summarise(
        value = phaseFun(Phase, na.rm = T),
        .groups = "drop"
      ) |> tidyr::drop_na() |>
      tidyr::pivot_wider(names_from = phaseBin, names_prefix = "DyadWTC_Phase_")
    df.agg = merge(df.agg, df.phase, all.x = T)
  }
  
  if (length(rsqLimits) > 2) {
    # aggregate Phase
    df.rsq = df |> filter(!is.na(rsqBin)) |>
      group_by(Dyad, Time, pseudoDyad, Feature, rsqBin) |>
      summarise(
        value = rsqFun(Rsq, na.rm = T),
        .groups = "drop"
      ) |> tidyr::drop_na() |>
      tidyr::pivot_wider(names_from = rsqBin, names_prefix = "DyadWTC_Rsq_")
    df.agg = merge(df.agg, df.rsq, all.x = T)
  }
  
  # save the data
  if (!is.null(rs.path)) saveRDS(df.agg, file = flnm)
  
  if (return) return(df.agg)
  
}