
#' Compute Windowed-Lagged Cross-Correlations (WLCC), potentially with Pseudo-Synchrony Comparison
#'
#' Computes empirical and pseudo windowed-lagged cross-correlations using backend tools from `rMEA`. 
#'
#' @note Requires package `rMEA` (>= 1.3.1) to support generating diagnostic heatmap distributions when `ABS = FALSE`.
#'
#' @param df Dataframe containing tracking data. Requires the variables `Dyad`, `Identifier`, `Frame` and the target numeric timeline (`colname`). 
#' @param rs.path Character. Path to destination directory. Files are written inside a generated directory nested under `[rs.path]/featWLCC[suffix]`.
#'        If empty (is.null(rs.path) == TRUE), then nothing is saved.
#' @param colname Character. Name of column vector inside `df` to isolate for analysis.
#' @param featname Character. Descriptive label for feature. Must not contain underscores.
#' @param win Numeric. Window size in seconds.
#' @param inc Numeric. Window increment step in seconds.
#' @param lag Numeric. Evaluated cross-correlation lag step in seconds.
#' @param fps Numeric. Sampling rate (frames per second).
#' @param suffix Character. Suffix string appended onto `rs.path`. Default is `""`.
#' @param parallel Logical. Enables multi-core clusters inside the parent `MEAccf` routine. Default is `TRUE`.
#' @param pseudoDyad Logical. Flags whether to generate pseudo-WLCC benchmarks using dyad shuffling. Default is `TRUE`.
#' @param nDyad Numeric. Total number of synthetic dyad simulations to execute when `pseudoDyad = TRUE`. Default is `100`.
#' @param bfThreshold Numeric. Credibility threshold for evaluation of log-transformed Bayes Factors. Default is `log(3)`.
#' @param pseudoShuffling Logical. Flags whether to generate pseudo-WLCC benchmarks using within-series random shuffles. Default is `FALSE`.
#' @param shuffleMethod Character. Assigns segment shuffling (`"Seg"`) or cell shuffling (`"Data"`). Default is `"Seg"`.
#' @param nShuffle Numeric. Shuffling iterations for `pseudoShuffling` per dyad. Default is `100`.
#' @param pseudoPass Numeric. Passing threshold to determine credible lags based on different pseudo-WLCC methods. Default is `1`.
#' @param credibleThreshold Numeric. Percentile cutoff for credibility based on permutations. Default is `90`.
#' @param method Character. Perform either peak extraction (`"peak"`) or window averages (`"mean"`) for observed WLCC. Default is `"peak"`.
#' @param r2Z Logical. Applies Fisher's Z transformation over raw cross-correlation indexes. Default is `TRUE`.
#' @param ABS Logical. Maps vectors using absolute values, filtering out directionality of WlCC. Default is `TRUE`.
#' @param seed Universal flag supporting reproducibility. Takes an integer seed or `"random"`. Default is `"random"`.
#' @param plot Logical. Generates PDF evaluation outputs in the output folder. Default is `TRUE`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return Summary dataframe arrays when `return = TRUE`. Creates a structured directory bundle featuring 
#'  `_ccf.rds`, `_df-agg.csv`, `_df.csv` datasets for observed and pseudo-WLCC, 
#'  statistical tables (`_pseudo-comp.csv`), and structural visualization plots (`.pdf`). This bundle
#'  is only created if `!is.null(rs.path)`.
#' 
#' @import dplyr
#' @import rMEA
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 
featWLCC = function(df, rs.path, colname, featname,
                    win, inc, lag, fps, suffix = '', parallel = T,
                    pseudoDyad = T, nDyad = 100, bfThreshold = log(3),
                    pseudoShuffling = F, shuffleMethod = 'Seg', nShuffle = 100,
                    pseudoPass = 1, credibleThreshold = 90, 
                    method = "peak", r2Z = T, ABS = T, seed = 'random', plot = T, 
                    verbose = T, recompute = F, return = T) {
  
  # check if featname contains an underscore - problem with rMEA
  if (grepl("_", featname)) stop("featname cannot include and underscore")
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    wlcc.path = fl.wlcc = flnm = flcsv = ''
  } else {
    # create filename
    fl.wlcc = sprintf("featWLCC-%s-w%d-l%d-i%d", featname, win, lag, inc)
    # create the directory
    wlcc.path = file.path(rs.path, sprintf("featWLCC%s", suffix))
    dir.create(wlcc.path, showWarnings = F, recursive = T)
    # combine
    flcsv = file.path(wlcc.path, sprintf("%s_df-agg.csv", fl.wlcc))
    flnm  = file.path(wlcc.path, fl.wlcc)
  }
  
  # get a random seed
  if (!is.numeric(seed)) {
    seed = sample(1000:9999, 1)
  }
  
  # initialise a list
  ls.mea = c()
  
  # if no recompute and the CSV with aggregated observed WLCC exists, load it
  if (!recompute & file.exists(flcsv)) {
    df.out = readr::read_csv(flcsv, show_col_types = F)
  } else {
    # if no recompute and the RDS file exists, it is simply loaded - observed WLCC
    if (!recompute & file.exists(paste0(flnm, "_ccf.rds"))) {
      ls.ccf = readRDS(paste0(flnm, "_ccf.rds"))
    } else {
      
      # ensure that the dataframe is properly arranged
      df = df |> 
        arrange(Dyad, Identifier, Frame)
      
      # check if the data frame contains the actors, otherwise create the column
      if ("Actor" %in% colnames(df)) {
        if (length(symdiff(unique(df$Actor), c("actor0", "actor1"))) > 0) {
          # rewrite the column but keep the original for later
          df = df |>
            group_by(Dyad, Identifier) |>
            mutate(
              Actor_Original = Actor,
              Actor = if_else(gsub("(.+)-.*", "\\1", Dyad) == Identifier,
                              "actor0", "actor1")
            ) |> ungroup()
        } 
      } else {
        # create the column
        df = df |>
          group_by(Dyad, Identifier) |>
          mutate(
            Actor = if_else(gsub("(.+)-.*", "\\1", Dyad) == Identifier,
                            "actor0", "actor1")
          ) |> ungroup()
      }

      # extract dyads and number of dyads      
      dyads = unique(df$Dyad)
      nod = length(dyads)
      
      if (verbose) cat("----------- Computing", featname, "WLCC for", nod, "dyads -----------\n")
      
      # loop through all dyads
      for (i in 1:nod){
        
        d = dyads[i]
          
        # create fake MEA object
        mea = MEAfake(df |> filter(Dyad == d & Actor == "actor0") |>
                        select(matches(colname)) |> pull(), 
                      df |> filter(Dyad == d & Actor == "actor1") |>
                        select(matches(colname)) |> pull(), 
                      fps, group = featname, id = d,
                      s1Name = "actor0", s2Name = "actor1") 
        
        # add object together into a list
        if (length(ls.mea) == 0) {
          ls.mea = mea
        } else {
          ls.mea = c(ls.mea, mea)
        }
        
      }
      
      # time lagged windowed cross-correlations
      ls.ccf = MEAccf(ls.mea, lagSec = lag, winSec = win, incSec = inc, 
                      r2Z = r2Z, ABS = ABS, cores = parallel) 
      
      if (!is.null(rs.path)) saveRDS(ls.ccf, paste0(flnm, "_ccf.rds"))
      
      # plotting with heatmaps
      if (plot & !is.null(rs.path)) {
        pdf(file.path(wlcc.path, sprintf("%s.pdf", fl.wlcc)))
        for (j in 1:length(ls.ccf)) {
          # configure heatmap
          par(col.main='white')                  # set plot title to white
          if (ABS) {
            heatmap = MEAheatmap(ls.ccf[[j]])
          } else {
            heatmap = MEAheatmap(ls.ccf[[j]],
                                 colors = c("#0048FF", "#86E89E", "#F5FBFF", "#FFF83F","#FF3700"), 
                                 mirror = F)
          }
          par(col.main='black')                  # set plot title back to black
          title(main = names(ls.ccf)[j])  # alternative title
        }
      }
      # dev off for plots
      if (plot & !is.null(rs.path)) dev.off()
      
    }
    
    # check whether pseudoWLCC based on Dyad shuffling should be performed
    if (pseudoDyad) {
      fl.dyad = sprintf("%s_pseudoDyad_seed%04d-n%04d", flnm, seed, nDyad)
      # if no recompute and the RDS file exists, it is simply loaded
      if (!recompute & file.exists(sprintf("%s_df-agg.rds", fl.dyad))) {
        df.pseudoDyad = readRDS(sprintf("%s_df.rds", fl.dyad))
        df.pseudoDyad.agg = readRDS(sprintf("%s_df-agg.rds", fl.dyad))
      } else {
        if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Computing pseudoDyad-WLCC\n")
        # set the random seed
        set.seed(seed)
        # shuffle the dyads while keeping Conf and Part separate
        ls.dyad = shuffle(ls.ccf, size = nDyad, keepRoles = T)
        # calculate the pseudosynchrony
        ls.dyad = MEAccf(ls.dyad, lagSec = lag, winSec = win, incSec = inc, 
                         r2Z = r2Z, ABS = ABS, cores = parallel) 
        # save the outcome
        if (!is.null(rs.path)) saveRDS(ls.dyad, sprintf("%s_ccf.rds", fl.dyad))
        
        # convert to dataframe
        if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Converting pseudoDyad-WLCC to dataframe\n")
        df.pseudoDyad = getCCF(ls.dyad, type = "fullMatrix") |> 
          bind_rows(.id = "name") |>
          tibble::rownames_to_column(var = "window") |> #
          mutate(
            window = gsub("\\..*", "", window)
          ) |> 
          tidyr::pivot_longer(cols = starts_with("lag"), names_to = "lag", values_to = "pseudo") |>
          mutate(
            lag  = as.numeric(substr(lag, 4, nchar(lag))),
            Feature = gsub("^(.+)_.*_1_.*_\\|_.*", "\\1", name)
          )
        
        # aggregate for overall comparison
        if (method == "peak") {
          df.pseudoDyad.agg = df.pseudoDyad |>
            group_by(window, name, Feature) |>
            drop_na() |> 
            # summarise by finding the peak
            summarise(pseudo = max(pseudo)) |>
            group_by(name, Feature) |>
            summarise(pseudo = mean(pseudo))
        } else {
          df.pseudoDyad.agg = df.pseudoDyad |>
            group_by(window, name, Dyad, Feature) |>
            drop_na() |> 
            # summarise using the mean
            summarise(pseudo = mean(pseudo)) |>
            group_by(name, Feature) |>
            summarise(pseudo = mean(pseudo))
        }
        
        # save the dataframe
        if (!is.null(rs.path)) saveRDS(df.pseudoDyad.agg |> ungroup(), sprintf("%s_df-agg.rds", fl.dyad))
        
        # aggregate for the lag comparison 
        df.pseudoDyad = df.pseudoDyad |> 
          group_by(name, Feature, lag) |>
          # drop NAs 
          filter(!is.na(pseudo)) |>
          summarise(
            mean = mean(pseudo, na.rm = T),
            peak = max(pseudo, na.rm = T)
          ) |> 
          tidyr::pivot_longer(names_to = "Stat", cols = c(mean, peak), values_to = "pseudo") |> 
          ungroup() |> 
          filter(!is.na(pseudo)) |> mutate(Method = 'Dyad')
        
        # save the dataframe
        if (!is.null(rs.path)) saveRDS(df.pseudoDyad |> ungroup(), sprintf("%s_df.rds", fl.dyad))
      }
    }
    
    # check whether pseudoWLCC based on Data of Segment shuffling should be performed
    if (pseudoShuffling) {
      # create the filename 
      fl.shuffle = sprintf("%s_pseudo%s_seed%04d-n%04d", flnm, shuffleMethod, seed, nShuffle)
      # if no recompute and the RDS file exists, it is simply loaded
      if (!recompute & 
          file.exists(paste0(fl.shuffle, "_df.rds"))) {
        df.pseudoShuff = readRDS(paste0(fl.shuffle, "_df.rds"))
        df.pseudoShuff.agg = readRDS(paste0(fl.shuffle, "_df-agg.rds"))
      } else {
        if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Computing pseudoShuffling-WLCC\n")
        # use the function to compute pseudoWLCC with shuffling
        df.pseudoShuff = pseudoWLCC(shuffleMethod, ls.ccf, wlcc.path, fl.wlcc,
                                    fps = fps, win = win, inc = inc, lag = lag, 
                                    recompute = recompute, verbose = verbose, 
                                    n = nShuffle, seed = seed)
      
        # aggregate for overall comparison
        if (method == "peak") {
          df.pseudoShuff.agg = df.pseudoShuff |>
            group_by(window, name, sim) |>
            drop_na() |> 
            # summarise by finding the peak
            summarise(pseudo = max(pseudo)) |>
            group_by(name, sim) |>
            summarise(pseudo = mean(pseudo)) |>
            separate(name, into = c("Feature", "Dyad", "Session"), sep = "_") |>
            select(-Session)
        } else {
          df.pseudoShuff.agg = df.pseudoShuff |>
            group_by(window, name, sim) |>
            drop_na() |> 
            # summarise using the mean
            summarise(pseudo = mean(pseudo)) |>
            group_by(name, sim) |>
            summarise(pseudo = mean(pseudo)) |>
            separate(name, into = c("Feature", "Dyad", "Session"), sep = "_") |>
            select(-Session)
        }
        
        # save the dataframe
        if (!is.null(rs.path)) saveRDS(df.pseudoShuff.agg, paste0(fl.shuffle, "_df-agg.rds"))
        
        # aggregate for lag comparison
        df.pseudoShuff = df.pseudoShuff |>
          group_by(name, sim, lag) |>
          # drop NAs 
          filter(!is.na(pseudo)) |>
          summarise(
            mean = mean(pseudo, na.rm = T),
            peak = max(pseudo, na.rm = T)
          ) |> ungroup() |>
          tidyr::pivot_longer(names_to = "Stat", cols = c(mean, peak), values_to = "pseudo") |> 
          filter(!is.na(pseudo)) |> mutate(Method = shuffleMethod) |>
          mutate(
            Dyad    = gsub(".*_(.+)_.*", "\\1", name),
            Feature = gsub("^(.+)_.*_.*", "\\1", name),
            # reconstruct name with simulation number
            name = sprintf("%s_%s_%03d", Feature, Dyad, sim)
          ) |> select(-sim)
        
        # save the dataframe
        if (!is.null(rs.path)) saveRDS(df.pseudoShuff, paste0(fl.shuffle, "_df.rds"))
      }
    }
    
    # convert to dataframe
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Converting WLCC values to dataframe\n")
    # convert from mea list to dataframe
    df.ccf = getCCF(ls.ccf, type = "fullMatrix") |> 
      bind_rows(.id = "name") |>
      tibble::rownames_to_column(var = "window") |> #
      mutate(
        window = gsub("\\..*", "", window)
      ) |> separate(col = name, into = c("Feature", "Dyad", "Session"), sep = "_") |>
      select(-Session) |>
      tidyr::pivot_longer(cols = starts_with("lag"), names_to = "lag", values_to = "WLCC") |>
      mutate(
        lag  = as.numeric(substr(lag, 4, nchar(lag))),
        Type = case_when(
          lag > 0 ~ "actor0",  # positive lag is s1 lead = actor0
          lag < 0 ~ "actor1",  # negative lag is s2 lead = actor1
          lag == 0 ~ "simultaneous")
      )
    # save the resulting dataframe
    if (!is.null(rs.path)) saveRDS(df.ccf |> ungroup(), paste0(fl.shuffle, "_df.rds"))
    
    # choose which pseudoWLCC dataframes to use for lag comparison
    if (all(c("df.pseudoDyad", "df.pseudoShuff") %in% ls())) {
      df.pseudo = rbind(
        df.pseudoShuff |> ungroup() |>
          select(Feature, Method, lag, Stat, pseudo),
        df.pseudoDyad |>
          select(Feature, Method, lag, Stat, pseudo))
    } else if ("df.pseudoDyad" %in% ls()) {
      df.pseudo = df.pseudoDyad |>
        select(Feature, Method, lag, Stat, pseudo)
    } else if ("df.pseudoShuff" %in% ls()) (
      df.pseudo = df.pseudoShuff |> ungroup() |>
        select(Feature, Method, lag, Stat, pseudo)
    )
    # merge it and compare it to the observed WLCC
    if ("df.pseudo" %in% ls()) {
      df.comp = merge(
        # aggregate observed WLCC
        df.ccf |> 
          group_by(Dyad, Feature, lag) |>
          # drop NAs 
          filter(!is.na(WLCC)) |>
          summarise(
            mean = mean(WLCC, na.rm = T),
            peak = max(WLCC, na.rm = T)
          ) |>
          tidyr::pivot_longer(names_to = "Stat", cols = c(mean, peak)) |>
          group_by(Feature, lag, Stat) |>
          summarise(
            observed = mean(value, na.rm = T),
            observed.sd = sd(value, na.rm = T)
          ),
        df.pseudo
      ) |>
        group_by(Method, Feature, lag, Stat) |>
        summarise(
          prob = mean(observed > pseudo)*100,
          observed.sd = mean(observed.sd, na.rm = T),
          observed = mean(observed, na.rm = T),
          pseudo.sd = sd(pseudo, na.rm = T),
          pseudo = mean(pseudo, na.rm = T)
        ) |> group_by(Feature, lag, Stat) |>
        mutate(
          count = sum(prob > credibleThreshold)
        ) |> ungroup() |> 
        mutate(
          credible = if_else(prob > credibleThreshold, "credible lags", "not credible")
        ) |> ungroup()
      if (!is.null(rs.path)) readr::write_csv(df.comp |> ungroup(), paste0(flnm, "_pseudo-comp.csv"))
    }
    
    # summarise overall WLCC regardless of who is leading / following
    if (method == "peak") {
      df.ccf.dyad = df.ccf |> 
        group_by(Dyad, Feature, window) |> 
        drop_na() |> summarise(WLCC = max(WLCC)) |>
        group_by(Dyad, Feature) |>
        summarise(DyadWLCC = mean(WLCC))
    } else {
      df.ccf.dyad = df.ccf |> 
        group_by(Dyad, Feature, window) |> 
        drop_na() |> summarise(WLCC = mean(WLCC)) |>
        group_by(Dyad, Feature) |>
        summarise(DyadWLCC = mean(WLCC))
    }
    
    # check whether overall observed WLCC higher than pseudo WLCC
    shuffle = "not tested"
    dyad    = "not tested"
    if ("df.pseudoShuff.agg" %in% ls()) {
      t.shuff = BayesFactor::ttestBF(df.ccf.dyad$DyadWLCC, df.pseudoShuff.agg$pseudo)
      if ((t.shuff@bayesFactor$bf > bfThreshold) & 
          (mean(df.ccf.dyad$DyadWLCC) > mean(df.pseudoShuff.agg$pseudo))) {
        shuffle = "credible"
      } else {
        shuffle = "not credible"
      }
    }
    if ("df.pseudoDyad.agg" %in% ls()) {
      t.dyad = BayesFactor::ttestBF(df.ccf.dyad$DyadWLCC, df.pseudoDyad.agg$pseudo)
      if ((t.dyad@bayesFactor$bf > bfThreshold) & 
          (mean(df.ccf.dyad$DyadWLCC) > mean(df.pseudoDyad.agg$pseudo))) {
        dyad = "credible"
      } else {
        dyad = "not credible"
      }
    }
    
    # add to the dataframe
    df.ccf.dyad = df.ccf.dyad |>
      mutate(
        testShuffle = shuffle,
        testDyad = dyad
      )
    
    # aggregate all lags over the windows: first, either peak of each window or average
    if (method == "peak") {
      df.ccf.agg = df.ccf |> 
        group_by(Dyad, Type, Feature, window) |> 
        drop_na() |> summarise(WLCC = max(WLCC))
    } else {
      df.ccf.agg = df.ccf |> 
        group_by(Dyad, Type, Feature, window) |> 
        drop_na() |> summarise(WLCC = mean(WLCC))
    }
      
    # then, average these values for each Dyad and Identifier
    df.out = merge(
      # get the leading of the two interaction partners
      df.ccf.agg |> 
        filter(Type != "simultaneous") |>
        group_by(Dyad, Type, Feature) |> 
        summarise(WLCC = mean(WLCC)) |> 
        mutate(
          Identifier = case_when(
            Type == "actor0" ~ gsub("(.+)-.*", "\\1", Dyad), 
            Type == "actor1" ~ gsub(".*-(.+)", "\\1", Dyad))
        ),
      # merge with the dyad WLCC
      df.ccf.dyad, all = T
    ) |> ungroup() |> select(-Type) |>
      relocate(Dyad, Identifier, Feature)
    
    # save feature WLCC dataframe
    if (!is.null(rs.path)) readr::write_csv(df.out, flcsv)
  }
  
  # return aggregated dataframe
  if (return) return(df.out)
  
}

#' Construct Synthetic MEA Object
#'
#' Mock constructor function forming an empty `MEA` structured class object.
#'
#' @param s1 Numeric vector. Metric distribution values mapped for actor track 1.
#' @param s2 Numeric vector. Metric distribution values mapped for actor track 2.
#' @param fps Numeric. Frame sampling frequency rate observed per second.
#' @param s1Name Character. Descriptive label assigned to track `s1`. Default is `"s1Name"`.
#' @param s2Name Character. Descriptive label assigned to track `s2`. Default is `"s2Name"`.
#' @param group Character. Class grouping metadata key attached straight onto the returned object. Default is `"all"`.
#' @param id Character. Unique identity structural key tracking references. Default is `"ID"`.
#' @param session Numeric. Reference indexing counter tracing tracking epochs. Default is `1`.
#'
#' @return A mock class object formatted to match structure schemas native to `rMEA`, labeled as `[group]_[id]_[session]`.
#' 
#' @note s1 and s2 must be two corresponding timecourses, i.e., of the same length and arranged by time.
#' 
#' @import rMEA
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export

MEAfake = function(s1, s2, fps, s1Name = "s1Name", s2Name = "s2Name", 
                   group = "all", id = "ID", session = 1) {
  data = list(s1Name = s1, s2Name = s2)
  names(data) = c(s1Name, s2Name)
  all  = list(all_01_01 = structure(list(MEA = structure(data, row.names = c(NA, -length(s1)), class = "data.frame"), 
                                         ccf = NULL, ccfRes = NULL), id = id, session = session, group = group, sampRate = fps, 
                                    filter = "raw", ccf = "", s1Name = s1Name, s2Name = s2Name, uid = sprintf("%s_%s_%d", group, id, session), 
                                    class = c("MEA","list")))
  names(all) = sprintf("%s_%s_%d", group, id, session)
  mea = structure(all, class = "MEAlist", nId = 1L, n = 1L, groups = group, sampRate = fps, 
                  filter = "raw", s1Name = s1Name, s2Name = s2Name, ccf = "")
  return(mea)
}


#' Permutation-Based Pseudosynchrony for Windowed-Lagged Cross-Correlation
#'
#' Evaluates baseline pseudosynchrony indicators leveraging `rMEA::CCF`. 
#' Based on methods documented in Moulder et al. (2018).
#' 
#' * Assigning `"Data"` tests time independence: shuffles all data points in one of the two sequences.
#' * Assigning `"Seg"` tests regional independence: slices one sequence into windows of length `win` and shuffles these segments.
#'
#' @param shuffleMethod Character string. Specifies the shuffling paradigm; accept options are either `"Data"` or `"Seg"`.
#' @param mea.orig List containing `MEA` objects (can be created via [MEAfake()]).
#' @param rs.path Character. Path to destination directory where the output files will be saved.
#'        If empty (is.null(rs.path) == TRUE), then nothing is saved.
#' @param fl.wlcc Character. Prefix for files saved to disk.
#' @param win Numeric. Window size in seconds.
#' @param inc Numeric. Window increment step in seconds.
#' @param lag Numeric. Evaluated cross-correlation lag step in seconds.
#' @param fps Numeric. Sampling rate (frames per second).
#' @param n Numeric. Total quantity of individual permutations for each `MEA` object. Default is `100`.
#' @param r2Z Logical. Applies Fisher's Z transformation over raw cross-correlation indexes. Default is `TRUE`.
#' @param ABS Logical. Maps vectors using absolute values, filtering out directionality Default is `TRUE`.
#' @param seed Universal flag supporting reproducibility. Takes an integer seed or `"random"`. Default is `"random"`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns aggregated pseudo-WLCC values. Automatically saves `ls.psync` to disk if `rs.path` is provided.
#' 
#' @references Moulder, R. G., et al. (2018). Psychol Methods.
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' @import rMEA
#' @import dplyr
#' 

pseudoWLCC = function(shuffleMethod, mea.orig, rs.path, fl.wlcc, 
                      win = win, inc = inc, lag = lag, fps = fps, 
                      n = 100, r2Z = T, ABS = T, seed = 'random',
                      verbose = T, recompute = F, return = T) {
  
  # get a random seed if none was provided
  if (!is.numeric(seed)) {
    seed = sample(1000:9999, 1)
  }
  set.seed(seed)
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    savePath = flcsv = ''
  } else {
    # add shuffleMethod and seed to savePath
    savePath = file.path(rs.path, sprintf('%s_pseudo%s_seed%04d-n%04d', 
                                          fl.wlcc, shuffleMethod, seed, n))
  }
  
  # check if recompute and if exist
  if (file.exists(paste0(savePath, '_df.rds')) & !recompute) {
    # load the existing list of mea objects
    if (verbose) cat("----------- Loading", basename(savePath), "-----------\n")
    df.pseudo = readRDS(paste0(savePath, '_df.rds'))
    if (return) return(df.pseudo)
  }
  
  # check if the ls exists
  if (!file.exists(paste0(savePath, "_ccf.rds")) | recompute) {
    # start with empty ls.psync
    ls.psync = list()
    
    # check if there are elements in the list
    n.dyads = length(mea.orig)
    if (n.dyads < 1) stop('There are no elements in the list provided!')
    
    # log the progress
    prg = sprintf("----------- Computing pseudoShuff WLCC for %i dyads: %s -----------\n", 
                  n.dyads, gsub("(.+)_.*_1", "\\1", names(mea.orig)[1]))
    
    # if n is not dividable by 2, adjust and add info
    if (n %% 2 > 0) {
      prg = sprintf("%sNumber of simulations was adjusted from %d to %d", prg, n, n+1)
      n = n + 1
    }
    
    counter = 0
    
    # shuffle each side n/2 times
    for (j in 1:(n/2)) {
      
      # go through both sides
      for (s in 1:2) {
        
        # print progress
        counter = counter + 1
        prg = sprintf("%s\n%s : Starting %i of %i simulations", prg, 
                      format(Sys.time(), "%X %Z"),
                      counter, n)
        if (verbose) cat(prg)
        
        # start again with the original MEA data
        mea = mea.orig
        
        # go through all list elements and shuffle them
        for (i in 1:length(mea)) {    
          
          # shuffling one person based on s
          if (shuffleMethod == "Data") {
            mea[[i]][["MEA"]][,s] = sample(mea[[i]][["MEA"]][,s])
          } else if (shuffleMethod == "Seg") {
            if (s == 1) {
              mea[[i]][["MEA"]][,1] = unlist(sample(split(mea[[i]][["MEA"]][[1]], 
                                                          floor(seq_along(mea[[i]][["MEA"]][[1]])/(win*fps)))))
              mea[[i]][["MEA"]][,2] = unlist(split(mea[[i]][["MEA"]][[2]], 
                                                   floor(seq_along(mea[[i]][["MEA"]][[2]])/(win*fps))))
            } else {
              mea[[i]][["MEA"]][,1] = unlist(split(mea[[i]][["MEA"]][[1]], 
                                                   floor(seq_along(mea[[i]][["MEA"]][[1]])/(win*fps))))
              mea[[i]][["MEA"]][,2] = unlist(sample(split(mea[[i]][["MEA"]][[2]], 
                                                          floor(seq_along(mea[[i]][["MEA"]][[2]])/(win*fps)))))
            }
          } else {
            stop("Only Data or Seg can be used as shuffleMethod.")
          }
          
        }
        
        # compute synchrony based on shuffled elements
        mea = MEAccf(mea, lagSec = lag, winSec = win, incSec = inc, 
                     r2Z = r2Z, ABS = ABS, cores = parallel)
        
        # add to the to be returned list
        if (length(ls.psync) == 0) {
          ls.psync = mea
        } else {
          ls.psync = c(ls.psync, mea)
        }
        
      }
      
    }
    
    # save the list for this AU
    if (!is.null(rs.path)) saveRDS(ls.psync, paste0(savePath, "_ccf.rds"))
  } else {
    ls.psync = readRDS(paste0(savePath, "_ccf.rds"))
  }
  
  # convert info to dataframe
  if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Converting pseudo-WLCC to dataframe\n")
  df.pseudo = getCCF(ls.psync, type = "fullMatrix") |> 
    bind_rows(.id = "name") |>
    tibble::rownames_to_column(var = "window") |> #
    mutate(
      window = gsub("\\..*", "", window)
    ) |>
    # add the simulation number
    group_by(name, window) |>
    mutate(sim = row_number()) |> ungroup() |>
    # convert the columns which are the lags into rows
    tidyr::pivot_longer(cols = starts_with("lag"), names_to = "lag", values_to = "pseudo") |> 
    # extract the timings of the lags
    mutate(
      lag  = as.numeric(substr(lag, 4, nchar(lag)))
      )
  
  if (return) return(df.pseudo |> ungroup())
  
}
