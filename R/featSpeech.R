
#' Aggregate Acoustic Information Extracted via Praat
#'
#' This function aggregates all of the speech information. There two options: 
#' 1. Information extracted from `featSpeech.praat`, which includes the aggregated 
#'    scores from the uhm-o-meter (de Jong et al., 2021). This option is chosen by
#'    providing `praat.path` and `praat.prefix`.
#' 2. Information based on VERSE audio tracking. If `is.null(praat.path) == TRUE`,
#'    the columns `Speaking` and `Listening` from VERSE are used to compute available
#'    features. 
#'
#' @param df.speak Dataframe containing all information about the sounding instances,
#'   typically created using [convertGrid()] for option 1 or a VERSE dataframe for option 2
#'   with columns `Dyad`, `Time`, `Identifier`, `Frame`, `Timestamp` and `Speaking`.
#' @param praat.path Character. Path to the directory containing the Praat output files.
#'   Needs to contain a file of the name `[praat.prefix]_pitchIntensity.csv` or be empty
#'   such that `is.null(praat.path) == TRUE`.
#' @param praat.prefix Character. Prefix used in the Praat script for the output files.
#' @param rs.path Character. Path to the directory where the output files will be saved.
#'   If empty (`is.null(rs.path) == TRUE`), nothing is saved to disk. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `FALSE`.
#'
#' @return If `return = TRUE`, returns the processed dataframe. Saves `featSpeech[suffix].csv` to disk if `rs.path` is provided.
#' 
#' @import dplyr
#' @references de Jong & Wempe (2009). Behavior Research Methods.
#' @references de Jong, Pacilly & Heeren (2021). Assessment in Education: Principles, Policy and Practice.
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export

featSpeech = function(df.speak, praat.path, praat.prefix, rs.path = c(), suffix = '',
                      verbose = T, recompute = F, return = F) {
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename
    flnm = file.path(rs.path, sprintf("featSpeech%s.csv", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading speech features\n")
    df.out = readr::read_csv(flnm, show_col_types = F)
  } else {
    
    if (!is.null(praat.path)) {
      # OPTION 1: ASSUMING PRAAT DATAFRAME
      
      checkDF(df.speak, c("Dyad", "Identifier", "Turn", "Start", "End", "Duration", "nSyll"))
      
      if (!file.exists(file.path(praat.path, paste0(praat.prefix, "_pitchIntensity.csv")))) {
        stop("Specified praat path and prefix do not lead to file ", paste0(praat.prefix, "_pitchIntensity.csv"))
      }
      
      # give some info
      if (verbose) cat("----------- Extracting and aggregating Speech features -----------\n")
      
      # read in the praat output capturing pitch and intensity
      df.pint = readr::read_csv(file.path(praat.path, paste0(praat.prefix, "_pitchIntensity.csv")),
                                show_col_types = F) |>
        tidyr::separate(Name, into = c("tmp1", "Dyad", "Identifier", "tmp2"), sep = "_") |>
        select(-tmp1, -tmp2)
      
      # remove all speaking instances that do not have syllables detected - these
      # are most likely just breathing sounds mistaken for speech
      df.speak = df.speak |> 
        # focus on speaking where there is at least one syllable
        filter(nSyll > 0)
      
    }
    else {
      # OPTION 2: ASSUMING VERSE DATAFRAME
      
      cols = c("Dyad", "Identifier", "Time", "Frame", "Timestamp", "Speaking")
      checkDF(df.speak, cols)
      
      # focus on the relevant columns
      df.speak = df.speak |> select(all_of(cols)) |>
        arrange(Dyad, Identifier, Time, Timestamp) |>
        # number consecutive Frames belonging to one sounding or silence instances
        group_by(Dyad, Identifier, Time) |>
        mutate(run_id = consecutive_id(Speaking),
               Exp.Start = min(Timestamp), 
               Exp.End   = max(Timestamp),
               Exp.Duration = as.numeric(Exp.End - Exp.Start, unit = "secs"),
               Timestamp = Timestamp - Exp.Start) |>
        # only keep the sounding instances
        filter(Speaking) |>
        # reset the numbering to only capture sounding instances
        mutate(Turn = consecutive_id(run_id)) |>
        # aggregate the information for each turn
        group_by(Dyad, Identifier, Time, Turn, Exp.Duration) |>
        summarise(
          Start = as.numeric(min(Timestamp), unit = "secs"),
          End   = as.numeric(max(Timestamp), unit = "secs"),
          Duration = End - Start,
          nSyll = NA,
          .groups = "drop"
        ) |> ungroup()
      
      df.pint = df.speak |>
        select(Dyad, Identifier, Time, Exp.Duration) |>
        distinct() |> rename(Duration = Exp.Duration)
      
    }
    
    # summarise the articulation rate (number of syllables / phonation duration)  
    # and the silence-to-turn ratio (level of the dyad)
    df = df.speak |>
      group_by(Dyad, Identifier) |>
      summarise(
        nSyll = sum(nSyll),
        PhonationDuration = sum(Duration),
        ArticulationRate = nSyll/PhonationDuration,
        .groups = "drop"
      ) |>
      full_join(df.pint, by = c("Dyad", "Identifier")) |>
      group_by(Dyad) |>
      mutate(
        # compute silence-to-turn ratio: higher means more silence
        DyadSPCH_SilenceToTurn = (mean(Duration) - sum(PhonationDuration))/sum(PhonationDuration)
      ) |> 
      select(Dyad, Identifier, any_of(c('PitchSD', 'IntensitySD')), ArticulationRate, PhonationDuration, DyadSPCH_SilenceToTurn) |>
      rename_with(~ paste0("SPCH_", .x), .cols = any_of(c('PitchSD', 'IntensitySD', 'ArticulationRate', 'PhonationDuration')))
    
    # extract the PhonationBalance for each participant
    df = df |>
      full_join(
        df |> select(Dyad, Identifier, SPCH_PhonationDuration) |>
          mutate(
            actor = if_else(gsub("(.+)-.*", "\\1", Dyad) == Identifier, "left", "right")
          ) |> select(-Identifier) |>
          tidyr::pivot_wider(names_from = actor, values_from = SPCH_PhonationDuration) |>
          mutate(
            tmp   = left/right,
            right = right/left
          ) |> select(-left) |> rename(left = tmp) |>
          tidyr::pivot_longer(cols = c(right, left), names_to = 'Identifier', values_to = 'SPCH_PhonationBalance') |>
          mutate(
            Identifier = if_else(Identifier == "right", 
                                 gsub(".*-(.+)", "\\1", Dyad), 
                                 gsub("(.+)-.*", "\\1", Dyad))
          ),
        by = c("Dyad", "Identifier")
      )
    
    # use the function to detect turns
    df.turns = detectTurns(df.speak, rs.path = rs.path, suffix = suffix,
                           verbose = verbose, recompute = T, return = T)
    
    # aggregate and merge all the information
    if (verbose) cat(format(Sys.time(), "%X %Z"), ": Aggregate and save features\n")
    df.out = df.turns |> 
      group_by(Dyad, Identifier) |> 
      summarise(SPCH_TurnGapsMedian = median(TTG, na.rm = T),
                SPCH_TurnGapsSD     = sd(TTG, na.rm = T),
                .groups = "drop") |>
      full_join(df.turns |> group_by(Dyad) |> summarise(DyadSPCH_nTurns = max(Turn),
                                                        .groups = "drop"),
                by = 'Dyad') |>
      full_join(df, by = c('Dyad', 'Identifier')) |>
      select(where(~ any(!is.na(.x))))
    
    # save the features
    if (!is.null(rs.path)) readr::write_csv(df.out, file = flnm)
    
  }
  
  if (return) return(df.out |> ungroup())
  
}

#' Detect Turns based on Sounding Instances
#'
#' This function takes a datframe containing sounding instances and detects turns -
#' continuous sounding instances by one speaker.  Specifically, all sounding 
#' instances that are completely engulfed by the counterpart's sounding  
#' instance are disregarded. Then, a turn goes from the start of the first 
#' until the end of the last consecutive sounding instance of one speaker.
#'
#' @param df.speak Dataframe containing all information about the sounding instances,
#'   can be created using [convertGrid()] and needs to contain the columns
#'   `Dyad`, `Identifier`, `Turn`, `Start`, `End` and `Duration`.
#' @param rs.path Character. Path to the directory where the output files will be saved.
#'   If empty (`is.null(rs.path) == TRUE`), nothing is saved to disk. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `FALSE`.
#'
#' @return If `return = TRUE`, returns the processed dataframe. Saves `dataTurn[suffix].rds` to disk if `rs.path` is provided.
#' 
#' @import dplyr
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 

detectTurns = function(df.speak, rs.path = c(), suffix = '',
                       verbose = T, recompute = F, return = F) {
  
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename
    flnm = file.path(rs.path, sprintf("dataTurn%s.rds", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading turns\n")
      df.turns = readRDS(flnm)
    }
  } else {
  
    checkDF(df.speak, c("Dyad", "Identifier", "Turn", "Start", "End", "Duration"))
    
    # ensure that the dataframe is properly arranged
    df.speak = df.speak |>
      arrange(Dyad, Start, End) |>
      mutate(
        engulfed = F
      )
    
    # we need to get rid of all sounding instance that are completely engulfed in another
    if (verbose) cat(format(Sys.time(), "%X %Z"), ": Remove engulfed sounds\n")
    for (i in 2:nrow(df.speak)) {
      if (sum((df.speak$Start[i] >= df.speak[(df.speak$Dyad == df.speak$Dyad[i]),]$Start) &  
              (df.speak$End[i]   <= df.speak[(df.speak$Dyad == df.speak$Dyad[i]),]$End)) > 1 ) { 
        df.speak$engulfed[i] = T
      } 
    }
    df.speak = df.speak |> filter(engulfed == F) |> select(-c(engulfed))
    
    # identify turns: here, turns are defined as starting with the first sounding
    # instance of a person until the end of the last sounding instance of this 
    # person before a non-engulfed sounding instance of another person
    if (verbose) cat(format(Sys.time(), "%X %Z"), ": Detect turns\n")
    df.turns = df.speak |>
      ungroup() |>
      mutate(rown = row_number()) |>              # add row number
      group_by(Dyad, Identifier) |>               # group by the person speaking
      mutate(
        tn = cumsum(c(TRUE, diff(rown) > 1))      # always keep the lowest row number of this turn as turn number
      ) |>
      ungroup() |>
      mutate(
        Turn = paste0(Identifier, "_", tn)        # add this turn number to the person speaking
      ) |>
      group_by(Dyad, Identifier, Turn) |>         # summarise by Dyad, Identifier and Turn
      summarise(
        StartTurn = min(Start, na.rm = T),        # take the start of the first sounding instance
        EndTurn   = max(End, na.rm = T),          # take the end of the last sounding instance
        DurTurn   = EndTurn - StartTurn,          # compute duration of the turn
        .groups = "drop") |> 
      arrange(Dyad, StartTurn) |>
      group_by(Dyad) |>
      mutate(
        Turn = row_number()
      ) |>
      group_by(Dyad) |>
      mutate(
        TTG = StartTurn - lag(EndTurn)
      )
    
    # save the features
    if (!is.null(rs.path)) saveRDS(df.turns, file = flnm)
    
  }
  
  if (return) return(df.turns |> ungroup())
  
}


#' Convert TextGrid Outputs into Dataframes
#'
#' This function converts the TextGrid output produced by the uhm-o-meter into a 
#' dataframe containing all sounding instances and the total number of syllables.
#'
#' @param ls.files Character vector. Paths for the TextGrid files, including filename and extension.
#' @param rs.path Character. Path to the directory where the output files will be saved.
#'   If empty (`is.null(rs.path) == TRUE`), nothing is saved to disk. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param prefix Character. Prefix for the filename allowing extraction of Dyad and Identifier from filenames. Default is `""`.
#' @param extract Logical. Whether Dyad and Identifier variables are extracted from the file names.
#'   Requires the filename structure: `"[prefix][Dyad]_[Identifier]_*"`. Default is `TRUE`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns the processed dataframe. Saves `dataUhm[suffix].rds` to disk if `rs.path` is provided.
#' 

#' @references de Jong & Wempe (2009). Behavior Research Methods.
#' @references de Jong, Pacilly & Heeren (2021). Assessment in Education: Principles, Policy and Practice.
#' 
#' @import dplyr
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 
convertGrid = function(ls.files, rs.path = c(), suffix = '', prefix = '', extract = T, 
                       verbose = T, recompute = F, return = F) {
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename
    flnm = file.path(rs.path, sprintf("dataUhm%s.rds", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading sounding instances\n")
      df.speak = readRDS(flnm)
    }
  } else {
  
    # give some info
    if (verbose) cat("----------- Extracting speak df from", length(ls.files), "TextGrids -----------\n")
    
    # initialise a dataframe
    df.speak = data.frame()
    
    # loop through the paths
    for (path in ls.files) {
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Converting", basename(path), "\n")
      
      # read in the TextGrid file
      txt = scan(path, what = "", sep = "\n", quiet = T)
      
      # create a dataframe with the turns
      idx = which(grepl("text = \"[0-9]+\"", txt))
      df.tmp = data.frame(Turn = 1:length(idx),
                          Start = as.numeric(gsub(".*xmin = (.+)", "\\1", txt[idx-2])), 
                          End   = as.numeric(gsub(".*xmax = (.+)", "\\1", txt[idx-1]))) |>
        mutate(
          Duration = End - Start, 
          Path = path
        )
      
      # create a dataframe with syllables
      idx = which(grepl("mark = \"[0-9]+\"", txt))
      df.syll = data.frame(x = as.numeric(gsub(".*number = (.+)", "\\1", txt[idx-1])))
      
      # add the syllables to the turns
      df.tmp$nSyll = NA
      for (i in 1:nrow(df.tmp)) {
        df.tmp$nSyll[i] = nrow(df.syll |> filter(x <= df.tmp$End[i] & x >= df.tmp$Start[i]))
      }
      
      # add to the dataframe
      df.speak = rbind(df.speak, df.tmp)
    }
    
    # extract the Dyad and the Identifier
    if (extract) {
      df.speak = df.speak |>
        mutate(
          Dyad = gsub(sprintf("^%s(.+)_.*_.*", prefix), "\\1", basename(Path)),
          Identifier = gsub(sprintf("^%s.*_(.+)_.*", prefix), "\\1", basename(Path)),
        ) |> select(-Path) |> relocate(Dyad, Identifier)
    }
  }
  
  # potentially save to disk
  if (!is.null(rs.path)) saveRDS(df.speak, flnm)
  
  if (return) return(df.speak |> ungroup())
  
}

#' Add Conversational States Based on Praat Output
#'
#' This function adds the `Listening`, `Speaking`, and `Communication` state columns
#' based on speech analysis performed in Praat using the uhm-o-meter developed by 
#' de Jong et al. (2021). If these columns already exist in the dataset, they will 
#' be dynamically renamed with the suffix `"_Original"`.
#'
#' @param df Dataframe containing the tracked data. Must contain the columns `Dyad`,
#'   `Time`, `Frame` and `Timestamp` (in POSIX format).
#' @param df.speak Dataframe containing all information about the sounding instances,
#'   typically created using [convertGrid()]. 
#'   Must contain the columns `Dyad`, `Identifier`, `Start`, and `End` (both in seconds).
#' @param rs.path Character. Path to the directory where the output files will be saved.
#'   If empty (`is.null(rs.path) == TRUE`), nothing is saved to disk. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns the processed dataframe. Saves `data[suffix].rds` to disk if `rs.path` is provided.
#' 
#' @references de Jong & Wempe (2009). Behavior Research Methods.
#' @references de Jong, Pacilly & Heeren (2021). Assessment in Education: Principles, Policy and Practice.
#' 
#' @import dplyr
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 
addCommunication = function(df, df.speak, rs.path = c(), suffix = '',
                            verbose = T, recompute = F, return = T) { 
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename
    flnm = file.path(rs.path, sprintf("data%s.rds", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading uhm-adjusted data\n")
      df = readRDS(flnm)
    }
  } else {
    # give some info
    if (verbose) cat("----------- Adding speaking info from uhm-o-meter to df -----------\n")
    
    # check for columns
    checkDF(df, c("Dyad", "Time", "Frame", "Timestamp"))
    checkDF(df.speak, c("Dyad", "Identifier", "Start", "End"))
    
    # rename the Listening, Speaking, Communication columns if they exist
    if (any(c("Listening", "Speaking", "Communication") %in% colnames(df))) {
      df = df |>
        rename_with(~ paste0(.x, "_Original"), 
                    any_of(c("Listening", "Speaking", "Communication")))
    }
    
    # add a Timepoint
    df = df |> group_by(Dyad, Identifier, Time) |> 
      mutate(Timepoint = Timestamp - min(Timestamp)) |> 
      ungroup()
    
    # Ensure dataframes are data.tables
    data.table::setDT(df)
    data.table::setDT(df.speak)
    
    # enable multi-threading for data.table (speeds up operations automatically)
    data.table::setDTthreads(0) # uses all available cores
    
    # Add a temporary row index to track rows after joining
    data.table::set(df, j = "row_id", value = seq_len(nrow(df)))
    
    # check if Identifier is speaking
    idxSpeaking = df.speak[
      df, 
      on = .(Dyad = Dyad, Identifier = Identifier, Start <= Timepoint, End >= Timepoint), 
      nomatch = NULL, 
      unique(i.row_id)
    ]
    
    # assign the speaking to the dataframe
    data.table::set(df, j = "Speaking", value = FALSE)
    data.table::set(df, i = which(df$row_id %in% idxSpeaking), j = "Speaking", value = TRUE)

    # additional preprocessing    
    df = df |>
      # add the Listening - only two rows per Dyad and Timepoint
      group_by(Dyad, Time, Timepoint) |>
      mutate(Listening = rev(Speaking)) |>
      ungroup() |>
      # merge speaking and listening columns
      mutate(
        Communication = case_when(
          Speaking & Listening ~ "Both",
          Speaking ~ "Speaking",
          Listening ~ "Listening",
          T ~ "None"
        )
      ) |> select(-row_id)
  }
  
  # save to disk
  if (!is.null(rs.path)) saveRDS(df, flnm)
  
  if (return) return(df |> ungroup())
  
}

#' Filter and Rewrite TextGrids for Manual Auditory Inspection
#'
#' The output of the uhm-o-meter requires visual and manual auditory inspection. 
#' Most misclassifications affect sounding instances with few syllables (`nSyll`), 
#' capturing loud breathing rather than authentic speaking. This utility function exports 
#' targeted TextGrids containing only specific syllable-count ranges to ease down-stream validation.
#'
#' @param ls.files Character vector. Full system paths and filenames pointing to target TextGrid files.
#' @param minSyll Numeric. Minimum syllables a valid sounding instance must contain. Default is `1`.
#' @param maxSyll Numeric. Maximum syllables a valid sounding instance can contain. Default is `4`.
#' @param max.nos Numeric. Maximum total number of sounding instances to keep. Default is `Inf`.
#' @param min.dist Numeric. Minimum distance in seconds required between separate sounding instances. Default is `0`.
#'
#' @return Saves reconstructed TextGrid files to disk in the same folder as original file named as `[filename]_check-[nos].TextGrid`.
#' 
#' @references de Jong & Wempe (2009). Behavior Research Methods.
#' @references de Jong, Pacilly & Heeren (2021). Assessment in Education: Principles, Policy and Practice.
#' 
#' @import dplyr
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export

rewriteGrid = function(ls.files, rs.path, minSyll = 1, maxSyll = 4,
                       max.nos = Inf, min.dist = 0) {
  
  for (path in ls.files) {

    # extract the turns and nSyll, then filter to those that should be removed
    df.speak = convertGrid(path, c(), verbose = F, 
                           return = T, extract = F) |>
      select(-Path) |>
      # exclude all turns that do not fit the nSyll range
      mutate(exclude = nSyll < minSyll | nSyll > maxSyll) |>
      # check whether distance to the neighbouring instances exceeds max.dist
      arrange(Start) |>
      mutate(DistanceStart = Start - lag(End, default = -Inf), # never exclude first based on distance to start
             DistanceEnd   = lead(Start, default = Inf) - End,     # never exclude last based on distance to end
             exclude = if_else(exclude == F & (DistanceStart < min.dist | DistanceEnd < min.dist),
                               T, exclude)) |>
      # additionally exclude all that exceed the maximum number of speaking instances
      # arrange by nSyll to keep the ones with the least number of Syllables
      arrange(nSyll) |> group_by(exclude) |>
      mutate(row = row_number(),
             exclude = if_else(!exclude & row > max.nos, T, exclude)) |>
      filter(exclude) |> arrange(Turn)
    
    # read in the TextGrid file
    txt = scan(path, what = "", sep = "\n", quiet = T)
    
    # extract the indices of all sounding instances
    idx.txt = which(grepl("text = \"[0-9]+\"", txt))
    
    # only keep the ones that are supposed to be deleted
    idx.rel = idx.txt[df.speak$Turn] 
    
    # remove these indices from the txt
    idx.keep = setdiff(1:length(txt),
                         # all of the indices for this sounding
                       c(idx.rel, idx.rel-1, idx.rel-2, idx.rel-3,
                         # and the following silence
                         idx.rel+1, idx.rel+2, idx.rel+3, idx.rel+4))
    txt.new = txt[idx.keep]
    
    # adjust interval size
    noi = length(which(grepl("text = ", txt.new)))
    idx.intro = which(grepl("intervals: size = ", txt.new))
    txt.new[idx.intro] = gsub("size = .*", sprintf("size = %d ", noi), txt.new[idx.intro])
    
    # save the new TextGrid with number of sounding instances in the name
    nos = length(which(grepl("text = \"[0-9]+\"", txt.new)))
    fl = file(gsub(".TextGrid", sprintf("_check-%d.TextGrid", nos), path))
    writeLines(txt.new, fl)
    close(fl)
    
  }

}
