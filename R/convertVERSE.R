# This R script contains helper functions to read in and extract data from the 
# virtual reality software VERSE. 
# (c) Irene Sophia Plank, 10planki@gmail.com

# list of packages to be installed but not loaded
ls = c("pacman", "lubridate", "data.table")
for (package in ls) {
  if(!(package %in% installed.packages()[,"Package"])) install.packages(package)
}
# then use pacman for all packages that should be loaded - if not installed, installed now
pacman::p_load(tidyverse, jsonlite)

#' Load tracked data from VERSE environment
#' 
#' Utility function to read in VERSE data streams from the TrackingLog.csv
#' 
#' @param df.info Dataframe with one row per file to be read in. 
#' Must contain the column Filename (full path to the CSV file), 
#' can additionally contain columns to describe start and end point of the data to be read in 
#' (start.use : first Timestamp to be used, 
#' end.use : last Timestamp to be used, 
#' frame.use : number of frames, starting either at start.use or at the beginning)
#' @param rs.path Directory path to where the rds file will be saved, if empty nothing is saved
#' @param timezone Timezone in which the data collection was conducted
#' @param suffix Suffix to be added to the file saved to disk (default: "")
#' @param resetTime Switch to toggle whether time shall be reset to 0 for anonymisation (default: FALSE)
#' @param verbose switch to toggle whether output shall be printed to the console (default: TRUE)
#' @param recompute switch to toggle whether existing data should be recomputed and overwritten (default: FALSE)
#' @param return switch to toggle whether a dataframe is returned (default: TRUE)
#' @returns dataVERSE[suffix].rds saved to disk in the rs.path [Optional], with a dataframe of the same data potentially returned [Optional]
#' @seealso [readCSVs()]
#' @export
#' 
extractData = function(df.info, rs.path, timezone, suffix = '',
                       resetTime = F, nos = 2,
                       verbose = T, recompute = F, return = T) {
  
  # check whether the data should be saved
  if (is_empty(rs.path)) {
    save = F
    rs.path = ''
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(file.path(rs.path, sprintf("dataVERSE%s.rds", suffix)))) {
    # give some info
    if (verbose) cat("----------- Loading data from VERSE experiments -----------\n")
    df = readRDS(file.path(rs.path, sprintf("dataVERSE%s.rds", suffix)))
    # check if the file content corresponds to the df.info file
    if (length(c(setdiff(unique(df$Time), df.info$Time), setdiff(df.info$Time, unique(df$Time))) > 0)) {
      stop("Input df.info file differs from loaded file. Please set recompute = T to recompute the file.")
    }
  } else {
    
    # give some info
    if (verbose) cat("----------- Extracting data from", nrow(df.info), "VERSE experiments -----------\n")
    
    # check whether dyadic or individual data
    if (nos == 2 ) {
      
      # extract the header from on of the files
      header = as.character(
        read_delim(df.info$Filename[1], delim = ";", col_select = c(1:186), col_names = F,
                   n_max = 1, show_col_types = F))
      
      # add a column for the conversation Partner and for Listening
      header = c(header, "Partner", "Listening")
      
      # read in the left and the right participant separately
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Reading data from left participants\n")
      df0 = readCSVs(df.info$Filename, c(1:188))
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Reading data from right participants\n")
      df1 = readCSVs(df.info$Filename, c(1, 187:371, 2:3))

      # set the colnames
      colnames(df0) = c("Filename", header)
      colnames(df1) = c("Filename", header)
      
      # now we can combine this data
      df = rbind(df0 %>% mutate(Actor = "actor0"), df1 %>% mutate(Actor = "actor1")) %>%
        mutate(
          # merge speaking and listening columns
          Communication = case_when(
            Speaking & Listening ~ "Both",
            Speaking ~ "Speaking",
            Listening ~ "Listening",
            T ~ "None"
          ), 
          # force the correct timezone
          Timestamp = lubridate::force_tz(Timestamp, tzone = timezone)
        )
      
    } else {
      # if it's just one person, you can simply load the data
      df = readCSVs(df.info$Filename)
    }
    
    # any other information in the df.info is added - needed for cutting data
    df = df %>% 
      merge(., df.info)
    
    # if the df.info contains start and end times, then cut out everything in-between
    if (sum(c("start.use", "end.use") %in% colnames(df.info)) == 2) {
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Cutting out relevant time window\n")
      df = df %>%
        filter(Timestamp >= start.use & Timestamp <= end.use) %>%
        select(-start.use, -end.use)
    } else if ("start.use" %in% colnames(df.info)) {
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Cutting out based on starting point\n")
      df = df %>%
        filter(Timestamp >= start.use) %>%
        select(-start.use)
    } else if ("end.use" %in% colnames(df.info)) {
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Cutting out based on end point\n")
      df = df %>%
        filter(Timestamp <= end.use) %>%
        select(-end.use)
    }
    
    # add Frames and size of avatar
    df = df %>%
      group_by(Time, Identifier) %>%
      mutate(
        # add a frame number
        Frame = row_number(),
        # add the avatar size for movement standardisation
        Size  = mean(NeckY) - mean(HipsY)
      )
    
    # if the df.info contains a number of frames, then only keep that many
    if ("frame.use" %in% colnames(df.info)) {
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Cutting out number of frames\n")
      df = df %>%
        filter(Frame <= frame.use) %>%
        select(-frame.use)
    }
    
    # now let's do some data wrangling
    if (verbose) cat(format(Sys.time(), "%X %Z"), ": Adding information\n")
    df = df %>%
      group_by(Identifier) %>% 
      mutate(
        # extract the starting time of this experiment from the Filename
        Time = as.POSIXct(gsub(".*/(.+)/TrackingDataLog.*", "\\1", Filename), 
                          format = "%Y-%m-%d_%H-%M-%S", tz = timezone),
        # cluster the gaze and eye targets, focusing on the one with the shortest distance
        # and replacing the IDs with Self versus Other
        AOI.gaze  = gsub("Name: (.+) Distance.*", "\\1", Gaze_Targets),
        AOI.gaze  = gsub(sprintf("%s's", Identifier[1]), "Self", AOI.gaze),
        AOI.gaze  = if_else(!grepl(".*'s .*", AOI.gaze), AOI.gaze, 
                            gsub(".*'s ", "Other ", AOI.gaze)),
        AOI.left  = gsub("Name: (.+) Distance.*", "\\1", LeftEye_Targets),
        AOI.left  = gsub(sprintf("%s's", Identifier[1]), "Self", AOI.left),
        AOI.left  = if_else(!grepl(".*'s .*", AOI.left), AOI.left, 
                            gsub(".*'s ", "Other ", AOI.left)),
        AOI.right  = gsub("Name: (.+) Distance.*", "\\1", RightEye_Targets),
        AOI.right  = gsub(sprintf("%s's", Identifier[1]), "Self", AOI.right),
        AOI.right  = if_else(!grepl(".*'s .*", AOI.right), AOI.right, 
                             gsub(".*'s ", "Other ", AOI.right))
      )
    
    if (nos == 2) {
      df = df %>%
        arrange(Actor, Identifier, Timestamp) %>%
        group_by(Time, Timestamp) %>%
        mutate(
          # add a Dyad identifier: [Identifier actor0]-[Identifier actor1]
          Dyad     = paste(Identifier, collapse = "-")
        ) %>% group_by(Dyad, Identifier, Time) %>%
        mutate(
          # get the duration between two frames
          Duration = Timestamp - lag(Timestamp)
        ) %>%
        select(-Filename)
    } else {
      df = df %>%
        arrange(Identifier, Timestamp) %>%
        group_by(Identifier, Time) %>%
        mutate(
          # get the duration between two frames
          Duration = Timestamp - lag(Timestamp)
        ) %>%
        select(-Filename)
    }
    
    # if the df.info contains avatars, then add this info, otherwise set to NA
    if (sum(c("avatar0", "avatar1") %in% colnames(df.info)) == 2) {
      df = df %>%
        mutate(Avatar = if_else(Actor == "actor0", avatar0, avatar1)) %>%
        select(-avatar0, -avatar1)
    } else if (sum("avatar0" %in% colnames(df.info)) == 1) {
      df = df %>%
        rename(Avatar = avatar0)
    } else {
      df = df %>% mutate(Avatar = NA)
    }
    
    # arrange the order of columns
    df = df %>%
      relocate(any_of('Dyad'), Time, Identifier, Avatar, Frame, Timestamp, Duration)
    
    # potentially reset the time to anonymise the data
    if (resetTime) {
      df = df %>% 
        mutate(Timestamp = as.POSIXct(as.numeric(Timestamp) - as.numeric(Time)),
               Time      = as.POSIXct(0))
    }
    
    if (verbose) cat(format(Sys.time(), "%X %Z"), ": Saving the data\n")
    
    # save the data frame
    if (save) saveRDS(df %>% ungroup(), file.path(rs.path, sprintf("dataVERSE%s.rds", suffix)))
  }
  
  # return the ungrouped dataframe
  if (return) return(df %>% ungroup())
  
}

#' Read in multiple files with fread
#' 
#' This function reads in a list of files using the data.table::fread function
#' 
#' @param Filename List of filepaths + filenames pointing to TrackingDataLog.csv files
#' @param cols Which columns to read in, default is 1:186 corresponding to actor0 in VERSE
#' @return Returns the df with the data read in from the files in Filename
#' @seealso [data.table::fread()]
#' @export
#' 
readCSVs = function(Filename, cols = 1:186) {
  
  # add the names to the Filename variable so that these can be used as the ID column
  names(Filename) = Filename
  
  # map over paths and bind them
  df = data.table::rbindlist(
    lapply(Filename, function(f) {
      data.table::fread(
        file = f,
        sep = ";",
        select = cols,
        header = FALSE,
        skip = 1
      )
    }),
    idcol = "Filename"
  )
  
  # return as a data.frame
  return(df)
  
}

#' Extract information from event files created by VERSE environment
#' 
#' Utility function to read in VERSE event files from EventLog.txt
#' @param fl.ls Vector of paths including Filenames for EventLog.txt files
#' @param timezone Timezone in which the data collection was conducted
#' @param type String to choose which type of variable should be created, either 
#' "list" (creates a list containing everthing from the Event file) or 
#' "df" (creates a dataframe with important information and Events)
#' @return list or dataframe
#' @export
#' 

extractEvents = function(fl.ls, timezone, type = "list") {
  
  # return all the information as a list
  if (type == "list") {
    # allocate an empty list
    ls.info = vector("list", length(fl.ls))
    # loop through the file paths
    for (i in 1:length(fl.ls)) {
      # read in the Event file
      txt = scan(fl.ls[i], what = "", sep = "\n")
      # extract the json part from it
      ls.info[[i]] = parse_json(txt[3:grep("^}", txt)])
      # grab the VERSE version
      ls.info[[i]]$Version = gsub(".*: (.+)", "\\1", txt[1])
      # include an Events dataframe
      ls.info[[i]]$Events = data.frame(Event = txt[(grep("^}", txt)+1):length(txt)]) %>% 
        mutate(Timestamp = as.POSIXct(substr(Event, 1, 23), 
                                      format="%Y-%m-%d %X",
                                      tz = timezone),
               Event = gsub(".*?: (.+)", "\\1", Event))
      # name the list entry after the directory
      names(ls.info)[i] = nth(strsplit(fl.ls[i], "/")[[1]], -2)
    }
    # return the list
    return(ls.info)
  } else if (type == "df") {
    # allocate an empty dataframe
    df.info = data.frame()
    # loop through the file paths
    for (i in 1:length(fl.ls)) {
      # read in the Event file
      txt = scan(fl.ls[i], what = "", sep = "\n")
      # create the general dataframe
      df.txt = data.frame(Event = txt[(grep("^}", txt)+1):length(txt)]) %>% 
        mutate(
          # extract the Timestamp from the Event
          Timestamp = as.POSIXct(substr(Event, 1, 23), 
                                 format="%Y-%m-%d %X",
                                 tz = timezone),
          # extract what happened
          Event = gsub(".*?: (.+)", "\\1", Event),
          # add additional information from the Event file
          Version = gsub(".*: (.+)", "\\1", txt[1]),
          Environment = gsub(".*\": \"(.+)\",", "\\1", 
                             txt[grep("\"Environment\":", txt)]),
          Time = as.POSIXct(gsub(".*/(.+)/EventLog.txt", "\\1", fl.ls[i]), 
                            format = "%Y-%m-%d_%H-%M-%S", tz = timezone),
          Filename = fl.ls[i]
        )
      # add the actors / participants, depending on their number
      if (length(grep("\"DefaultAvatar\":", txt)) == 1) {
        df.txt = df.txt %>% 
          mutate(
            actor0  = gsub(".*\"(.+)\"", "\\1", txt[grep("\"Participants\":", txt)+1]),
            avatar0 = gsub(".*\": \"(.+)\",", "\\1", txt[grep("\"DefaultAvatar\":", txt)[1]]), 
            actor1  = NA,
            avatar1 = NA
          )
      } else {
        df.txt = df.txt %>% 
          mutate(
            actor0  = gsub(".*\"(.+)\",", "\\1", txt[grep("\"Participants\":", txt)+1]),
            actor1  = gsub(".*\"(.+)\"", "\\1", txt[grep("\"Participants\":", txt)+2]),
            avatar0 = gsub(".*\": \"(.+)\",", "\\1", txt[grep("\"DefaultAvatar\":", txt)[1]]),
            avatar1 = gsub(".*\": \"(.+)\",", "\\1", txt[grep("\"DefaultAvatar\":", txt)[2]])
          )
      }
      df.info = rbind(df.info, df.txt)
    }
    # return the dataframe
    return(df.info)
  } else {
    stop("Type must be df or list.")
  }
}
