#' Load Tracked Data from VERSE Environment
#' 
#' Utility function to read in VERSE data streams from TrackingLog.csv file(s).
#' 
#' @param df.info Dataframe with one row per file to be read in. Must contain the 
#'   column `Filename` (full path to the CSV file) and `Time` which is the name
#'   of the folder automatically by VERSE. All TrackingLogData.csv to which the 
#'   Filenames point have to contain the same number of Social Actors.
#'   Can optionally contain parameters defining windows: 
#'   `start.use` (first Timestamp), `end.use` (last Timestamp), 
#'   or `frame.use` (number of frames starting at `start.use` or sequence origin).
#' @param timezone Character. Timezone in which the data collection was conducted. Default is `"UTC"`.
#' @param rs.path Path to the directory where the output file will be saved, if empty (is.null(rs.path) == TRUE), 
#'   then nothing is saved. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param anonymise Logical. Switch to toggle whether Identifiers should be anonymised and Time should be 
#'   reset to 0 for anonymisation. Default is `FALSE`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns the imported VERSE dataframe. Saves `dataVERSE[suffix].rds` to `rs.path` if provided.
#' 
#' @seealso [readCSVs()]
#' @import dplyr
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export

extractDataVERSE = function(df.info, timezone = "UTC", 
                            rs.path = c(), suffix = '', anonymise = F, 
                            verbose = T, recompute = F, return = T) {
  
  if (verbose) cat("-------------- Extracting data from VERSE experiments --------------\n")
  
  # check whether the data should be saved
  if (is.null(rs.path)) {
    save = F
    rs.path = ''
  } else {
    save = T
  }
  
  flnm = file.path(rs.path, sprintf("dataVERSE%s.arrow", suffix))
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    # give some info
    if (verbose) cat(format(Sys.time(), "%X %Z"), ": Loading VERSE experiments\n")
    df = arrow::read_feather(flnm)
    # check if the file content corresponds to the df.info file
    if (length(c(setdiff(unique(df$Time), df.info$Time), setdiff(df.info$Time, unique(df$Time))) > 0)) {
      stop("Input df.info file differs from loaded file. Please set recompute = T to recompute the file.")
    }
  } else {
    
    # give some info
    if (verbose) cat(format(Sys.time(), "%X %Z"), ": Extracting data from", nrow(df.info), "VERSE experiments\n")

    checkDF(df.info, c("Filename", "Time"))
    
    # extract the header from one of the files
    header = as.character(
      readr::read_delim(df.info$Filename[1], delim = ";", col_names = F,
                        n_max = 1, show_col_types = F))
    # check whether dyadic or individual data
    if (length(header) > 200) nos = 2 else nos = 1
    
    # only use the header for one person
    header = header[1:186]
    
    # read in dependent on whether nos = 2 or nos = 1
    if (nos == 2 ) {
      
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
      df = rbind(df0 |> mutate(Actor = "actor0"), df1 |> mutate(Actor = "actor1")) |>
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
      # set the colnames
      colnames(df) = c("Filename", header)
      # add information 
      df = df |> mutate(
        # Actor column for consistency, even though it is only one
        Actor = "actor0",
        # force the correct timezone
        Timestamp = lubridate::force_tz(Timestamp, tzone = timezone)
        )
    }
    
    # any other information in the df.info is added - needed for cutting data
    df = df |> 
      # disregard any timestamps that may be in df.info
      merge(df.info |> select(-any_of("Timestamp")))
    
    # if the df.info contains start and end times, then cut out everything in-between
    if (sum(c("start.use", "end.use") %in% colnames(df.info)) == 2) {
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Cutting out relevant time window\n")
      df = df |>
        filter(Timestamp >= start.use & Timestamp <= end.use) |>
        select(-start.use, -end.use)
    } else if ("start.use" %in% colnames(df.info)) {
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Cutting out based on starting point\n")
      df = df |>
        filter(Timestamp >= start.use) |>
        select(-start.use)
    } else if ("end.use" %in% colnames(df.info)) {
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Cutting out based on end point\n")
      df = df |>
        filter(Timestamp <= end.use) |>
        select(-end.use)
    }
    
    # add Frames and size of avatar
    df = df |>
      group_by(Time, Identifier) |>
      mutate(
        # add a frame number
        Frame = row_number(),
        # add the avatar size for movement standardisation
        Size  = mean(NeckY) - mean(HipsY)
      )
    
    # if the df.info contains a number of frames, then only keep that many
    if ("frame.use" %in% colnames(df.info)) {
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Cutting out number of frames\n")
      df = df |>
        filter(Frame <= frame.use) |>
        select(-frame.use)
    }
    
    # now let's do some data wrangling
    if (verbose) cat(format(Sys.time(), "%X %Z"), ": Adding information\n")
    df = df |>
      group_by(Identifier) |> 
      mutate(
        # extract the starting time of this experiment from the Filename
        Time = as.POSIXct(gsub(".*/(.+)/TrackingDataLog.*", "\\1", Filename), 
                          format = "%Y-%m-%d_%H-%M-%S", tz = timezone),
        # focus on the target with the shortest distance
        AOI.gaze  = sub("^Name: (.*?) Distance.*$", "\\1", Gaze_Targets, perl = TRUE),
        AOI.left  = sub("^Name: (.*?) Distance.*$", "\\1", LeftEye_Targets, perl = TRUE),
        AOI.right  = sub("^Name: (.*?) Distance.*$", "\\1", RightEye_Targets, perl = TRUE),
        # replace the Identifier with Self
        AOI.gaze  = gsub(sprintf("%s's", Identifier[1]), "Self", AOI.gaze),
        AOI.left  = gsub(sprintf("%s's", Identifier[1]), "Self", AOI.left),
        AOI.right  = gsub(sprintf("%s's", Identifier[1]), "Self", AOI.right),
        # replace all other phrases with possessive with Other
        AOI.gaze  = if_else(!grepl(".*'s .*", AOI.gaze), AOI.gaze, 
                            gsub(".*'s ", "Other ", AOI.gaze)),
        AOI.left  = if_else(!grepl(".*'s .*", AOI.left), AOI.left, 
                            gsub(".*'s ", "Other ", AOI.left)),
        AOI.right  = if_else(!grepl(".*'s .*", AOI.right), AOI.right, 
                             gsub(".*'s ", "Other ", AOI.right))
      )
    
    if (nos == 2) {
      df = df |>
        arrange(Actor, Identifier, Timestamp) |>
        group_by(Time, Timestamp) |>
        mutate(
          # add a Dyad identifier: [Identifier actor0]-[Identifier actor1]
          Dyad     = paste(Identifier, collapse = "-")
        ) |> group_by(Dyad, Identifier, Time) |>
        mutate(
          # get the duration between two frames
          Duration = Timestamp - lag(Timestamp)
        ) |>
        select(-Filename)
    } else {
      df = df |>
        arrange(Identifier, Timestamp) |>
        group_by(Identifier, Time) |>
        mutate(
          # add a Dyad which only contains this Identifier + a suffix
          Dyad = paste0(Identifier, "-solo"),
          # get the duration between two frames
          Duration = Timestamp - lag(Timestamp)
        ) |>
        select(-Filename)
    }
    
    # if the df.info contains avatars, then add this info, otherwise set to NA
    if (sum(c("avatar0", "avatar1") %in% colnames(df.info)) == 2) {
      df = df |>
        mutate(Avatar = if_else(Actor == "actor0", avatar0, avatar1)) |>
        select(-avatar0, -avatar1)
    } else if (sum("avatar0" %in% colnames(df.info)) == 1) {
      df = df |>
        rename(Avatar = avatar0)
    } else {
      df = df |> mutate(Avatar = NA)
    }
    
    # arrange the order of columns
    df = df |>
      relocate(Dyad, Time, Identifier, Avatar, Frame, Timestamp, Duration)
    
    # potentially anonymise the data
    if (anonymise) {
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Anonymising the data\n")
      # reset the time to 0
      df = df |> 
        mutate(Timestamp = as.POSIXct(as.numeric(Timestamp) - as.numeric(Time)),
               Time      = as.POSIXct(0))
      # get all Identifier IDs and rename them with sub-[number]
      nos = length(unique(df$Identifier))
      fmt = paste0("sub-%0", nchar(as.character(nos)), "d")
      recode = sprintf(fmt, 1:nos)
      names(recode) = unique(df$Identifier)
      df$Identifier = recode[df$Identifier]
    }
    
    # save the data frame
    if (save) {
      if (verbose) cat(format(Sys.time(), "%X %Z"), ": Saving the data\n")
      arrow::write_feather(df |> ungroup(), flnm, compression = "zstd")
      }
  }
  
  # return the ungrouped dataframe
  if (return) return(df |> ungroup())
  
  if (verbose) cat(format(Sys.time(), "%X %Z"), ": Done\n")
  
}

#' Read Multiple Tracking Log Files In Parallel
#' 
#' Fast pipeline tool to parse collections of filepaths leveraging `data.table::fread`.
#' 
#' @param Filename Character vector. Structured file paths and locations pointing to individual `TrackingDataLog.csv` files.
#' @param cols Numeric vector. Index maps of specific column dimensions to fetch. Default is `1:186` (corresponds to actor0 inside VERSE configurations).
#'
#' @return A consolidated dataframe representing data parsed from file targets declared in `Filename`.
#' 
#' @seealso [data.table::fread()]
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export

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


#' Extract Information from VERSE Environment Event Files
#' 
#' Utility function to read in raw VERSE environmental logs from an underlying `EventLog.txt` structure.
#' 
#' @param fl.ls Character vector. File system paths pointing to target `EventLog.txt` documents.
#' @param timezone Character. Timezone in which the data collection was conducted. Default is `"UTC"`.
#' @param type Character string. Target format of parsing outcome. Use `"list"` to compile all 
#'   unprocessed elements dynamically, or `"df"` to extract structured tabular fields. Default is `"df"`.
#'
#' @return A composite `list` or integrated `dataframe` depending on selection assigned to `type`.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @import dplyr
#' @export
#' 

extractEventsVERSE = function(fl.ls, timezone = "UTC", type = "df") {
  
  # return all the information as a list
  if (type == "list") {
    # allocate an empty list
    ls.info = vector("list", length(fl.ls))
    # loop through the file paths
    for (i in 1:length(fl.ls)) {
      # read in the Event file
      txt = scan(fl.ls[i], what = "", sep = "\n")
      # extract the json part from it
      ls.info[[i]] = jsonlite::parse_json(txt[3:grep("^}", txt)])
      # grab the VERSE version
      ls.info[[i]]$Version = gsub(".*: (.+)", "\\1", txt[1])
      # include an Events dataframe
      ls.info[[i]]$Events = data.frame(Event = txt[(grep("^}", txt)+1):length(txt)]) |> 
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
      df.txt = data.frame(Event = txt[(grep("^}", txt)+1):length(txt)]) |> 
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
          Filename = gsub("EventLog.txt", "TrackingDataLog.csv", fl.ls[i])
        )
      # add the actors / participants, depending on their number
      if (length(grep("\"DefaultAvatar\":", txt)) == 1) {
        df.txt = df.txt |> 
          mutate(
            actor0  = gsub(".*\"(.+)\"", "\\1", txt[grep("\"Participants\":", txt)+1]),
            avatar0 = gsub(".*\": \"(.+)\",", "\\1", txt[grep("\"DefaultAvatar\":", txt)[1]]), 
            actor1  = "",
            avatar1 = ""
          )
      } else {
        df.txt = df.txt |> 
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
