#' Dwell time extraction from gaze patterns
#'
#' This function extracts dwell times and joint attention from gaze hits. The function
#' assumes a dataframe with one row per frame. It assumes that AOIs are captured
#' in two columns, one for the left eye (AOI.left) and one for the right eye (AOI.right). 
#' AOI identifiers can be provided which are then used to classify the fixation. 
#' The AOI columns should be "None" if no AOI was hit in this frame. 
#' 
#' @param df Dataframe containing the tracked data. Must contain the columns Dyad, 
#' Identifier and Frame as well as AOI.left and AOI.right. 
#' @param ls.AOI List of character vectors. If not empty, these descriptions will be used to classify the fixation, 
#' ignoring all other possible targets by setting them to "None". 
#' @param rs.path Path to the directory where the output csv will be saved, if empty (is_empty(rs.path) == TRUE), then nothing is saved
#' @param suffix Suffix added to the file saved to disk (default: "")
#' @param verbose Whether output is printed to the console (BOOLEAN, default: TRUE)
#' @param recompute Whether existing data is recomputed and overwritten (BOOLEAN, default: FALSE)
#' @param return Whether the dataframe is returned (BOOLEAN, default: TRUE)
#' @return csv containing aggregated dwell times saved to disk [Optional] or returned [Optional]
#' @import tidyverse
#' @export
#' 

featDwell = function(df, ls.AOI, rs.path, suffix = "", 
                     verbose = T, recompute = F, return = T, save = T) {
  
  # check rs.path
  if (is_empty(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename
    flnm = file.path(rs.path, sprintf("featDwell%s.csv", suffix))
  }
  
  # combine the AOI list into a pattern
  if (!is_empty(ls.AOI)) pattern = paste(ls.AOI, collapse = "|")
  
  # create an Actor column containing actor0 and actor1
  df = df %>%
    group_by(Dyad, Identifier) %>%
    mutate(
      Actor = if_else(gsub("(.+)-.*", "\\1", Dyad) == Identifier,
                      "actor0", "actor1")
    ) %>% ungroup()
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading dwell times\n")
      df.out = read_csv(flnm, show_col_types = F)
    }
  } else {
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Preprocessing dwell times\n")
    # if ls.AOI is given, classify according to this
    if (!is_empty(ls.AOI)) {
      df = df %>% 
        mutate(
          AOI.left = coalesce(str_extract(AOI.left, pattern), "None"),
          AOI.right = coalesce(str_extract(AOI.right, pattern), "None")
          )
    }
    
    # combine the two eyes into one gaze fixation
    df.dwell = df %>% 
      mutate(
        AOI = case_when(
          AOI.left == AOI.right ~ AOI.left, 
          grepl("None", AOI.left) ~ AOI.right,
          grepl("None", AOI.right) ~ AOI.left,
          T ~ AOI.left
          )
      )
    
    # aggregate the dwell times
    df.dwell.agg = df.dwell %>% 
      group_by(Dyad, Identifier) %>%
      mutate(
        # get the total number of frames
        Frames.total = n()
      ) %>%
      group_by(Dyad, AOI, Identifier, Frames.total) %>%
      summarise(
        AOI.frames = n()
      ) %>% ungroup() %>%
      mutate(
        Dwell = AOI.frames * 100 / Frames.total
      ) %>% select(-AOI.frames, -Frames.total) %>%
      pivot_wider(names_from = AOI, values_from = Dwell,
                  names_glue = "{.value}_{AOI}_Total")
    # potentially add the values depending on Communication
    if ("Communication" %in% colnames(df)) {
      df.dwell.agg = merge(
        df.dwell.agg, 
        df.dwell %>% 
          group_by(Dyad, Identifier) %>%
          mutate(
            # get the total number of frames
            Frames.total = n()
          ) %>%
          group_by(Dyad, AOI, Identifier, Communication, Frames.total) %>%
          summarise(
            AOI.frames = n()
          ) %>% ungroup() %>%
          mutate(
            Dwell = AOI.frames * 100 / Frames.total
          ) %>% select(-AOI.frames, -Frames.total) %>%
          pivot_wider(names_from = c(AOI, Communication), values_from = Dwell,
                      names_glue = "{.value}_{AOI}_{Communication}")
      )
    }
    
    # joint attention 
    df.dwell.joint = df.dwell %>%
      select(Dyad, Actor, Frame, AOI) %>% filter(AOI != "None") %>%
      pivot_wider(names_from = Actor, values_from = AOI) %>%
      filter(actor0 == actor1) %>%
      rename(AOI = actor0) %>%
      group_by(Dyad, AOI) %>%
      summarise(
        value = n()*100/f.total
      ) %>% pivot_wider(names_from = AOI,
                        names_glue = "DyadDwell_{AOI}_Total")
    
    df.out = merge(df.dwell.agg, df.dwell.joint, all.x = T) %>% 
      replace(is.na(.), 0)
    
    # save speech dwell dataframe
    if (!is_empty(rs.path)) write_csv(df.out, flnm)
    
  }
  
  # return speech dwell dataframe
  if (return) return(df.out)

}
