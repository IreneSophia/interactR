#' Dwell Time Extraction From Gaze Vectors
#'
#' Extracts descriptive dwell times and indexes joint-attention allocations from spatial gaze data. 
#' The function assumes a dataframe schema with one row per sample frame. 
#' AOIs can be captured separately for eyes in `AOI.left` and `AOI.right` or in one column `AOI`.
#'
#' @param df Dataframe containing tracking data streams. Must explicitly feature columns `Dyad`, 
#'   `Identifier`, `Frame`, either `AOI.left` and `AOI.right` or `AOI`.
#' @param ls.AOI List of character vectors. When specified, values isolate targets for AOI classification, 
#'   automatically re-coding undeclared targets to `"None"`. 
#'   If empty (`is_empty(ls.AOI) == TRUE`), existing classification is used.
#' @param rs.path Character. Path to the directory where the output files will be saved.
#'   If empty (`is_empty(rs.path) == TRUE`), nothing is saved to disk.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns the dataframe or saves consolidated summary CSV to `rs.path` if provided.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @import tidyverse
#' @export
#' 

featDwell = function(df, ls.AOI, rs.path, suffix = "", 
                     verbose = T, recompute = F, return = T) {
  
  # check rs.path
  if (is_empty(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename
    flnm = file.path(rs.path, sprintf("featDwell%s.csv", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (return) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Loading dwell times\n")
      df.out = read_csv(flnm, show_col_types = F)
    }
  } else {
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Preprocessing dwell times\n")
    # combine the AOI list into a pattern
    if (!is_empty(ls.AOI)) pattern = paste(ls.AOI, collapse = "|")
    
    # create an Actor column containing actor0 and actor1
    df = df %>%
      group_by(Dyad, Identifier) %>%
      mutate(
        Actor = if_else(gsub("(.+)-.*", "\\1", Dyad) == Identifier,
                        "actor0", "actor1")
      ) %>% ungroup()
    
    # if ls.AOI is given, classify according to this
    if (!is_empty(ls.AOI)) {
      if ("AOI" %in% colnames(df)) {
        df = df %>% 
          mutate(
            AOI = coalesce(str_extract(AOI, pattern), "None")
          )
      } else {
        df = df %>% 
          mutate(
            AOI.left = coalesce(str_extract(AOI.left, pattern), "None"),
            AOI.right = coalesce(str_extract(AOI.right, pattern), "None")
          )
      }
    }
    
    # if necessary, combine the two eyes into one gaze fixation
    if ("AOI" %in% colnames(df)) {
      df.dwell = df
    } else {
      df.dwell = df %>% 
        mutate(
          AOI = case_when(
            AOI.left == AOI.right ~ AOI.left, 
            grepl("None", AOI.left) ~ AOI.right,
            grepl("None", AOI.right) ~ AOI.left,
            T ~ AOI.left
          )
        )
      }
    
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
