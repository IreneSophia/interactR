#' Dwell Time Extraction From Gaze Vectors
#'
#' Extracts descriptive dwell times and indexes joint-attention allocations from spatial gaze data. 
#' The function assumes a dataframe schema with one row per sample frame. 
#' AOIs can be captured separately for eyes in `AOI.left` and `AOI.right` or in one column `AOI`.
#'
#' @param df Dataframe containing tracking data streams. Must explicitly feature columns `Dyad`, 
#'   `Identifier`, `Frame`, `Time`, either `AOI.left` and `AOI.right` or `AOI`.
#' @param ls.AOI List of character vectors. When specified, values isolate targets for AOI classification, 
#'   ignoring undeclared targets. All but alphabet characters will be removed, both in the AOI columns 
#'   and in this list. The order of the AOIs matters: if more than one were to fit, then the first
#'   AOI is chosen. E.g., if ls.AOI = c("Self", "Laptop"), then "Self Laptop" is classified as "Self".
#'   If empty (`is.null(ls.AOI) == TRUE`), existing classification is used.
#' @param rs.path Character. Path to the directory where the output files will be saved.
#'   If empty (`is.null(rs.path) == TRUE`), nothing is saved to disk. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns the dataframe or saves consolidated summary CSV to `rs.path` if provided.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @import dplyr
#' @export
#' 

featDwell = function(df, ls.AOI, rs.path = c(), suffix = "", 
                     verbose = T, recompute = F, return = T) {
  
  if (verbose) cat("------------------ Extracting dwell time features ------------------\n")
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = ''
  } else {
    # create filename
    flnm = file.path(rs.path, sprintf("featDwell%s.csv", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flnm)) {
    if (verbose) cat(format(Sys.time(), "%X"), ": Loading dwell times\n")
    df.out = readr::read_csv(flnm, show_col_types = F)
  } else {
    if (verbose) cat(format(Sys.time(), "%X"), ": Preprocessing dwell times\n")
    
    # check columns
    checkDF(df, c("Dyad", "Identifier", "Frame", "Time"))
    if (!("AOI" %in% colnames(df)) & !all(c("AOI.left", "AOI.right") %in% colnames(df))) {
      stop("Dataframe df must contain either column AOI or columns AOI.left & AOI.right")
    }
    
    # combine the AOI list into a pattern
    if (!is.null(ls.AOI)) pattern = paste(gsub("[^a-zA-Z]", "", ls.AOI), collapse = "|")
    
    # create an Actor column containing actor0 and actor1
    df = df |>
      group_by(Dyad, Identifier, Time) |>
      mutate(
        Actor = if_else(gsub("(.+)-.*", "\\1", Dyad) == Identifier,
                        "actor0", "actor1")
      ) |> ungroup()
    
    # if ls.AOI is given, classify according to this
    if (!is.null(ls.AOI)) {
      if ("AOI" %in% colnames(df)) {
        df = df |> 
          mutate(
            AOI = coalesce(stringr::str_extract(gsub("[^a-zA-Z]", "", AOI), pattern), "noAOI")
          )
      } else {
        df = df |> 
          mutate(
            AOI.left = coalesce(stringr::str_extract(gsub("[^a-zA-Z]", "", AOI.left), pattern), "noAOI"),
            AOI.right = coalesce(stringr::str_extract(gsub("[^a-zA-Z]", "", AOI.right), pattern), "noAOI")
          )
      }
    }
    
    # if necessary, combine the two eyes into one gaze fixation
    if ("AOI" %in% colnames(df)) {
      df.dwell = df
    } else {
      df.dwell = df |> 
        mutate(
          AOI = case_when(
            AOI.left == AOI.right ~ AOI.left, 
            grepl("noAOI", AOI.left) ~ AOI.right,
            grepl("noAOI", AOI.right) ~ AOI.left,
            T ~ AOI.left
          )
        )
    }
    
    # add total number of frames
    df.dwell = df.dwell |>
      group_by(Dyad, Identifier, Time) |>
      mutate(
        Frames.total = n()
      ) |> ungroup() |>
      # remove any where there was no fixation on an AOI
      filter(AOI != "noAOI")
    
    # aggregate the dwell times
    df.dwell.agg = df.dwell |>
      group_by(Dyad, Time, AOI, Identifier, Frames.total) |>
      summarise(
        AOI.frames = n(),
        .groups = "drop"
      ) |> ungroup() |>
      mutate(
        Dwell = AOI.frames * 100 / Frames.total
      ) |> select(-AOI.frames, -Frames.total) |>
      tidyr::pivot_wider(names_from = AOI, values_from = Dwell,
                  names_glue = "{.value}_{AOI}_Total")
    # potentially add the values depending on Communication
    if ("Communication" %in% colnames(df)) {
      df.dwell.agg = merge(
        df.dwell.agg, 
        df.dwell |> 
          group_by(Dyad, Time, AOI, Identifier, Communication, Frames.total) |>
          summarise(
            AOI.frames = n(),
            .groups = "drop"
          ) |> ungroup() |>
          mutate(
            Dwell = AOI.frames * 100 / Frames.total
          ) |> select(-AOI.frames, -Frames.total) |>
          tidyr::pivot_wider(names_from = c(AOI, Communication), values_from = Dwell,
                      names_glue = "{.value}_{AOI}_{Communication}")
      )
    }
    
    # joint attention 
    df.dwell.joint = df.dwell |>
      select(Dyad, Time, Actor, Frame, AOI, Frames.total) |> filter(AOI != "noAOI") |>
      tidyr::pivot_wider(names_from = Actor, values_from = AOI) |>
      filter(actor0 == actor1) |>
      rename(AOI = actor0) |>
      group_by(Dyad, Time, AOI, Frames.total) |>
      summarise(
        value = n()*100,
        .groups = "drop"
      ) |> mutate(value = value/Frames.total) |>
      tidyr::pivot_wider(names_from = AOI,
                         names_glue = "DyadDwell_{AOI}_Total") |>
      ungroup() |> select(-Frames.total)
    
    df.out = merge(df.dwell.agg, df.dwell.joint, all.x = T) |> 
      mutate(across(where(is.numeric), \(x) coalesce(x, 0)))
    
    # save speech dwell dataframe
    if (!is.null(rs.path)) {
      if (verbose) cat(format(Sys.time(), "%X"), ": Saving the Dwell feature csv\n")
      readr::write_csv(df.out, flnm)
      }
    
  }
  
  if (verbose) cat(format(Sys.time(), "%X"), ": Done\n")
  
  # return feature dwell dataframe
  if (return) return(df.out)

}
