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
#' @param fps Numeric. Frame processing rate frequency profile (frames per second) of the dataset.
#' @param rs.path Character. Path to the directory where the output files will be saved.
#'   If empty (`is.null(rs.path) == TRUE`), nothing is saved to disk. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `TRUE`.
#'
#' @return If `return = TRUE`, returns the dataframe or saves periods of fixations and a consolidated summary to file in `rs.path` if provided.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @import dplyr
#' @export
#' 

featDwell = function(df, ls.AOI, fps, rs.path = c(), suffix = "", 
                     verbose = T, recompute = F, return = T) {
  
  if (verbose) cat("------------------ Extracting dwell time features ------------------\n")
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flcsv = ''
  } else {
    # create filename
    flcsv = file.path(rs.path, sprintf("featDwell%s.csv", suffix))
    fldat = file.path(rs.path, sprintf("dataDwell%s.arrow", suffix))
  }
  
  # if no recompute and the file exists, it is simply loaded
  if (!recompute & file.exists(flcsv)) {
    if (verbose) cat(format(Sys.time(), "%X"), ": Loading dwell times\n")
    df.out = readr::read_csv(flcsv, show_col_types = F)
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
    
    # if ls.AOI is given, classify according to this [!CHECK: None or noAOI here????]
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
    
    # smooth the AOI extraction
    df.dwell = df.dwell |>
      select(Dyad, Time, Identifier, Actor, Frame, Timestamp, AOI, Communication) |>
      #  iterative smoothing to get rid of one-frame AOIs encased by the same other AOI
      arrange(Dyad, Time, Identifier, Frame) |>
      mutate(
        AOI_smooth = iterSmoothing(AOI)
      )
    
    # extract initiations and solo attention to the other's head
    df.ini = df.dwell %>%
      group_by(Dyad) %>%
      group_modify(~ classifyInitiation(.x)) |>
      select(-Actor) |>
      tidyr::pivot_wider(names_from = Classification, values_from = count,
                         names_prefix = "Dwell_")
    
    # extract periods of fixation
    df.blocks = df.dwell |>
      arrange(Dyad, Time, Identifier, Frame) |>
      mutate(block = consecutive_id(AOI_smooth)) |>
      group_by(Dyad, Time, Identifier, block) |>
      filter(n() > 1) |> # get rid of everything that is just one sample
      ungroup() |>
      mutate(block = consecutive_id(AOI_smooth)) |>
      group_by(Dyad, Time, Identifier, block, AOI_smooth) |>
      summarise(
        minFrame = min(Frame),
        maxFrame = max(Frame),
        minTime  = min(Timestamp),
        maxTime  = max(Timestamp),
        .groups = 'drop'
      ) |>
      select(-block) |>
      filter(AOI_smooth != "noAOI") |> 
      mutate(Duration = (maxFrame - minFrame)/fps)
    
    arrow::write_feather(df.blocks, fldat, compression = "zstd")
    
    # aggregate the block durations
    df.blocks = df.blocks |>
      group_by(Dyad, Time, Identifier, AOI_smooth) |>
      summarise(
        AVG = mean(Duration, na.rm = T),
        SD  = sd(Duration, na.rm = T),
        MED = median(Duration, na.rm = T),
        .groups = "drop"
      ) |>
      tidyr::pivot_wider(names_from = AOI_smooth, values_from = c(AVG, SD, MED),
                         names_glue = "DwellBlocks_{AOI_smooth}_{.value}")
    
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
                         names_glue = "{.value}_{AOI}_Total") |>
      left_join(df.blocks, by = c("Dyad", "Time", "Identifier")) |>
      left_join(df.ini, by = c("Dyad", "Time", "Identifier"))
    
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
      readr::write_csv(df.out, flcsv)
    }
    
  }
  
  if (verbose) cat(format(Sys.time(), "%X"), ": Done\n")
  
  # return feature dwell dataframe
  if (return) return(df.out)
  
}

#' Iteratively Smooth Categorical Data
#'
#' Removes single- and two-frame transient noise or "stranglers" from a categorical  
#' vector (such as eye-tracking AOIs) by iteratively replacing isolated values 
#' that are sandwiched between identical preceding and succeeding values.
#'
#' @param x Character. A character or factor vector ordered by frame or time, i.e. AOI classifications
#' @param niter Numeric. An integer specifying the maximum number of smoothing passes to perform 
#'   (default is `10`). The loop terminates early if no further changes are detected.
#'
#' @return A character or factor vector of the same length as `x`, with isolated 
#'   one or two-frame discrepancies smoothed out.
#'
#' @details 
#' The function operates in an iterative loop, performing two sequential checks per pass:
#' 1. **Single-frame gaps (width 3):** Evaluates overlapping three-frame windows. If the left and right 
#'    elements match each other (`left == right`) but differ from the middle element (`mid != left`), 
#'    the middle element is overwritten with the outer value.
#' 2. **Two-frame gaps (width 6):** Evaluates six-frame windows to check if a block of two identical 
#'    values (`m1 == m2`) is encased by matching outer blocks of two frames each (`l1 == l2`, `r1 == r2`, 
#'    and `l1 == r1`), where the middle differs from the outer (`m1 != l1`). If matched, both middle 
#'    frames are overwritten with the outer value.
#'    
#' Single-frame gaps are cleared first, followed immediately by two-frame gaps, ensuring cascading 
#' artifacts are fully resolved across multiple iterations until the vector stabilises.
#'
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#'
#' @examples
#' \dontrun{
#' 
#' # Example sequence with isolated "Desk" glitches inside "Head" looks
#' AOI = c("Head", "Head", "Desk", "Head", "Desk", "Head", "Head")
#' 
#' iterSmoothing(AOI)
#' # Returns: c("Head" "Head" "Head" "Head" "Head" "Head" "Head")
#' 
#' }
#' 
#' @export
#' 
iterSmoothing = function(x, niter = 4) {
  
  for (i in seq_len(niter)) {
    prev_x = x
    n = length(x)
    
    # handle single-frame gaps (width 3)
    if (n >= 3) {
      left  = x[1:(n - 2)]
      mid   = x[2:(n - 1)]
      right = x[3:n]
      
      mid[(left == right) & (mid != left)] = left[(left == right) & (mid != left)]
      x[2:(n - 1)] = mid
    }
    
    # handle two-frame gaps (width 6: 2 left, 2 middle, 2 right)
    n = length(x)
    if (n >= 6) {
      l1 = x[1:(n - 5)]
      l2 = x[2:(n - 4)]
      m1 = x[3:(n - 3)]
      m2 = x[4:(n - 2)]
      r1 = x[5:(n - 1)]
      r2 = x[6:n]
      
      to_change_2 = (l1 == l2) & (r1 == r2) & (l1 == r1) & (m1 == m2) & (m1 != l1)
      
      if (any(to_change_2)) {
        replacement_val = l1[to_change_2]
        
        x[which(to_change_2) + 2] = replacement_val
        x[which(to_change_2) + 3] = replacement_val
      }
    }
    
    # stop early if no changes occurred in this pass
    if (identical(x, prev_x)) break
    
  }
  return(x)
}

#' Classify AOI Dwell Periods
#'
#' Takes a dataframe containing the smoothed AOI data of one Dyad distinguished 
#' by Identifiers and the Actor column. The function focuses on the target AOI
#' and extracts periods of continuous fixation considering both Actors. For each
#' of these periods, it is determined if only one Actor fixated on the target
#' AOI (SoloAttention) or both did in which case it is also determined who 
#' initiated the period (Initiation).
#'
#' @param df.dwell A dataframe containing the smoothed AOI data of one Dyad. Must
#'   contain the columns `AOI_smooth`, `Actor`, `Identifier`, `Time` and `Frame`.
#' @param target Character. Character describing the target AOI. Must be found 
#'   in `AOI_smooth`. Default is `OtherHead`.
#'
#' @return An aggregated dataframe with the columns `Identifier`, `Actor`,
#'   `Classification`, `count` and `Time`.
#'
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' 
#' @export
#' 
classifyInitiation = function(df.dwell, target = "OtherHead") {
  
  # focus on the target AOI and extract start and end frames
  df.dyad = df.dwell |>
    filter(AOI_smooth == target) |>
    group_by(Actor) |>
    mutate(grp = cumsum(c(0, diff(Frame) != 1))) |> group_by(Actor, grp) |>
    summarise(start = min(Frame), end = max(Frame), .groups = "drop") |>
    select(Actor, start, end) |> 
    arrange(start) |>
    mutate(
      BlockID = 0
    )
  
  # add a block ID: as long as overlapping Frames of any actor, same block
  cend = -Inf   # current end of the block
  cid  = 0      # current id number
  
  # loop through and check 
  for (i in seq_len(nrow(df.dyad))) {
    # is this start after the current end point? always T for first block
    if (df.dyad$start[i] > cend) {
      cid = cid + 1         # increase the current block id
      cend = df.dyad$end[i] # set the new current end point to the end of this row
    } else {
      # use the later end point: current maximum or the end of this row
      cend = max(cend, df.dyad$end[i])
    }
    # add the current block ID for this row
    df.dyad$BlockID[i] = cid
  }
  
  # summarise into individual face attention and eye contact initiation
  df.dyad |>
    group_by(BlockID) |>
    # get the lowest starting point for each block for each actor
    summarise(
      Actor0 = if (any(Actor == "actor0")) min(start[Actor == "actor0"]) else NA_integer_,
      Actor1 = if (any(Actor == "actor1")) min(start[Actor == "actor1"]) else NA_integer_,
      .groups = "drop"
    ) |>
    # classify into eye contact initiation and solo face attention
    mutate(
      tmp = case_when(
        Actor0 < Actor1 ~ "actor0_Initiation",
        Actor0 > Actor1 ~ "actor1_Initiation",
        is.na(Actor0)   ~ "actor1_SoloAttention",
        is.na(Actor1)   ~ "actor0_SoloAttention"
      )
    ) |> filter(!is.na(tmp)) |>
    count(tmp, name = "count") |>
    tidyr::separate(col = tmp, into = c("Actor", "Classification")) |>
    mutate(Classification = paste0(target, "_", Classification)) |>
    # combine with the dyad information
    left_join(df.dwell |> select(Identifier, Time, Actor) |>
                distinct(), by = "Actor")
}
