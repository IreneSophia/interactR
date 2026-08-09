
#' Combine ARKit52 Facial Expressions into Emotional Facial Expressions (EFE)
#'
#' This function combines individual facial expressions into emotional facial expressions
#' (EFE) either based on the classification proposed in Aldenhoven et al. (2026) or
#' by using a provided list of EFEs and associated columns. 
#' 
#' #' @details 
#' If `catEFE` is set to Aldenhoven, then emotions are aggregates as the mean of the following Facial Expressions:
#' 1. **Anger:** `"Facial_BrowDown.*"`, `"Facial_.*Eye_Squint"`, `"Facial_.*Eye_Wide"`, 
#'   `"Facial_MouthPucker"`
#' 2. **Disgust:** `"Facial_BrowDown.*"`, `"Facial_MouthFrown.*"`, `"Facial_MouthLowerDown.*"`, 
#'   `"Facial_NoseSneer.*"`
#' 3. **Joy:** `"Facial_CheekSquint.*"`, `"Facial_MouthSmile.*"`
#' 4. **Fear:** `"Facial_BrowDown.*"`, `"Facial_BrowInnerUp"`, `"Facial_BrowOuterUp.*"`, 
#'   `"Facial_.*Eye_Squint"`, `"Facial_.*Eye_Wide"`, `"Facial_JawOpen"`, `"Facial_MouthStretch.*"`
#' 5. **Sadness:** `"Facial_BrowDown.*"`, `"Facial_BrowInnerUp"`, `"Facial_MouthFrown.*"`
#' 6. **Surprise:** `"Facial_BrowInnerUp"`, `"Facial_BrowOuterUp.*"`, `"Facial_.*Eye_Wide"`, 
#'   `"Facial_JawOpen"`
#'
#' @param df Dataframe containing the facial data. Must contain the columns `Dyad`, 
#'   `Identifier`, `Time` and `Frame`, as well as columns detailed in `catEFE`.
#'   If the dataframe contains a `Communication` column, EFEs are also
#'   aggregated based on Speaking, Listening, Both and None.
#' @param catEFE List or character. If it is character, then it must be "Aldenhoven" to use
#'   the predefined categorisation of EFEs. If it is a list, then each list entry
#'   should have the name of the EFE and contain a character vector with the associated column names. 
#' @param rescaleVERSE Logical. Switch whether to correct for composite scores in VERSE Blendshapes. Default is `TRUE`.
#' @param rs.path Character. Path to the directory where the output files will be saved.
#'   If empty (`is.null(rs.path) == TRUE`), nothing is saved to disk. Default is `c()`.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `FALSE`.
#'
#' @return If `return = TRUE`, returns a dataframe with aggregated results (one row per participant). 
#'   Otherwise, returns `NULL` invisibly. Saves `dataEFE[suffix].rds` and `featEFE[suffix].csv` to disk if `rs.path` is provided.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @import dplyr
#' @references Aldenhoven et al. (2026). Sensors.
#' @export

featEFE = function(df, catEFE = "Aldenhoven2026", rescaleVERSE = T, 
                   rs.path = c(), suffix = "", verbose = T,
                   recompute = F, return = T) {
  
  # check rs.path
  if (is.null(rs.path)) {
    # create empty filename because nothing will be saved
    flcsv = flrds = ''
  } else {
    # create filenames
    flcsv = file.path(rs.path, sprintf("featEFE%s.csv", suffix))
    flrds = file.path(rs.path, sprintf("dataEFE%s.rds", suffix))
  }
  
  # if no recompute and the CSV file exists, it is simply loaded
  if (!recompute & file.exists(flcsv)) {
    df.out = readr::read_csv(flcsv, show_col_types = F)
  } else {
    # no recompute and the RDS file exists, it is simply loaded
    if (!recompute & file.exists(flrds)) {
      df.face = readRDS(flrds)
    } else {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Preprocessing facial expressions\n")
      
      # check the general columns
      checkDF(df, c("Dyad", "Identifier", "Frame", "Time"))
      
      # check if specified catEFE or Aldenhoven
      if (!is.list(catEFE)) {
        if (catEFE == "Aldenhoven2026") {
          # use the Aldenhoven et al. (2026, Sensors) categorisation
          catEFE = list(
            "Anger" = c("Facial_BrowDownRight", "Facial_BrowDownLeft", "Facial_RightEye_Squint",
                        "Facial_LeftEye_Squint", "Facial_LeftEye_Wide", "Facial_RightEye_Wide",
                        "Facial_MouthPucker"),
            "Disgust" = c("Facial_BrowDownRight",       "Facial_BrowDownLeft",        "Facial_MouthFrownLeft",      "Facial_MouthFrownRight",    
                          "Facial_MouthLowerDownLeft",  "Facial_MouthLowerDownRight", "Facial_NoseSneerRight",      "Facial_NoseSneerLeft" ),
            "Joy" = c("Facial_CheekSquintLeft",  "Facial_CheekSquintRight", "Facial_MouthSmileRight",  "Facial_MouthSmileLeft"),
            "Fear" = c("Facial_BrowDownRight",   "Facial_BrowDownLeft",    "Facial_BrowInnerUp",     "Facial_BrowOuterUpLeft",
                       "Facial_BrowOuterUpRight","Facial_RightEye_Squint", "Facial_LeftEye_Squint",  "Facial_LeftEye_Wide",   
                       "Facial_RightEye_Wide",   "Facial_JawOpen",         "Facial_MouthStretchLeft","Facial_MouthStretchRight"),
            "Sadness" = c("Facial_BrowDownRight", "Facial_BrowDownLeft",  "Facial_BrowInnerUp",   "Facial_MouthFrownLeft","Facial_MouthFrownRight"),
            "Surprise" = c("Facial_BrowInnerUp",    "Facial_BrowOuterUpLeft","Facial_BrowOuterUpRight", "Facial_LeftEye_Wide",  
                           "Facial_RightEye_Wide",  "Facial_JawOpen")
          )
        }
      }
      
      # check whether all columns named in catEFE are in the dataframe
      if (length(setdiff(unlist(catEFE), names(df))) > 0) stop("The following columns do not exist in `df`: ", 
                                                               paste(setdiff(unlist(catEFE), names(df)), collapse = ", "))
      
      # rescale some of the facial expressions from VERSE if they are composites
      # of multiple shapes
      if (rescaleVERSE) {
        df = df |>
          mutate(
            Facial_BrowInnerUp      = Facial_BrowInnerUp / 2, # combines InnerBrowRaiserR and InnerBrowRaiserL
            Facial_BrowOuterUpLeft  = Facial_BrowOuterUpLeft / 1.5, # takes 0.015 instead of 0.010
            Facial_BrowOuterUpRight = Facial_BrowOuterUpRight / 1.5, # takes 0.015 instead of 0.010
            Facial_CheekPuff        = Facial_CheekPuff / 2, # combines CheekPuffL and CheekPuffR
            Facial_CheekSquintLeft  = Facial_CheekSquintLeft / 2, # combines CheekRaiserL and LipCornerPullerL
            Facial_CheekSquintRight = Facial_CheekSquintRight / 2, 
            Facial_LeftEye_Squint   = Facial_LeftEye_Squint / 1.25, # LidTightenerL and 0.0025 of LipCornerPullerL
            Facial_MouthFunnel      = Facial_MouthFunnel / 1.5, # combines 0.0025 of 6 Blendshapes
            Facial_MouthPucker      = Facial_MouthPucker / 1.3, # combines LipPuckerL and R with each 0.0065
            Facial_RightEye_Squint  = Facial_RightEye_Squint / 1.25 # same as Left
          )
      }
      
      # preprocessing facial expressions
      df.face = df |> 
        (\(.data) {
          colsEFE = purrr::imap(catEFE, ~ rowMeans(.data[.x], na.rm = TRUE))
          names(colsEFE) = paste0("EFE_", names(colsEFE))
          # 3. Bind the prefixed new columns to the dataframe
          bind_cols(.data, as_tibble(colsEFE))
        })()
      # get a list of columns that don't contain any data
      ls.cols = df.face |> 
        summarise(across(where(is.numeric), sum),
                  .groups = "drop") |>
        tidyr::pivot_longer(cols = everything()) |> 
        filter(value == 0) |> pull(name)
      # remove them from the dataframe
      df.face = df.face |>
        select(-all_of(ls.cols)) |>
        mutate(
          EFE_All = rowMeans(select(., starts_with("EFE_")))
        )
      # save the preprocessed facial data
      if (!is.null(rs.path)) saveRDS(df.face, flrds)
    }
    # convert to long
    df.face = df.face |> select(-starts_with("Facial_")) |> 
      tidyr::pivot_longer(cols = starts_with("EFE_"), 
                          names_to = "Emotion")
    # aggregate the emotional expressions
    df.out = merge(
      # individually for the emotions
      df.face |> 
        group_by(Dyad, Identifier, Time, across(any_of('Partner')), Emotion) |> 
        summarise(value = mean(value),
                  .groups = "drop") |> 
        tidyr::pivot_wider(names_from = Emotion, 
                           names_glue = "{Emotion}_Total"),
      # over all emotions
      df.face |> 
        group_by(Dyad, Identifier, Time, across(any_of('Partner'))) |> 
        summarise(EFE_Total = mean(value),
                  .groups = "drop"))
    # potentially adding values depending on Communication
    if ("Communication" %in% colnames(df.face)) {
      df.out = merge(
        df.out, 
        df.face |> 
          group_by(Dyad, Identifier, Time, across(any_of('Partner')), Communication, Emotion) |> 
          summarise(value = mean(value),
                    .groups = "drop") |> 
          tidyr::pivot_wider(names_from = c(Emotion, Communication), 
                             names_glue = "{Emotion}_{Communication}")) |>
        merge(
          df.face |> 
            group_by(Dyad, Identifier, Time, across(any_of('Partner')), Communication) |> 
            summarise(value = mean(value),
                      .groups = "drop") |> 
            tidyr::pivot_wider(names_from = Communication, 
                               names_glue = "EFE_{Communication}")
        )
    }
    # replace any NAs with 0 - for example, if one person is not speaking at all
    df.out = df.out |> 
      mutate(across(everything(), ~ tidyr::replace_na(.x, 0)))
    
    # save feature face dataframe
    if (!is.null(rs.path)) readr::write_csv(df.out, flcsv)
  }
  
  # return aggregated dataframe
  if (return) return(df.out)
  
}