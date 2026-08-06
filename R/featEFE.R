
#' Combine ARKit52 Facial Expressions into Emotional Facial Expressions (EFE)
#'
#' This function combines individual facial expressions into emotional facial expressions
#' (EFE) based on the classification proposed in Aldenhoven et al. (2026).
#' While these are based on the ARKit52 system collected in the VERSE environment, 
#' similar facial expressions can also be extracted using computer vision (e.g., OpenFace),
#' provided the column names are adjusted to fit the system described below.
#' 
#' #' @details 
#' Emotions are aggregates as the mean of the following Facial Expressions:
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
#' 7. **Contempt:** Difference between sides in `"Facial_MouthDimple.*"`, `"Facial_MouthSmile.*"`
#'
#' @param df Dataframe containing the facial data. Must contain the columns `Dyad`, 
#'   `Identifier`, `Time` and `Frame`, as well as columns with ARKit52 values. 
#'   Specifically, the algorithm uses the following patterns (`.*` denotes sides):
#'   `"Facial_BrowDown.*"`, `"Facial_.*Eye_Squint"`, `"Facial_.*Eye_Wide"`, 
#'   `"Facial_MouthPucker"`, `"Facial_MouthFrown.*"`, `"Facial_MouthLowerDown.*"`, 
#'   `"Facial_NoseSneer.*"`, `"Facial_CheekSquint.*"`, `"Facial_MouthSmile.*"`, 
#'   `"Facial_BrowInnerUp"`, `"Facial_BrowOuterUp.*"`, `"Facial_JawOpen"`, 
#'   `"Facial_MouthStretch.*"` and `"Facial_MouthDimple.*"`.
#'   If the dataframe contains a `Communication` column, EFEs are also
#'   aggregated based on Speaking, Listening, Both and None.
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

featEFE = function(df, rs.path = c(), suffix = "", verbose = T,
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
      
      # check the columns
      checkDF(df, c("Dyad", "Identifier", "Frame", "Time"))
      
      # check the facial columns
      if (length(colnames(df |> 
                          select(matches(
                            c("Facial_BrowDown.*", "Facial_.*Eye_Squint", "Facial_.*Eye_Wide", 
                              "Facial_MouthPucker", "Facial_MouthFrown.*", "Facial_MouthLowerDown.*", 
                              "Facial_NoseSneer.*", "Facial_CheekSquint.*", "Facial_MouthSmile.*", 
                              "Facial_BrowInnerUp", "Facial_BrowOuterUp.*", "Facial_JawOpen", 
                              "Facial_MouthStretch.*", "Facial_MouthDimple.*"))))) < 25) {
        warning("Dataframe df is containing less columns for extracting the emotions than expected.\n",
                "If you aggregated across both sides of the face, this can be expected and is potentially no need to worry.\n",
                "To be sure, please check the documentation of the function to assess whether you have\n",
                "enough information to meaningfully interpret each emotion.")
      }
      
      # rescale some of the facial expressions from VERSE if they are composites
      # of multiple shapes
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
      
      # preprocessing facial expressions
      df.face = df |> 
        select(Dyad, Identifier, Time, Frame, any_of(c('Time', 'Timestamp', 'Partner', 'Actor', 'Communication')),
               starts_with("Facial_")) |>
        # code as emotions following Aldenhoven et al. (2026, Sensors)
        mutate(
          Anger = rowMeans(across(contains(c("BrowDown", "Eye_Squint", 
                                             "Eye_Wide", "MouthPucker")))),
          Disgust = rowMeans(across(contains(c("BrowDown", "MouthFrown", 
                                               "MouthLowerDown", "NoseSneer")))),
          Joy  = rowMeans(across(contains(c("CheekSquint", "MouthSmile")))),
          Fear = rowMeans(across(contains(c("BrowDown", "BrowInnerUp", "BrowOuterUp", 
                                            "Eye_Squint", "Eye_Wide", 
                                            "JawOpen", "MouthStretch")))),
          Sadness = rowMeans(across(contains(c("BrowDown", "BrowInnerUp", "MouthFrown")))),
          Surprise = rowMeans(across(contains(c("BrowInnerUp", "BrowOuterUp", 
                                                "Eye_Wide", "JawOpen")))),
          # for contempt, Aldenhoven state only one side is active
          # thus, we use the difference score between both sides
          Contempt = abs(rowMeans(across(contains(c("MouthDimpleLeft", "MouthSmileLeft")))) - 
                           rowMeans(across(contains(c("MouthDimpleRight", "MouthSmileRight")))))
        )
      # get a list of columns that don't contain any data
      ls.cols = df.face |> 
        summarise(across(where(is.numeric), sum),
                  .groups = "drop") |>
        tidyr::pivot_longer(cols = everything()) |> 
        filter(value == 0) |> pull(name)
      # remove them from the dataframe
      df.face = df.face |>
        select(-all_of(ls.cols))
      # save the preprocessed facial data
      if (!is.null(rs.path)) saveRDS(df.face, flrds)
    }
    # convert to long
    df.face = df.face |> select(-starts_with("Facial_")) |> 
      tidyr::pivot_longer(cols = c(Anger, Disgust, Joy, Fear, Sadness, Surprise, Contempt), 
                          names_to = "Emotion", names_prefix = "EFE.")
    # aggregate the emotional expressions
    df.out = merge(
      # individually for the emotions
      df.face |> 
        group_by(Dyad, Identifier, Time, across(any_of('Partner')), Emotion) |> 
        summarise(value = mean(value),
                  .groups = "drop") |> 
        tidyr::pivot_wider(names_from = Emotion, 
                    names_glue = "EFE_{Emotion}_Total"),
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
                      names_glue = "EFE_{Emotion}_{Communication}")) |>
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