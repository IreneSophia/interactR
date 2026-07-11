
#' Combine ARKit52 Facial Expressions into Emotional Facial Expressions (EFE)
#'
#' This function combines individual facial expressions into emotional facial expressions
#' (EFE) based on the classification proposed in Aldenhoven et al. (2026).
#' While these are based on the ARKit52 system collected in the VERSE environment, 
#' similar facial expressions can also be extracted using computer vision (e.g., OpenFace),
#' provided the column names are adjusted to fit the system described below.
#'
#' @param df Dataframe containing the facial data. Must contain the columns `Dyad`, 
#'   `Identifier` and `Frame`, as well as columns with ARKit52 values. 
#'   Specifically, the algorithm uses the following patterns (`.*` denotes sides):
#'   `"Facial_BrowDown.*"`, `"Facial_.*Eye_Squint"`, `"Facial_.*Eye_Wide"`, 
#'   `"Facial_MouthPucker"`, `"Facial_MouthFrown.*"`, `"Facial_MouthLowerDown.*"`, 
#'   `"Facial_NoseSneer.*"`, `"Facial_CheekSquint.*"`, `"Facial_MouthSmile.*"`, 
#'   `"Facial_BrowInnerUp"`, `"Facial_BrowOuterUp.*"`, `"Facial_JawOpen"`, 
#'   `"Facial_MouthStretch.*"` and `"Facial_MouthDimple.*"`.
#'   If the dataframe contains a `Communication` column, EFEs are also
#'   aggregated based on Speaking, Listening, Both and None.
#' @param rs.path Character. Path to the directory where the output files will be saved.
#'   If empty (`is.null(rs.path) == TRUE`), nothing is saved to disk.
#' @param suffix Character. Suffix to be added to the files saved to disk. Default is `""`.
#' @param verbose Logical. Whether progress and output are printed to the console. Default is `TRUE`.
#' @param recompute Logical. Whether existing data on disk should be recomputed and overwritten. Default is `FALSE`.
#' @param return Logical. Whether the processed dataframe should be returned by the function. Default is `FALSE`.
#'
#' @return If `return = TRUE`, returns a dataframe with aggregated results (one row per participant). 
#'   Otherwise, returns `NULL` invisibly. Saves `dataEFE[suffix].rds` and `featEFE[suffix].csv` to disk if `rs.path` is provided.
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @import tidyverse
#' @references Aldenhoven et al. (2026). Sensors.
#' @export

featEFE = function(df, rs.path, suffix = "", verbose = T,
                    recompute = F, return = F) {
  
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
    df.out = read_csv(flcsv, show_col_types = F)
  } else {
    # no recompute and the RDS file exists, it is simply loaded
    if (!recompute & file.exists(flrds)) {
      df.face = readRDS(flrds)
    } else {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Preprocessing facial expressions\n")
      # preprocessing facial expressions
      df.face = df %>% 
        select(Dyad, Identifier, Frame, any_of(c('Time', 'Timestamp', 'Partner', 'Actor', 'Communication')),
               starts_with("Facial_")) %>%
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
      ls.cols = df.face %>% 
        summarise(across(where(is.numeric), sum)) %>%
        pivot_longer(cols = everything()) %>% 
        filter(value == 0) %>% pull(name)
      # remove them from the dataframe
      df.face = df.face %>%
        select(-all_of(ls.cols))
      # save the preprocessed facial data
      if (!is.null(rs.path)) saveRDS(df.face, flrds)
    }
    # aggregate the emotional expressions
    df.out = df.face %>% select(-starts_with("Facial_")) %>% 
      pivot_longer(cols = c(Anger, Disgust, Joy, Fear, Sadness, Surprise, Contempt), 
                   names_to = "Emotion", names_prefix = "EFE.") %>% 
      group_by(Dyad, Identifier, Time, Partner, Emotion) %>% 
      summarise(value = mean(value)) %>% 
      pivot_wider(names_from = Emotion, 
                  names_glue = "EFE_{Emotion}_Total")
    # potentially adding values depending on Communication
    if ("Communication" %in% colnames(df.face)) {
      df.out = merge(
        df.out, 
        df.face %>% select(-starts_with("Facial_")) %>% 
          pivot_longer(cols = c(Anger, Disgust, Joy, Fear, Sadness, Surprise, Contempt), 
                       names_to = "Emotion", names_prefix = "EFE.") %>% 
          group_by(Dyad, Identifier, Time, Partner, Communication, Emotion) %>% 
          summarise(value = mean(value)) %>% 
          pivot_wider(names_from = c(Emotion, Communication), 
                      names_glue = "EFE_{Emotion}_{Communication}"))
    }
    # save feature face dataframe
    if (!is.null(rs.path)) write_csv(df.out, flcsv)
  }
  
  # return aggregated dataframe
  if (return) return(df.out)

}