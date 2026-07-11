# Functions to extract features from Facial Expressions. 
# (c) Irene Sophia Plank, 10planki@gmail.com

# if packman is not installed yet, install it
if(!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
pacman::p_load(tidyverse)

# This function combines facial expressions to emotional facial expressions
#  (EFE) based on the classification proposed in Aldenhoven et al. (2026).
# While these are based on the ARKit52 system, similar facial expressions can 
# also be extracted using computer vision, e.g., OpenFace. Then, the column
# names have to be adjusted to fit the system described below. 
# (c) Irene Sophia Plank, 10planki@gmail.com
#
# Inputs: 
#   * df : dataframe containing the Facial data. Must contain the columns Dyad, 
#               Identifier, Frame as well as columns with ARKit52 values, 
#               specifically the algorithm uses (.* for Sides): 
#               "Facial_BrowDown.*", "Facial_.*Eye_Squint", "Facial_.*Eye_Wide", 
#               "Facial_MouthPucker", "Facial_MouthFrown.*", "Facial_MouthLowerDown.*", 
#               "Facial_NoseSneer.*", "Facial_CheekSquint.*", "Facial_MouthSmile.*", 
#               "Facial_BrowInnerUp", "Facial_BrowOuterUp.*", "Facial_JawOpen", 
#               "Facial_MouthStretch.*", "Facial_MouthDimple.*"
#               If the data frame contains a Communication column, EFEs are also
#               aggregated based on Speaking, Listening, Both and None. 
#   * rs.path : path to the directory where the output csv will be saved, if 
#               empty (is_empty(rs.path) == TRUE), then nothing is saved
#   * suffix : suffix to be added to the file saved to disk (default: '')
#   * verbose : boolean, whether output is printed to the console (default: TRUE)
#   * recompute : boolean, whether existing data is recomputed and overwritten (default: FALSE)
#   * return : boolean, whether the dataframe is returned (default: FALSE)
# 
# Output:
#   * returns data frame with aggregated results (one row per participant) [Optional]
#   * dataEFE[suffix].rds and featEFE[suffix].csv saved to disk [Optional]
# 
featEFE = function(df, rs.path, suffix = "", verbose = T,
                    recompute = F, return = F) {
  
  # check rs.path
  if (is_empty(rs.path)) {
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
      if (!is_empty(rs.path)) saveRDS(df.face, flrds)
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
    if (!is_empty(rs.path)) write_csv(df.out, flcsv)
  }
  
  # return aggregated dataframe
  if (return) return(df.out)

}