# Functions to extract features from head movement.
# (c) Irene Sophia Plank, 10planki@gmail.com

# if packman is not installed yet, install it
if(!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
pacman::p_load(tidyverse)

# Helper function to correct for circularity in one rotational data column.
# Caution: x has to be in the correct order, i.e., by Frame / Timecourse. 
# Inspired by Hale et al. (2020)
fixCirc = function(x, th = 270) {
  
  # do nothing, if x only has one entry
  if (length(x) <= 1) return(x)
  
  dz = diff(x)
  
  # initialise a vector of zeros for the shifts
  shifts = rep(0, length(x))
  
  # identify jumps and assign the +360 or -360 correction factor
  up   = which(dz > th)
  down = which(dz < -th)
  
  # *subsequent* indices get shifted
  shifts[up   + 1] = -360
  shifts[down + 1] = +360
  
  # accumulate the shifts forward using cumsum
  return(x + cumsum(shifts))
}

# Helper function to correct for circularity when computing the difference 
# between two timecourses, e.g. head and neck. 
rotDiff = function(x, y) {
  return((((x - y) + 180) %% 360) - 180)
}

# This function preprocesses head motion data. Rotational columns are corrected 
# for circularity or adjusted by a "baseline" body part. Translational columns 
# are detrended.
# 
# Inputs: 
#   * df : dataframe containing all head movement data. Must also contain the 
#               columns Dyad, Identifier and Frame
#   * rs.path : path to the directory where the output csv and rds is saved, if 
#               empty (is_empty(rs.path) == TRUE), then nothing is saved
#   * rotnames : names of columns containing rotational head movement 
#   * tranames : names of the columns containing translational head movement
#   * suffix : suffix to be added to files saved to disk (default: '')
#   * performFixCirc : boolean, whether to fix the circularity (default: TRUE)
#   * cornames : list of column names of the same length of rotnames, these are
#               used to perform a "baseline correction", i.e., the spine for 
#               head rotation independent of the body. Column names have to be 
#               in the same order as the rotnames (default: c(), not performed)
#   * verbose : boolean, whether output is printed to the console (default: TRUE)
#   * recompute : boolean, whether existing data is recomputed and overwritten (default: FALSE)
#   * return : boolean, whether the dataframe is returned (default: FALSE)
# 
# Output:
#   * dataframe with adjusted head movement [Optional]
#   * dataHead[suffix].rds saved to disk in rs.path [Optional]
preproHead = function(df, rs.path, rotnames, tranames, suffix = '',
                      performFixCirc = T, cornames = c(),
                      correct = '', verbose = T, recompute = F, return = F) {
  
  # check rs.path
  if (is_empty(rs.path)) {
    # create empty filename because nothing will be saved
    flnm = flrds = ''
  } else {
    # create filenames
    flnm = file.path(rs.path, sprintf("dataHead%s.rds", suffix))
  }
  
  # give some info
  if (verbose) cat("----------- Preprocess head motion data -----------\n")
  if (file.exists(flnm) & !recompute) {
    df = readRDS(flnm)
  } else {
    # focus on relevant columns
    df = df %>% 
      select(Dyad, Identifier, Frame, 
             any_of(c("Speaking", "Listening", "Communication")),
             any_of(c(rotnames, tranames, cornames))) %>%
      arrange(Dyad, Identifier, Frame)
    
    if (performFixCirc) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Fixing circularity\n")
      df = df %>%
        group_by(Dyad, Identifier) %>%
        # fix circularity based on the algorithm of Hale et al. (2020),
        # default threshold is 270
        mutate(across(all_of(rotnames), fixCirc, .names = "{.col}_fixCirc")) %>%
        ungroup()
    }
    if (length(cornames) == length(rotnames)) {
      if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Rotational difference with to other body part\n")
      # compute rotational difference to a different body part
      for (i in 1:length(cornames)) {
        new_name = paste0(rotnames[i], "_rotDiff")
        df[[new_name]] = rotDiff(df[[rotnames[i]]], df[[cornames[i]]])
      }
    }
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Detrend translational columns\n")
    fixnames = paste0(rotnames, "_fixCirc")
    df = df %>%
      # de-trended translational and fixCirc values by subtracting mean value
      group_by(Dyad, Identifier) %>%
      mutate(across(all_of(c(tranames, fixnames)), ~ .x - mean(.x), .names = "{.col}_detrended")) %>%
      ungroup() %>% arrange(Dyad, Identifier, Frame)
    
    if (verbose) cat(format(Sys.time(), "%x %X %Z"), ": Save data\n")
    if (!is_empty(rs.path)) saveRDS(df, file = flnm)
    
  }
  
  if (return) return(df)

}
