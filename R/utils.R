
#' Bin numeric vector by counts
#'
#' This function bins a numeric vector into equally spaced intervals and returns 
#' the count of values in each bin by computing a histogram without rendering the plot.
#'
#' @param x Numeric. Data to be binned.
#' @param Bins Numeric. Number of bins. 
#' @param minLimit Numeric. Lower limit of Bins, needs to be <= min(x). Default is `-pi`.
#' @param maxLimit Numeric. Upper limit of Bins, needs to be >= max(x). Default is `-pi`.
#'
#' @return Returns a numeric vector of length Bins containing the counts. 
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 
binVector = function(x, Bins, minLimit = -pi, maxLimit = pi) {
  if (all(is.na(x))) return(rep(NA_integer_, Bins))
  if ((max(x) > maxLimit) | (min(x) < minLimit)) stop("minLimit and maxLimit need to be at least ", min(x), " and ", max(x), ".")
  graphics::hist(x, breaks = seq(minLimit, maxLimit, length.out = Bins + 1), plot = FALSE)$counts
}


#' Checks column names of dataframe
#'
#' Assesses whether each column named in `colnames` is a column of the dataframe `df`.
#'
#' @param df Dataframe. The dataset to be checked.
#' @param colnames Character vector. The exact name or names of the column(s) in \code{df} 
#' which must be in `df`.
#'
#' @return Throws an error if any of the `colnames` are missing in `df`. 
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 
checkDF = function(df, colnames) {
  colDiff = paste(setdiff(colnames, colnames(df)), collapse = ", ")
  if (nchar(colDiff) > 0) stop("Missing columns in df: ", colDiff)
}

#' Shuffle dyads to create pseudo-dyads
#'
#' Takes a dataframe of dyads and returns a specific number of pseudo-dyads of
#' interaction partners who have not interacted with each other. Identifiers are
#' not paired with themselves, existing dyads are not re-paired, regardless of time,
#' and, if necessary, identifiers from the same time are not paired with each other. 
#' Assignment (left versus right identifier in dyad ID) can be kept or ignored to
#' allow for the inclusion of conditions.
#'
#' @param df Dataframe. Must contain the columns `Dyad`, `Time` and `Identifier`
#'   to unequivocally identify one session. Each Dyad has one row per Identifier. 
#'   Dyad ID must consist of `[Identifier]-[Identifier]`. 
#' @param seed Numeric or other. Seed for reproducibility. If seed is not 
#'   numeric, then a random seed is chosen. Default is `NULL`.
#' @param considerTime Logical. Flags whether to exclude pseudo dyads with same
#'   Time value. Default is `TRUE` corresponding to the exclusion. 
#' @param inOrder Logical. Flags whether to keep assignment of Identifiers. Default
#'   is `FALSE`, allowing left Identifiers to be assigned to right Identifier
#'   in the pseudo-dyads.
#' @param nsim Numeric. Maximum number of pseudo-dyads to be returned. Default is `100`.
#'
#' @return Returns a dataframe with pseudo-dyads, including the information of
#'   `Dyad`, now consisting of both original dyad IDs, as well as each two 
#'   columns for the left and the right `Identifer` and `Time`, thus, one row
#'   per pseudo-dyad.
#' 
#' @import dplyr
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 
shuffleDyads = function(df, seed = NULL, considerTime = T, inOrder = F, nsim = 100) {
  
  checkDF(df, c("Dyad", "Time", "Identifier"))

  # extract the timezone
  timezone = attr(df$Time[1],"tzone")
  
  # set the seed
  if (is.numeric(seed)) set.seed(seed) else set.seed(sample(1000:9999, 1))
  
  # get the original dyads - order here is irrelevant
  ls.dyads = df |> pull(Dyad)
  ls.dyads = unique(
    c(ls.dyads, 
      stringr::str_replace(ls.dyads, "([^-]+)-([^-]+)", "\\2-\\1")))
  
  if (!inOrder) {
    # combine all info into one value
    ls.sess = df |> mutate(sess = paste(Time, Dyad, Identifier, sep = "_")) |> 
      pull(sess)
    # generate combinations
    df.out = as.data.frame(t(combn(ls.sess, 2)))
  } else {
    # combine left versus right info into one value each
    df.sess = df |> 
      mutate(side = if_else(gsub("-.*", "", Dyad) == Identifier, "left", "right"),
             sess = paste(Time, Dyad, Identifier, sep = "_"),
             tmp  = as.numeric(as.factor(paste(Time, Dyad)))) |> 
      select(tmp, Dyad, side, sess) |>
      tidyr::pivot_wider(names_from = side, values_from = sess) |> select(-tmp)
    # generate combinations
    df.out = expand.grid(V1 = df.sess$left, V2 = df.sess$right)
  }
  
  # generate combinations
  df.out = df.out |>
    # extract information again
    tidyr::separate(V1, into = c("left_Time", "left_Dyad", "left_Identifier"), sep = "_") |>
    tidyr::separate(V2, into = c("right_Time", "right_Dyad", "right_Identifier"), sep = "_") |>
    # filter forbidden combinations
    filter(
      left_Dyad       != right_Dyad,        # remove same Identifier
      left_Identifier != right_Identifier   # remove same Dyad
    )
  
  if (nrow(df.out) == 0) stop("No possible new combinations.")
  
  # potentially filter out pseudoDyads with the same Time 
  if (considerTime) {
    df.out = df.out |>
      filter(left_Time != right_Time)   # remove same Time
  }
  
  df.out = df.out |> 
    # create the new Dyad information
    mutate(
      Dyad = paste0(left_Identifier, "-", right_Identifier)
    ) |> 
    # filter out original dyads
    filter(!(Dyad %in% ls.dyads)) |>
    mutate(Dyad = paste0(left_Dyad, "|", right_Dyad)) |>
    select(-left_Dyad, -right_Dyad) |>
    mutate(across(ends_with("Time"), ~ as.POSIXct(.x, tz = timezone)))
  
  # select nsim random dyads
  if (nsim < nrow(df.out)) df.out = df.out[sample.int(nrow(df.out), nsim, replace = F),]
  
  return(df.out)
  
}

#' Shuffle one identifier of a dyad to create pseudo-dyads
#'
#' Takes a dataframe of dyads and returns a specific number of pseudo-dyads of
#' interaction partners who have not interacted with each other by shuffling  either
#' the left or the right Identifier.
#'
#' @param df Dataframe. Must contain the columns `Dyad`, `Time` and `Identifier`
#'   to unequivocally identify one session. Each Dyad has one row per Identifier. 
#'   Dyad ID must consist of `[Identifier]-[Identifier]`. 
#' @param seed Numeric or other. Seed for reproducibility. If seed is not 
#'   numeric, then a random seed is chosen. Default is `NULL`.
#' @param side Character. Which side to shuffle. Default is `right`.
#' @param nsim Numeric. Number of pseudo-dyads to be returned. Default is `100`.
#'
#' @return Returns a dataframe with pseudo-dyads, including the information of
#'   `Dyad`, now consisting of both original dyad IDs, as well as each two 
#'   columns for the left and the right `Identifer` and `Time`, thus, one row
#'   per pseudo-dyad.
#' 
#' @import dplyr
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @export
#' 
shuffleIdentifier = function(df, seed = NULL, side = "right", nsim = 100) {
  
  checkDF(df, c("Dyad", "Time", "Identifier"))
  
  # extract the timezone
  timezone = attr(df$Time[1],"tzone")
  
  # set the seed
  if (is.numeric(seed)) set.seed(seed) else set.seed(sample(1000:9999, 1))
  
  # preprocess the dataframe
  df = df |>
    # extract the side
    mutate(
      Side = if_else(gsub("-.*", "", Dyad) == Identifier, "left", "right"),
      # get a temporary unique ID for each Dyad, in case they interact more often
      tmp = as.numeric(as.factor(paste0(Dyad, "_", Time)))
    ) |>
    # unite dyad, time and identifier into one session information
    tidyr::unite("Session", c(Dyad, Time, Identifier)) |>
    tidyr::pivot_wider(values_from = Session, names_from = Side) |>
    select(-tmp)
  
  df.out = data.frame()
  x = 0
  
  while (nrow(df.out) < nsim) {
    
    x = x + 1
    
    # get the original preprocessed dataframe
    df.tmp = df
    
    # shuffle either the right or the left ensuring no entry is in its original place
    shuffled = sample(df.tmp[[side]])
    while (any(shuffled == df.tmp[[side]])) {
      shuffled = sample(df.tmp[[side]])
    }
    
    # replace original column
    df.tmp[[side]] = shuffled
    
    # extract needed columns
    df.tmp = df.tmp |>
      mutate(pseudoDyad = row_number()) |>
      tidyr::separate(right, into = c("right_Dyad", "right_Time", "right_Identifier"), sep = "_") |>
      tidyr::separate(left,  into = c("left_Dyad",  "left_Time",  "left_Identifier"), sep = "_") |>
      mutate(
        Dyad = paste0(left_Dyad, "|", right_Dyad)
      ) |> select(-left_Dyad, -right_Dyad) |>
      mutate(across(ends_with("Time"), ~ as.POSIXct(.x, tz = timezone)))
    
    # add to df.out and only keep unique rows
    df.out = rbind(df.out, df.tmp |> select(-pseudoDyad)) |> distinct()
    
  }
  
  # potentially select nsim random dyads
  if (nsim < nrow(df.out)) df.out = df.out[sample.int(nrow(df.out), nsim, replace = F),]
  
  return(df.out)
  
}