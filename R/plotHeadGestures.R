#' Plot Head Gestures Extracted from Time Series Data
#'
#' Input should be data of one Identifier preprocessed using \code{\link{preproHead}}. This 
#' data is then used to compute and plot the results from \code{\link{featHeadGestures}}. 
#' One subplot for the nodding and one for the head shaking. These gestures are extracted
#' by detecting Zero Crossing using \code{\link{extractZCrossing}}. Gestures are classified
#' exclusively as nodding, head shaking or nothing for each Frame.
#'
#' @param df Dataframe. The dataset containing the variables to be processed, potentially created by \code{\link{preproHead}}. 
#'   Must explicitly feature columns `Identifier`, `Frame` and the columns `colNodding` and `colShaking`. This dataframe
#'   will be processed using \code{\link{featHeadGestures}} to extract nodding and head shaking, assuming that only one can happen at a time. 
#' @param colNodding Character. The exact name of the column in \code{df} from which
#'   to extract, and then plot, nodding. 
#' @param colShaking Character. The exact name of the column in \code{df} from which
#'   to extract, and then plot, head shaking. 
#' @param fps Numeric. Frame processing rate frequency profile (frames per second) of the dataset.
#' @param minDegree Numeric. How many degree of rotational difference are needed for the movement to be considered relevant. 
#' Depends on the fps and the specific movement. Setting to negative number leads to no thresholding based on degrees. 
#' @param minFrame Numeric or NULL. First Frame to be plotted. If `NULL` minimum available Frame is used. Default is `NULL`.
#' @param maxFrame Numeric or NULL. Last Frame to be plotted. If `NULL` maximum available Frame is used. Default is `NULL`.
#' @param win Numeric. Window duration scale evaluated in seconds for the moving frequency summary. Default is \code{2}.
#' @param minFreq Numeric. The lower cutoff boundary of the targeted frequency band in Hz. Default is \code{1.5}.
#' @param maxFreq Numeric. The upper cutoff boundary of the targeted frequency band in Hz. Default is \code{7}.
#' @param winCentre Numeric. Seconds for detrending before zero crossings are extracted. 
#' Default is `NULL` translating to same size as `win`. Setting it to \code{0} translates to no centring.
#' @param winSmooth Numeric. Seconds for smoothing to remove noise based on majority presence - only if head gesture is 
#' present in more than half of this time window, this translates to true. Default is \code{0} translating to no smoothing.
#' @param ID.cols Character vector of hex colours. If there are two Identifiers, then two colours must be provided. Default is colourblind-friendly blue and dark green.
#' @param legend Boolean. Switch for the legend. Default is `TRUE`.
#'
#' @return ggplot element
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @seealso \code{\link{featHeadGestures}} \code{\link{preproHead}}
#' @import dplyr
#' @import ggplot2
#' @export
#' 
plotHeadGestures = function(df, colNodding, colShaking, fps, minDegree, 
                            minFrame = NULL, maxFrame = NULL, 
                            minFreq = 1.5, maxFreq = 6.5, win = 2, 
                            winCentre = NULL, winSmooth = 0, 
                            ID.cols = c("#1E88E5", "#004D40"),
                            legend = T) {

  checkDF(df, c("Identifier", "Frame", colNodding, colShaking))
  
  # if none is provided, get the Frame range
  if (is.null(maxFrame)) maxFrame = max(df$Frame)
  if (is.null(minFrame)) minFrame = min(df$Frame)
  
  # ensure that only one dataset
  if (max(df |> count(Frame) |> pull(n)) > 1) stop("This function is for plotting one dataset of one person")
  
  settings = sprintf("Settings: min %.2f degree, %.0f to %.0fHz, %.0fs centering, %.0fs smoothing",
                     minDegree, minFreq, maxFreq, winCentre, winSmooth)
  
  # process the dataframe to extract head gestures
  df = featHeadGestures(df, colNodding, colShaking, fps, minDegree, 
                        win = win, minFreq = minFreq, maxFreq = maxFreq, 
                        winCentre = winCentre, winSmooth = winSmooth, verbose = F)
  
  
  # check if there is a Timecourse in the data, if not create it
  if (!("Timecourse" %in% colnames(df))) {
    df = df |>
      mutate(
        Timecourse = as.numeric(difftime(Timestamp, min(Timestamp), units = "secs"))
      )
  }
  
  # preprocess the dataframe
  df = df |>
    # extract the frames
    filter(Frame >= minFrame & Frame <= maxFrame) |>
    # focus on the relevant columns
    select(Dyad, Identifier, Frame, Timecourse,
           starts_with(c("Nodding", "Shaking")))
  
  # if smoothing, then take the smoothing column
  if (winSmooth > 0) {
    df = df |>
      select(!ends_with("rel")) |>
      rename(
        "Shaking_rel" = "Shaking_rel_smooth",
        "Nodding_rel" = "Nodding_rel_smooth"
      )
  }
  
  # finish preprocessing
  df = df |>
    # wrangle to long format
    tidyr::pivot_longer(cols = starts_with(c("Nodding", "Shaking"))) |>
    mutate(
      name = if_else(!grepl("_", name), paste0(name, "_signal"), name)
    ) |>
    tidyr::separate(name, into = c("Gesture", "name"), sep = "_") |>
    # back to mixed format - Gestures as rows but type of signal as columns
    tidyr::pivot_wider(names_prefix = "V_") |>
    rename("V" = "V_signal")
  
  # get the shift for the raw data / scale for the frequency
  scaleFactor = ceiling(max(df$V)/(maxFreq*1.2))
  
  p = df |> 
    ggplot(aes(x = Timecourse, fill = Gesture)) + 
    geom_hline(yintercept = 0, linewidth = 0.5) + 
    geom_col(aes(y = V_sum, alpha = "All"), na.rm = T, width = 1/fps, 
             position = "identity") + 
    geom_col(data = df |> filter(V_rel == 1), na.rm = T,
             aes(y = V_sum, alpha = "Within Frequency Band"), width = 1/fps, 
             position = "identity") + 
    geom_line(aes(y = V/scaleFactor, colour = Gesture, linetype = "Centered"), 
              linewidth = 1) + 
    geom_line(aes(y = V/scaleFactor, colour = Gesture, linetype = "Input"), 
              linewidth = 1) + 
    geom_vline(data = df |> filter(V_zcT == 1), 
               aes(xintercept = Timecourse, linetype = "Zero Crossing"), alpha = 0.3) + 
    scale_x_continuous(
      breaks = scales::breaks_width(win),
      expand = c(0.02, 0.02)
    ) + 
    scale_fill_manual(values = ID.cols) + 
    scale_colour_manual(values = ID.cols) + 
    scale_alpha_manual(
      name   = "Sums of Zero Crossings",
      values = c("All" = 0.2, "Within Frequency Band" = 0.6)
    ) +
    scale_linetype_manual(
      name   = "Signal",
      values = c("Centered" = "solid", "Input" = "dotted", "Zero Crossing" = "dashed")
    ) + 
    scale_y_continuous(
      name = "Hz", limits = c(-maxFreq*1.2, maxFreq*1.2), 
      sec.axis = sec_axis(transform = ~ . * scaleFactor, name = "Signal")
    ) +
    xlab(sprintf("Seconds (window size %d s)", win)) + 
    facet_wrap(. ~ Gesture, nrow = 2) + 
    theme_bw() + labs(title = sprintf("Head gestures: %s", df$Identifier[1]),
                      subtitle = settings) + 
    theme(legend.position = "bottom", 
          legend.direction = "vertical")
  
  # remove legend if necessary
  if (!legend) p = p + theme(legend.position = "none")
  
  suppressWarnings(return(p))
  
}
