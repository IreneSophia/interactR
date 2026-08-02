#' Plot Zero-Crossings Extracted from Time Series Data
#'
#' Plot the results from \code{\link{featZCrossing}} for one Dyad or one Identifier. 
#'
#' @param df Dataframe. The dataset containing the variables to be processed, potentially created by \code{\link{preproHead}}. 
#'   Must explicitly feature columns `Identifier`, `Frame` and the column `colname`. This dataframe
#'   will be processed using \code{\link{featZCrossing}} to extract relevant zero crossings. 
#'   If `Communication` is a column, Speaking and Listening information is highlighted. If
#'   no `Communication` column is provided, the plot focuses on the different steps to extract
#'   relevant Zero Crossings from the signal. 
#' @param colname Character. The exact name of the column in \code{df} from which
#'   to extract and plot zero-crossing features. 
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
#' @param winSmooth Numeric. Seconds for smoothing to remove noise based on majority presence - only if relevant ZCrossings are 
#' present in more than half of this time window, this translates to true. Default is \code{0} translating to no smoothing.
#' @param ID.cols Character vector of hex colours. If there are two Identifiers, then two colours must be provided. Default is colourblind-friendly blue and dark green.
#' @param legend Boolean. Switch for the legend. Default is `TRUE`.
#'
#' @return ggplot element
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @seealso \code{\link{featZCrossing}} \code{\link{preproHead}}
#' @import dplyr
#' @import ggplot2
#' @export
#' 
plotZCrossings = function(df, colname, fps, minDegree, 
                          minFrame = NULL, maxFrame = NULL, 
                          minFreq = 1.5, maxFreq = 6.5, win = 2, 
                          winCentre = NULL, winSmooth = 0, 
                          ID.cols = c("#1E88E5", "#004D40"),
                          legend = T) {
  
  # if none is provided, get the Frame range
  if (is.null(maxFrame)) maxFrame = max(df$Frame)
  if (is.null(minFrame)) minFrame = min(df$Frame)
  
  # check if Dyad column exists, if not, create it
  if (!("Dyad") %in% colnames(df)) df = df |> mutate(Dyad = "tmp1-dyad")
  
  # process the dataframe to extract Zero Crossings
  df = featZCrossing(df, c(), colname, fps, minDegree,
                     win = win, minFreq = minFreq, maxFreq = maxFreq, 
                     winCentre = winCentre, winSmooth = winSmooth, verbose = F)
  
  # rename the columns depending on whether centred or not
  if (winCentre > 0) {
    df = df |>
      rename_with(~ gsub(paste0(colname, "_centred"), "V", .x), .cols = matches(colname)) |>
      filter(Frame >= minFrame & Frame <= maxFrame)
  } else {
    df = df |>
      rename_with(~ gsub(colname, "V", .x), .cols = matches(colname)) |>
      filter(Frame >= minFrame & Frame <= maxFrame)
  }
  
  # get the shift for the raw data / scale for the frequency
  shift = ceiling(abs(min(df$V)) + maxFreq + 1)
  shift_max = shift + max(df$V)
  scaleFactor = ceiling(max(df$V)/(maxFreq*1.2))
  
  # check whether dyad or solo
  IDs = unique(df$Identifier)
  if (length(IDs) == 2) dyad = TRUE else dyad = FALSE
  if (length(IDs) > 2 | length(IDs) < 1) stop("Function works with solo (one Identifier) or dyad (two Identifiers) data.")
  
  if (dyad & ('Communication' %in% colnames(df))) {
    
    # compute the relevant zero crossings
    df = df |>
      mutate(
        # capture which detected crossing is relevant based on frequency band
        relevant = na_if(V_rel, 0),
        # need to be adjusted to correspond to the frequency - but how???
        zc_win  = (V_sum / win),
        zc_rel = relevant * zc_win
      )
    
    # Communication dataframe for colouring
    df.speak = df |> 
      filter(Identifier == IDs[1]) |>
      arrange(Frame) |>
      mutate(
        Communication = case_when(
          Communication == "Speaking"  ~ sprintf("%s speaking", IDs[1]),
          Communication == "Listening" ~ sprintf("%s speaking", IDs[2]),
          T ~ "None"
        ),
        run_id = data.table::rleid(Communication)) |> 
      group_by(run_id, Communication) |> 
      summarise(
        xmin = min(Frame),
        xmax = max(Frame),
        ymin = 0,
        ymax = shift_max
      ) |> ungroup() |> filter(Communication != "None")
    df.speak = rbind(df.speak |> mutate(Identifier = IDs[1]),
                     df.speak |> mutate(Identifier = IDs[2]))
    
    # loop through both participants
    p = df |>
      ggplot(aes(x = Frame)) +
      # highlight Communication
      geom_rect(data = df.speak, 
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                    fill = Communication),
                alpha = 0.2, inherit.aes = FALSE) + 
      # plot the zero crossing signal
      geom_line(aes(y = zc_win), color = "black", linetype = "dotted") +
      # plot the frequency limits
      geom_hline(yintercept = maxFreq, color = "black", linetype = "dashed") +
      geom_hline(yintercept = minFreq, color = "black", linetype = "dashed")  + 
      # plot the actual signal shifted above
      geom_line(aes(y = V + shift, color = Identifier), linetype = "solid") +
      # highlight zero crossings that fit the frequency band
      geom_point(aes(y = zc_rel, color = Identifier), shape = 20, na.rm = TRUE) +
      xlab("Seconds") + 
      facet_grid(rows = vars(Identifier)) + 
      scale_x_continuous(
        breaks = seq(minFrame, maxFrame, by = win * fps),
        labels = round(seq(minFrame / fps, maxFrame / fps, by = win)),
        limits = c(minFrame, maxFrame),
        expand = c(0.02, 0.02)
      ) +
      scale_y_continuous(
        breaks = c(0, 2, 4, 6, 8, shift),
        labels = c("0", "2", "4", "6", "8", "signal"),
        limits = c(0, shift_max),
        expand = c(0, 0)
      ) + 
      # set the colours
      scale_colour_manual(values = ID.cols) + 
      scale_fill_manual(values = ID.cols) + 
      theme_bw() + labs(title = colname) + 
      theme(legend.position = "bottom", legend.title = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.title.y = element_blank())
  } else if (!dyad & ('Communication' %in% colnames(df))) {
    # get one colour
    ID.cols = ID.cols[1]
    
    # compute the relevant zero crossings
    df = df |>
      mutate(
        # capture which detected crossing is relevant based on frequency band
        relevant = na_if(V_rel, 0),
        # need to be adjusted to correspond to the frequency - but how???
        zc_win  = (V_sum / win),
        zc_rel = relevant * zc_win
      )
    
    # Communication dataframe for colouring
    df.speak = df |> 
      arrange(Frame) |>
      mutate(
        Communication = case_when(
          Communication == "Speaking"  ~ sprintf("%s speaking", IDs),
          Communication == "Listening" ~ sprintf("%s listening", IDs),
          T ~ "None"
        ),
        run_id = data.table::rleid(Communication)) |> 
      group_by(run_id, Communication) |> 
      summarise(
        xmin = min(Frame),
        xmax = max(Frame),
        ymin = 0,
        ymax = shift_max
      ) |> ungroup() |> filter(Communication != "None") |>
      mutate(Identifier = IDs)
    
    # loop through both participants
    p = df |>
      ggplot(aes(x = Frame)) +
      # highlight Communication
      geom_rect(data = df.speak, 
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                    fill = Communication),
                alpha = 0.2, inherit.aes = FALSE) + 
      # plot the zero crossing signal
      geom_line(aes(y = zc_win), color = "black", linetype = "dotted") +
      # plot the frequency limits
      geom_hline(yintercept = maxFreq, color = "black", linetype = "dashed") +
      geom_hline(yintercept = minFreq, color = "black", linetype = "dashed")  + 
      # plot the actual signal shifted above
      geom_line(aes(y = V + shift, color = Identifier), linetype = "solid") +
      # highlight zero crossings that fit the frequency band
      geom_point(aes(y = zc_rel, color = Identifier), shape = 20, na.rm = TRUE) +
      xlab("Seconds") + 
      scale_x_continuous(
        breaks = seq(minFrame, maxFrame, by = win * fps),
        labels = round(seq(minFrame / fps, maxFrame / fps, by = win)),
        limits = c(minFrame, maxFrame),
        expand = c(0.02, 0.02)
      ) +
      scale_y_continuous(
        breaks = c(0, 2, 4, 6, 8, shift),
        labels = c("0", "2", "4", "6", "8", "signal"),
        limits = c(0, shift_max),
        expand = c(0, 0)
      ) + 
      # set the colours
      scale_colour_manual(values = ID.cols) + 
      scale_fill_manual(values = ID.cols) + 
      theme_bw() + labs(title = colname) + 
      theme(legend.position = "bottom", legend.title = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.title.y = element_blank())
  } else {
    if (!dyad) ID.cols = ID.cols[1]
    p = df |> 
      ggplot(aes(x = Frame, fill = Identifier)) + 
      geom_hline(yintercept = 0, linewidth = 0.5) + 
      geom_col(aes(y = V_sum, alpha = "All"), na.rm = T, width = 1) + 
      geom_col(data = df |> filter(V_rel == 1), na.rm = T,
               aes(y = V_sum, alpha = "Within Frequency Band"), width = 1) + 
      geom_line(aes(y = V/scaleFactor, colour = Identifier, linetype = "Centered"), 
                linewidth = 1) + 
      geom_line(aes(y = V/scaleFactor, colour = Identifier, linetype = "Input"), 
                linewidth = 1) + 
      geom_vline(data = df |> filter(V_zc == 1), 
                 aes(xintercept = Frame, linetype = "Zero Crossing"), alpha = 0.3) + 
      scale_x_continuous(
        breaks = seq(minFrame, maxFrame, by = win*fps),
        labels = round(seq(minFrame / fps, maxFrame / fps, by = win)),
        limits = c(minFrame, maxFrame),
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
      theme_bw() + labs(title = colname) + 
      theme(legend.position = "bottom", 
            legend.direction = "vertical")
    if (dyad) p = p + facet_wrap(. ~ Identifier, nrow = length(IDs))
  }
  
  # remove legend if necessary
  if (!legend) p = p + theme(legend.position = "none")
  
  suppressWarnings(return(p))
}
