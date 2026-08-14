#' Plot Dwell Time periods
#'
#' Input should be data of one dyad preprocessed using \code{\link{featDwell}}. 
#'
#' @param df Dataframe. The dataset containing the variables to be processed, created by \code{\link{featDwell}}. 
#' @param fps Numeric. Frame processing rate frequency profile (frames per second) of the dataset.
#' @param noSecs Numeric. How many seconds to plot in one facet row. Default is `10`.
#' @param minFrame Numeric or NULL. First Frame to be plotted. If `NULL` minimum available Frame is used. Default is `NULL`.
#' @param maxFrame Numeric or NULL. Last Frame to be plotted. If `NULL` maximum available Frame is used. Default is `NULL`.
#'
#' @return ggplot element
#' 
#' @author Irene Sophia Plank (\email{10planki@@gmail.com})
#' @seealso \code{\link{featDwell}}
#' @import dplyr
#' @import ggplot2
#' @export
#' 
plotETperiods = function(df, fps, noSecs = 10,
                         minFrame = NULL, maxFrame = NULL) {
  
  # if none is provided, get the Frame range
  if (is.null(maxFrame)) maxFrame = max(df$Frame)
  if (is.null(minFrame)) minFrame = min(df$Frame)
  
  # ensure that only one dataset
  if (df |> select(Identifier) |> unique() |> count() |> pull(n) != 2) stop("This function is for plotting the dataset of one dyad")
  
  # check if there is a Timecourse in the data, if not create it
  if (!(all(c("minTimecourse", "maxTimecourse") %in% colnames(df)))) {
    df = df |>
      mutate(
        minTimecourse = minFrame/fps,
        maxTimecourse = maxFrame/fps
      )
  }
  
  # height of the boxes
  boxHeight = 0.8
  
  # preprocess dataframe
  df = df |>
    mutate(
      yCentre = as.numeric(as.factor(Identifier)),
      Chunk = floor(minTimecourse / noSecs) * noSecs,
      maxChunk = Chunk + noSecs - 1/fps,
      ymin = yCentre - boxHeight/2,
      ymax = yCentre + boxHeight/2
      )
  
  # find periods that go across the border
  df.border = df[df$maxTimecourse > df$maxChunk,] |>
    # adjust the minimum timecourse and the chunk to the next one
    mutate(
      Chunk = Chunk + noSecs,
      minTimecourse = Chunk
    )
  
  # adjust the ones in df and merge together
  df.adjusted = df |>
    mutate(maxTimecourse = if_else(maxTimecourse > maxChunk, maxChunk, maxTimecourse)) |>
    rbind(df.border)
  
  # plot this adjusted dataframe
  ggplot(df.adjusted) +
    geom_rect(aes(xmin = minTimecourse, xmax = maxTimecourse, 
                  ymin = ymin, ymax = ymax, 
                  fill = AOI), color = "black", alpha = 0.8) +
    facet_wrap(~ Chunk, ncol = 1, scales = "free_x") +
    scale_y_continuous(
      breaks = df |> select(Identifier, yCentre) |> distinct() |> pull(yCentre),
      labels = df |> select(Identifier, yCentre) |> distinct() |> pull(Identifier),
      limits = c(0.5, max(df |> select(Identifier, yCentre) |> distinct() |> pull(yCentre)) + 0.5)
    ) +
    labs(title = "Dwell time periods", x = "Timecourse (s)") +
    theme_bw() + 
    theme(legend.position = "bottom", legend.title = element_blank(),
          axis.title.y = element_blank())
  
}
