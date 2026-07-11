# Plotting functions for the interactR
# (c) Irene Sophia Plank, 10planki@gmail.com

# if packman is not installed yet, install it
if(!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
pacman::p_load(tidyverse, ggpattern)

# Function to plot the comparison of credible lags. This function takes as input
# the csv files which are created by using the featWLCC pipeline. It creates 
# line plots comparing observed and pseudoWlCC with lags on the x-axis and WLCC
# on the y-axis. Credible lags are highlighted. Facet_wraps is used for the 
# different features in the dataframe, and, if present, for the different Stat
# (either mean or peak). 
#
# Input: 
#  * df : ".*_pseudo-comp.csv" written by featWLCC
#  * ncol : number of columns for facet_wrap
#  * cols.wlcc : named character vector with colours for WLCC and different 
#           types of pseudo-WLCC
#  * cols.cred : named character vector with colours for credbile and not lags
# 
# Output: returns a ggplot object
#
plotWLCCcomp = function(df, ncol = 3,
                        cols.wlcc = c("WLCC" = "#1B9E77", "pseudoSegment-WLCC" = "#FFB6C1", 
                                      "pseudoDyad-WLCC" = "#CAE1FF", "pseudoData-WLCC" = "#C9A0DC"), 
                        cols.cred = c("credible" = "#FFFC80", "not credible" = "white")) {
  
  # get the peak lags
  df.lag = df %>%
    group_by(Stat, Feature) %>%
    filter(observed == max(observed)) %>%
    mutate(
      label = sprintf("Peak: %.2fs", lag),
      ypos  = pseudo - pseudo.sd/2
    )
  
  if (length(unique(df$Stat)) > 1) {
    
    return(df %>%
             pivot_longer(cols = starts_with(c("observed", "pseudo"))) %>%
             mutate(temp = if_else(grepl(".*sd", name), "SD", "mean"),
                    Type = case_when(grepl("observ.*", name) ~ "WLCC", 
                                     Method == "Dyad" ~ "pseudoDyad-WLCC", 
                                     Method == "Seg"  ~ "pseudoSegment-WLCC",
                                     Method == "Data" ~ "pseudoData-WLCC")) %>%
             select(-name) %>%
             pivot_wider(names_from = temp) %>%
             group_by(Type, Feature) %>%
             mutate(yTile = mean(mean, na.rm = T)) %>%
             group_by(Feature, Stat, lag) %>% 
             mutate(credible = case_when(sum(credible == "credible lags") > 0 ~ "credible",
                                         T ~ "not credible")) %>%
             ungroup() %>% 
             ggplot(.) +  
             geom_tile(aes(fill = credible, x = lag, y = yTile, height = Inf), alpha = 0.66) +
             geom_ribbon(aes(x = lag, y = mean, group = Type, fill = Type, 
                             ymin = mean - SD, ymax = mean + SD), 
                         alpha = 0.5) +
             geom_line(aes(x = lag, y = mean, group = Type, colour = Type), linewidth = 1) +
             geom_vline(data = df.lag, aes(xintercept = lag), colour = "black") + 
             geom_label(data = df.lag, aes(x = lag, y = ypos, label = label),  size = 2.5) + 
             facet_wrap(Stat ~ Feature, scale = "free_y", ncol = ncol) +  
             scale_fill_manual(values = c(cols.cred, 
                                          cols.wlcc), 
                               breaks = c("WLCC", "pseudoDyad-WLCC", "pseudoSegment-WLCC", "pseudoData-WLCC", 
                                          "credible")) + 
             scale_colour_manual(values = cols.wlcc, 
                                 guide = "none") +
             theme_bw() + 
             theme(legend.position = "bottom", legend.title = element_blank(),
                   axis.title.y = element_blank())
    )
  } else {
    
    return(df %>%
             pivot_longer(cols = starts_with(c("observed", "pseudo"))) %>%
             mutate(temp = if_else(grepl(".*sd", name), "SD", "mean"),
                    Type = case_when(grepl("observ.*", name) ~ "WLCC", 
                                     Method == "Dyad" ~ "pseudoDyad-WLCC", 
                                     Method == "Seg"  ~ "pseudoSegment-WLCC",
                                     Method == "Data" ~ "pseudoData-WLCC")) %>%
             select(-name) %>%
             pivot_wider(names_from = temp) %>%
             group_by(Type, Feature) %>%
             mutate(yTile = mean(mean, na.rm = T)) %>%
             group_by(Feature, Stat, lag) %>% 
             mutate(credible = case_when(sum(credible == "credible lags") > 0 ~ "credible",
                                         T ~ "not credible")) %>%
             ungroup() %>% 
             ggplot(.) +  
             geom_tile(aes(fill = credible, x = lag, y = yTile, height = Inf), alpha = 0.66) +
             geom_ribbon(aes(x = lag, y = mean, group = Type, fill = Type, 
                             ymin = mean - SD, ymax = mean + SD), 
                         alpha = 0.5) +
             geom_line(aes(x = lag, y = mean, group = Type, colour = Type), linewidth = 1) +
             geom_vline(data = df.lag, aes(xintercept = lag), colour = "black") + 
             geom_label(data = df.lag, aes(x = lag, y = ypos, label = label),  size = 2.5) + 
             facet_wrap(. ~ Feature, scale = "free_y", ncol = ncol) +  
             scale_fill_manual(values = c(cols.cred, 
                                          cols.wlcc), 
                               breaks = c("WLCC", "pseudoDyad-WLCC", "pseudoSegment-WLCC", "pseudoData-WLCC", 
                                          "credible")) + 
             scale_colour_manual(values = cols.wlcc, 
                                 guide = "none") +
             theme_bw() + 
             theme(legend.position = "bottom", legend.title = element_blank(),
                   axis.title.y = element_blank())
    )
  }
  
}