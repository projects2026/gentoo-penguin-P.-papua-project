# VISUALIZATION of individual and population niche with ellipses area #

library(ggplot2)
#Load data
DF=read.csv("ellipses_data.csv", header=TRUE, sep=";", dec=",")

#DF structure is: 	
#iso1 = d13C
#iso2 = d15N
#group = 0 for population values (that is, all individual values); 1:17 for individual values
#community = same as for group

#Create a variable to facilitate aesthetic mapping:
DF$ellipse_group <- "Others"
DF$ellipse_group[DF$group == 0]  <- "Population"
DF$ellipse_group[DF$group == 5]  <- "ID 5"
DF$ellipse_group[DF$group == 16] <- "ID 16"
DF$ellipse_group[DF$group == 17] <- "ID 17"

DF$ellipse_group <- factor(
  DF$ellipse_group,
  levels = c("Population", "ID 5", "ID 16", "ID 17", "Others")
)

# Figure
# Note that stat_ellipse is used multiple times to control the drawing order.

ellipse_plot <- ggplot(DF, aes(x = iso1, y = iso2)) +
  ## Points 
  geom_jitter(
    aes(color = ellipse_group),
    size = 3.5
  ) +
  
  ## 1) Others: IDs drawn in the background
  stat_ellipse(
    data = subset(DF, ellipse_group == "Others"),
    aes(group = group),
    type = "norm",
    geom = "polygon",
    level = 0.40,
    color = "lightgray",
    fill  = "lightgray",
    alpha = 0.1,
    linewidth  = 1
  ) +
  
  ## 2) ID 16 (darkmagenta)
  stat_ellipse(
    data = subset(DF, ellipse_group == "ID 16"),
    aes(group = group, color = ellipse_group, fill = ellipse_group),
    type = "norm",
    geom = "polygon",
    level = 0.40,
    alpha = 0.2,
    linewidth = 1
  ) +
  
  ## 3) ID 17 (darkcyan)
  stat_ellipse(
    data = subset(DF, ellipse_group == "ID 17"),
    aes(group = group, color = ellipse_group, fill = ellipse_group),
    type = "norm",
    geom = "polygon",
    level = 0.40,
    alpha = 0.2,
    linewidth = 1
  ) +
  
  ## 4) Population (black dashed outline)
  stat_ellipse(
    data = subset(DF, ellipse_group == "Population"),
    aes(group = group, color = ellipse_group, fill = ellipse_group),
    type = "norm",
    geom = "polygon",
    level = 0.40,
    alpha = 0.3,
    linewidth = 1,
    linetype = "dashed"
  ) +
  
  ## 5) ID 5 (forestgreen)
  stat_ellipse(
    data = subset(DF, ellipse_group == "ID 5"),
    aes(group = group, color = ellipse_group, fill = ellipse_group),
    type = "norm",
    geom = "polygon",
    level = 0.40,
    alpha = 0.2,
    linewidth = 1
  ) +
  
  ## Scales
  scale_color_manual(
    values = c(
      "Population" = "black",
      "ID 5"      = "forestgreen",
      "ID 16"     = "darkmagenta",
      "ID 17"     = "darkcyan",
      "Others"     = "lightgray"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Population" = "black",
      "ID 5"      = "forestgreen",
      "ID 16"     = "darkmagenta",
      "ID 17"     = "darkcyan",
      "Others"     = "lightgray"
    )
  ) +
  
  ## Theme
  labs(
    x = expression({delta}^13*C~'\u2030'),
    y = expression({delta}^15*N~'\u2030')
  ) +
  scale_y_continuous(
    breaks = seq(6, 9, by = 0.5),
    labels = function(x) sprintf("%.1f", x)
  ) +
  coord_cartesian(ylim = c(6, 9), xlim=c(-25.35,-23.0))+
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.title.x = element_text(size = 28, colour = "black"),
    axis.title.y = element_text(size = 28, colour = "black"),
    axis.text.x  = element_text(size = 26, colour = "black"),
    axis.text.y  = element_text(size = 26, colour = "black"),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size=20))

ellipse_plot
