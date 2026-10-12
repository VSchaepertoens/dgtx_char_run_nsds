
color_mapping_condition <- c(
  "CT" = "#E6641E",
  "TS" = "#4B288C"
)


theme_tidy <- function() {
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
    plot.background = element_rect(fill = NA, colour = NA),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA),
    # panel.border = element_rect(
    #   fill = NA,
    #   colour = "black",
    #   linewidth = 0.5
    # ),
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    text = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )
}
  