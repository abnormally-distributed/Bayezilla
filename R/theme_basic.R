theme_basic = function (base_size = 16,
                        base_family = "") {
  theme_foundation = function (base_size = 16,
                               base_family = "")
  {
    theme_foundation() + theme(
      line = element_line(
        colour = "black",
        lineend = "round",
        linetype = "solid"
      ),
      rect = element_rect(
        fill = "white",
        colour = "black",
        linetype = "solid"
      ),
      text = element_text(
        colour = "black",
        face = "plain",
        family = base_family,
        size = base_size,
        vjust = 0.5,
        hjust = 0.5,
        lineheight = 1
      ),
      panel.grid = element_blank(),
      strip.background = element_rect(colour = NA),
      legend.key = element_rect(colour = NA),
      title = element_text(size = rel(1)),
      plot.title = element_text(size = rel(1.2),
                                face = "bold"),
      strip.text = element_text(),
      axis.ticks.length = unit(0.5,
                               "lines")
    )
  }



  theme_foundation() + theme(
    line = element_line(
      colour = "black",
      lineend = "round",
      linetype = "solid"
    ),
    rect = element_rect(
      fill = "white",
      colour = "black",
      linetype = "solid"
    ),
    text = element_text(
      colour = "black",
      face = "plain",
      family = base_family,
      size = base_size,
      vjust = 0.5,
      hjust = 0.5,
      lineheight = 1
    ),
    panel.grid = element_blank(),
    strip.background = element_rect(colour = NA),
    legend.key = element_rect(colour = NA),
    title = element_text(size = rel(1)),
    plot.title = element_text(size = rel(1.2),
                              face = "bold"),
    strip.text = element_text(),
    axis.ticks.length = unit(0.5,
                             "lines")
  )



}
