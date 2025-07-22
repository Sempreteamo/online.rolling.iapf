  p <- ggplot() +
    geom_density(
      data = empirical_particles_df,
      aes(x = value, group = time, color = "Standardized Particles"),
      alpha = 1,
      linewidth = 0.1,
      key_glyph = "path" 
    ) +
    stat_function(
      fun = dnorm,
      args = list(mean = 0, sd = 1),
      aes(color = "Standard Normal"),
      linewidth = 0.5,
      key_glyph = "path" 
    ) +
    scale_color_manual(
      name = "Distribution Type",
      values = c("Standardized Particles" = "deepskyblue",
                 "Standard Normal" = "red"),
    ) +
    labs(
      title = paste("d =", d_),
      x = NULL,
      y = NULL
    ) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none", 
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.2),
      panel.grid.minor = element_blank()
    )
  
  plot_list[[as.character(d_)]] <- p
}


library(dplyr)
library(ggplot2)
library(patchwork) 


y_axis_label_plot <- ggplot() +
  labs(y = "Density") + 
  theme_void() + 
  theme(
    axis.title.y = element_text(size = 16, angle = 90, face = "bold", margin = margin(r = 15)), 
    plot.margin = margin(0,0,0,0) 
  )



final_plot <- y_axis_label_plot + wrap_plots(plot_list, ncol = 3) +
  plot_layout(
    widths = c(0.05, 1), 
    guides = "collect"   
  ) +
  plot_annotation(
 
    caption = "Standardized Value", 
    theme = theme(
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10), 
      plot.caption = element_text(size = 16, hjust = 0.5, face = "bold", margin = margin(t = 15)) 
    )
  ) & 
  theme(
    legend.position = "bottom" 
  )


print(final_plot)
