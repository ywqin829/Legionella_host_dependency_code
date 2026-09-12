# Shared publication theme - LARGE fonts (user requirement)
library(ggplot2)
theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      axis.title = element_text(size = base_size + 2, face = "bold"),
      axis.text = element_text(size = base_size + 1, colour = "black"),
      legend.title = element_text(size = base_size + 1, face = "bold"),
      legend.text = element_text(size = base_size),
      strip.text = element_text(size = base_size + 1, face = "bold"),
      plot.title = element_text(size = base_size + 3, face = "bold", hjust = 0.5),
      panel.grid = element_blank()
    )
}
save_pub_pdf <- function(p, file, w = 7.2, h = 5.5) {
  ggsave(paste0(file, ".pdf"), p, width = w, height = h, device = cairo_pdf)
  message("wrote ", file, ".pdf")
}
