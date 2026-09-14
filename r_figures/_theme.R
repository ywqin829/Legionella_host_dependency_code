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
# Panel labels (A, B, C, ...) at the top-left corner of each panel.
# The panel letter is the first character of the panel title, so left-aligning
# the title anchors it to the top-left corner of the panel instead of centring
# it. patchwork aligns the panels themselves, so letters in one column share a
# single x position and letters in one row share a single y position.
# Apply to the assembled figure, e.g.
#   fig <- (pA + pB + pC + pD + plot_layout(ncol = 2)) & panel_labels_topleft()
#   fig <- (pA | pB) & panel_labels_topleft()
panel_labels_topleft <- function() {
  theme(plot.title = element_text(hjust = 0))
}
