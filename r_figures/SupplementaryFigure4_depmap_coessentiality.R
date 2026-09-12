# Supplementary Figure 4. DepMap pan-cancer co-essentiality (exploratory)
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
suppressMessages({library(dplyr); library(readr); library(tidyr)})

DEP <- file.path(ROOT, "01_data_raw/depmap/CRISPRGeneEffect.csv")
machinery <- c("RAB1A", "RAB1B", "RAB10", "RABIF", "SEC22B", "SAR1A", "SEC23A", "SEC24C",
               "SEC31A", "COPA", "COPB1", "COPE", "TRAPPC3", "TRAPPC10", "STX5", "GOSR1",
               "GOLGA2", "VAPA", "RAB8A", "RAB8B")

# read header once to map "SYMBOL (entrez)" -> SYMBOL
hdr <- readLines(DEP, n = 1)
cols <- strsplit(hdr, ",")[[1]]
symmap <- setNames(cols, sub(" \\(.*", "", cols))
use <- intersect(machinery, names(symmap))
# read only the cell-line column + the machinery gene columns
library(data.table)
dt <- fread(DEP, select = c(1, match(symmap[use], cols)), header = TRUE)
dt <- as.data.frame(dt)
colnames(dt)[1] <- "cell_line"
g <- colnames(dt)[-1]
cormat <- cor(dt[g], use = "pairwise.complete.obs")
# order by family groups
C <- as.data.frame(cormat) %>% mutate(gene1 = rownames(cormat)) %>%
  pivot_longer(-gene1, names_to = "gene2", values_to = "r")
C$gene1 <- factor(C$gene1, levels = g)
C$gene2 <- factor(C$gene2, levels = rev(g))

# The correlation matrix is symmetric. Retain one off-diagonal triangle only;
# the diagonal (self-correlation, r = 1) is not informative for co-essentiality.
C <- C %>%
  mutate(i = match(as.character(gene1), g),
         j = match(as.character(gene2), rev(g))) %>%
  filter(i + j < length(g) + 1) %>%
  mutate(label = ifelse(abs(r) >= 0.20, sprintf("%.2f", r), ""))

p <- ggplot(C, aes(gene1, gene2, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 2.6, colour = "black", na.rm = TRUE) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C1121F",
                       midpoint = 0, limits = c(-0.7, 0.7), oob = scales::squish,
                       name = "Pearson r") +
  labs(x = "Rab/ER-Golgi gene", y = "Rab/ER-Golgi gene",
       title = "DepMap co-essentiality of Rab/ER-Golgi machinery (exploratory)",
       subtitle = "One symmetric triangle shown; labels indicate |r| >= 0.20") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9))
save_pub_pdf(p, "SupplementaryFigure4_depmap_coessentiality", w = 8.5, h = 7.5)
