# Supplementary Figure 3. Species composition of evidence for leading candidates
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
SRC <- file.path(ROOT, "13_postfreeze_interpretation")
suppressMessages({library(dplyr); library(tidyr); library(readr); library(patchwork)})

B <- read_tsv(file.path(SRC, "b2_species_dependency_scores.tsv"))
# leading candidates: RAB1A, RAB10 + top tier2 from fragility analysis
leads <- c("RAB1A", "RAB10", "GNA13", "HNRNPA2B1", "HSPA8", "PRRC2C", "ACBD3", "CAP1",
           "EHD4", "ELOVL5", "GPR84", "LDLR", "ACTR3", "VAPA", "RAB7A", "RAB8A", "RAB8B")
B <- B %>% filter(human_symbol %in% leads)

# Panel A: stacked human/mouse expression studies
Bp <- B %>% select(human_symbol, n_human, n_mouse) %>%
  pivot_longer(c(n_human, n_mouse), names_to = "species", values_to = "n") %>%
  mutate(species = ifelse(species == "n_human", "Human", "Mouse"))
pA <- ggplot(Bp, aes(human_symbol, n, fill = species)) +
  geom_col(position = "stack", colour = "black", linewidth = 0.2, width = 0.7) +
  scale_fill_manual(values = c("Human" = "#1B4965", "Mouse" = "#C9A227")) +
  coord_flip() +
  labs(x = NULL, y = "Expression studies", fill = NULL, title = "A. Expression by species") +
  theme_pub() + theme(axis.text.y = element_text(size = 9))

# Panel B: LCV source presence (human WT / ankB / mouse)
# Retain the full table: Panel A uses n_human/n_mouse and Panel B uses lcv_state.
hm <- function(s) grepl("HUMAN", s); ha <- function(s) grepl("ANKB", s); mm <- function(s) grepl("MOUSE", s)
L <- data.frame(gene = B$human_symbol,
                human_WT = as.integer(hm(B$lcv_state) & !ha(B$lcv_state) | grepl("WT_AND_ANKB", B$lcv_state)),
                human_ankB = as.integer(ha(B$lcv_state)),
                mouse = as.integer(mm(B$lcv_state)))
L_long <- L %>%
  pivot_longer(-gene, names_to = "source", values_to = "val") %>%
  mutate(source = factor(source,
                         levels = c("human_WT", "human_ankB", "mouse"),
                         labels = c("Human\nWT", "Human\nankB mutant", "Mouse")))
pB <- ggplot(L_long, aes(source, gene, fill = factor(val))) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#00897B"), name = "Detected") +
  labs(x = "LCV proteomics source", y = NULL, title = "B. LCV-associated source") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5, lineheight = 0.9),
        axis.text.y = element_text(size = 9))

fig <- pA | pB
save_pub_pdf(fig, "SupplementaryFigure3_species_composition", w = 9, h = 6.5)
