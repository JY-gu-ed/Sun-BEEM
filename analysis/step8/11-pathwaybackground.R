if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step8/"))

infile  <- file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step8/results/2_olink_BEEM_lm_trend_results_FDR.csv")
out_dir <- file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step8/new_pic/Figure3A")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- readr::read_csv(infile, show_col_types = FALSE)

stopifnot(all(c("biomarker", "estimate", "p_adj") %in% names(dat)))

dat2 <- dat %>%
  mutate(
    neglog10_padj = -log10(p_adj),
    sig = !is.na(p_adj) & p_adj < 0.01,
    plot_group = case_when(
      sig & estimate > 0 ~ "Upregulated",
      sig & estimate < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

n_up   <- dat2 %>% filter(sig, estimate > 0) %>% nrow()
n_down <- dat2 %>% filter(sig, estimate < 0) %>% nrow()

sig_up_all <- dat2 %>%
  filter(sig, estimate > 0) %>%
  arrange(p_adj)

sig_down_all <- dat2 %>%
  filter(sig, estimate < 0) %>%
  arrange(p_adj)

sig_all <- bind_rows(sig_up_all, sig_down_all) %>%
  arrange(p_adj)

readr::write_csv(sig_up_all,   file.path(out_dir, "Figure3A_sig_up_all_FDR0.01.csv"))
readr::write_csv(sig_down_all, file.path(out_dir, "Figure3A_sig_down_all_FDR0.01.csv"))
readr::write_csv(sig_all,      file.path(out_dir, "Figure3A_sig_all_FDR0.01_up_down_combined.csv"))

sig_summary <- tibble(
  FDR_threshold = 0.01,
  n_up = n_up,
  n_down = n_down,
  n_total = n_up + n_down
)
readr::write_csv(sig_summary, file.path(out_dir, "Figure3A_sig_summary_counts_FDR0.01.csv"))

top_up <- sig_up_all %>% slice_head(n = 10)
top_down <- sig_down_all %>% slice_head(n = 10)
top_labels <- bind_rows(top_up, top_down)

p_volcano <- ggplot(dat2, aes(x = estimate, y = neglog10_padj)) +
  geom_point(aes(color = plot_group), alpha = 0.8, size = 1.6) +
  scale_color_manual(values = c(
    "Upregulated" = "#D73027",
    "Downregulated" = "#4575B4",
    "Not significant" = "grey80"
  )) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  ggrepel::geom_text_repel(
    data = top_labels,
    aes(label = biomarker),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.2,
    min.segment.length = 0
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("FDR < 0.01\nUp: ", n_up, "\nDown: ", n_down),
    hjust = 1.1, vjust = 1.1, size = 3.5
  ) +
  labs(
    x = "Effect size (β per one-category increase in Sun-BEEM)",
    y = expression(-log[10]("adjusted P value")),
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

print(p_volcano)

ggsave(file.path(out_dir, "Figure3A_volcano_top10up_top10down.png"),
       p_volcano, width = 7, height = 5, dpi = 300)
ggsave(file.path(out_dir, "Figure3A_volcano_top10up_top10down.pdf"),
       p_volcano, width = 7, height = 5)

sym2entrez <- function(symbols) {
  if (length(symbols) == 0) return(tibble(SYMBOL = character(), ENTREZID = character()))
  clusterProfiler::bitr(
    symbols,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  ) %>% distinct(SYMBOL, ENTREZID)
}

sig_up_symbols <- dat2 %>%
  filter(sig, estimate > 0) %>%
  mutate(SYMBOL = toupper(biomarker)) %>%
  distinct(SYMBOL) %>%
  pull(SYMBOL) %>%
  sort()

sig_down_symbols <- dat2 %>%
  filter(sig, estimate < 0) %>%
  mutate(SYMBOL = toupper(biomarker)) %>%
  distinct(SYMBOL) %>%
  pull(SYMBOL) %>%
  sort()

bg_symbols <- dat2 %>%
  mutate(SYMBOL = toupper(biomarker)) %>%
  distinct(SYMBOL) %>%
  pull(SYMBOL) %>%
  sort()

message("Figure 3A significant proteins (SYMBOL):")
message("  Up:   ", length(sig_up_symbols))
message("  Down: ", length(sig_down_symbols))
message("Background (ALL measured SYMBOL): ", length(bg_symbols))

map_up   <- sym2entrez(sig_up_symbols)
map_down <- sym2entrez(sig_down_symbols)
map_bg   <- sym2entrez(bg_symbols)

ent_up   <- unique(map_up$ENTREZID)
ent_down <- unique(map_down$ENTREZID)
ent_bg   <- unique(map_bg$ENTREZID)

message("Mapping to ENTREZ:")
message("  Up:   ", length(ent_up),   " / ", length(sig_up_symbols))
message("  Down: ", length(ent_down), " / ", length(sig_down_symbols))
message("  BG:   ", length(ent_bg),   " / ", length(bg_symbols))

readr::write_csv(tibble(unmapped_symbol = setdiff(sig_up_symbols, map_up$SYMBOL)),
                 file.path(out_dir, "Figure3A_up_unmapped_symbols.csv"))
readr::write_csv(tibble(unmapped_symbol = setdiff(sig_down_symbols, map_down$SYMBOL)),
                 file.path(out_dir, "Figure3A_down_unmapped_symbols.csv"))
readr::write_csv(tibble(unmapped_symbol = setdiff(bg_symbols, map_bg$SYMBOL)),
                 file.path(out_dir, "Figure3A_background_unmapped_symbols.csv"))

run_gobp <- function(entrez_vec, universe_vec) {
  if (length(entrez_vec) < 5) return(NULL)
  enrichGO(
    gene          = entrez_vec,
    universe      = universe_vec,     # <<< KEY: background = all measured
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )
}

run_kegg <- function(entrez_vec, universe_vec) {
  if (length(entrez_vec) < 5) return(NULL)
  enrichKEGG(
    gene          = entrez_vec,
    universe      = universe_vec,     # <<< KEY: background = all measured
    organism      = "hsa",
    keyType       = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20
  )
}

ego_up   <- run_gobp(ent_up,   ent_bg)
ego_down <- run_gobp(ent_down, ent_bg)
ekk_up   <- run_kegg(ent_up,   ent_bg)
ekk_down <- run_kegg(ent_down, ent_bg)

save_enrich <- function(obj, prefix) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) {
    message(prefix, ": no significant terms (or too few genes).")
    return(invisible(NULL))
  }
  readr::write_csv(as.data.frame(obj), file.path(out_dir, paste0(prefix, ".csv")))
}

save_enrich(ego_up,   "Figure3A_UP_GOBP_universe_all_measured")
save_enrich(ego_down, "Figure3A_DOWN_GOBP_universe_all_measured")
save_enrich(ekk_up,   "Figure3A_UP_KEGG_universe_all_measured")
save_enrich(ekk_down, "Figure3A_DOWN_KEGG_universe_all_measured")

plot_dot <- function(obj, title, prefix, n_show = 15) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) return(invisible(NULL))
  p <- dotplot(obj, showCategory = n_show) +
    ggtitle(title) +
    theme_bw(base_size = 12)
  
  ggsave(file.path(out_dir, paste0(prefix, ".png")), p, width = 9, height = 6, dpi = 300)
  ggsave(file.path(out_dir, paste0(prefix, ".pdf")), p, width = 9, height = 6)
}

plot_dot(ego_up,   "GO BP (Upregulated proteins, Sun-BEEM)\nUniverse: all measured",
         "Figure3A_UP_GOBP_dotplot_universe_all_measured")
plot_dot(ego_down, "GO BP (Downregulated proteins, Sun-BEEM)\nUniverse: all measured",
         "Figure3A_DOWN_GOBP_dotplot_universe_all_measured")
plot_dot(ekk_up,   "KEGG (Upregulated proteins, Sun-BEEM)\nUniverse: all measured",
         "Figure3A_UP_KEGG_dotplot_universe_all_measured")
plot_dot(ekk_down, "KEGG (Downregulated proteins, Sun-BEEM)\nUniverse: all measured",
         "Figure3A_DOWN_KEGG_dotplot_universe_all_measured")

message("All Figure 3A outputs saved to: ", out_dir)
