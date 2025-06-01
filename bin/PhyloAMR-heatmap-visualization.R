#!/usr/bin/env Rscript

# plot_tree_with_heatmap_with_mutation_rate.R
# ---------------------------------------------------------------------------------
# Reads:
#   1) "amr_matrix.csv"   — a 20×10 binary matrix (0/1),
#                          rownames = Sample_1..Sample_20,
#                          column names will be forced to "Gene_1".."Gene_10"
#   2) "sample_tree.nwk"  — a Newick tree whose tip labels exactly match Sample_1..Sample_20
#   3) (Optional) "mutation_rates.csv" — a two-column file: tip label and numeric mutation rate
#
# It then:
#   • Installs / loads required CRAN & Bioconductor packages (ggplot2, ggtree, cowplot, etc.)
#   • Builds a rectangular ggtree on the left (with tip labels and optional mutation rate text)
#   • Attaches a perfectly aligned heatmap (white = 0/Absent, tomato = 1/Present)
#     on the right using gheatmap(), and explicitly suppresses any legend.
#   • Removes the horizontal x-axis scale (so no axis line or ticks appear under the heatmap).
#   • Pushes the column names ("Gene_1" … "Gene_10") up higher (colnames_offset_y = 0.04)
#     and shrinks their font (font.size = 2.0), ensuring they do NOT overlap Sample_16.
#   • (Optional) Reads "mutation_rates.csv" and adds rate values as text to tree tips.
#   • Creates a separate “legend only” canvas (0 vs. 1) using cowplot, and saves it as its own file.
#
# Outputs:
#   • tree_plus_amr_heatmap_withRates_noLegend.{pdf,png}
#   • amr_legend_only.{pdf,png}
# ---------------------------------------------------------------------------------

# (0) Ensure that BiocManager exists (for Bioconductor packages):
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# (1) Define required CRAN + Bioconductor packages
cran_pkgs <- c("ggplot2", "grid", "gridExtra", "cowplot")
bioc_pkgs <- c("ggtree", "treeio", "ape")

# (2) Install any missing CRAN packages
installed_cran <- rownames(installed.packages())
for (pkg in cran_pkgs) {
  if (!pkg %in% installed_cran) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
}

# (3) Install any missing Bioconductor packages
for (pkg in bioc_pkgs) {
  if (!suppressWarnings(requireNamespace(pkg, quietly = TRUE))) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# (4) Load all libraries (in correct order to ensure ggtitle() is available)
suppressPackageStartupMessages({
  library(ggplot2)     # must load before ggtree() to ensure ggtitle() works
  library(ggtree)      # for ggtree() and gheatmap()
  library(treeio)      # dependency of ggtree
  library(ape)         # for read.tree()
  library(grid)
  library(gridExtra)
  library(cowplot)     # for get_legend()
})

# (5) Read the 20×10 AMR presence/absence matrix
amr_csv <- "amr_matrix.csv"
if (!file.exists(amr_csv)) {
  stop("Error: 'amr_matrix.csv' not found in the working directory.")
}
amr_df <- read.csv(
  amr_csv,
  row.names       = 1,
  check.names     = FALSE,
  stringsAsFactors = FALSE
)

# Convert to numeric 0/1, and verify
amr_mat <- as.matrix(amr_df)
mode(amr_mat) <- "numeric"
if (!all(amr_mat %in% c(0, 1))) {
  stop("Error: 'amr_matrix.csv' must contain only 0 and 1 values.")
}

# (6) Read the Newick tree
tree_file <- "sample_tree.nwk"
if (!file.exists(tree_file)) {
  stop("Error: 'sample_tree.nwk' not found in the working directory.")
}
tree <- read.tree(tree_file)

# (7) (Optional) Read mutation rates if provided
mutation_file <- "mutation_rates.csv"
mutation_df <- NULL
if (file.exists(mutation_file)) {
  mutation_df <- read.csv(
    mutation_file,
    header          = TRUE,
    stringsAsFactors = FALSE
  )
  # Expect two columns: tip label and numeric rate
  if (!all(c("sample", "rate") %in% colnames(mutation_df))) {
    warning("'mutation_rates.csv' must contain columns 'sample' and 'rate'. Skipping rates.")
    mutation_df <- NULL
  } else {
    rownames(mutation_df) <- mutation_df$sample
  }
}

# (8) Ensure that rownames(amr_mat) exactly match tree tip labels
tips    <- tree$tip.label
samples <- rownames(amr_mat)
if (!identical(sort(samples), sort(tips))) {
  stop(
    "Error: The rownames in 'amr_matrix.csv' do not match the tree's tip labels exactly.\n",
    "       Check: sort(rownames(amr_mat)) vs sort(tree$tip.label)."
  )
}

# (9) Reorder the AMR matrix so its rows follow the tree’s tip order
amr_mat_reordered <- amr_mat[tips, , drop = FALSE]

# (10) Force‐rename columns to "Gene_1" … "Gene_10"
#     (If your CSV already has these names, this simply reassigns them)
n_genes <- ncol(amr_mat_reordered)
colnames(amr_mat_reordered) <- paste0("Gene_", seq_len(n_genes))

# Convert to a data.frame for gheatmap() (rownames must be tip labels)
amr_for_gheat <- as.data.frame(amr_mat_reordered)

# (11) Build the base ggtree with tip labels (with optional mutation-rate text)
p_tree <- ggtree(tree, layout = "rectangular", branch.length = "branch.length") +
  geom_tiplab(
    size     = 2.5,
    align    = TRUE,
    linetype = "dotted",
    linesize = 0.3,
    offset   = 0.02
  )

if (!is.null(mutation_df)) {
  # Add mutation rate next to tip labels
  p_tree <- p_tree +
    geom_treescale() +
    geom_text2(
      aes(label = sprintf("%.3f", mutation_df[.data$label, "rate"])),
      hjust = -0.1,
      size  = 2.5
    )
}

p_tree <- p_tree +
  theme_tree2() +
  ggtitle("Phylogenetic (Clustering) Tree with AMR Heatmap (20 Samples)") +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "none"   # remove any tree-side legend
  )

# (12) Attach the heatmap via gheatmap(), suppressing any legend:
#     • offset            = 0.15   → leave a small gap so the first column doesn’t collide with tip labels
#     • width             = 0.50   → half of the remaining horizontal space is used by the heatmap
#     • colnames_position = "top", colnames_angle = 90 → rotate gene labels vertically
#     • colnames_offset_y = 0.04   → push column names well above the top-row tiles
#     • font.size         = 2.0    → shrink column-name font to avoid overlap
#     • color             = "grey80" → draw a light grid-line around each tile
#     • low               = "white", high = "tomato"  → 0=white (Absent), 1=tomato (Present)
p_heat_noLegend <- gheatmap(
  p_tree,
  amr_for_gheat,
  offset            = 0.15,
  width             = 0.50,
  colnames_position = "top",
  colnames_angle    = 90,
  colnames_offset_y = 0.04,
  font.size         = 2.0,
  color             = "grey80",
  low               = "white",
  high              = "tomato"
) +
  # (13) Remove the horizontal x-axis scale line / ticks under the heatmap area
  theme(
    legend.position      = "none",          # no legend in this combined plot
    axis.title.y         = element_blank(),   # drop any y-axis title
    axis.title.x         = element_blank(),
    axis.text.x          = element_blank(),
    axis.ticks.x         = element_blank(),
    axis.line.x          = element_blank()
  )

# (14) Save the combined tree + heatmap WITHOUT a legend
pdf("tree_plus_amr_heatmap_withRates_noLegend.pdf", width = 10, height = 7)
print(p_heat_noLegend)
dev.off()

png("tree_plus_amr_heatmap_withRates_noLegend.png", width = 1000, height = 700, res = 150)
print(p_heat_noLegend)
dev.off()

# (15) Build a small “dummy” ggplot purely to create a standalone legend
#       We create two tiles: one labeled “0 (Absent)” (white), one “1 (Present)” (tomato).
dummy_df <- data.frame(
  status = factor(c("0 (Absent)", "1 (Present)"),
                  levels = c("0 (Absent)", "1 (Present)")),
  y      = 1
)

legend_plot <- ggplot(dummy_df, aes(x = status, y = y, fill = status)) +
  geom_tile(show.legend = TRUE) +
  scale_fill_manual(
    name   = "AMR Status",
    values = c("0 (Absent)"  = "white",
               "1 (Present)" = "tomato")
  ) +
  guides(fill = guide_legend(
    direction      = "horizontal",
    title.position = "top",
    title.hjust    = 0.5,
    label.position = "right",
    nrow           = 1
  )) +
  theme_void() +
  theme(
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.6, "cm"),
    legend.margin   = margin(t = 5, r = 5, b = 5, l = 5)
  )

# Extract that legend as a grob
legend_grob <- get_legend(legend_plot)

# (16) Save the “legend only” as PDF + PNG
pdf("amr_legend_only.pdf", width = 3, height = 2)
grid.newpage()
grid.draw(legend_grob)
dev.off()

png("amr_legend_only.png", width = 300, height = 200, res = 150)
grid.newpage()
grid.draw(legend_grob)
dev.off()

message(
  "\n✅ Done! Four files generated:\n",
  "   • tree_plus_amr_heatmap_withRates_noLegend.pdf\n",
  "   • tree_plus_amr_heatmap_withRates_noLegend.png\n",
  "   • amr_legend_only.pdf\n",
  "   • amr_legend_only.png\n\n"
)
