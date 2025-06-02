# circular_phylo_tree.R 

# Install required packages if not already installed
packages <- c("ggtree", "ggtreeExtra", "ggplot2", "treeio", "tidyr", "dplyr", "ggnewscale", "RColorBrewer", "viridis")
installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!(pkg %in% installed)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Load necessary libraries
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(treeio)
library(tidyr)
library(dplyr)
library(ggnewscale)
library(RColorBrewer)
library(viridis)

# Load tree
tree <- read.tree("Mtb_tree_50tips.nwk")

# Load metadata
metadata <- read.delim("Mtb_metadata_50tips.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure the label column exists and is character
colnames(metadata)[1] <- "label"
metadata$label <- as.character(metadata$label)

# Match metadata to tree tips
tip_labels <- tree$tip.label
metadata <- metadata %>% filter(label %in% tip_labels)

# If any tips are missing metadata, add rows with default values
missing_tips <- setdiff(tip_labels, metadata$label)
if (length(missing_tips) > 0) {
  metadata <- bind_rows(
    metadata,
    data.frame(
      label = missing_tips,
      Location = rep(NA, length(missing_tips)),
      Year = rep(NA, length(missing_tips)),
      Mtb_Lineage = rep(NA, length(missing_tips)),
      stringsAsFactors = FALSE
    )
  )
}

# Ensure columns exist before coercing to factor
if (!"Location" %in% names(metadata)) metadata$Location <- NA
if (!"Year" %in% names(metadata)) metadata$Year <- NA
if (!"Mtb_Lineage" %in% names(metadata)) metadata$Mtb_Lineage <- NA

# Ensure correct factor levels
metadata$Location <- factor(metadata$Location)
metadata$Mtb_Lineage <- factor(metadata$Mtb_Lineage)

# Filter out unused levels
metadata <- droplevels(metadata)

# Define color palettes
n_location <- length(levels(metadata$Location))
n_lineage <- length(levels(metadata$Mtb_Lineage))
location_palette <- colorRampPalette(brewer.pal(8, "Dark2"))(n_location)
year_palette <- brewer.pal(n = 4, name = "Blues")
lineage_palette <- viridis(n_lineage, option = "D", direction = -1)

# Standardized ring spacing and width
ring_width <- 0.36 * 3
ring_spacing <- ring_width + 0.08

# Build the tree plot
p <- ggtree(tree, layout = "circular", branch.length = "none") %<+% metadata +
  geom_tiplab(size = 2, align = TRUE, linesize = .3, linetype = "solid", offset = ring_spacing * 3, fontface = "bold") +

  # Ring 1: Location
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = Location),
    width = ring_width,
    offset = ring_spacing * 0
  ) +
  scale_fill_manual(
    name = "Location (inner ring)",
    values = location_palette,
    guide = guide_legend(title.theme = element_text(face = "bold"))
  ) +

  # Ring 2: Year Group
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = cut(Year, breaks = 4)),
    width = ring_width,
    offset = ring_spacing * 0.2
  ) +
  scale_fill_manual(
    name = "Year Group (middle ring)",
    values = year_palette,
    guide = guide_legend(title.theme = element_text(face = "bold"))
  ) +

  # Ring 3: MTB Lineage
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = Mtb_Lineage),
    width = ring_width,
    offset = ring_spacing * 0.2
  ) +
  scale_fill_manual(
    name = "Mtb Lineage (outer ring)",
    values = lineage_palette,
    guide = guide_legend(title.theme = element_text(face = "bold"))
  )

# Save output
ggsave("circular_tree_plot_annotated.pdf", plot = p, width = 10, height = 10, dpi = 300)
message("\u2705 Plot saved as circular_tree_plot_annotated.pdf")
