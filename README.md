# PhyloAMR-heatmap-visualization
PhyloAMR-heatmap-visualization presents R scripts &amp; data to generate a combined phylogenetic tree and AMR presence/absence heatmap for bacterial isolates. It outputs publication-quality PDF/PNG figures &amp; a standalone legend, automating package installation, reproducible visualization, version control compatibility &amp; scalable workflow support. This repository contains an R-script for visualizing antimicrobial resistance (AMR) patterns across samples alongside their phylogenetic relationships. It reads a Newick-format tree & a corresponding binary AMR matrix (0 = absent, 1 = present), producing a composite figure combining a phylogenetic tree & heatmap.

## Input files

* `sample_tree.nwk`: Newick tree with 20 tip labels 
* `amr_matrix.csv`: 20x10 binary matrix of AMR gene presence/absence 

## Features

* Visualizes tree using `ggtree`
* Attaches gene presence/absence heatmap using `gheatmap`


## Required R packages

* CRAN: `ggplot2`, `grid`, `gridExtra`, `cowplot`
* Bioconductor: `ggtree`, `treeio`, `ape`

## Outputs

* `tree_plus_amr_heatmap_noLegend.pdf/png`: Final tree + heatmap
* `amr_legend_only.pdf/png`: AMR legend (standalone)


## Example
![Phylogenetic_tree_with_heatmap.png]

## Sample execution

Run the script in your *R* environment or via terminal:

```bash

Rscript PhyloAMR-heatmap-visualization.R

```

## Citation

If you use this script in your work, please cite:

* Xu et al., *Ggtree: A serialized data object for visualization of a phylogenetic tree and annotation data*, iMeta 2022


---

**Repo Purpose:** Visualize phylogenetic relationships & gene presence for AMR surveillance.
