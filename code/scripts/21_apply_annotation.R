#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
source("scripts/lib_common.R")
theme_set(theme_nature_contract())

args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[[1]] else "data/processed/scrna/scrna_annotation_evidence.rds"
map_path <- if (length(args) >= 2) args[[2]] else "config/annotation_cluster_map.tsv"
processed_root <- if (length(args) >= 3) args[[3]] else "data/processed"
results_root <- if (length(args) >= 4) args[[4]] else "results"

object <- readRDS(input)
mapping <- utils::read.delim(map_path, check.names = FALSE)
observed <- sort(unique(as.character(object$seurat_clusters)))
validate_annotation_map(mapping, observed)
lookup <- setNames(seq_len(nrow(mapping)), as.character(mapping$cluster))
indices <- lookup[as.character(object$seurat_clusters)]
object$cell_type <- mapping$cell_type[indices]
object$cell_state_label <- mapping$cell_state_label[indices]
object$annotation_confidence <- mapping$confidence[indices]

DefaultAssay(object) <- "SoupX"
Idents(object) <- "cell_type"
saveRDS(object, file.path(processed_root, "scrna", "scrna_annotated.rds"), compress = FALSE)

metadata <- object[[]]
metadata$cell_id <- rownames(metadata)
write_tsv(metadata, file.path(results_root, "tables", "annotation", "annotated_cell_metadata.tsv.gz"))

counts <- composition_from_cells(object@meta.data)
write_tsv(counts, file.path(results_root, "tables", "composition", "cell_type_counts_and_proportions.tsv"))

marker_panel <- unique(unlist(strsplit(paste(mapping$canonical_marker_evidence, collapse = ","), ",")))
marker_panel <- intersect(marker_panel, rownames(object))
cell_groups <- split(seq_len(ncol(object)), object$cell_type)
log_matrix <- GetAssayData(object, assay = "SoupX", layer = "data")
count_matrix <- GetAssayData(object, assay = "SoupX", layer = "counts")
average <- sapply(cell_groups, function(index) Matrix::rowMeans(log_matrix[marker_panel, index, drop = FALSE]))
percent <- sapply(cell_groups, function(index) Matrix::rowMeans(count_matrix[marker_panel, index, drop = FALSE] > 0))
specificity_rows <- expand.grid(gene = marker_panel, cell_type = colnames(average), stringsAsFactors = FALSE)
specificity_rows$average_log_expression <- mapply(function(gene, type) average[gene, type], specificity_rows$gene, specificity_rows$cell_type)
specificity_rows$percent_expressing <- mapply(function(gene, type) percent[gene, type], specificity_rows$gene, specificity_rows$cell_type)
write_tsv(specificity_rows, file.path(results_root, "tables", "annotation", "marker_specificity_by_cell_type.tsv"))

umap <- as.data.frame(Embeddings(object, "umap_unintegrated"))
umap$cell_id <- rownames(umap)
umap$cell_type <- object$cell_type
umap$patient_id <- object$patient_id
umap$side <- object$side
write_tsv(umap, file.path(results_root, "tables", "annotation", "annotated_umap.tsv"))

cell_type_levels <- c(
  "fibro-adipogenic progenitor", "muscle stem/progenitor cell", "macrophage/monocyte",
  "T/NK cell", "B/plasma cell", "endothelial cell", "pericyte",
  "smooth-muscle-like cell", "Schwann/perineurial cell"
)
type_palette <- setNames(c("#D8A65D", "#7B61A8", "#C84A4A", "#4D78B5", "#70A6D2",
                           "#33A69A", "#E28E2C", "#8C6D31", "#8A8A8A"), cell_type_levels)

p_umap <- ggplot(umap, aes(x = umapunintegrated_1, y = umapunintegrated_2, colour = cell_type)) +
  geom_point(size = 0.12, alpha = 0.75) + coord_equal() +
  scale_colour_manual(values = type_palette) +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  labs(x = "UMAP 1", y = "UMAP 2", colour = NULL, title = "Conservative cell-type annotation")

counts$side <- factor(counts$side, levels = c("convex", "concave"))
p_comp <- ggplot(counts, aes(x = sample_id, y = proportion, fill = cell_type)) +
  geom_col(width = 0.78) + scale_fill_manual(values = type_palette) +
  guides(fill = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Cell proportion", fill = NULL, title = "Per-sample composition") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

p_pair <- ggplot(counts, aes(x = side, y = proportion, group = patient_id, colour = patient_id)) +
  geom_line(linewidth = 0.35) + geom_point(size = 1.1) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 3) +
  scale_x_discrete(limits = c("convex", "concave")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = NULL, y = "Cell proportion", colour = "Patient", title = "Paired composition by cell type") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

figure2 <- p_umap + p_comp + p_pair +
  plot_layout(design = "AAB\nCCC", heights = c(1, 1.15), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"), legend.position = "right")
save_pub_r(figure2, file.path(results_root, "figures", "Figure_2_cell_atlas_composition"), 183, 150)
message("Annotation applied and Figure 2 exported")
