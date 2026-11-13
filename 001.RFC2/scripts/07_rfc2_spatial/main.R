# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
if (TRUE) {
    library(scRNAtoolVis)
    library(tidyverse)
    library(ggplot2)
    library(SpaCET)
    library(Seurat)
    library(ggsci)
    library(dplyr)
    library(scop)
    gene = "RFC2"
    invisible(NULL)
}

sc <- prad_input("input_049")

plot_sc_UMAP_cell(scRNA = sc, plotname = "CellType2", reduction = "umap", label = "CellType2", alphas = 1)

invisible(NULL)

Idents(sc) <- "CellType2"

df.all <- FindAllMarkers(object = sc, only.pos = TRUE, logfc.threshold = 0.585)

df <- df.all %>% dplyr::mutate(diff = pct.1 - pct.2) %>% dplyr::group_by(cluster) %>% dplyr::arrange(desc(diff)) %>% 
    dplyr::slice_head(n = 5)

genes <- df$gene %>% as.character() %>% unique()

mycol <- c(pal_npg()(10), (ggsci::pal_lancet())(10), (ggsci::pal_aaas())(10))

mycol <- mycol[1:length(unique(sc@meta.data$CellType2))]

names(mycol) <- unique(sc@meta.data$CellType2)

invisible(NULL)

averageHeatmap(object = sc, markerGene = genes, annoCol = TRUE, myanCol = mycol)

invisible(NULL)

scop::FeatureStatPlot(srt = sc, stat.by = gene, group.by = "CellType2", add_box = TRUE, legend.position = "none", 
    xlab = "")

invisible(NULL)

plot_sc_UMAP_Feature_3(scRNA = sc, marker = gene)

sc <- prad_input("input_050")

sc[["RNA"]] <- as(sc[["RNA"]], Class = "Assay")

plot_sc_Markers_Simplified(scRNA = sc)

sc@meta.data$CellType2 <- as.character(sc@meta.data$CellType)

sc@meta.data$CellType2[sc@meta.data$seurat_clusters == 20] <- "Prolif"

sc@meta.data$CellType2[sc@meta.data$seurat_clusters == 26] <- "Prolif"

sc@meta.data$CellType2 <- factor(sc@meta.data$CellType2, levels = c("Epithelial_cell", "T_cell", "Macrophage", 
    "Fibroblast", "Endothelial_cell", "Prolif", "Mast_cell", "MSCs", "B_cell"))

table(sc@meta.data$CellType2)

plot_sc_UMAP_cell(scRNA = sc, plotname = "CellType2", reduction = "umap", label = "CellType2", alphas = 1)

invisible(NULL)

Idents(sc) <- "CellType2"

df.all <- FindAllMarkers(object = sc, only.pos = TRUE, logfc.threshold = 0.585, group.by = "CellType2")

df <- df.all %>% dplyr::mutate(diff = pct.1 - pct.2) %>% dplyr::group_by(cluster) %>% dplyr::arrange(desc(diff)) %>% 
    dplyr::slice_head(n = 5)

genes <- df$gene %>% as.character() %>% unique()

mycol <- c(pal_npg()(10), (ggsci::pal_lancet())(10), (ggsci::pal_aaas())(10))

mycol <- mycol[1:length(unique(sc@meta.data$CellType2))]

names(mycol) <- unique(sc@meta.data$CellType2)

invisible(NULL)

ht <- plot_sc_markers_heatmap(object = sc, markerGene = genes)

ComplexHeatmap::draw(ht)

invisible(NULL)

scop::FeatureStatPlot(srt = sc, stat.by = gene, group.by = "CellType2", add_box = TRUE, legend.position = "none", 
    xlab = "")

invisible(NULL)

plot_sc_UMAP_Feature_3(scRNA = sc, marker = gene)

sc <- prad_input("input_051")

sc[["RNA"]] <- as(sc[["RNA"]], Class = "Assay")

scop::FeatureStatPlot(srt = sc, stat.by = gene, group.by = "CellType2", add_box = TRUE, legend.position = "none", 
    xlab = "")

invisible(NULL)

plot_sc_UMAP_Feature_3(scRNA = sc, marker = gene)

invisible(NULL)

SpaCET_obj <- prad_input("input_052")

SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes = 100)

SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType = "PRAD", coreNo = 20)

SpaCET_obj <- SpaCET.CCI.colocalization(SpaCET_obj)

SpaCET_obj <- SpaCET.CCI.LRNetworkScore(SpaCET_obj, coreNo = 20)

SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair = c("", ""))

SpaCET_obj <- SpaCET.identify.interface(SpaCET_obj)

SpaCET_obj <- SpaCET.combine.interface(SpaCET_obj, cellTypePair = c("", ""))

SpaCET_obj <- SpaCET.deconvolution.malignant(SpaCET_obj, coreNo = 20)

SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "Hallmark")

SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "CancerCellState")

SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "TLS")

gmt1 <- list(Tcell = c("CD2", "CD3E", "CD3D"), Myeloid = c("SPI1", "FCER1G", "CSF1R"))

SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets = gmt1)

SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", colors = c("#CCCCCC", "#00599F", "#D01910"), 
    spatialFeatures = "Interface", pointSize = 1.5)

invisible(NULL)

SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "GeneExpression", spatialFeatures = gene, pointSize = 1.5)

invisible(NULL)

dat_exp <- SpaCET_obj@input$counts %>% as.matrix() %>% as.data.frame()

dat_use <- dat_exp[gene, , drop = FALSE] %>% t() %>% as.data.frame()

dat_pos <- SpaCET_obj@results$CCI$interface %>% t() %>% as.data.frame()

com_pos <- intersect(rownames(dat_use), rownames(dat_pos))

dat_pos <- dat_pos[com_pos, , drop = FALSE]

dat_use <- dat_use[com_pos, , drop = FALSE]

dat_use <- cbind(dat_use, dat_pos)

plot_expr_diff_Bar(data = dat_use, Group = "Interface", variable = gene, levels = c("Tumor", "Interface", "Stroma"), 
    alphas = 1)

invisible(NULL)

