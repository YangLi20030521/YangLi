# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
if (TRUE) {
    library(ComplexHeatmap)
    library(patchwork)
    library(tidyverse)
    library(dorothea)
    library(circlize)
    library(CellChat)
    library(pheatmap)
    library(ggplot2)
    library(Seurat)
    library(irGSEA)
    library(dplyr)
    library(scop)
    gene = "RFC2"
    invisible(NULL)
    sc <- prad_input("input_062")
}

sc_prolif <- sc[, sc@meta.data$CellType2 == "Prolif"]

sc_prolif <- do_sc_p_n_gene(sc = sc_prolif, gene = gene, expression = 1)

plot_sc_CellType_Cortplot(sc = sc_prolif, idents = "status", version = "v4")

sc_non <- sc[, sc@meta.data$CellType2 != "Prolif" & sc@meta.data$CellType2 == "Epithelials"]

sc_non@meta.data$CellType3 <- sc_non@meta.data$CellType2 %>% as.character()

sc_prolif@meta.data$CellType3 <- sc_prolif@meta.data$status %>% as.character()

sc_prolif@meta.data$CellType3[sc_prolif@meta.data$CellType3 == "RFC2+"] <- "H_RFC_Prolif"

sc_prolif@meta.data$CellType3[sc_prolif@meta.data$CellType3 == "RFC2-"] <- "L_RFC_Prolif"

sc <- merge(sc_prolif, sc_non)

table(sc@meta.data$orig.ident)

sc <- do_sc_Pipeline(scRNA = sc, resolution = 0.8, integration = "css")

invisible(NULL)

invisible(NULL)

sc <- prad_input("input_063")

table(sc@meta.data$scMalignantFinder_prediction)

plot_sc_UMAP_cell(scRNA = sc, plotname = "scMalignantFinder_prediction", reduction = "umap", label = "scMalignantFinder_prediction", 
    alphas = 1)

invisible(NULL)

sc <- do_sc_GIS(sc = sc, version = "v4")

scop::FeatureStatPlot(srt = sc, stat.by = "GIS", legend.position = "none", xlab = "", comparisons = list(c("Malignant", 
    "Normal")), add_box = TRUE, group.by = "scMalignantFinder_prediction", palette = "npg")

invisible(NULL)

sc_mal <- sc[, sc@meta.data$scMalignantFinder_prediction == "Malignant"]

sc_nonmal <- sc[, sc@meta.data$scMalignantFinder_prediction == "Normal"]

table(sc_mal@meta.data$orig.ident, sc_mal@meta.data$CellType3)

sc_mal <- sc_mal[, sc_mal@meta.data$orig.ident != "GSM4089153"]

sc_mal@meta.data$CellType3 <- as.character(sc_mal@meta.data$CellType3)

sc_mal@meta.data$CellType3[sc_mal@meta.data$CellType3 == "Epithelials"] <- "Malignants"

sc_mal <- do_sc_Pipeline(scRNA = sc_mal, resolution = 0.8, integration = "css")

plot_sc_UMAP_cell(scRNA = sc_mal, plotname = "CellType3", reduction = "umap", label = "CellType3", alphas = 1)

invisible(NULL)

plot_sc_CellType_Cortplot(sc = sc_mal, idents = "CellType3", version = "v4")

res_cyto <- do_sc_CytoTRACE2(sc = sc_mal, species = "human", version = "v4")

table(res_cyto@meta.data$preKNN_CytoTRACE2_Potency, res_cyto@meta.data$CellType3)

scop::CellDimPlot(srt = res_cyto, group.by = "CellType3")

scop::FeatureDimPlot(srt = res_cyto, features = "CytoTRACE2_Score")

res_cyto@meta.data$CellType2 <- as.character(res_cyto@meta.data$CellType2)

res_cyto@meta.data$CellType2[res_cyto@meta.data$CellType2 == "Epithelials"] <- "Malignant"

res_cyto@meta.data$CellType2 <- factor(res_cyto@meta.data$CellType2, levels = c("Prolif", "Malignant"))

scop::FeatureStatPlot(srt = res_cyto, stat.by = "CytoTRACE2_Score", group.by = "CellType2", add_box = TRUE, comparisons = list(c("Prolif", 
    "Malignant")), legend.position = "none", palette = "npg", xlab = "")

invisible(NULL)

res_cyto@meta.data$CellType3 <- as.character(res_cyto@meta.data$CellType3)

res_cyto@meta.data$CellType3 <- factor(res_cyto@meta.data$CellType3, levels = c("H_RFC_Prolif", "L_RFC_Prolif", 
    "Malignant"))

scop::FeatureStatPlot(srt = res_cyto, stat.by = "CytoTRACE2_Score", group.by = "CellType3", add_box = TRUE, comparisons = list(c("H_RFC_Prolif", 
    "L_RFC_Prolif"), c("H_RFC_Prolif", "Malignant")), legend.position = "none", palette = "npg", xlab = "")

invisible(NULL)

plot_sc_UMAP_cell(scRNA = res_cyto, plotname = "CytoTRACE2_Potency", reduction = "umap", label = "CytoTRACE2_Potency", 
    alphas = 1)

invisible(NULL)

invisible(NULL)

res_irGSEA <- irGSEA.score(object = sc_mal, assay = "RNA", slot = "data", seeds = 123, ncores = 25, min.cells = 3, 
    min.feature = 0, custom = F, geneset = NULL, msigdb = T, species = "Homo sapiens", category = "H", subcategory = NULL, 
    geneid = "symbol", method = c("AUCell"), aucell.MaxRank = NULL, ucell.MaxRank = NULL, kcdf = "Gaussian")

dat_hallmark <- res_irGSEA@assays$AUCell$data %>% as.matrix() %>% as.data.frame()

res_irGSEA <- irGSEA.score(object = sc_mal, assay = "RNA", slot = "data", seeds = 123, ncores = 25, min.cells = 3, 
    min.feature = 0, custom = F, geneset = NULL, msigdb = T, species = "Homo sapiens", category = "C5", subcategory = "GO:BP", 
    geneid = "symbol", method = c("AUCell"), aucell.MaxRank = NULL, ucell.MaxRank = NULL, kcdf = "Gaussian")

dat_gobp <- res_irGSEA@assays$AUCell$data %>% as.matrix() %>% as.data.frame()

res_irGSEA <- irGSEA.score(object = sc_mal, assay = "RNA", slot = "data", seeds = 123, ncores = 25, min.cells = 3, 
    min.feature = 0, custom = F, geneset = NULL, msigdb = T, species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY", 
    geneid = "symbol", method = c("AUCell"), aucell.MaxRank = NULL, ucell.MaxRank = NULL, kcdf = "Gaussian")

dat_kegg <- res_irGSEA@assays$AUCell$data %>% as.matrix() %>% as.data.frame()

dat_kegg <- as.data.frame(t(dat_kegg))

identical(rownames(dat_kegg), rownames(sc_mal@meta.data))

dat_kegg <- dat_kegg %>% dplyr::mutate(CellType = sc_mal@meta.data$CellType3)

dat_kegg$CellType <- factor(dat_kegg$CellType, levels = c("H_RFC_Prolif", "L_RFC_Prolif", "Malignant"))

celltype_means <- dat_kegg %>% dplyr::group_by(CellType) %>% dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
    column_to_rownames("CellType") %>% t() %>% as.data.frame()

other_avg <- rowMeans(celltype_means[, c("L_RFC_Prolif", "Malignant")])

significant_high <- celltype_means[celltype_means$H_RFC_Prolif > other_avg * 1.2, ]

significant_high$Other_Avg <- other_avg[rownames(significant_high)]

significant_high <- significant_high[order(significant_high$H_RFC_Prolif, decreasing = TRUE), ]

significant_high$Other_Avg <- NULL

str(significant_high)

pheatmap(mat = as.matrix(significant_high), scale = "row", color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), 
    cluster_rows = TRUE, cluster_cols = FALSE, clustering_distance_rows = "euclidean", clustering_method = "complete", 
    show_rownames = TRUE, show_colnames = TRUE, fontsize = 11, fontsize_row = 10, fontsize_col = 11, main = "", 
    legend = FALSE, border_color = "black", cellwidth = 25, cellheight = 20, treeheight_row = 0, treeheight_col = 0, 
    display_numbers = FALSE, silent = FALSE)

dat_hallmark <- as.data.frame(t(dat_hallmark))

identical(rownames(dat_hallmark), rownames(sc_mal@meta.data))

dat_hallmark <- dat_hallmark %>% dplyr::mutate(CellType = sc_mal@meta.data$CellType3)

dat_hallmark$CellType <- factor(dat_hallmark$CellType, levels = c("H_RFC_Prolif", "L_RFC_Prolif", "Malignant"))

celltype_means <- dat_hallmark %>% dplyr::group_by(CellType) %>% dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
    column_to_rownames("CellType") %>% t() %>% as.data.frame()

other_avg <- rowMeans(celltype_means[, c("L_RFC_Prolif", "Malignant")])

significant_high <- celltype_means[celltype_means$H_RFC_Prolif > other_avg * 1.011, ]

significant_high$Other_Avg <- other_avg[rownames(significant_high)]

significant_high <- significant_high[order(significant_high$H_RFC_Prolif, decreasing = TRUE), ]

significant_high$Other_Avg <- NULL

str(significant_high)

pheatmap(mat = as.matrix(significant_high), scale = "row", color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), 
    cluster_rows = TRUE, cluster_cols = FALSE, clustering_distance_rows = "euclidean", clustering_method = "complete", 
    show_rownames = TRUE, show_colnames = TRUE, fontsize = 11, fontsize_row = 10, fontsize_col = 11, main = "", 
    legend = FALSE, border_color = "black", cellwidth = 25, cellheight = 20, treeheight_row = 0, treeheight_col = 0, 
    display_numbers = FALSE, silent = FALSE)

dat_gobp <- as.data.frame(t(dat_gobp))

identical(rownames(dat_gobp), rownames(sc_mal@meta.data))

dat_gobp <- dat_gobp %>% dplyr::mutate(CellType = sc_mal@meta.data$CellType3)

dat_gobp$CellType <- factor(dat_gobp$CellType, levels = c("H_RFC_Prolif", "L_RFC_Prolif", "Malignant"))

celltype_means <- dat_gobp %>% dplyr::group_by(CellType) %>% dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
    column_to_rownames("CellType") %>% t() %>% as.data.frame()

other_avg <- rowMeans(celltype_means[, c("L_RFC_Prolif", "Malignant")])

significant_high <- celltype_means[celltype_means$H_RFC_Prolif > other_avg * 2.55, ]

significant_high$Other_Avg <- other_avg[rownames(significant_high)]

significant_high <- significant_high[order(significant_high$H_RFC_Prolif, decreasing = TRUE), ]

significant_high$Other_Avg <- NULL

str(significant_high)

pheatmap(mat = as.matrix(significant_high), scale = "row", color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), 
    cluster_rows = TRUE, cluster_cols = FALSE, clustering_distance_rows = "euclidean", clustering_method = "complete", 
    show_rownames = TRUE, show_colnames = TRUE, fontsize = 11, fontsize_row = 10, fontsize_col = 11, main = "", 
    legend = FALSE, border_color = "black", cellwidth = 25, cellheight = 20, treeheight_row = 0, treeheight_col = 0, 
    display_numbers = FALSE, silent = FALSE)

sc_mal@meta.data$CellType3 <- factor(sc_mal@meta.data$CellType3, levels = c("H_RFC_Prolif", "L_RFC_Prolif", "Malignant"))

scop::FeatureStatPlot(srt = sc_mal, stat.by = "STMN1", group.by = "CellType3", add_box = TRUE, comparisons = list(c("H_RFC_Prolif", 
    "L_RFC_Prolif"), c("H_RFC_Prolif", "Malignant")), legend.position = "none", palette = "npg", xlab = "")

invisible(NULL)

scop::FeatureStatPlot(srt = sc_mal, stat.by = "FEN1", group.by = "CellType3", add_box = TRUE, comparisons = list(c("H_RFC_Prolif", 
    "L_RFC_Prolif"), c("H_RFC_Prolif", "Malignant")), legend.position = "none", palette = "npg", xlab = "")

invisible(NULL)

scop::FeatureStatPlot(srt = sc_mal, stat.by = "BAX", group.by = "CellType3", add_box = TRUE, comparisons = list(c("H_RFC_Prolif", 
    "L_RFC_Prolif"), c("H_RFC_Prolif", "Malignant")), legend.position = "none", palette = "npg", xlab = "")

invisible(NULL)

sc <- prad_input("input_064")

sc_nonmal@meta.data$CellType3 <- "Nomal_Epi"

sc_non <- sc[, !sc@meta.data$CellType2 %in% c("Prolif", "Epithelials")]

sc_non@meta.data$CellType3 <- sc_non@meta.data$CellType2

sc_mal@meta.data$CellType3 <- as.character(sc_mal@meta.data$CellType3)

sc_mal@meta.data$CellType3[sc_mal@meta.data$CellType3 == "H_RFC_Prolif"] <- "H_RFC2_Prolif"

sc_mal@meta.data$CellType3[sc_mal@meta.data$CellType3 == "L_RFC_Prolif"] <- "L_RFC2_Prolif"

sc_list <- list(sc_non = sc_non, sc_mal = sc_mal, sc_nonmal = sc_nonmal)

sc <- merge(sc_list[[1]], sc_list[-1])

sc <- do_sc_Pipeline(scRNA = sc, resolution = 0.8, integration = "harmony")

invisible(NULL)

cellchat <- createCellChat(object = sc@assays$RNA@data, meta = sc@meta.data, group.by = "CellType3")

cellchat <- setIdent(cellchat, ident.use = "CellType3")

groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human

CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"), key = "annotation")

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)

cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

netVisual_heatmap(cellchat, measure = "weight")

netVisual_heatmap(cellchat, measure = "count")

dat_plot1 <- cellchat@net$weight %>% as.data.frame()

dat_plot_long <- as.data.frame(dat_plot1) %>% rownames_to_column("row_names")

dat_plot_long <- dat_plot_long %>% pivot_longer(cols = -row_names, names_to = "col_names", values_to = "number")

dat_plot2 <- dat_plot_long %>% dplyr::mutate(indexa = row_names, indexb = col_names, total_number = NA) %>% dplyr::select(indexa, 
    indexb, number, total_number)

dat_plot <- dat_plot2

dat_plot$total_number <- round(dat_plot$number, digits = 2)

dat_plot <- statdf %>% dplyr::filter(indexa %in% c("B_cells", "Endothelials", "Mast_cells", "Malignant", "L_RFC2_Prolif", 
    "H_RFC2_Prolif", "Fibroblasts", "Myeloids", "NK_cells", "Normal_Epi", "Plasmas", "SMCs", "T_cells"))

dat_plot <- dat_plot %>% dplyr::filter(indexb %in% c("B_cells", "Endothelials", "Mast_cells", "Malignant", "L_RFC2_Prolif", 
    "H_RFC2_Prolif", "Fibroblasts", "Myeloids", "NK_cells", "Normal_Epi", "Plasmas", "SMCs", "T_cells"))

limits.max <- max(dat_plot$total_number)

color.palette <- c("#4393C3", "#ffdbba", "#B2182B")

size.of.text <- c(18, 14, 14, 12)

indexa_levels <- c("B_cells", "Endothelials", "Mast_cells", "Malignant", "L_RFC2_Prolif", "H_RFC2_Prolif", "Fibroblasts", 
    "Myeloids", "NK_cells", "Normal_Epi", "Plasmas", "SMCs", "T_cells")

indexb_levels <- c("B_cells", "Endothelials", "Mast_cells", "Malignant", "L_RFC2_Prolif", "H_RFC2_Prolif", "Fibroblasts", 
    "Myeloids", "NK_cells", "Normal_Epi", "Plasmas", "SMCs", "T_cells")

dat_plot$indexb <- factor(dat_plot$indexb, levels = indexb_levels)

dat_plot$indexa <- factor(dat_plot$indexa, levels = indexa_levels)

p1 <- ggplot(dat_plot, aes(x = indexa, y = indexb, fill = total_number)) + geom_tile(color = "white") + geom_text(aes(label = total_number), 
    color = "white", size = 5) + scale_fill_gradientn(colours = color.palette, limits = c(0, limits.max), name = "Total Number") + 
    scale_x_discrete("Sender") + scale_y_discrete("Receivor") + theme_minimal() + theme(axis.title = element_text(size = size.of.text[1]), 
    axis.text.x.bottom = element_text(hjust = 1, angle = 45, size = size.of.text[2], color = "black"), axis.text.y.left = element_text(size = size.of.text[2], 
        color = "black"), legend.title = element_text(size = size.of.text[3]), legend.text = element_text(size = size.of.text[4]), 
    panel.grid = element_blank(), legend.position = "none", plot.margin = margin(r = 0))

summary_dat_plot <- dat_plot %>% dplyr::group_by(indexb) %>% dplyr::summarise(total_sum = sum(total_number))

summary_dat_plot$indexb <- factor(summary_dat_plot$indexb, levels = indexb_levels)

max_sum <- max(summary_dat_plot$total_sum)

p2 <- ggplot(summary_dat_plot, aes(x = total_sum, y = indexb)) + geom_col(fill = "#0072B5FF") + scale_x_continuous("") + 
    scale_y_discrete("") + geom_text(aes(label = total_sum), hjust = -0.1, size = 5) + coord_cartesian(xlim = c(0, 
    max_sum * 1.1)) + theme_minimal() + theme(axis.title.x = element_text(size = size.of.text[1]), axis.title.y = element_blank(), 
    axis.text.x = element_text(size = size.of.text[2], color = "black"), axis.text.y = element_blank(), panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(), panel.grid.major.x = element_line(linetype = "dotted"), panel.grid.minor.x = element_blank(), 
    plot.margin = margin(l = 1))

p1 + p2 + plot_layout(widths = c(3, 1))

invisible(NULL)

counts <- as.matrix(sc@assays$RNA@data)

invisible(NULL)

meta_data <- cbind(rownames(sc@meta.data), sc@meta.data[, "CellType3", drop = FALSE]) %>% as.matrix()

colnames(meta_data) <- c("barcode", "celltype")

invisible(NULL)

invisible(NULL)

pvalues = prad_input("input_065")

pvalues = pvalues[, 12:dim(pvalues)[2]]

statdf = as.data.frame(colSums(pvalues < 0.05))

colnames(statdf) = c("number")

statdf$indexb = stringr::str_replace(rownames(statdf), "^.*\\|", "")

statdf$indexa = stringr::str_replace(rownames(statdf), "\\|.*$", "")

statdf$total_number = 0

for (i in 1:dim(statdf)[1]) {
    tmp_indexb = statdf[i, "indexb"]
    tmp_indexa = statdf[i, "indexa"]
    if (tmp_indexa == tmp_indexb) {
        statdf[i, "total_number"] = statdf[i, "number"]
    }
    else {
        statdf[i, "total_number"] = statdf[statdf$indexb == tmp_indexb & statdf$indexa == tmp_indexa, "number"] + 
            statdf[statdf$indexa == tmp_indexb & statdf$indexb == tmp_indexa, "number"]
    }
}

rankname = sort(unique(statdf$indexa))

statdf$indexa = factor(statdf$indexa, levels = rankname)

statdf$indexb = factor(statdf$indexb, levels = rankname)

dat_plot <- statdf %>% dplyr::filter(indexa %in% c("B_cells", "Endothelials", "Mast_cells", "Malignant", "L_RFC2_Prolif", 
    "H_RFC2_Prolif", "Fibroblasts", "Myeloids", "NK_cells", "Normal_Epi", "Plasmas", "SMCs", "T_cells"))

dat_plot <- dat_plot %>% dplyr::filter(indexb %in% c("B_cells", "Endothelials", "Mast_cells", "Malignant", "L_RFC2_Prolif", 
    "H_RFC2_Prolif", "Fibroblasts", "Myeloids", "NK_cells", "Normal_Epi", "Plasmas", "SMCs", "T_cells"))

limits.max <- max(dat_plot$total_number)

color.palette <- c("#4393C3", "#ffdbba", "#B2182B")

size.of.text <- c(18, 14, 14, 12)

indexa_levels <- c("B_cells", "Endothelials", "Mast_cells", "Malignant", "L_RFC2_Prolif", "H_RFC2_Prolif", "Fibroblasts", 
    "Myeloids", "NK_cells", "Normal_Epi", "Plasmas", "SMCs", "T_cells")

indexb_levels <- c("B_cells", "Endothelials", "Mast_cells", "Malignant", "L_RFC2_Prolif", "H_RFC2_Prolif", "Fibroblasts", 
    "Myeloids", "NK_cells", "Normal_Epi", "Plasmas", "SMCs", "T_cells")

dat_plot$indexb <- factor(dat_plot$indexb, levels = indexb_levels)

dat_plot$indexa <- factor(dat_plot$indexa, levels = indexa_levels)

p1 <- ggplot(dat_plot, aes(x = indexa, y = indexb, fill = total_number)) + geom_tile(color = "white") + geom_text(aes(label = total_number), 
    color = "white", size = 5) + scale_fill_gradientn(colours = color.palette, limits = c(0, limits.max), name = "Total Number") + 
    scale_x_discrete("Sender") + scale_y_discrete("Receivor") + theme_minimal() + theme(axis.title = element_text(size = size.of.text[1]), 
    axis.text.x.bottom = element_text(hjust = 1, angle = 45, size = size.of.text[2], color = "black"), axis.text.y.left = element_text(size = size.of.text[2], 
        color = "black"), legend.title = element_text(size = size.of.text[3]), legend.text = element_text(size = size.of.text[4]), 
    panel.grid = element_blank(), legend.position = "none", plot.margin = margin(r = 0))

summary_dat_plot <- dat_plot %>% dplyr::group_by(indexb) %>% dplyr::summarise(total_sum = sum(total_number))

summary_dat_plot$indexb <- factor(summary_dat_plot$indexb, levels = indexb_levels)

max_sum <- max(summary_dat_plot$total_sum)

p2 <- ggplot(summary_dat_plot, aes(x = total_sum, y = indexb)) + geom_col(fill = "#0072B5FF") + scale_x_continuous("") + 
    scale_y_discrete("") + geom_text(aes(label = total_sum), hjust = -0.1, size = 5) + coord_cartesian(xlim = c(0, 
    max_sum * 1.1)) + theme_minimal() + theme(axis.title.x = element_text(size = size.of.text[1]), axis.title.y = element_blank(), 
    axis.text.x = element_text(size = size.of.text[2], color = "black"), axis.text.y = element_blank(), panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(), panel.grid.major.x = element_line(linetype = "dotted"), panel.grid.minor.x = element_blank(), 
    plot.margin = margin(l = 1))

p1 + p2 + plot_layout(widths = c(3, 1))

invisible(NULL)

for_scenic <- sc[, sc@meta.data$CellType3 %in% c("H_RFC2_Prolif", "L_RFC2_Prolif", "Malignant", "Normal_Epi")]

for_scenic <- FindVariableFeatures(object = for_scenic, nfeatures = 3000)

for_scenic <- for_scenic[VariableFeatures(for_scenic), ]

for_scenic[["RNA"]] <- as(for_scenic[["RNA"]], Class = "Assay")

invisible(NULL)

for_scenic <- sc[, sc@meta.data$CellType3 %in% c("H_RFC2_Prolif", "L_RFC2_Prolif", "Malignant", "Normal_Epi")]

net <- decoupleR::get_collectri(split_complexes = T, organism = "human")

dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

regulon <- dorothea_regulon_human %>% dplyr::filter(confidence %in% c("A", "B", "C"))

res_doro <- dorothea::run_viper(input = for_scenic, regulons = regulon, options = list(method = "scale", minsize = 3, 
    eset.filter = FALSE, cores = 1, verbose = FALSE))

DefaultAssay(res_doro) <- "dorothea"

Idents(res_doro) <- "CellType3"

df.all <- FindAllMarkers(res_doro, only.pos = TRUE, logfc.threshold = 0.585)

df <- df.all %>% dplyr::group_by(cluster) %>% dplyr::mutate(diff = pct.1 - pct.2) %>% dplyr::arrange(desc(diff)) %>% 
    dplyr::filter(cluster == "H_RFC2_Prolif")

mycolor = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF", 
    "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", 
    "#B09C85FF", "#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF", "#008280FF", "#BB0021FF", "#5F559BFF", "#A20056FF", 
    "#808180FF", "#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "#ADB6B6FF", 
    "#374E55FF", "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF", "#1f77b4", "#ff7f0e", 
    "#279e68", "#d62728", "#aa40fc", "#8c564b", "#e377c2", "#b5bd61", "#17becf", "#aec7e8")

res_doro@meta.data$CellType3 <- factor(res_doro@meta.data$CellType3, levels = c("H_RFC2_Prolif", "L_RFC2_Prolif", 
    "Malignant", "Normal_Epi"))

mycol <- mycolor[1:length(unique(res_doro@meta.data$CellType3))]

names(mycol) <- unique(res_doro@meta.data$CellType3)

invisible(NULL)

plot_sc_markers_heatmap(object = res_doro, markerGene = df$gene, annoCol = TRUE, myanCol = mycol, cluster_columns = TRUE, 
    cluster_rows = TRUE)

invisible(NULL)

DefaultAssay(res_doro) <- "dorothea"

scop::FeatureStatPlot(srt = res_doro, stat.by = "HSF1", group.by = "CellType3", add_box = TRUE, legend.position = "none", 
    xlab = "", comparisons = list(c("H_RFC2_Prolif", "L_RFC2_Prolif"), c("H_RFC2_Prolif", "Malignant")), palette = "npg")

invisible(NULL)

scop::FeatureStatPlot(srt = res_doro, stat.by = "YY1", group.by = "CellType3", add_box = TRUE, legend.position = "none", 
    xlab = "", comparisons = list(c("H_RFC2_Prolif", "L_RFC2_Prolif"), c("H_RFC2_Prolif", "Malignant")), palette = "npg")

invisible(NULL)

scop::FeatureStatPlot(srt = res_doro, stat.by = "E2F2", group.by = "CellType3", add_box = TRUE, legend.position = "none", 
    xlab = "", comparisons = list(c("H_RFC2_Prolif", "L_RFC2_Prolif"), c("H_RFC2_Prolif", "Malignant")), palette = "npg")

invisible(NULL)

scop::FeatureStatPlot(srt = res_doro, stat.by = "MYCN", group.by = "CellType3", add_box = TRUE, legend.position = "none", 
    xlab = "", comparisons = list(c("H_RFC2_Prolif", "L_RFC2_Prolif"), c("H_RFC2_Prolif", "Malignant")), palette = "npg")

invisible(NULL)

scop::FeatureStatPlot(srt = res_doro, stat.by = "SMAD5", group.by = "CellType3", add_box = TRUE, legend.position = "none", 
    xlab = "", comparisons = list(c("H_RFC2_Prolif", "L_RFC2_Prolif"), c("H_RFC2_Prolif", "Malignant")), palette = "npg")

invisible(NULL)

DefaultAssay(res_doro) <- "RNA"

scop::FeatureStatPlot(srt = res_doro, stat.by = "HSF1", group.by = "CellType3", add_box = TRUE, legend.position = "none", 
    xlab = "", comparisons = list(c("H_RFC2_Prolif", "L_RFC2_Prolif"), c("H_RFC2_Prolif", "Malignant")), palette = "npg")

invisible(NULL)

scop::FeatureStatPlot(srt = res_doro, stat.by = "YY1", group.by = "CellType3", add_box = TRUE, legend.position = "none", 
    xlab = "", comparisons = list(c("H_RFC2_Prolif", "L_RFC2_Prolif"), c("H_RFC2_Prolif", "Malignant")), palette = "npg")

invisible(NULL)

scop::FeatureStatPlot(srt = res_doro, stat.by = "E2F2", group.by = "CellType3", add_box = TRUE, legend.position = "none", 
    xlab = "", comparisons = list(c("H_RFC2_Prolif", "L_RFC2_Prolif"), c("H_RFC2_Prolif", "Malignant")), palette = "npg")

invisible(NULL)

scop::FeatureStatPlot(srt = res_doro, stat.by = "MYCN", group.by = "CellType3", add_box = TRUE, legend.position = "none", 
    xlab = "", comparisons = list(c("H_RFC2_Prolif", "L_RFC2_Prolif"), c("H_RFC2_Prolif", "Malignant")), palette = "npg")

invisible(NULL)

scop::FeatureStatPlot(srt = res_doro, stat.by = "SMAD5", group.by = "CellType3", add_box = TRUE, legend.position = "none", 
    xlab = "", comparisons = list(c("H_RFC2_Prolif", "L_RFC2_Prolif"), c("H_RFC2_Prolif", "Malignant")), palette = "npg")

invisible(NULL)

invisible(NULL)

