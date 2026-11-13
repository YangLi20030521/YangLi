# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
if (TRUE) {
    library(scRNAtoolVis)
    library(org.Hs.eg.db)
    library(BioEnricher)
    library(ClusterGVis)
    library(factoextra)
    library(tidyverse)
    library(openxlsx)
    library(ggthemes)
    library(ggplot2)
    library(GseaVis)
    library(DESeq2)
    library(ggpubr)
    library(Seurat)
    library(broom)
    library(ggsci)
    library(IOBR)
    library(scop)
    invisible(NULL)
    prad_restore("input_009")
}

sc <- prad_input("input_010")

table(sc@meta.data$CellType1)

sc@meta.data$CellType2 <- sc@meta.data$CellType1 %>% as.character()

sc@meta.data$CellType2[sc@meta.data$seurat_clusters == 23] <- "Prolif"

sc@meta.data$CellType2[sc@meta.data$CellType2 == "Cycling_T"] <- "Prolif"

sc@meta.data$CellType2[sc@meta.data$CellType2 == "TNKs"] <- "CTLs"

table(sc@meta.data$CellType2)

sc@meta.data$Group1 <- sc@meta.data$Grade %>% as.character()

sc@meta.data$Group1[sc@meta.data$Group1 == "high"] <- "High_Grade"

sc@meta.data$Group1[sc@meta.data$Group1 == "low"] <- "Low_Grade"

sc@meta.data$ISUP <- case_when(sc@meta.data$Gleason %in% c("3+3=6", "3+3=6, tertiary 4") ~ "Group1", sc@meta.data$Gleason == 
    "3+4=7" ~ "Group2", sc@meta.data$Gleason == "4+3=7" ~ "Group3", sc@meta.data$Gleason %in% c("3+5=8", "4+4", 
    "4+4=8", "5+3=8") ~ "Group4", sc@meta.data$Gleason %in% c("4+5=9", "5+5=10") ~ "Group5")

plot_sc_UMAP_cell(scRNA = sc, plotname = "CellType2", reduction = "umap", label = "CellType2", alphas = 1)

invisible(NULL)

invisible(NULL)

do_sc_ROE(sc = sc, col_celltype = "CellType2", col_group = "ISUP", col_sample = "sample", levels = c("Group1", 
    "Group2", "Group3", "Group4", "Group5"))

invisible(NULL)

invisible(NULL)

do_sc_ROE(sc = sc, col_celltype = "CellType2", col_group = "Group2", col_sample = "sample", levels = c("Tumor", 
    "Normal"))

invisible(NULL)

invisible(NULL)

do_sc_ROE(sc = sc, col_celltype = "CellType2", col_group = "Group1", col_sample = "sample", levels = c("High_Grade", 
    "Low_Grade"))

invisible(NULL)

sc[["RNA"]] <- as(sc[["RNA"]], Class = "Assay")

Idents(sc) <- "CellType2"

df.all <- FindAllMarkers(object = sc, group.by = "CellType2", only.pos = TRUE, logfc.threshold = 0.585)

df.all <- df.all %>% dplyr::mutate(diff = pct.1 - pct.2)

gene <- df.all %>% dplyr::group_by(cluster) %>% dplyr::arrange(desc(diff)) %>% dplyr::slice_head(n = 5) %>% dplyr::pull(gene) %>% 
    unique()

mycol <- c(pal_npg()(10), (ggsci::pal_lancet())(10), (ggsci::pal_aaas())(10))

mycol <- mycol[1:length(unique(sc@meta.data$CellType2))]

names(mycol) <- unique(sc@meta.data$CellType2)

invisible(NULL)

averageHeatmap(object = sc, markerGene = gene, annoCol = TRUE, myanCol = mycol)

invisible(NULL)

gene <- df.all %>% dplyr::group_by(cluster) %>% dplyr::arrange(desc(diff)) %>% dplyr::filter(avg_log2FC > 0.585 & 
    p_val_adj < 0.001) %>% dplyr::slice_head(n = 5)

dat_use <- prepareDataFromscRNA(object = sc, diffData = gene, showAverage = TRUE)

dat_enrich <- enrichCluster(object = dat_use, OrgDb = org.Hs.eg.db, type = "KEGG", organism = "hsa", pvalueCutoff = 0.05, 
    topn = 3, seed = 5201314)

markGenes <- unique(gene$gene)

invisible(NULL)

visCluster(object = dat_use, plot.type = "both", column_names_rot = 45, show_row_dend = FALSE, markGenes = markGenes, 
    markGenes.side = "left", annoTerm.data = dat_enrich, line.side = "left", cluster.order = c(1:16), kegg.col = rep(jjAnno::useMyCol("stallion", 
        n = 16), each = 5), add.bar = TRUE)

invisible(NULL)

df <- df.all %>% dplyr::mutate(diff = pct.1 - pct.2) %>% dplyr::filter(cluster == "Prolif") %>% dplyr::arrange(desc(diff))

gene1 <- df %>% dplyr::filter(diff > 0.2) %>% dplyr::pull(gene) %>% as.character()

dat_plot <- df %>% dplyr::select(c(avg_log2FC, diff)) %>% as.data.frame()

ggplot(dat_plot, aes(x = diff, y = avg_log2FC)) + geom_point(color = "#4A90E2", size = 2.5, alpha = 0.6, shape = 16) + 
    theme_minimal() + theme(plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white", 
    color = NA), panel.grid.major = element_line(color = "white", size = 0.3), panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "#333333", size = 0.8), axis.ticks = element_line(color = "#333333", size = 0.5), 
    axis.ticks.length = unit(0.15, "cm"), axis.title = element_text(size = 13, color = "#333333", face = "bold", 
        family = "sans"), axis.text = element_text(size = 11, color = "#555555", family = "sans"), plot.title = element_text(size = 15, 
        face = "bold", color = "#333333", hjust = 0.5, family = "sans"), plot.margin = margin(20, 20, 20, 20)) + 
    labs(x = "Difference (pct.1 - pct.2)", y = "Average Log2 Fold Change (avg_log2FC)", title = "") + scale_x_continuous(expand = expansion(mult = c(0.05, 
    0.05))) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))

invisible(NULL)

res_enrich <- lzq_ORA.integrated_1(genes = gene1)

dat_plot <- res_enrich$GO@result %>% dplyr::slice_head(n = 15)

plot_ORA_Bar(dat = dat_plot, color = "#00599F")

invisible(NULL)

table(sc@meta.data$Group2)

plot_sc_UMAP_cell(scRNA = sc, plotname = "Group2", reduction = "umap", label = "Group2", alphas = 1)

invisible(NULL)

invisible(NULL)

dat_use <- prad_input("input_011")

dat_exp <- dat_use$dat_exp %>% as.data.frame()

dat_cli <- dat_use$dat_cli %>% as.data.frame()

dat_cli$Group <- ifelse(dat_cli$Treatment == "post-ADT", "post", "pre")

dat_cli <- subset(dat_cli, dat_cli$Group == "post")

dat_cli$Group <- ifelse(dat_cli$Recurrence == 0, "noRecurrence", "Recurrence")

dat_cli <- dat_cli %>% dplyr::select(Group)

dat_exp <- dat_exp[, rownames(dat_cli)] %>% as.data.frame()

dds <- DESeqDataSetFromMatrix(countData = as.matrix(dat_exp), colData = dat_cli, design = ~Group)

dds <- DESeq(dds)

dds$Group <- factor(dds$Group, levels = c("noRecurrence", "Recurrence"))

cat("use:", resultsNames(dds)[2])

resLFC <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm")

Diff <- na.omit(as.data.frame(resLFC))

Diff <- Diff %>% dplyr::mutate(status = ifelse(pvalue < 0.05 & log2FoldChange > 0.585, "Up", ifelse(pvalue < 0.05 & 
    log2FoldChange < -0.585, "Down", "Not")))

Diff$SYMBOL <- rownames(Diff)

grlist <- Diff$log2FoldChange

names(grlist) <- Diff$SYMBOL

grlist <- sort(grlist, decreasing = TRUE)

table(Diff$status)

my_label <- paste0("FC:1.5 ; Pvalue:0.05 ; ", "Up:", table(Diff$status)[3], " ; ", "Down:", table(Diff$status)[1])

Diff$logPvalue <- -log10(Diff$pvalue)

ggscatter(Diff, x = "log2FoldChange", y = "logPvalue", color = "status", size = 2, main = paste0("ADT Recurrence -vs- ADT noRecurrence"), 
    xlab = "log2FoldChange", ylab = "-log10(P.value)", palette = c("#00599F", "#CCCCCC", "#D01910"), ylim = c(-0, 
        2), xlim = c(-8, 8)) + theme_base() + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#222222") + 
    geom_vline(xintercept = 0.58, linetype = "dashed", color = "#222222") + geom_vline(xintercept = -0.58, linetype = "dashed", 
    color = "#222222") + labs(subtitle = my_label) + theme(plot.background = element_blank())

invisible(NULL)

enrich_up <- lzq_ORA.integrated_1(genes = rownames(subset(Diff, Diff$status == "Up")), organism = "Human")

enrich_down <- lzq_ORA.integrated_1(genes = rownames(subset(Diff, Diff$status == "Down")), organism = "Human", 
    perform.Cancer.Gene.Network = FALSE)

enrich_gsea <- lzq_GSEA.integrated_1(genes = grlist, organism = "Human", pvalue.cutoff = 1)

gseaNb(object = enrich_gsea$WikiPathways, geneSetID = "WP2023", newGsea = T, addPval = T, pvalX = 0.75, pvalY = 0.8, 
    pCol = "black", pHjust = 0)

invisible(NULL)

gseaNb(object = enrich_gsea$GO, geneSetID = "GO:0043950", newGsea = T, addPval = T, pvalX = 0.75, pvalY = 0.8, 
    pCol = "black", pHjust = 0)

invisible(NULL)

dat_plot <- enrich_up$GO@result %>% slice_head(n = 20)

plot_ORA_Bar(dat = dat_plot, color = "#D01910")

invisible(NULL)

dat_plot <- enrich_up$KEGG@result %>% slice_head(n = 20)

plot_ORA_Bar(dat = dat_plot, color = "#D01910")

invisible(NULL)

gene2 <- subset(Diff, Diff$status == "Up") %>% dplyr::pull(SYMBOL) %>% as.character()

sc <- do_sc_AUCell(scRNA = sc, markers = gene2, cores_use = 25, name_use = "ADT Recurrence")

scop::FeatureStatPlot(srt = sc, stat.by = "ADT Recurrence", group.by = "CellType2", add_box = TRUE, legend.position = "none", 
    xlab = "")

invisible(NULL)

candidata <- intersect(gene1, gene2)

dat_plot <- list(Prolif = gene1, Recurrence = gene2)

plot_venn_1(data.list = dat_plot)

invisible(NULL)

dat_use <- prad_input("input_012")

dat_exp <- dat_use$dat_exp %>% as.data.frame()

dat_cli <- dat_use$dat_cli %>% as.data.frame()

dat_exp <- as.data.frame(t(dat_exp))

dat_exp <- IOBR::count2tpm(countMat = dat_exp, idType = "Symbol", org = "hsa")

res_enrich <- do_bulk_gsva_ssgsea(dat_exp = dat_exp, geneset = candidata, id = "RGDPs")

res_gsva <- res_enrich$res_gsva %>% as.data.frame() %>% cbind(., dat_cli)

plot_expr_diff_Violin(data = res_gsva, Group = "Group", variable = "RGDPs", colors = c(Tumor = "#D01910", Normal = "#00599F"), 
    method = "t.test")

invisible(NULL)

dat_use <- prad_input("input_013")

dat_exp <- dat_use$dat_exp

rownames(dat_exp) <- gsub(" ", "", rownames(dat_exp))

dat_cli <- dat_use$dat_cli

res_enrich <- do_bulk_gsva_ssgsea(dat_exp = dat_exp, geneset = candidata, id = "RGDPs")

res_gsva <- res_enrich$res_gsva %>% as.data.frame() %>% cbind(., dat_cli) %>% dplyr::filter(Group %in% c("Normal", 
    "Stroma", "Low Tumor", "Tumor"))

plot_expr_diff_Ridge(data = res_gsva, Group = "Group", variable = "RGDPs", levels = c("Normal", "Stroma", "Low Tumor", 
    "Tumor"), test.methods = "kruskal.test")

invisible(NULL)

dat_exp <- total_expr_list$TCGA_PRAD %>% as.data.frame()

dat_cli <- total_clin_list$TCGA_PRAD %>% dplyr::select(c("PFS", "PFS.time")) %>% na.omit()

com_sam <- intersect(rownames(dat_cli), colnames(dat_exp))

dat_cli <- dat_cli[com_sam, ]

dat_exp <- dat_exp[, com_sam]

res_enrich <- do_bulk_gsva_ssgsea(dat_exp = dat_exp, geneset = candidata, id = "RGDPs")

res_gsva <- res_enrich$res_gsva %>% as.data.frame() %>% cbind(., dat_cli)

plot_survival_KM_1(data = res_gsva, minprop = 0.3, GeneName = "RGDPs", time = "PFS.time", event = "PFS", Timeunit = "y")

dat_exp <- total_expr_list$GSE54460 %>% as.data.frame()

dat_cli <- total_clin_list$GSE54460 %>% as.data.frame()

com_sam <- intersect(rownames(dat_cli), colnames(dat_exp))

dat_cli <- dat_cli[com_sam, ]

dat_exp <- dat_exp[, com_sam]

res_enrich <- do_bulk_gsva_ssgsea(dat_exp = dat_exp, geneset = candidata, id = "RGDPs")

res_gsva <- res_enrich$res_gsva %>% as.data.frame() %>% cbind(., dat_cli)

plot_survival_KM_1(data = res_gsva, minprop = 0.3, GeneName = "RGDPs", time = "RFS.time", event = "RFS", Timeunit = "y")

dat_exp <- total_expr_list$TCGA_PRAD %>% as.data.frame()

dat_cli <- total_clin_list$TCGA_PRAD %>% dplyr::select(c("Age", "Race", "T", "N", "Gleason", "PSA", "Radiotherapy"))

com_sam <- intersect(rownames(dat_cli), colnames(dat_exp))

dat_cli <- dat_cli[com_sam, ]

dat_exp <- dat_exp[, com_sam]

res_enrich <- do_bulk_gsva_ssgsea(dat_exp = dat_exp, geneset = candidata, id = "RGDPs")

res_gsva <- res_enrich$res_gsva %>% as.data.frame() %>% cbind(., dat_cli)

res_cli <- clinical_heatmap_plot(data = res_gsva, score_var = "RGDPs", output_file = "./output/plot_73f94288.pdf")

dat_signature_hallmark <- IOBR::calculate_sig_score(pdata = NULL, eset = as.matrix(dat_exp), signature = hallmark, 
    method = "zscore", mini_gene_count = 2)

dat_signature_hallmark <- dat_signature_hallmark %>% column_to_rownames("ID") %>% as.data.frame() %>% dplyr::mutate(RGDPs = res_enrich$res_gsva$RGDPs)

res_cor <- do_cor_test_col(data = dat_signature_hallmark, variable = "RGDPs", methods = "pearson")

plot_expr_cor_Dot_plus(dat = res_cor, col_cor = "cor", col_pvalue = "pvalue")

invisible(NULL)

dat_signature_kegg <- IOBR::calculate_sig_score(pdata = NULL, eset = as.matrix(dat_exp), signature = kegg, method = "zscore", 
    mini_gene_count = 2)

dat_signature_kegg <- dat_signature_kegg %>% column_to_rownames("ID") %>% as.data.frame() %>% dplyr::mutate(RGDPs = res_enrich$res_gsva$RGDPs)

res_cor <- do_cor_test_col(data = dat_signature_kegg, variable = "RGDPs", methods = "pearson")

res_cor <- res_cor %>% dplyr::arrange(desc(abs(cor))) %>% dplyr::slice_head(n = 50)

plot_expr_cor_Dot_plus(dat = res_cor, col_cor = "cor", col_pvalue = "pvalue")

invisible(NULL)

dat_tip <- prad_input("input_014")

dat_tip <- dat_tip$PRAD %>% t() %>% as.data.frame()

dat_use <- res_enrich$res_gsva

rownames(dat_use) <- gsub("-", ".", rownames(dat_use))

rownames(dat_use) <- substr(rownames(dat_use), 1, 15)

com_sam <- intersect(rownames(dat_use), rownames(dat_tip))

dat_use <- dat_use[com_sam, , drop = FALSE]

dat_tip <- dat_tip[com_sam, ]

dat_use <- cbind(dat_use, dat_tip) %>% dplyr::mutate(Group = ifelse(RGDPs > median(RGDPs), "High RGDPs", "Low RGDPs")) %>% 
    dplyr::select(-RGDPs)

dat_use <- reshape2::melt(dat_use)

colnames(dat_use)

plot_expr_multidiff_Box(data = dat_use, variable = "variable", Group = "Group", Expression = "value", color_use = c("#D01910", 
    "#00599F"))

invisible(NULL)

dat_use <- prad_input("input_015")

dat_exp <- dat_use$dat_exp %>% as.data.frame()

dat_exp <- IOBR::count2tpm(countMat = dat_exp, idType = "Symbol", org = "hsa")

dat_cli <- dat_use$dat_cli %>% as.data.frame()

dat_cli$Group <- ifelse(dat_cli$Treatment == "post-ADT", "post", "pre")

dat_cli <- subset(dat_cli, dat_cli$Group == "post")

dat_cli$Group <- ifelse(dat_cli$Recurrence == 0, "ADT-noRecurrence", "ADT-Recurrence")

dat_exp <- dat_exp[, rownames(dat_cli)] %>% as.data.frame()

res_enrich <- do_bulk_gsva_ssgsea(dat_exp = dat_exp, geneset = candidata, id = "RGDPs")

res_gsva <- res_enrich$res_gsva %>% as.data.frame() %>% cbind(., dat_cli)

plot_expr_diff_Violin(data = res_gsva, Group = "Group", variable = "RGDPs", colors = c(`ADT-Recurrence` = "#D01910", 
    `ADT-noRecurrence` = "#00599F"), method = "t.test")

invisible(NULL)

CRPC <- paste("GSM", 1817707:1817719, sep = "")

Normal <- paste0("GSM", 1817833:1817905, seq = "")

dat_exp <- total_expr_list$GSE70768 %>% as.data.frame()

dat_cli <- total_clin_list$GSE70768 %>% as.data.frame() %>% dplyr::mutate(Group = case_when(ID %in% CRPC ~ "CRPC", 
    ID %in% Normal ~ "Normal", TRUE ~ "Tumor"))

res_enrich <- do_bulk_gsva_ssgsea(dat_exp = dat_exp, geneset = candidata, id = "RGDPs")

res_gsva <- res_enrich$res_gsva %>% as.data.frame() %>% cbind(., dat_cli)

plot_expr_diff_Ridge(data = res_gsva, Group = "Group", variable = "RGDPs", levels = c("Normal", "CRPC", "Tumor"), 
    test.methods = "kruskal.test")

invisible(NULL)

invisible(NULL)

