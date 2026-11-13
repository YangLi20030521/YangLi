# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
if (TRUE) {
    library(data.table)
    library(tidyverse)
    library(pheatmap)
    library(dplyr)
    invisible(NULL)
    gene = c("CCNB1", "CDKN3", "CLSPN", "FERMT3", "NUF2", "PTTG1", "RFC2", "RGS19", "WDR54")
    cancer = "PRAD"
}

dat1 <- prad_input("input_038") %>% as.data.frame()

rownames(dat1) <- dat1[, 1]

dat1 <- dat1[, -1]

colnames(dat1) <- gsub(" \\(.*\\)", "", colnames(dat1))

print(table(duplicated(colnames(dat1))))

dat2 <- prad_input("input_039") %>% as.data.frame()

rownames(dat2) <- dat2[, 1]

dat2 <- dat2[, -1]

colnames(dat2) <- gsub(" \\(.*\\)", "", colnames(dat2))

print(table(duplicated(colnames(dat2))))

model <- prad_input("input_040") %>% as.data.frame() %>% column_to_rownames("ModelID")

dat_exp <- prad_input("input_041") %>% as.data.frame() %>% dplyr::filter(IsDefaultEntryForModel == "Yes" & IsDefaultEntryForMC == 
    "Yes") %>% column_to_rownames("ModelID") %>% dplyr::select(-c(1:5))

colnames(dat_exp) <- gsub(" \\(.*\\)", "", colnames(dat_exp))

print(table(duplicated(colnames(dat_exp))))

com_gene <- intersect(colnames(dat1), gene)

dat_use1 <- dat1[, com_gene, drop = FALSE] %>% as.data.frame()

dat_use2 <- model %>% dplyr::filter(DepmapModelType == cancer)

com_cell <- intersect(rownames(dat_use1), rownames(dat_use2))

dat_use1 <- dat_use1[com_cell, , drop = FALSE]

dat_use2 <- dat_use2[com_cell, , drop = FALSE] %>% dplyr::select(StrippedCellLineName) %>% dplyr::rename(CellLine = "StrippedCellLineName")

print(identical(rownames(dat_use1), rownames(dat_use2)))

dat_use <- cbind(dat_use1, dat_use2) %>% tibble::remove_rownames()

mat <- t(as.matrix(dat_use[, 1:ncol(dat_use) - 1]))

colnames(mat) <- dat_use$CellLine

biorender_colors <- colorRampPalette(c("#F7F7F7", "#FEE0D2", "#FC9272", "#DE2D26", "#A50F15"))(100)

pheatmap(mat, color = biorender_colors, border_color = "black", cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, 
    show_colnames = TRUE, fontsize = 11, fontsize_row = 11, fontsize_col = 8, cellwidth = 25, cellheight = 25, 
    main = "", legend_breaks = seq(0, 1, 0.2), legend_labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"), legend = TRUE, 
    treeheight_row = 0, treeheight_col = 0)

com_gene <- intersect(colnames(dat2), gene)

dat_use1 <- dat2[, com_gene, drop = FALSE] %>% as.data.frame()

dat_use2 <- model %>% dplyr::filter(DepmapModelType == cancer)

com_cell <- intersect(rownames(dat_use1), rownames(dat_use2))

dat_use1 <- dat_use1[com_cell, , drop = FALSE]

dat_use2 <- dat_use2[com_cell, , drop = FALSE] %>% dplyr::select(StrippedCellLineName) %>% dplyr::rename(CellLine = "StrippedCellLineName")

print(identical(rownames(dat_use1), rownames(dat_use2)))

dat_use <- cbind(dat_use1, dat_use2) %>% tibble::remove_rownames()

mat <- t(as.matrix(dat_use[, 1:ncol(dat_use) - 1]))

colnames(mat) <- dat_use$CellLine

biorender_colors <- colorRampPalette(c("#2166AC", "#67A9CF", "#F7F7F7", "#F4A582", "#B2182B"))(100)

pheatmap(mat, color = biorender_colors, border_color = "black", cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, 
    show_colnames = TRUE, fontsize = 11, fontsize_row = 11, fontsize_col = 8, cellwidth = 25, cellheight = 25, 
    main = "", legend_breaks = seq(-2, 2, 0.4), legend_labels = c("-2.0", "-1.6", "-1.2", "-0.8", "-0.4", "0.0", 
        "0.4", "0.8", "1.2", "1.6", "2.0"), legend = TRUE, treeheight_row = 0, treeheight_col = 0)

for (i in 1:length(gene)) {
    gg <- gene[i]
    com_gene <- Reduce(intersect, list(colnames(dat_exp), colnames(dat2), gg))
    dat_use1 <- dat_exp[, com_gene, drop = FALSE]
    colnames(dat_use1) <- paste0(colnames(dat_use1), "_Expression")
    dat_use2 <- dat2[, com_gene, drop = FALSE]
    colnames(dat_use2) <- paste0(colnames(dat_use2), "_GeneEffect")
    dat_use3 <- model %>% dplyr::filter(DepmapModelType == cancer)
    com_cell <- Reduce(intersect, list(rownames(dat_use1), rownames(dat_use2), rownames(dat_use3)))
    dat_use1 <- dat_use1[com_cell, , drop = FALSE]
    dat_use2 <- dat_use2[com_cell, , drop = FALSE]
    dat_use3 <- dat_use3[com_cell, , drop = FALSE] %>% dplyr::select(StrippedCellLineName) %>% dplyr::rename(CellLine = "StrippedCellLineName")
    dat_use <- cbind(dat_use1, dat_use2, dat_use3) %>% tibble::remove_rownames()
    plot_expr_cor_Scatter_1(data = dat_use, var1 = paste0(gg, "_Expression"), var2 = paste0(gg, "_GeneEffect"))
}

for (i in 1:length(gene)) {
    gg <- gene[i]
    com_gene <- Reduce(intersect, list(colnames(dat_exp), colnames(dat1), gg))
    dat_use1 <- dat_exp[, com_gene, drop = FALSE]
    colnames(dat_use1) <- paste0(colnames(dat_use1), "_Expression")
    dat_use2 <- dat1[, com_gene, drop = FALSE]
    colnames(dat_use2) <- paste0(colnames(dat_use2), "_GeneDependency")
    dat_use3 <- model %>% dplyr::filter(DepmapModelType == cancer)
    com_cell <- Reduce(intersect, list(rownames(dat_use1), rownames(dat_use2), rownames(dat_use3)))
    dat_use1 <- dat_use1[com_cell, , drop = FALSE]
    dat_use2 <- dat_use2[com_cell, , drop = FALSE]
    dat_use3 <- dat_use3[com_cell, , drop = FALSE] %>% dplyr::select(StrippedCellLineName) %>% dplyr::rename(CellLine = "StrippedCellLineName")
    dat_use <- cbind(dat_use1, dat_use2, dat_use3) %>% tibble::remove_rownames()
    plot_expr_cor_Scatter_1(data = dat_use, var1 = paste0(gg, "_Expression"), var2 = paste0(gg, "_GeneDependency"))
}

dat_use <- prad_input("input_042")

dat_plot <- dat_use$dat_coef %>% as.data.frame() %>% tibble::remove_rownames()

dat_plot <- dat_plot %>% dplyr::arrange(coefficient) %>% dplyr::mutate(gene = factor(gene, levels = gene))

ggplot(dat_plot, aes(x = gene, y = coefficient, fill = coefficient > 0)) + geom_bar(stat = "identity", width = 0.7, 
    color = "black", size = 0.5) + geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.8) + 
    scale_fill_manual(values = c(`TRUE` = "#D73027", `FALSE` = "#4575B4"), guide = "none") + labs(x = "Gene", y = "Coefficient", 
    title = "Model Gene Coefficients") + theme_classic(base_size = 12) + theme(axis.line = element_line(color = "black", 
    size = 0.8), axis.ticks = element_line(color = "black", size = 0.6), axis.ticks.length = unit(0.15, "cm"), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 11, face = "italic"), 
    axis.text.y = element_text(color = "black", size = 11), axis.title = element_text(color = "black", size = 13, 
        face = "bold"), plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) + 
    coord_cartesian(clip = "off")

invisible(NULL)

invisible(NULL)

