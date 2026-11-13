# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
if (TRUE) {
    library(tidyverse)
    library(ggplot2)
    library(dplyr)
    gene = "RFC2"
    invisible(NULL)
    prad_restore("input_045")
}

dat_use <- prad_input("input_046")

dat_exp <- dat_use$dat_exp[, gene, drop = FALSE] %>% as.data.frame() %>% log1p()

dat_plot <- dat_exp %>% cbind(., dat_use$dat_cli)

t.test(dat_plot$RFC2 ~ dat_plot$Group)

dat_use <- prad_input("input_047")

dat_exp <- dat_use$dat_exp[gene, , drop = FALSE] %>% t() %>% as.data.frame() %>% na.omit()

dat_cli <- dat_use$dat_cli %>% as.data.frame()

com_sam <- intersect(rownames(dat_cli), rownames(dat_exp))

dat_exp <- dat_exp[com_sam, , drop = FALSE]

dat_cli <- dat_cli[com_sam, , drop = FALSE]

dat_plot <- cbind(dat_exp, dat_cli)

t.test(dat_plot$RFC2 ~ dat_plot$Group)

plot_expr_diff_Violin(data = dat_plot, Group = "Group", variable = gene, colors = c(Adjacent = "#00599F", Tumor = "#D01910"), 
    method = "t.test")

invisible(NULL)

dat_exp <- total_expr_list$TCGA_PRAD[gene, , drop = FALSE] %>% t() %>% as.data.frame()

dat_cli <- total_clin_list$TCGA_PRAD %>% as.data.frame() %>% dplyr::select(-"ID")

com_sam <- intersect(rownames(dat_cli), rownames(dat_exp))

dat_exp <- dat_exp[com_sam, , drop = FALSE]

dat_cli <- dat_cli[com_sam, , drop = FALSE] %>% dplyr::rename(RFS = "PFS", RFS.time = "PFS.time")

dat_plot <- cbind(dat_exp, dat_cli)

plot_survival_KM_1(data = dat_plot, minprop = 0.3, GeneName = gene, time = "RFS.time", event = "RFS", Timeunit = "y")

dat_plot <- cbind(dat_exp, dat_cli) %>% dplyr::select(c(gene, "OS", "OS.time"))

plot_survival_KM_1(data = dat_plot, minprop = 0.3, GeneName = gene, time = "OS.time", event = "OS", Timeunit = "y")

dat_plot <- cbind(dat_exp, dat_cli) %>% dplyr::select(c(gene, "DFS", "DFS.time"))

plot_survival_KM_1(data = dat_plot, minprop = 0.3, GeneName = gene, time = "DFS.time", event = "DFS", Timeunit = "y")

dat_exp <- total_expr_list$DKFZ2018[gene, , drop = FALSE] %>% t() %>% as.data.frame() %>% na.omit()

dat_cli <- total_clin_list$DKFZ2018 %>% as.data.frame() %>% dplyr::select(-"ID")

com_sam <- intersect(rownames(dat_cli), rownames(dat_exp))

dat_exp <- dat_exp[com_sam, , drop = FALSE]

dat_cli <- dat_cli[com_sam, , drop = FALSE]

dat_plot <- cbind(dat_exp, dat_cli)

plot_survival_KM_1(data = dat_plot, minprop = 0.3, GeneName = gene, time = "RFS.time", event = "RFS", Timeunit = "y")

dat_exp <- total_expr_list$GSE54460[gene, , drop = FALSE] %>% t() %>% as.data.frame() %>% na.omit()

dat_cli <- total_clin_list$GSE54460 %>% as.data.frame() %>% dplyr::select(-"ID")

com_sam <- intersect(rownames(dat_cli), rownames(dat_exp))

dat_exp <- dat_exp[com_sam, , drop = FALSE]

dat_cli <- dat_cli[com_sam, , drop = FALSE]

dat_plot <- cbind(dat_exp, dat_cli)

plot_survival_KM_1(data = dat_plot, minprop = 0.3, GeneName = gene, time = "RFS.time", event = "RFS", Timeunit = "y")

dat_exp <- total_expr_list$GSE21034[gene, , drop = FALSE] %>% t() %>% as.data.frame() %>% na.omit()

dat_cli <- total_clin_list$GSE21034 %>% as.data.frame() %>% dplyr::select(-"ID")

com_sam <- intersect(rownames(dat_cli), rownames(dat_exp))

dat_exp <- dat_exp[com_sam, , drop = FALSE]

dat_cli <- dat_cli[com_sam, , drop = FALSE]

dat_plot <- cbind(dat_exp, dat_cli)

plot_survival_KM_1(data = dat_plot, minprop = 0.3, GeneName = gene, time = "RFS.time", event = "RFS", Timeunit = "y")

dat_exp <- total_expr_list$GSE116918[gene, , drop = FALSE] %>% t() %>% as.data.frame() %>% na.omit()

dat_cli <- total_clin_list$GSE116918 %>% as.data.frame() %>% dplyr::select(-"ID")

com_sam <- intersect(rownames(dat_cli), rownames(dat_exp))

dat_exp <- dat_exp[com_sam, , drop = FALSE]

dat_cli <- dat_cli[com_sam, , drop = FALSE]

dat_plot <- cbind(dat_exp, dat_cli)

plot_survival_KM_1(data = dat_plot, minprop = 0.3, GeneName = gene, time = "RFS.time", event = "RFS", Timeunit = "y")

dat_exp <- total_expr_list$GSE107299[gene, , drop = FALSE] %>% t() %>% as.data.frame() %>% na.omit()

dat_cli <- total_clin_list$GSE107299 %>% as.data.frame() %>% dplyr::select(-"ID")

com_sam <- intersect(rownames(dat_cli), rownames(dat_exp))

dat_exp <- dat_exp[com_sam, , drop = FALSE]

dat_cli <- dat_cli[com_sam, , drop = FALSE]

dat_plot <- cbind(dat_exp, dat_cli)

plot_survival_KM_1(data = dat_plot, minprop = 0.3, GeneName = gene, time = "RFS.time", event = "RFS", Timeunit = "y")

dat_exp <- total_expr_list$GSE70768[gene, , drop = FALSE] %>% t() %>% as.data.frame() %>% na.omit()

dat_cli <- total_clin_list$GSE70768 %>% as.data.frame() %>% dplyr::select(-"ID")

com_sam <- intersect(rownames(dat_cli), rownames(dat_exp))

dat_exp <- dat_exp[com_sam, , drop = FALSE]

dat_cli <- dat_cli[com_sam, , drop = FALSE]

dat_plot <- cbind(dat_exp, dat_cli)

plot_survival_KM_1(data = dat_plot, minprop = 0.3, GeneName = gene, time = "RFS.time", event = "RFS", Timeunit = "y")

dat_exp <- total_expr_list$GSE70769[gene, , drop = FALSE] %>% t() %>% as.data.frame() %>% na.omit()

dat_cli <- total_clin_list$GSE70769 %>% as.data.frame() %>% dplyr::select(-"ID")

com_sam <- intersect(rownames(dat_cli), rownames(dat_exp))

dat_exp <- dat_exp[com_sam, , drop = FALSE]

dat_cli <- dat_cli[com_sam, , drop = FALSE]

dat_plot <- cbind(dat_exp, dat_cli)

plot_survival_KM_1(data = dat_plot, minprop = 0.3, GeneName = gene, time = "RFS.time", event = "RFS", Timeunit = "y")

dat_use <- prad_input("input_048")

dat_exp <- dat_use$dat_exp[gene, , drop = FALSE] %>% t() %>% as.data.frame() %>% na.omit()

dat_cli <- dat_use$dat_cli %>% as.data.frame()

com_sam <- intersect(rownames(dat_cli), rownames(dat_exp))

dat_exp <- dat_exp[com_sam, , drop = FALSE]

dat_cli <- dat_cli[com_sam, , drop = FALSE]

dat_plot <- cbind(dat_exp, dat_cli) %>% dplyr::select(c(gene, "RFS", "RFS.time")) %>% na.omit()

plot_survival_KM_1(data = dat_plot, minprop = 0.3, GeneName = gene, time = "RFS.time", event = "RFS", Timeunit = "y")

dat_exp <- total_expr_list$TCGA_PRAD[gene, , drop = FALSE] %>% t() %>% as.data.frame()

dat_cli <- total_clin_list$TCGA_PRAD %>% as.data.frame() %>% dplyr::select(-"ID")

com_sam <- intersect(rownames(dat_cli), rownames(dat_exp))

dat_exp <- dat_exp[com_sam, , drop = FALSE]

dat_cli <- dat_cli[com_sam, , drop = FALSE] %>% dplyr::rename(RFS = "PFS", RFS.time = "PFS.time")

dat_plot <- cbind(dat_exp, dat_cli) %>% dplyr::select(gene, "Gleason") %>% na.omit()

plot_expr_diff_Violin_plus(data = dat_plot, Group = "Gleason", variable = gene, comparisons = list(c("≤6", "7"), 
    c("7", "8-10"), c("≤6", "8-10")), method = "wilcox.test", colors = c(`≤6` = "#00599F", `7` = "#FFA500", 
    `8-10` = "#D01910"))

invisible(NULL)

dat_plot <- cbind(dat_exp, dat_cli) %>% dplyr::select(gene, "PSA") %>% na.omit()

plot_expr_diff_Violin(data = dat_plot, Group = "PSA", variable = gene, colors = c(Normal = "#00599F", Abnormal = "#D01910"), 
    method = "t.test")

invisible(NULL)

dat_plot <- cbind(dat_exp, dat_cli) %>% dplyr::select(gene, "N") %>% na.omit()

plot_expr_diff_Violin(data = dat_plot, Group = "N", variable = gene, colors = c(N0 = "#00599F", N1 = "#D01910"), 
    method = "t.test")

invisible(NULL)

dat_exp <- total_expr_list$TCGA_PRAD %>% as.data.frame()

dat_plot <- do_sc_bulk_Progeny(dat_exp = dat_exp, sc = NULL, dat_type = "Bulk", org = "Human")

dat_use <- total_expr_list$TCGA_PRAD[gene, , drop = FALSE] %>% t() %>% as.data.frame()

dat_plot <- cbind(dat_use, dat_plot)

pathway_data <- transform_pathway_correlation(dat_plot, gene_col = "RFC2")

plot_tcga_correlation(data = pathway_data, x_var = "gene_expression", y_var = "pathway_score", cancer_col = "pathway", 
    ncol = 7, point_size = 1.5, point_alpha = 0.6)

invisible(NULL)

dat_exp <- total_expr_list$TCGA_PRAD %>% as.data.frame()

dat_plot <- data.frame(row.names = colnames(dat_exp))

for (i in 1:length(cancer_sea_list)) {
    gene_i <- cancer_sea_list[[i]]
    name <- names(cancer_sea_list[i])
    dat_plot_i <- do_bulk_gsva_ssgsea(dat_exp = dat_exp, geneset = gene_i, id = name)
    dat_plot_i <- dat_plot_i$res_gsva %>% as.data.frame()
    dat_plot <- cbind(dat_plot, dat_plot_i)
}

dat_use <- total_expr_list$TCGA_PRAD[gene, , drop = FALSE] %>% t() %>% as.data.frame()

dat_plot <- cbind(dat_use, dat_plot)

pathway_data <- transform_pathway_correlation(dat_plot, gene_col = "RFC2")

plot_tcga_correlation(data = pathway_data, x_var = "gene_expression", y_var = "pathway_score", cancer_col = "pathway", 
    ncol = 7, point_size = 1.5, point_alpha = 0.6)

invisible(NULL)

