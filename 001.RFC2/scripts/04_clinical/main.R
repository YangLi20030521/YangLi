# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
if (TRUE) {
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(tidyverse)
    library(circlize)
    library(TCGAplot)
    library(ggplot2)
    library(dplyr)
    library(IOBR)
    invisible(NULL)
    prad_restore("input_027")
    prad_restore("input_028")
    dat_use <- prad_input("input_029")
    gene <- prad_input("input_030")$candidata %>% as.character()
}

calculate_RS_score <- function(exp_data, coef_data) {
    genes_in_coef <- coef_data$gene
    coefficients <- coef_data$coefficient
    matching_genes <- intersect(genes_in_coef, colnames(exp_data))
    if (length(matching_genes) != length(genes_in_coef)) {
        missing_genes <- setdiff(genes_in_coef, matching_genes)
        cat("警告: 以下基因在表达数据中未找到:\n")
        cat(paste(missing_genes, collapse = ", "), "\n")
    }
    exp_subset <- exp_data[, matching_genes, drop = FALSE]
    coef_matching <- coef_data[match(colnames(exp_subset), coef_data$gene), ]
    coefficients_matching <- coef_matching$coefficient
    RS_scores <- as.matrix(exp_subset) %*% coefficients_matching
    RS_scores <- as.numeric(RS_scores)
    RS_scores_normalized <- (RS_scores - mean(RS_scores))/sd(RS_scores)
    result <- data.frame(sample = rownames(exp_data), RS_score_raw = RS_scores, RS_score_normalized = RS_scores_normalized)
    return(result)
}

add_significance <- function(p_values) {
    sig_levels <- c(`***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1)
    cut_points <- c(0, 0.001, 0.01, 0.05, 1)
    significance <- cut(p_values, breaks = cut_points, labels = c("****", "***", "**", "*"), right = FALSE)
    return(as.character(significance))
}

normfun <- function(data, ymin = -2, ymax = 2) {
    xmax = max(data)
    xmin = min(data)
    y = (ymax - ymin) * (data - xmin)/(xmax - xmin) + ymin
    return(y)
}

t_test_filter_genes <- function(data, p_threshold = 0.05) {
    if (!"Group" %in% colnames(data)) {
        stop("数据中必须包含名为'Group'的分组列")
    }
    gene_cols <- colnames(data)[colnames(data) != "Group"]
    results <- data.frame(Gene = character(), t_statistic = numeric(), p_value = numeric(), mean_group1 = numeric(), 
        mean_group2 = numeric(), stringsAsFactors = FALSE)
    for (gene in gene_cols) {
        expr_data <- data[[gene]]
        group_data <- data$Group
        group1_data <- expr_data[group_data == levels(as.factor(group_data))[1]]
        group2_data <- expr_data[group_data == levels(as.factor(group_data))[2]]
        t_test_result <- t.test(group1_data, group2_data, var.equal = TRUE)
        mean1 <- mean(group1_data, na.rm = TRUE)
        mean2 <- mean(group2_data, na.rm = TRUE)
        results <- rbind(results, data.frame(Gene = gene, t_statistic = t_test_result$statistic, p_value = t_test_result$p.value, 
            mean_group1 = mean1, mean_group2 = mean2, stringsAsFactors = FALSE))
    }
    significant_genes <- results[results$p_value < p_threshold, ]
    significant_genes <- significant_genes[order(significant_genes$p_value), ]
    return(significant_genes)
}

dat_exp <- total_expr_list$TCGA_PRAD %>% as.data.frame()

res_enrich <- do_bulk_gsva_ssgsea(dat_exp = dat_exp, geneset = gene, id = "RGDPs")

res_gsva <- res_enrich$res_gsva %>% as.data.frame()

dat_risk <- dat_use$dat_risk$TCGA_PRAD %>% as.data.frame() %>% column_to_rownames("ID")

com_sam <- intersect(rownames(dat_risk), rownames(res_gsva))

dat_risk <- dat_risk[com_sam, , drop = FALSE]

res_gsva <- res_gsva[com_sam, , drop = FALSE]

dat_risk <- cbind(dat_risk, res_gsva)

invisible(NULL)

plot_expr_cor_Scatter_1(data = dat_risk, var1 = "RS", var2 = "RGDPs", method = "pearson")

invisible(NULL)

dat_list <- prad_input("input_031")

dat_exp <- dat_list$dat_exp %>% as.data.frame()

dat_exp <- IOBR::count2tpm(countMat = dat_exp, idType = "Symbol", org = "hsa")

dat_cli <- dat_list$dat_cli %>% as.data.frame()

dat_cli$Group <- ifelse(dat_cli$Treatment == "post-ADT", "post", "pre")

dat_cli <- subset(dat_cli, dat_cli$Group == "post")

dat_cli$Group <- ifelse(dat_cli$Recurrence == 0, "ADT-noRecurrence", "ADT-Recurrence")

dat_exp <- dat_exp[, rownames(dat_cli)] %>% as.data.frame()

coef <- dat_use$dat_coef %>% tibble::remove_rownames() %>% as.data.frame()

dat_risk <- calculate_RS_score(exp_data = t(dat_exp), coef_data = coef)

print(identical(dat_risk$sample, rownames(dat_cli)))

dat_risk <- cbind(dat_risk, dat_cli) %>% dplyr::select(c(Group, RS_score_normalized)) %>% dplyr::rename(RS = "RS_score_normalized") %>% 
    tibble::remove_rownames()

t.test(dat_risk$RS ~ dat_risk$Group)

p = 1

stop("AUDIT: outcome-dependent subsampling disabled; repeat this comparison using a prespecified eligible cohort.")

plot_expr_diff_Violin(data = dat_sample, Group = "Group", variable = "RS", colors = c(`ADT-Recurrence` = "#D01910", 
    `ADT-noRecurrence` = "#00599F"))

invisible(NULL)

dat_tip <- prad_input("input_032")

dat_tip <- dat_tip$PRAD %>% t() %>% as.data.frame()

dat_risk <- dat_use$dat_risk$TCGA_PRAD %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID") %>% 
    dplyr::select(RS)

rownames(dat_risk) <- gsub("-", ".", rownames(dat_risk))

rownames(dat_risk) <- substr(rownames(dat_risk), 1, 15)

com_sam <- intersect(rownames(dat_risk), rownames(dat_tip))

dat_risk <- dat_risk[com_sam, , drop = FALSE]

dat_tip <- dat_tip[com_sam, ]

dat_risk <- cbind(dat_risk, dat_tip) %>% dplyr::mutate(Group = ifelse(RS > quantile(RS, 0.7), "High RS", ifelse(RS < 
    quantile(RS, 0.3), "Low RS", "Medium"))) %>% dplyr::select(-RS) %>% dplyr::filter(Group %in% c("High RS", "Low RS"))

dat_risk <- reshape2::melt(dat_risk)

colnames(dat_risk)

plot_expr_multidiff_Box(data = dat_risk, variable = "variable", Group = "Group", Expression = "value", color_use = c("#D01910", 
    "#00599F"))

invisible(NULL)

prad_restore("input_033")

dat_risk <- calculate_RS_score(exp_data = t(exp), coef_data = coef)

print(identical(dat_risk$sample, rownames(pdata)))

dat_risk <- dat_risk %>% dplyr::mutate(Response = pdata$binaryResponse) %>% na.omit() %>% tibble::remove_rownames() %>% 
    column_to_rownames("sample") %>% dplyr::select(RS_score_normalized, Response) %>% dplyr::rename(RS = "RS_score_normalized")

plot_expr_diff_Violin(data = dat_risk, Group = "Response", variable = "RS", colors = c(`CR/PR` = "#D01910", `SD/PD` = "#00599F"))

invisible(NULL)

prad_restore("input_034")

dat_risk <- calculate_RS_score(exp_data = t(exp), coef_data = coef)

print(identical(dat_risk$sample, rownames(pdata)))

dat_risk <- dat_risk %>% dplyr::mutate(Response = pdata$response) %>% na.omit() %>% tibble::remove_rownames() %>% 
    column_to_rownames("sample") %>% dplyr::select(RS_score_normalized, Response) %>% dplyr::rename(RS = "RS_score_normalized")

dat_risk <- dat_risk %>% dplyr::mutate(Response = case_when(Response == "CR" ~ "CR/PR", Response == "PR" ~ "CR/PR", 
    Response == "PD" ~ "SD/PD", Response == "SD" ~ "SD/PD", TRUE ~ Response))

plot_expr_diff_Violin(data = dat_risk, Group = "Response", variable = "RS", colors = c(`CR/PR` = "#D01910", `SD/PD` = "#00599F"))

p = 1

stop("AUDIT: outcome-dependent subsampling disabled; repeat this comparison using a prespecified eligible cohort.")

plot_expr_diff_Violin(data = dat_sample, Group = "Response", variable = "RS", colors = c(`CR/PR` = "#D01910", `SD/PD` = "#00599F"))

invisible(NULL)

prad_restore("input_035")

dat_risk <- calculate_RS_score(exp_data = t(exp), coef_data = coef)

print(identical(dat_risk$sample, rownames(pdata)))

dat_risk <- dat_risk %>% dplyr::mutate(Response = pdata$response) %>% na.omit() %>% tibble::remove_rownames() %>% 
    column_to_rownames("sample") %>% dplyr::select(RS_score_normalized, Response) %>% dplyr::rename(RS = "RS_score_normalized")

dat_risk <- dat_risk %>% dplyr::mutate(Response = case_when(Response == "CR" ~ "CR/PR", Response == "PR" ~ "CR/PR", 
    Response == "PD" ~ "SD/PD", Response == "SD" ~ "SD/PD", TRUE ~ Response))

plot_expr_diff_Violin(data = dat_risk, Group = "Response", variable = "RS", colors = c(`CR/PR` = "#D01910", `SD/PD` = "#00599F"))

p = 1

stop("AUDIT: outcome-dependent subsampling disabled; repeat this comparison using a prespecified eligible cohort.")

plot_expr_diff_Violin(data = dat_sample, Group = "Response", variable = "RS", colors = c(`CR/PR` = "#D01910", `SD/PD` = "#00599F"))

invisible(NULL)

dat_tmb <- TCGAplot::get_tmb()

dat_risk <- dat_use$dat_risk$TCGA_PRAD %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID") %>% 
    dplyr::select(RS)

com_sam <- intersect(rownames(dat_tmb), rownames(dat_risk))

dat_tmb <- dat_tmb[com_sam, , drop = FALSE]

dat_risk <- dat_risk[com_sam, , drop = FALSE]

dat_risk <- cbind(dat_risk, dat_tmb) %>% dplyr::mutate(Group = ifelse(TMB > median(TMB), "High TMB", "Low TMB"))

plot_expr_diff_Violin(data = dat_risk, Group = "Group", variable = "RS", colors = c(`High TMB` = "#D01910", `Low TMB` = "#00599F"))

invisible(NULL)

dat_immu <- cell_list$TCGA_PRAD %>% as.data.frame() %>% column_to_rownames("ID")

dat_risk <- dat_use$dat_risk$TCGA_PRAD %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID") %>% 
    dplyr::select(RS) %>% dplyr::mutate(Group = ifelse(RS > median(RS), "High_RS", "Low_RS"))

com_sam <- intersect(rownames(dat_immu), rownames(dat_risk))

dat_immu <- dat_immu[com_sam, ]

dat_risk <- dat_risk[com_sam, ]

res_test <- NULL

for (i in 1:length(dat_immu)) {
    cell = colnames(dat_immu)[i]
    x = dat_risk$Group
    y = dat_immu %>% dplyr::pull(cell) %>% as.numeric()
    res <- t.test(y ~ x)
    res_diff <- as.data.frame(res$estimate) %>% t() %>% as.data.frame() %>% tibble::remove_rownames() %>% dplyr::rename(High_RS_mean = "mean in group High_RS", 
        Low_RS_mean = "mean in group Low_RS")
    res_diff <- res_diff[, c("High_RS_mean", "Low_RS_mean")]
    res_test_i <- data.frame(cell = cell, pvalue = res$p.value)
    res_test_i <- cbind(res_test_i, res_diff)
    res_test <- rbind(res_test, res_test_i)
}

res_test <- res_test %>% dplyr::filter(pvalue < 0.05)

dat_immu <- dat_immu[, res_test$cell] %>% t() %>% as.data.frame()

identical(colnames(dat_immu), rownames(dat_risk))

res_test$sig <- add_significance(res_test$pvalue)

res_test$direction_sig <- ifelse(res_test$High_RS_mean > res_test$Low_RS_mean, paste0("High_RS ", res_test$sig), 
    paste0("Low_RS ", res_test$sig))

nc_colors <- list(low = "#00468BFF", mid = "white", high = "#ED0000FF")

mat <- as.matrix(dat_immu)

for (i in 1:ncol(mat)) {
    mat[, i] <- normfun(mat[, i], ymin = -4, ymax = 4)
}

feature_type <- sub(".*_([^_]+)$", "\\1", rownames(mat))

feature_type[feature_type == "estimate"] <- "ESTIMATE"

feature_type[feature_type == "ABS"] <- "CIBERSORT-ABS"

feature_type[feature_type == "quantiseq"] <- "Quantiseq"

feature_type_groups <- c("ESTIMATE", "CIBERSORT", "CIBERSORT-ABS", "TIMER", "Quantiseq", "MCPcounter", "EPIC", 
    "xCell")

feature_type_factor <- factor(feature_type, levels = feature_type_groups)

type_order <- order(feature_type_factor)

mat_sorted <- mat[type_order, ]

feature_type_sorted <- feature_type_factor[type_order]

new_rownames <- sub("^(.*)_([^_]+)$", "\\1", rownames(mat_sorted))

new_rownames <- sub("^(.*)_CIBERSORT$", "\\1", new_rownames)

rownames(mat_sorted) <- new_rownames

risk_groups <- factor(dat_risk$Group, levels = c("Low_RS", "High_RS"))

ha <- HeatmapAnnotation(Risk_Group = risk_groups, col = list(Risk_Group = c(Low_RS = "#1A5CA3", High_RS = "#EF3E2A")), 
    show_legend = FALSE, border = TRUE, annotation_name_side = "left")

row_ha_left <- rowAnnotation(Feature_Type = feature_type_sorted, col = list(Feature_Type = c(ESTIMATE = "#ADB6B6FF", 
    CIBERSORT = "#e08214", `CIBERSORT-ABS` = "#807dba", TIMER = "#b2182b", Quantiseq = "#2166ac", MCPcounter = "#4393c3", 
    EPIC = "#d6604d", xCell = "#925E9FFF")), show_annotation_name = FALSE, show_legend = FALSE, width = unit(0.4, 
    "cm"), border = TRUE)

max_val <- max(abs(mat_sorted), na.rm = TRUE)

col_fun <- colorRamp2(c(-max_val, 0, max_val), c(nc_colors$low, nc_colors$mid, nc_colors$high))

sig_annot <- res_test$direction_sig[type_order]

direction <- gsub("\\s\\*+$", "", sig_annot)

direction_colors <- ifelse(direction == "High_RS", "#EF3E2A", "#1A5CA3")

row_ha_right <- rowAnnotation(Significance = anno_text(sig_annot, gp = gpar(fontsize = 8, fontface = "bold"), width = max_text_width(sig_annot) + 
    unit(2, "mm"), show_name = FALSE), Direction = direction, col = list(Direction = c(High_RS = "#EF3E2A", Low_RS = "#1A5CA3")), 
    show_annotation_name = FALSE, show_legend = FALSE, border = TRUE)

ht <- Heatmap(mat_sorted, name = "z-score", col = col_fun, row_names_side = "left", row_names_gp = gpar(fontsize = 8, 
    fontface = "italic"), column_names_gp = gpar(fontsize = 6), show_column_names = FALSE, row_split = feature_type_sorted, 
    row_gap = unit(1.5, "mm"), row_title_rot = 0, row_title_gp = gpar(fontsize = 10, fontface = "bold"), column_split = risk_groups, 
    column_gap = unit(1.5, "mm"), column_title_gp = gpar(fontsize = 10, fontface = "bold"), cluster_columns = TRUE, 
    cluster_rows = FALSE, show_row_dend = FALSE, top_annotation = ha, left_annotation = row_ha_left, right_annotation = row_ha_right, 
    border = TRUE, rect_gp = gpar(col = "gray90", lwd = 0.5), show_heatmap_legend = FALSE, use_raster = TRUE, raster_quality = 2)

invisible(NULL)

ht_draw <- draw(ht, padding = unit(c(2, 2, 2, 2), "cm"))

n_row_slices <- length(unique(feature_type_sorted))

for (slice_index in 1:n_row_slices) {
    decorate_heatmap_body("z-score", {
        grid.rect(gp = gpar(fill = NA, col = "black", lwd = 2))
    }, slice = slice_index, column_slice = 1)
}

for (slice_index in 1:n_row_slices) {
    decorate_heatmap_body("z-score", {
        grid.rect(gp = gpar(fill = NA, col = "black", lwd = 2))
    }, slice = slice_index, column_slice = 2)
}

invisible(NULL)

ici <- prad_input("input_036")

ici_ac <- ici %>% dplyr::filter(Role %in% c("Active", "ACTIVE")) %>% dplyr::pull("Gene")

ici_in <- ici %>% dplyr::filter(Role == "Inhibition") %>% dplyr::pull("Gene")

dat_risk <- dat_use$dat_risk$TCGA_PRAD %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID") %>% 
    dplyr::select(RS)

dat_exp <- total_expr_list$TCGA_PRAD %>% as.data.frame()

com_gene <- intersect(rownames(dat_exp), ici_ac)

dat_plot <- dat_exp[com_gene, ] %>% t() %>% as.data.frame()

print(identical(rownames(dat_plot), rownames(dat_risk)))

dat_risk <- dat_risk %>% dplyr::mutate(Group = ifelse(RS > median(RS), "High RS", "Low RS")) %>% cbind(., dat_plot) %>% 
    dplyr::select(-RS)

sig_genes <- t_test_filter_genes(dat_risk, p_threshold = 0.001) %>% dplyr::pull("Gene") %>% as.character()

dat_risk <- dat_risk[, c("Group", all_of(sig_genes))]

dat_risk <- reshape2::melt(dat_risk)

colnames(dat_risk)

plot_expr_multidiff_Box(data = dat_risk, variable = "variable", Group = "Group", Expression = "value", color_use = c("#D01910", 
    "#00599F"))

invisible(NULL)

genelist <- IOBR::signature_collection

gene <- genelist$Immune_Checkpoint %>% as.character()

dat_exp <- total_expr_list$TCGA_PRAD %>% as.data.frame()

dat_risk <- dat_use$dat_risk$TCGA_PRAD %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID") %>% 
    dplyr::select(RS)

com_gene <- intersect(rownames(dat_exp), gene)

dat_plot <- dat_exp[com_gene, ] %>% t() %>% as.data.frame()

print(identical(rownames(dat_plot), rownames(dat_risk)))

dat_risk <- dat_risk %>% dplyr::mutate(Group = ifelse(RS > median(RS), "High RS", "Low RS")) %>% cbind(., dat_plot) %>% 
    dplyr::select(-RS)

dat_risk <- reshape2::melt(dat_risk)

colnames(dat_risk)

plot_expr_multidiff_Box(data = dat_risk, variable = "variable", Group = "Group", Expression = "value", color_use = c("#D01910", 
    "#00599F"))

invisible(NULL)

invisible(NULL)

dat_list <- prad_input("input_037")

dat_exp <- dat_list$dat_exp %>% as.data.frame()

dat_cli <- dat_list$dat_cli %>% as.data.frame()

dat_risk <- calculate_RS_score(exp_data = t(dat_exp), coef_data = coef)

print(identical(dat_risk$sample, rownames(dat_cli)))

dat_risk <- dat_risk %>% dplyr::mutate(Gleason = dat_cli$Gleason) %>% na.omit() %>% tibble::remove_rownames() %>% 
    column_to_rownames("sample") %>% dplyr::select("Gleason", "RS_score_normalized") %>% dplyr::rename(RS = "RS_score_normalized")

dat_summary <- dat_risk %>% dplyr::group_by("Gleason") %>% dplyr::summarise(., mean("RS"))

dat_summary <- dat_risk %>% dplyr::group_by(Gleason) %>% summarise(mean_RS = mean(RS, na.rm = TRUE), .groups = "drop")

