# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
pancan_exp_diff <- function(gene = "") {
    pan_TCGA_fpkm <- prad_input("input_044")
    res <- list()
    for (i in 1:length(pan_TCGA_fpkm)) {
        cancer_name <- names(pan_TCGA_fpkm)[i]
        exp <- pan_TCGA_fpkm[[i]] %>% as.data.frame()
        if (gene %in% colnames(exp)) {
            exp <- exp %>% dplyr::select(c("sample_type", gene)) %>% dplyr::mutate(Cancer = cancer_name) %>% dplyr::rename(Group = "sample_type")
            exp$Group[exp$Group == "TCGA_tumor"] <- "Tumor"
            exp$Group[exp$Group == "TCGA_normal"] <- "Normal"
        }
        else {
            print(paste0(gene, "不存在于", cancer_name))
        }
        if (length(table(exp$Group)) == 2) {
            exp_list <- list(a = exp)
            names(exp_list) <- cancer_name
            res <- append(res, exp_list)
        }
        else {
            print(paste0(cancer_name, "没有正常组织信息"))
        }
    }
    combine_res <- do.call(rbind, lapply(names(res), function(id) {
        df <- res[[id]]
        cbind(id = id, df)
    }))
    plot_expr_12 = function(data, tumor_colnames, group_use, var_use, freq_use = 3, color_use = c("#00599F", "#D01910")) {
        library(tidyverse)
        library(ggplot2)
        data$expr = data[, var_use]
        data$Type = data[, group_use]
        data$CancerType = data[, tumor_colnames]
        df_delete = data.frame(table(data$CancerType, data$Type))
        df_delete = df_delete[which(df_delete$Freq < freq_use), ]
        data = data[which(!data[, "CancerType"] %in% unique(df_delete$Var1)), ]
        pvalues = sapply(unique(data$CancerType), function(x) {
            res = wilcox.test(as.numeric(expr) ~ Type, data = subset(data, CancerType == x))
            res$p.value
        })
        pv = data.frame(cancer_type = unique(data$CancerType), pvalue = pvalues)
        pv$sigcode = cut(pv$pvalue, c(0, 1e-04, 0.001, 0.01, 0.05, 1), labels = c("****", "***", "**", "*", ""))
        count_N = data %>% group_by(CancerType, Type) %>% tally
        count_N$n = paste("n =", count_N$n)
        ylabname = var_use
        boxplot_res = ggplot(data, aes(x = CancerType, y = expr, color = Type, fill = Type)) + geom_point(shape = 21, 
            size = 2, position = position_jitterdodge(), alpha = 0.35) + geom_boxplot(notch = F, outlier.size = -1, 
            color = "black", lwd = 0.8, alpha = 0.7) + theme_classic() + scale_color_manual(values = color_use) + 
            scale_fill_manual(values = color_use) + theme(panel.border = element_rect(colour = "black", fill = NA, 
            size = 0.8)) + theme(axis.text.x = element_text(colour = "black", size = 11, angle = 45, hjust = 0.5, 
            vjust = 0.5), legend.position = "bottom") + geom_text(aes(x = cancer_type, y = max(data[, var_use]), 
            label = pv$sigcode), data = pv, inherit.aes = F) + ylab(ylabname) + geom_text(data = count_N, aes(label = n, 
            y = min(data[, var_use]) - 2, color = Type), position = position_dodge2(0.9), size = 3, angle = 90, 
            hjust = 0) + ylim(min(data[, var_use]) - 2, max(data[, var_use])) + theme(axis.text.x = element_text(angle = 45, 
            hjust = 1, vjust = 1), legend.position = "bottom", panel.border = element_rect(colour = "black", fill = NA, 
            size = 0.6))
        res_diff_pancancer = list(original_data = data, pvalue_res = pv, boxplot_res = boxplot_res)
    }
    combine_res$Group <- factor(combine_res$Group, levels = c("Normal", "Tumor"))
    p <- plot_expr_12(data = combine_res, tumor_colnames = "Cancer", group_use = "Group", var_use = gene, freq_use = 0)
    p
}

plot_expr_cor_Scatter_3 <- function(data, x_var = "", y_var = "", cancer_col = "", ncol = 7, point_size = 1.5, 
    point_alpha = 0.3) {
    library(tidyverse)
    library(cowplot)
    cancer_list <- unique(data[[cancer_col]])
    colors <- c("#4A90E2", "#E8743B", "#19A979", "#9B59B6", "#E74C3C", "#3498DB", "#1ABC9C", "#F39C12", "#E91E63", 
        "#9C27B0", "#00BCD4", "#4CAF50", "#FF9800", "#795548", "#607D8B", "#FF5722", "#673AB7", "#03A9F4", "#8BC34A", 
        "#FFC107", "#00897B", "#D32F2F", "#7B1FA2", "#0288D1", "#689F38", "#F57C00", "#5D4037", "#455A64", "#C62828", 
        "#6A1B9A", "#0277BD", "#558B2F")
    plot_list <- lapply(cancer_list, function(cancer) {
        cancer_data <- data %>% dplyr::filter(!!sym(cancer_col) == cancer)
        cor_val <- unique(cancer_data$correlation)
        p_val <- unique(cancer_data$p_value)
        if (p_val < 2.2e-16) {
            p_label <- "p < 2.2e-16"
        }
        else if (p_val < 0.001) {
            p_label <- sprintf("p = %.2e", p_val)
        }
        else {
            p_label <- sprintf("p = %.4f", p_val)
        }
        stat_label <- sprintf("R = %.2f, %s", cor_val, p_label)
        color_idx <- which(cancer_list == cancer)
        point_color <- colors[color_idx]
        p_main <- ggplot(cancer_data, aes(x = !!sym(x_var), y = !!sym(y_var))) + geom_point(color = point_color, 
            size = point_size, alpha = point_alpha) + geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", 
            linewidth = 1.2, linetype = "solid") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
            scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) + labs(x = NULL, y = NULL) + theme_classic(base_size = 10) + 
            theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6), axis.line = element_blank(), 
                axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = ggplot2::margin(0, 1, 
                  1, 1), panel.background = element_rect(fill = "white", color = NA))
        title_box <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) + annotate("rect", xmin = 0, xmax = 1, ymin = 0, 
            ymax = 1, fill = "#E8E8E8", color = NA) + annotate("text", x = 0.5, y = 0.5, label = paste0(cancer, 
            "\n", stat_label), size = 2.8, hjust = 0.5, vjust = 0.5, lineheight = 0.8) + scale_x_continuous(limits = c(0, 
            1), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + theme_void() + theme(panel.border = element_rect(fill = NA, 
            color = "black", linewidth = 0.8), plot.margin = ggplot2::margin(1, 1, 0, 1), panel.background = element_rect(fill = "#E8E8E8", 
            color = NA))
        combined <- plot_grid(title_box, p_main, ncol = 1, rel_heights = c(0.2, 1), align = "v", axis = "lr")
        return(combined)
    })
    combined_plot <- plot_grid(plotlist = plot_list, ncol = ncol, align = "hv", axis = "tblr")
    invisible(NULL)
}

transform_pathway_correlation <- function(data, gene_col = NULL, method = "pearson") {
    if (is.null(gene_col)) {
        gene_col <- names(data)[1]
    }
    if (!gene_col %in% names(data)) {
        stop(paste("基因列", gene_col, "不存在于数据中"))
    }
    pathway_cols <- setdiff(names(data), gene_col)
    if (length(pathway_cols) == 0) {
        stop("数据中没有通路列")
    }
    long_data <- data %>% tidyr::pivot_longer(cols = all_of(pathway_cols), names_to = "pathway", values_to = "pathway_score") %>% 
        dplyr::rename(gene_expression = !!sym(gene_col))
    correlation_results <- long_data %>% dplyr::group_by(pathway) %>% dplyr::summarise(correlation = cor(gene_expression, 
        pathway_score, method = method, use = "complete.obs"), p_value = cor.test(gene_expression, pathway_score, 
        method = method)$p.value, .groups = "drop")
    result_data <- long_data %>% dplyr::left_join(correlation_results, by = "pathway")
    result_data <- result_data %>% dplyr::select(gene_expression, pathway_score, pathway, correlation, p_value)
    return(result_data)
}

plot_expr_diff_Violin_plus <- function(data, Group = "Group", variable = NULL, colors = NULL, comparisons = NULL, 
    method = "wilcox.test", ncol = NULL, show_jitter = TRUE, jitter_size = 2.5, jitter_alpha = 0.8) {
    library(patchwork)
    library(ggplot2)
    library(ggpubr)
    library(dplyr)
    library(tidyr)
    if (!is.data.frame(data)) 
        stop("data必须为数据框")
    if (!Group %in% names(data)) 
        stop(paste(Group, "不存在于数据框中"))
    if (is.null(variable)) 
        stop("必须指定variable参数")
    missing_vars <- variable[!variable %in% names(data)]
    if (length(missing_vars) > 0) {
        stop(paste("以下变量不存在于数据框中:", paste(missing_vars, collapse = ", ")))
    }
    data[[Group]] <- factor(data[[Group]])
    group_levels <- levels(data[[Group]])
    n_groups <- length(group_levels)
    if (is.null(colors)) {
        if (n_groups == 2) {
            colors <- c("#00599F", "#D01910")
        }
        else if (n_groups == 3) {
            colors <- c("#00599F", "#FFA500", "#D01910")
        }
        else {
            colors <- (scales::hue_pal())(n_groups)
        }
        names(colors) <- group_levels
    }
    else {
        if (is.null(names(colors))) {
            if (length(colors) == n_groups) {
                names(colors) <- group_levels
            }
            else {
                stop("颜色参数长度必须与分组水平数量匹配")
            }
        }
        else {
            if (!all(group_levels %in% names(colors))) {
                stop("颜色名称需与分组水平完全匹配: ", paste(group_levels, collapse = ", "))
            }
            colors <- colors[group_levels]
        }
    }
    if (is.null(comparisons)) {
        if (n_groups == 2) {
            comparisons <- list(group_levels)
        }
        else {
            comparisons <- combn(group_levels, 2, simplify = FALSE)
        }
    }
    else {
        for (comp in comparisons) {
            if (!all(comp %in% group_levels)) {
                stop(paste("比较组", paste(comp, collapse = " vs "), "中存在不在分组中的水平"))
            }
        }
    }
    if (n_groups == 2) {
        if (!method %in% c("t.test", "wilcox.test")) {
            warning("两组比较建议使用t.test或wilcox.test,已自动使用wilcox.test")
            method <- "wilcox.test"
        }
    }
    else {
        if (method %in% c("t.test", "wilcox.test")) {
            warning("多组比较将对每对进行", method, "检验")
        }
    }
    create_single_plot <- function(var_name) {
        p <- ggplot(data, aes(x = .data[[Group]], y = .data[[var_name]], fill = .data[[Group]])) + geom_violin(trim = FALSE, 
            alpha = 0.8, width = 0.7, color = NA) + geom_boxplot(width = 0.2, color = "black", alpha = 0, outlier.shape = NA, 
            linewidth = 0.5)
        if (show_jitter) {
            p <- p + geom_jitter(size = jitter_size, shape = 21, color = "black", position = position_jitter(width = 0.1, 
                height = 0), alpha = jitter_alpha)
        }
        p <- p + scale_fill_manual(values = colors)
        p <- p + stat_compare_means(method = method, comparisons = comparisons, label = "p.signif", size = 6, tip.length = 0, 
            bracket.size = 0.6)
        if (n_groups > 2 && method %in% c("anova", "kruskal.test")) {
            p <- p + stat_compare_means(method = method, label.y.npc = 0.95, size = 3.5)
        }
        p <- p + labs(y = var_name, x = "", title = var_name) + theme_classic(base_size = 12) + theme(legend.position = "none", 
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), axis.text = element_text(color = "black", 
                size = 11), axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.y = element_text(size = 12, 
                margin = ggplot2::margin(r = 10)), axis.line = element_line(linewidth = 0.5), plot.title = element_text(hjust = 0.5, 
                size = 13, face = "bold"), plot.margin = unit(c(5, 15, 5, 5), "pt"))
        return(p)
    }
    if (length(variable) == 1) {
        p <- create_single_plot(variable)
        print(p)
        return(invisible(p))
    }
    plot_list <- lapply(variable, create_single_plot)
    names(plot_list) <- variable
    if (is.null(ncol)) {
        n_vars <- length(variable)
        if (n_vars <= 2) {
            ncol <- 2
        }
        else if (n_vars <= 4) {
            ncol <- 2
        }
        else if (n_vars <= 6) {
            ncol <- 3
        }
        else {
            ncol <- 3
        }
    }
    combined_plot <- wrap_plots(plot_list, ncol = ncol)
    print(combined_plot)
    return(invisible(list(combined = combined_plot, individual = plot_list)))
}

