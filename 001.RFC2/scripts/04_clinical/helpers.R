# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
plot_expr_diff_Violin_plus <- function(data, Group = "Group", variable = NULL, colors = NULL, comparisons = NULL, 
    method = "wilcox.test", ncol = NULL, show_jitter = TRUE, jitter_size = 2.5, jitter_alpha = 0.8) {
    library(ggplot2)
    library(ggpubr)
    library(dplyr)
    library(tidyr)
    library(patchwork)
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
        p <- p + stat_compare_means(method = method, comparisons = comparisons, label = "p.signif", size = 5, tip.length = 0, 
            bracket.size = 0.6)
        if (n_groups > 2 && method %in% c("anova", "kruskal.test")) {
            p <- p + stat_compare_means(method = method, label.y.npc = 0.95, size = 3.5)
        }
        p <- p + labs(y = var_name, x = "", title = var_name) + theme_classic(base_size = 12) + theme(legend.position = "none", 
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), axis.text = element_text(color = "black", 
                size = 11), axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.y = element_text(size = 12, 
                margin = margin(r = 10)), axis.line = element_line(linewidth = 0.5), plot.title = element_text(hjust = 0.5, 
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

transform_pathway_correlation <- function(data, gene_col = NULL, method = "pearson") {
    library(dplyr)
    library(tidyr)
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
    long_data <- data %>% pivot_longer(cols = all_of(pathway_cols), names_to = "pathway", values_to = "pathway_score") %>% 
        rename(gene_expression = !!sym(gene_col))
    correlation_results <- long_data %>% group_by(pathway) %>% summarise(correlation = cor(gene_expression, pathway_score, 
        method = method, use = "complete.obs"), p_value = cor.test(gene_expression, pathway_score, method = method)$p.value, 
        .groups = "drop")
    result_data <- long_data %>% left_join(correlation_results, by = "pathway")
    result_data <- result_data %>% select(gene_expression, pathway_score, pathway, correlation, p_value)
    return(result_data)
}

plot_tcga_correlation <- function(data, x_var = "MKI67_zscore", y_var = "CellCycle_combined_zscore", cancer_col = "cancer_type", 
    ncol = 11, point_size = 1.5, point_alpha = 0.6) {
    library(ggplot2)
    library(dplyr)
    library(cowplot)
    cancer_list <- unique(data[[cancer_col]])
    colors <- c("#4A90E2", "#E8743B", "#19A979", "#9B59B6", "#E74C3C", "#3498DB", "#1ABC9C", "#F39C12", "#E91E63", 
        "#9C27B0", "#00BCD4", "#4CAF50", "#FF9800", "#795548", "#607D8B", "#FF5722", "#673AB7", "#03A9F4", "#8BC34A", 
        "#FFC107", "#00897B", "#D32F2F", "#7B1FA2", "#0288D1", "#689F38", "#F57C00", "#5D4037", "#455A64", "#C62828", 
        "#6A1B9A", "#0277BD", "#558B2F", "#E65100")
    plot_list <- lapply(cancer_list, function(cancer) {
        cancer_data <- data %>% filter(!!sym(cancer_col) == cancer)
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
            size = point_size, alpha = point_alpha) + geom_smooth(method = "lm", se = FALSE, color = "#E74C3C", 
            linewidth = 1.2, linetype = "solid") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
            scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) + labs(x = NULL, y = NULL) + theme_classic(base_size = 10) + 
            theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), axis.line = element_blank(), 
                axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = margin(0, 1, 1, 1), panel.background = element_rect(fill = "white", 
                  color = NA))
        title_box <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) + annotate("rect", xmin = 0, xmax = 1, ymin = 0, 
            ymax = 1, fill = "#E8E8E8", color = NA) + annotate("text", x = 0.5, y = 0.5, label = paste0(cancer, 
            "\n", stat_label), fontface = "bold", size = 2.8, hjust = 0.5, vjust = 0.5, lineheight = 0.8) + scale_x_continuous(limits = c(0, 
            1), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + theme_void() + theme(panel.border = element_rect(fill = NA, 
            color = "black", linewidth = 0.8), plot.margin = margin(1, 1, 0, 1), panel.background = element_rect(fill = "#E8E8E8", 
            color = NA))
        combined <- plot_grid(title_box, p_main, ncol = 1, rel_heights = c(0.15, 1), align = "v", axis = "lr")
        return(combined)
    })
    combined_plot <- plot_grid(plotlist = plot_list, ncol = ncol, align = "hv", axis = "tblr")
    return(combined_plot)
}

