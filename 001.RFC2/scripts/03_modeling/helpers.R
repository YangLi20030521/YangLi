# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
plot_heatmap_1 <- function(hplot_input, bar_color_use = "steelblue", legend_heatmap_use, legend_topann_use, cellwidth_use = 1, 
    cellheight_use = 0.5, fontsize_use = 10, n = 50, topann_color_use = c("#BC3C29FF", "#0072B5FF", "#E18727FF", 
        "#20854EFF", "#7876B1FF"), heatmap_name_use = "heatmap_name_use") {
    library(tidyverse)
    library(ComplexHeatmap)
    hplot <- function(hplot_input, right_annotation_use, topann_color_use, bar_color_use, cellwidth, cellheight, 
        legend_heatmap_use, legend_topann_use, fontsize_use) {
        col_ha <- columnAnnotation(legend_topann_use = colnames(hplot_input), col = list(legend_topann_use = topann_color_use), 
            show_annotation_name = F)
        row_ha <- rowAnnotation(bar = anno_barplot(right_annotation_use, bar_width = 0.8, border = FALSE, gp = gpar(fill = bar_color_use, 
            col = NA), add_numbers = TRUE, numbers_offset = unit(-10, "mm"), axis_param = list(labels_rot = 0), 
            numbers_gp = gpar(fontsize = fontsize_use, col = "white"), width = unit(3, "cm")), show_annotation_name = FALSE)
        Heatmap(as.matrix(hplot_input), name = legend_heatmap_use, right_annotation = row_ha, top_annotation = col_ha, 
            col = c("#27879f", "#e7efc7", "#ef7740"), rect_gp = gpar(col = "black", lwd = 1), cluster_columns = FALSE, 
            cluster_rows = FALSE, show_column_names = FALSE, show_row_names = TRUE, row_names_side = "left", width = unit(cellwidth * 
                ncol(hplot_input) + 2, "cm"), height = unit(cellheight * nrow(hplot_input), "cm"), column_split = factor(colnames(hplot_input), 
                levels = colnames(hplot_input)), column_title = NULL, cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(label = format(hplot_input[i, j], digits = 2, nsmall = 2), x, y, gp = gpar(fontsize = 10))
            })
    }
    for (i in 1:ncol(hplot_input)) {
        hplot_input[, i] <- as.numeric(hplot_input[, i])
    }
    df <- hplot_input %>% mutate(mean = apply(hplot_input, 1, mean)) %>% arrange(., desc(mean)) %>% dplyr::slice_head(n = n)
    right_annotation_use <- round(df$mean, 2)
    hplot_input <- df[, -ncol(df)]
    topann_color_use <- setNames(topann_color_use[1:length(colnames(hplot_input))], colnames(hplot_input))
    hm <- hplot(hplot_input, right_annotation_use, topann_color_use, bar_color_use, cellwidth_use, cellheight_use, 
        legend_heatmap_use, legend_topann_use, fontsize_use)
    existing_color_mapping <- hm@top_annotation@anno_list[["legend_topann_use"]]@color_mapping
    colors <- existing_color_mapping@colors
    new_color_mapping <- ColorMapping(levels = colnames(hplot_input), colors = colors[colnames(hplot_input)], name = legend_topann_use)
    hm@top_annotation@anno_list[["legend_topann_use"]]@color_mapping <- new_color_mapping
    invisible(NULL)
    draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")
    invisible(invisible(NULL))
}

plot_Time_ROC <- function(dat, target, year = TRUE) {
    library(timeROC)
    library(ggplot2)
    if (year) {
        dat <- dat
    }
    else {
        stop("时间变量必须以年为单位")
    }
    if (all(c("RFS", "RFS.time") %in% colnames(dat))) {
        print("数据具有所分析变量，运行下面代码")
    }
    else {
        stop("数据框缺少RFS或RFS.time列，请检查")
    }
    marker <- dat[, target, drop = T]
    timeroc <- timeROC(T = dat$RFS.time, delta = dat$RFS, marker = marker, cause = 1, weighting = "cox", times = c(1, 
        3, 5), iid = F)
    dat1 <- timeroc$TP %>% as.data.frame()
    dat2 <- timeroc$FP %>% as.data.frame()
    t1 <- data.frame(tp = dat1$`t=1`, fp = dat2$`t=1`) %>% mutate(time = "1-Year")
    t2 <- data.frame(tp = dat1$`t=3`, fp = dat2$`t=3`) %>% mutate(time = "3-Year")
    t3 <- data.frame(tp = dat1$`t=5`, fp = dat2$`t=5`) %>% mutate(time = "5-Year")
    tmp <- rbind(t1, t2, t3)
    p <- ggplot(tmp, aes(fp, tp, color = time)) + geom_line(size = 1) + labs(x = "1-Specificity", y = "Sensitivity", 
        color = NULL) + theme_bw(base_rect_size = 1.5) + geom_abline(slope = 1, color = "grey70") + ggtitle("TCGA Cohort") + 
        theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 15, colour = "darkred", 
            face = "bold"), panel.grid = element_blank(), panel.grid.major = element_line(color = "white", linetype = "dashed"), 
            panel.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5, size = 18, 
                colour = "darkred", face = "bold"), legend.text = element_text(size = 14), legend.position = c(0.99, 
                0.013), legend.justification = c(1, 0)) + scale_color_manual(values = c("#2a9d8f", "#fcbf49", "#eb5e28")) + 
        scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01))
    print(p)
}

plot_Calibration_Curve <- function(dat, target, year = TRUE) {
    library(ggplot2)
    library(rms)
    tmp <- dat
    index <- which(colnames(tmp) == target)
    colnames(tmp)[index] <- "score"
    tmp$score <- as.numeric(tmp$score)
    dd <- rms::datadist(tmp)
    options(datadist = "dd")
    if (year) {
        tmp <- tmp
    }
    else {
        stop("时间变量必须以年为单位")
    }
    if (all(c("RFS", "RFS.time") %in% colnames(tmp))) {
        print("数据具有所分析变量，运行下面代码")
    }
    else {
        stop("数据框缺少RFS或RFS.time列，请检查")
    }
    cox1 <- cph(Surv(RFS.time, RFS) ~ score, data = tmp, x = TRUE, y = TRUE, surv = TRUE, time.inc = 1)
    cal1 <- tryCatch(calibrate(cox1, cmethod = "KM", method = "boot", u = 1, m = ceiling(nrow(tmp) * 0.95/5), B = 1000), 
        error = function(e) {
            NA
        })
    cox2 <- cph(Surv(RFS.time, RFS) ~ score, data = tmp, x = TRUE, y = TRUE, surv = TRUE, time.inc = 3)
    cal2 <- tryCatch(calibrate(cox2, cmethod = "KM", method = "boot", u = 3, m = ceiling(nrow(tmp) * 0.95/5), B = 1000), 
        error = function(e) {
            NA
        })
    cox3 <- cph(Surv(RFS.time, RFS) ~ score, data = tmp, x = TRUE, y = TRUE, surv = TRUE, time.inc = 5)
    cal3 <- tryCatch(calibrate(cox3, cmethod = "KM", method = "boot", u = 5, m = ceiling(nrow(tmp) * 0.95/5), B = 1000), 
        error = function(e) {
            NA
        })
    dd <- Reduce(rbind, list(tryCatch(cal1[, c("mean.predicted", "KM")], error = function(e) data.frame(mean.predicted = rep(NA, 
        5), KM = rep(NA, 5))), tryCatch(cal2[, c("mean.predicted", "KM")], error = function(e) data.frame(mean.predicted = rep(NA, 
        5), KM = rep(NA, 5))), tryCatch(cal3[, c("mean.predicted", "KM")], error = function(e) data.frame(mean.predicted = rep(NA, 
        5), KM = rep(NA, 5))))) %>% as.data.frame()
    time_cut <- c("1-Year AUC", "3-Year AUC", "5-Year AUC")
    dd$time <- paste0(rep(gsub(" AUC", "", time_cut), each = 5), " Calibration")
    dd <- na.omit(dd)
    colnames(dd)[1:2] <- c("Predicted", "Observed")
    p <- ggplot(dd, aes(Predicted, Observed)) + geom_abline(slope = 1, color = "grey70") + geom_line(size = 1, 
        alpha = 0.7, aes(color = time)) + theme_bw(base_rect_size = 1.5) + geom_point(size = 3, shape = 21, aes(fill = time), 
        stroke = 1) + labs(x = "Predicted survival (%)", y = "Observed survival (%)", title = "") + theme(axis.text.y = element_text(size = 12, 
        colour = "black"), axis.text.x = element_text(size = 12, colour = "black"), axis.title = element_text(size = 15, 
        colour = "darkred", face = "bold"), panel.grid = element_blank(), panel.grid.major = element_line(color = "white", 
        linetype = "dashed"), panel.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5, 
        size = 18, colour = "darkred", face = "bold"), legend.text = element_text(size = 14), legend.background = element_blank(), 
        legend.key = element_blank(), legend.title = element_blank(), legend.position = c(0.99, 0.013), legend.justification = c(1, 
            0)) + scale_color_manual(values = c("#2a9d8f", "#fcbf49", "#eb5e28")) + scale_fill_manual(values = c("#2a9d8f", 
        "#fcbf49", "#eb5e28")) + scale_x_continuous(expand = c(0.03, 0)) + scale_y_continuous(expand = c(0, 0.03))
    print(p)
}

plot_Time_ROC <- function(dat, target, year = TRUE) {
    library(timeROC)
    library(ggplot2)
    if (year) {
        dat <- dat
    }
    else {
        stop("时间变量必须以年为单位")
    }
    if (all(c("RFS", "RFS.time") %in% colnames(dat))) {
        print("数据具有所分析变量，运行下面代码")
    }
    else {
        stop("数据框缺少RFS或RFS.time列，请检查")
    }
    marker <- dat[, target, drop = T]
    timeroc <- timeROC(T = dat$RFS.time, delta = dat$RFS, marker = marker, cause = 1, weighting = "cox", times = c(3, 
        5, 8), iid = F)
    auc_3y <- round(timeroc$AUC[1], 3)
    auc_5y <- round(timeroc$AUC[2], 3)
    auc_8y <- round(timeroc$AUC[3], 3)
    dat1 <- timeroc$TP %>% as.data.frame()
    dat2 <- timeroc$FP %>% as.data.frame()
    t1 <- data.frame(tp = dat1$`t=3`, fp = dat2$`t=3`) %>% mutate(time = paste0("3-Year (AUC = ", auc_3y, ")"))
    t2 <- data.frame(tp = dat1$`t=5`, fp = dat2$`t=5`) %>% mutate(time = paste0("5-Year (AUC = ", auc_5y, ")"))
    t3 <- data.frame(tp = dat1$`t=8`, fp = dat2$`t=8`) %>% mutate(time = paste0("8-Year (AUC = ", auc_8y, ")"))
    tmp <- rbind(t1, t2, t3)
    p <- ggplot(tmp, aes(fp, tp, color = time)) + geom_line(size = 1) + labs(x = "1-Specificity", y = "Sensitivity", 
        color = NULL) + theme_bw(base_rect_size = 1.5) + geom_abline(slope = 1, color = "grey70") + ggtitle("TCGA Cohort") + 
        theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 15, colour = "darkred", 
            face = "bold"), panel.grid = element_blank(), panel.grid.major = element_line(color = "white", linetype = "dashed"), 
            panel.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5, size = 18, 
                colour = "darkred", face = "bold"), legend.text = element_text(size = 14), legend.position = c(0.99, 
                0.013), legend.justification = c(1, 0)) + scale_color_manual(values = c("#2a9d8f", "#fcbf49", "#eb5e28")) + 
        scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01))
    print(p)
    invisible(list(plot = p, AUC_3year = auc_3y, AUC_5year = auc_5y, AUC_8year = auc_8y))
}

