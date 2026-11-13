# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
if (TRUE) {
    library(tidyverse)
    library(survminer)
    library(openxlsx)
    library(survival)
    library(openxlsx)
    library(timeROC)
    library(ggplot2)
    library(Mime1)
    invisible(NULL)
}

prad_restore("input_025")

dat_list <- list()

for (i in 1:length(total_expr_list)) {
    id <- names(total_expr_list)[i]
    print(id)
    expi <- total_expr_list[[i]]
    clii <- total_clin_list[[i]]
    if ("ID" %in% colnames(clii)) {
        rownames(clii) <- NULL
        clii <- clii %>% column_to_rownames("ID")
    }
    tryCatch({
        if (!"RFS" %in% colnames(clii)) {
            stop(paste0(id, " 没有预后信息"))
        }
        else {
            clii <- clii %>% dplyr::select(RFS.time, RFS) %>% na.omit()
            com <- intersect(colnames(expi), rownames(clii))
            expi <- expi[, com]
            clii <- clii[com, ]
            expi <- as.data.frame(t(expi))
            if (!identical(rownames(clii), rownames(expi))) {
                stop("clii 和 expi 的行名不匹配")
            }
            dati <- cbind(clii, expi) %>% tibble::rownames_to_column("ID")
            dati <- list(a = dati)
            names(dati) <- id
            if (!exists("dat_list")) {
                dat_list <- list()
            }
            dat_list <- append(dat_list, dati)
        }
    }, error = function(e) {
        message("发生错误: ", e$message)
    }, finally = {
        message("处理完成或已跳过错误。")
    })
}

expi <- total_expr_list$TCGA_PRAD %>% as.data.frame()

clii <- total_clin_list$TCGA_PRAD %>% dplyr::select(PFS.time, PFS) %>% na.omit() %>% dplyr::rename(RFS.time = "PFS.time", 
    RFS = "PFS")

com <- intersect(colnames(expi), rownames(clii))

expi <- expi[, com]

clii <- clii[com, ]

expi <- as.data.frame(t(expi))

dati <- cbind(clii, expi) %>% tibble::rownames_to_column("ID")

dati <- list(TCGA_PRAD = dati)

dat_list <- append(dat_list, dati)

for (i in 1:length(dat_list)) {
    dati <- dat_list[[i]]$RFS.time
    print(range(dati))
}

for (i in 1:length(dat_list)) {
    dat_list[[i]]$RFS.time <- as.integer(round(dat_list[[i]]$RFS.time * 365))
}

for (i in 1:length(dat_list)) {
    dati <- dat_list[[i]]$RFS.time
    print(range(dati))
}

for (i in 1:length(dat_list)) {
    colnames(dat_list[[i]])[1:3] <- c("ID", "OS.time", "OS")
}

gene <- prad_input("input_026")

gene <- gene$candidata %>% as.character() %>% unique()

target_datasets <- c("DKFZ2018", "GSE70769", "GSE54460", "GSE21034", "TCGA_PRAD", "GSE107299", "GSE116918", "GSE70768")

omit_datasets <- setdiff(names(dat_list), target_datasets)

all_gene_sets_base <- append(list(gene), lapply(dat_list[target_datasets], colnames))

com_gene <- Reduce(intersect, all_gene_sets_base)

com_vector <- c("ID", "OS.time", "OS", all_of(com_gene))

dat_list <- dat_list[target_datasets]

for (i in 1:length(dat_list)) {
    dat_list[[i]] <- dat_list[[i]][, com_vector] %>% na.omit()
}

res_ML <- ML.Dev.Prog.Sig_modified(train_data = dat_list$TCGA_PRAD, list_train_vali_Data = dat_list, unicox.filter.for.candi = TRUE, 
    candidate_genes = com_gene, mode = "all", nodesize = 5, seed = 5201314)

dat_use <- res_ML$Cindex.res %>% as.data.frame()

wide_data <- dat_use %>% pivot_wider(names_from = Model, values_from = Cindex) %>% column_to_rownames("ID") %>% 
    t() %>% as.data.frame()

wide_data <- as.data.frame(apply(wide_data, 2, function(x) round(x, 2)))

col <- RColorBrewer::brewer.pal(n = length(colnames(wide_data)), name = "Set1")

plot_heatmap_1(hplot_input = wide_data, bar_color_use = "steelblue", legend_heatmap_use = "Heatmap Legend", legend_topann_use = "Top Annotation Legend", 
    cellwidth_use = 0.8, cellheight_use = 0.8, fontsize_use = 10, topann_color_use = col, n = 115, heatmap_name_use = "example_heatmap")

dat_use <- res_ML$Cindex.res %>% as.data.frame()

wide_data <- dat_use %>% pivot_wider(names_from = Model, values_from = Cindex) %>% column_to_rownames("ID") %>% 
    t() %>% as.data.frame() %>% dplyr::select(., c("DKFZ2018", "GSE70769", "GSE21034", "TCGA_PRAD", "GSE107299", 
    "GSE116918", "GSE70768"))

wide_data <- as.data.frame(apply(wide_data, 2, function(x) round(x, 3)))

col <- RColorBrewer::brewer.pal(n = length(colnames(wide_data)), name = "Set1")

plot_heatmap_1(hplot_input = wide_data, bar_color_use = "steelblue", legend_heatmap_use = "Heatmap Legend", legend_topann_use = "Top Annotation Legend", 
    cellwidth_use = 0.8, cellheight_use = 0.8, fontsize_use = 10, topann_color_use = col, n = 115, heatmap_name_use = "example_heatmap_plus")

model = "StepCox[forward] + Ridge"

dat_risk <- res_ML$riskscore[[model]]

cv_fit <- res_ML[["ml.res"]][["StepCox[forward] + Ridge"]]

coef_values <- beta_matrix[, cv_fit$index["min", "Lambda"]]

res_coef <- data.frame(gene = rownames(beta_matrix), coefficient = coef_values)

dat_cindex <- res_ML$Cindex.res %>% as.data.frame()

cv_data <- data.frame(log_lambda = log(cv_fit$lambda), cvm = cv_fit$cvm, cvlo = cv_fit$cvlo, cvup = cv_fit$cvup, 
    nzero = cv_fit$nzero)

ggplot(cv_data, aes(x = log_lambda, y = cvm)) + geom_point(color = "#D01910", size = 2) + geom_line(color = "#D01910") + 
    geom_errorbar(aes(ymin = cvlo, ymax = cvup), width = 0.1, alpha = 0.3, color = "gray50") + geom_vline(xintercept = log(cv_fit$lambda.min), 
    linetype = "dashed", color = "#D01910", size = 1) + geom_vline(xintercept = log(cv_fit$lambda.1se), linetype = "dashed", 
    color = "#087a0b", size = 1) + annotate("text", x = log(cv_fit$lambda.min), y = max(cv_data$cvup), label = "lambda.min", 
    color = "#00599F", vjust = -0.5) + annotate("text", x = log(cv_fit$lambda.1se), y = max(cv_data$cvup), label = "lambda.1se", 
    color = "#087a0b", vjust = -0.5) + labs(title = "Ridge Regression: Cross-Validation Curve", x = "Log(Lambda)", 
    y = "Partial Likelihood Deviance") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold", 
    size = 14))

invisible(NULL)

beta_matrix <- as.matrix(cv_fit$glmnet.fit$beta)

lambda_values <- cv_fit$glmnet.fit$lambda

coef_data <- data.frame(lambda = rep(lambda_values, each = nrow(beta_matrix)), log_lambda = rep(log(lambda_values), 
    each = nrow(beta_matrix)), gene = rep(rownames(beta_matrix), times = ncol(beta_matrix)), coefficient = as.vector(beta_matrix))

ggplot(coef_data, aes(x = log_lambda, y = coefficient, color = gene)) + geom_line(size = 1.2) + geom_vline(xintercept = log(cv_fit$lambda.min), 
    linetype = "dashed", color = "#D01910", size = 1) + geom_vline(xintercept = log(cv_fit$lambda.1se), linetype = "dashed", 
    color = "#087a0b", size = 1) + annotate("text", x = log(cv_fit$lambda.min), y = max(coef_data$coefficient), 
    label = "lambda.min", color = "#D01910", vjust = -0.5, size = 3.5) + annotate("text", x = log(cv_fit$lambda.1se), 
    y = max(coef_data$coefficient), label = "lambda.1se", color = "#087a0b", vjust = -0.5, size = 3.5) + labs(title = "Ridge Regression: Coefficient Path", 
    x = "Log(Lambda)", y = "Coefficients", color = "Gene") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, 
    face = "bold", size = 14), legend.position = "right", legend.text = element_text(size = 10))

invisible(NULL)

res_model <- list(dat_risk = dat_risk, dat_cindex = dat_cindex, dat_coef = res_coef)

invisible(NULL)

for (i in 1:length(dat_risk)) {
    colnames(dat_risk[[i]])[2:3] <- c("RFS.time", "RFS")
}

dat_use <- dat_risk$TCGA_PRAD %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID")

plot_survival_KM_1(data = dat_use, minprop = 0.3, GeneName = "RS", time = "RFS.time", event = "RFS", Timeunit = "d")

dat_use <- dat_risk$DKFZ2018 %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID")

plot_survival_KM_1(data = dat_use, minprop = 0.3, GeneName = "RS", time = "RFS.time", event = "RFS", Timeunit = "d")

dat_use <- dat_risk$GSE70769 %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID")

plot_survival_KM_1(data = dat_use, minprop = 0.3, GeneName = "RS", time = "RFS.time", event = "RFS", Timeunit = "d")

dat_use <- dat_risk$GSE54460 %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID")

plot_survival_KM_1(data = dat_use, minprop = 0.3, GeneName = "RS", time = "RFS.time", event = "RFS", Timeunit = "d")

dat_use <- dat_risk$GSE21034 %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID")

plot_survival_KM_1(data = dat_use, minprop = 0.3, GeneName = "RS", time = "RFS.time", event = "RFS", Timeunit = "d")

dat_use <- dat_risk$GSE107299 %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID")

plot_survival_KM_1(data = dat_use, minprop = 0.3, GeneName = "RS", time = "RFS.time", event = "RFS", Timeunit = "d")

dat_use <- dat_risk$GSE116918 %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID")

plot_survival_KM_1(data = dat_use, minprop = 0.3, GeneName = "RS", time = "RFS.time", event = "RFS", Timeunit = "d")

dat_use <- dat_risk$GSE70768 %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID")

plot_survival_KM_1(data = dat_use, minprop = 0.3, GeneName = "RS", time = "RFS.time", event = "RFS", Timeunit = "d")

dat_use <- dat_risk$TCGA_PRAD %>% as.data.frame() %>% column_to_rownames("ID")

dat_use$RFS.time <- dat_use$RFS.time/365

plot_Time_ROC(dat = dat_use, target = "RS", year = TRUE)

dat_use <- dat_risk$DKFZ2018 %>% as.data.frame() %>% column_to_rownames("ID")

dat_use$RFS.time <- dat_use$RFS.time/365

plot_Time_ROC(dat = dat_use, target = "RS", year = TRUE)

dat_use <- dat_risk$GSE70769 %>% as.data.frame() %>% column_to_rownames("ID")

dat_use$RFS.time <- dat_use$RFS.time/365

plot_Time_ROC(dat = dat_use, target = "RS", year = TRUE)

dat_use <- dat_risk$GSE21034 %>% as.data.frame() %>% column_to_rownames("ID")

dat_use$RFS.time <- dat_use$RFS.time/365

plot_Time_ROC(dat = dat_use, target = "RS", year = TRUE)

dat_use <- dat_risk$GSE70768 %>% as.data.frame() %>% column_to_rownames("ID")

dat_use$RFS.time <- dat_use$RFS.time/365

plot_Time_ROC(dat = dat_use, target = "RS", year = TRUE)

dat_use <- dat_risk$GSE116918 %>% as.data.frame() %>% column_to_rownames("ID")

dat_use$RFS.time <- dat_use$RFS.time/365

plot_Time_ROC(dat = dat_use, target = "RS", year = TRUE)

