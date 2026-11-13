# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
if (TRUE) {
    library(tidyverse)
    library(ggplot2)
    library(cowplot)
    library(dplyr)
    number = 20
    invisible(NULL)
    prad_restore("input_053")
    res_model <- prad_input("input_054")
}

dat_exp <- res_model$dat_risk$TCGA_PRAD %>% tibble::remove_rownames() %>% column_to_rownames("ID") %>% dplyr::select("RS")

dat_drug <- Drug_list$GDSC1$TCGA_PRAD %>% as.data.frame() %>% log1p()

ptypeOut <- dat_drug

data_meta <- dat_exp

top.pps <- data_meta[data_meta[, 1] >= quantile(data_meta[, 1], probs = seq(0, 1, 0.1))[10], , drop = F]

bot.pps <- data_meta[data_meta[, 1] <= quantile(data_meta[, 1], probs = seq(0, 1, 0.1))[2], , drop = F]

ccle.cor <- ccle.cor.p <- c()

com <- intersect(rownames(ptypeOut), rownames(data_meta))

data_meta <- data_meta[com, , drop = F]

ptypeOut <- ptypeOut[com, ]

for (i in 1:ncol(ptypeOut)) {
    d <- colnames(ptypeOut)[i]
    a <- as.numeric(ptypeOut[, d])
    b <- as.numeric(data_meta[, 1])
    r <- cor.test(a, b, method = "spearman")$estimate
    names(r) <- d
    p <- cor.test(a, b, method = "spearman")$p.value
    names(p) <- d
    ccle.cor <- c(ccle.cor, r)
    ccle.cor.p <- c(ccle.cor.p, p)
}

p.values <- as.numeric(ccle.cor.p)

names <- names(ccle.cor.p)

drug_data_1 <- data.frame(Drug_Name = names, p.value = p.values)

correlation <- as.numeric(ccle.cor)

names <- names(ccle.cor)

drug_data_2 <- data.frame(Drug_Name = names, correlation = correlation)

drug_data <- merge(drug_data_2, drug_data_1, by.x = 1, by.y = 1)

drug_data$state <- ifelse(drug_data$cor > 0, "positive", "negative")

drug <- drug_data %>% dplyr::group_by(state) %>% dplyr::filter(p.value < 0.05) %>% dplyr::arrange(desc(abs(correlation))) %>% 
    dplyr::slice_head(n = number) %>% dplyr::pull(Drug_Name)

darkblue <- "#0772B9"

lightblue <- "#48C8EF"

cor.data <- data.frame(drug = drug, r = ccle.cor[drug], p = -log10(ccle.cor.p[drug]))

cor.data$p[cor.data$p == Inf] <- 15

ggplot(data = cor.data, aes(r, forcats::fct_reorder(drug, r, .desc = T))) + geom_segment(aes(xend = 0, yend = drug), 
    linetype = 2) + geom_point(aes(size = p), col = darkblue) + scale_size_continuous(range = c(2, 8)) + scale_x_reverse(breaks = c(0, 
    -0.3, -0.5), expand = expansion(mult = c(0.01, 0.1))) + theme_classic() + labs(x = "Correlation coefficient", 
    y = "", size = bquote("-log"[10] ~ "(" ~ italic(P) ~ "-value)")) + theme(legend.position = "bottom", axis.line.y = element_blank())

invisible(NULL)

ccle.boxdata <- NULL

for (d in drug) {
    a <- as.numeric(ptypeOut[rownames(top.pps), d])
    b <- as.numeric(ptypeOut[rownames(bot.pps), d])
    p <- wilcox.test(a, b)$p.value
    s <- as.character(cut(p, c(0, 0.001, 0.01, 0.05, 1), labels = c("***", "**", "*", "")))
    ccle.boxdata <- rbind.data.frame(ccle.boxdata, data.frame(drug = d, auc = c(a, b), p = p, s = s, group = rep(c("High_signature", 
        "Low_signature"), c(nrow(top.pps), nrow(bot.pps))), stringsAsFactors = F), stringsAsFactors = F)
}

p2 <- ggplot(ccle.boxdata, aes(drug, auc, fill = group)) + geom_boxplot(aes(col = group), outlier.shape = NA) + 
    geom_text(aes(drug, y = max(auc)), label = ccle.boxdata$s, data = ccle.boxdata, inherit.aes = F) + scale_fill_manual(values = c(darkblue, 
    lightblue)) + scale_color_manual(values = c(darkblue, lightblue)) + xlab(NULL) + ylab("Estimated AUC value") + 
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 10), legend.position = "bottom", 
    legend.title = element_blank())

dat <- ggplot_build(p2)$data[[1]]

p2 + geom_segment(data = dat, aes(x = xmin, xend = xmax, y = middle, yend = middle), color = "white", inherit.aes = F)

invisible(NULL)

dat_exp <- res_model$dat_risk$TCGA_PRAD %>% tibble::remove_rownames() %>% column_to_rownames("ID") %>% dplyr::select("RS")

dat_drug <- Drug_list$GDSC2$TCGA_PRAD %>% as.data.frame() %>% log1p()

ptypeOut <- dat_drug

data_meta <- dat_exp

top.pps <- data_meta[data_meta[, 1] >= quantile(data_meta[, 1], probs = seq(0, 1, 0.1))[10], , drop = F]

bot.pps <- data_meta[data_meta[, 1] <= quantile(data_meta[, 1], probs = seq(0, 1, 0.1))[2], , drop = F]

ccle.cor <- ccle.cor.p <- c()

com <- intersect(rownames(ptypeOut), rownames(data_meta))

data_meta <- data_meta[com, , drop = F]

ptypeOut <- ptypeOut[com, ]

for (i in 1:ncol(ptypeOut)) {
    d <- colnames(ptypeOut)[i]
    a <- as.numeric(ptypeOut[, d])
    b <- as.numeric(data_meta[, 1])
    r <- cor.test(a, b, method = "spearman")$estimate
    names(r) <- d
    p <- cor.test(a, b, method = "spearman")$p.value
    names(p) <- d
    ccle.cor <- c(ccle.cor, r)
    ccle.cor.p <- c(ccle.cor.p, p)
}

p.values <- as.numeric(ccle.cor.p)

names <- names(ccle.cor.p)

drug_data_1 <- data.frame(Drug_Name = names, p.value = p.values)

correlation <- as.numeric(ccle.cor)

names <- names(ccle.cor)

drug_data_2 <- data.frame(Drug_Name = names, correlation = correlation)

drug_data <- merge(drug_data_2, drug_data_1, by.x = 1, by.y = 1)

drug_data$state <- ifelse(drug_data$cor > 0, "positive", "negative")

drug <- drug_data %>% dplyr::group_by(state) %>% dplyr::filter(p.value < 0.05) %>% dplyr::arrange(desc(abs(correlation))) %>% 
    dplyr::slice_head(n = number) %>% dplyr::pull(Drug_Name)

darkblue <- "#0772B9"

lightblue <- "#48C8EF"

cor.data <- data.frame(drug = drug, r = ccle.cor[drug], p = -log10(ccle.cor.p[drug]))

cor.data$p[cor.data$p == Inf] <- 15

ggplot(data = cor.data, aes(r, forcats::fct_reorder(drug, r, .desc = T))) + geom_segment(aes(xend = 0, yend = drug), 
    linetype = 2) + geom_point(aes(size = p), col = darkblue) + scale_size_continuous(range = c(2, 8)) + scale_x_reverse(breaks = c(0, 
    -0.3, -0.5), expand = expansion(mult = c(0.01, 0.1))) + theme_classic() + labs(x = "Correlation coefficient", 
    y = "", size = bquote("-log"[10] ~ "(" ~ italic(P) ~ "-value)")) + theme(legend.position = "bottom", axis.line.y = element_blank())

invisible(NULL)

ccle.boxdata <- NULL

for (d in drug) {
    a <- as.numeric(ptypeOut[rownames(top.pps), d])
    b <- as.numeric(ptypeOut[rownames(bot.pps), d])
    p <- wilcox.test(a, b)$p.value
    s <- as.character(cut(p, c(0, 0.001, 0.01, 0.05, 1), labels = c("***", "**", "*", "")))
    ccle.boxdata <- rbind.data.frame(ccle.boxdata, data.frame(drug = d, auc = c(a, b), p = p, s = s, group = rep(c("High_signature", 
        "Low_signature"), c(nrow(top.pps), nrow(bot.pps))), stringsAsFactors = F), stringsAsFactors = F)
}

p2 <- ggplot(ccle.boxdata, aes(drug, auc, fill = group)) + geom_boxplot(aes(col = group), outlier.shape = NA) + 
    geom_text(aes(drug, y = max(auc)), label = ccle.boxdata$s, data = ccle.boxdata, inherit.aes = F) + scale_fill_manual(values = c(darkblue, 
    lightblue)) + scale_color_manual(values = c(darkblue, lightblue)) + xlab(NULL) + ylab("Estimated AUC value") + 
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 10), legend.position = "bottom", 
    legend.title = element_blank())

dat <- ggplot_build(p2)$data[[1]]

p2 + geom_segment(data = dat, aes(x = xmin, xend = xmax, y = middle, yend = middle), color = "white", inherit.aes = F)

invisible(NULL)

dat_exp <- res_model$dat_risk$TCGA_PRAD %>% tibble::remove_rownames() %>% column_to_rownames("ID") %>% dplyr::select("RS")

dat_drug <- Drug_list$CTRP$TCGA_PRAD %>% as.data.frame() %>% log1p()

ptypeOut <- dat_drug

data_meta <- dat_exp

top.pps <- data_meta[data_meta[, 1] >= quantile(data_meta[, 1], probs = seq(0, 1, 0.1))[10], , drop = F]

bot.pps <- data_meta[data_meta[, 1] <= quantile(data_meta[, 1], probs = seq(0, 1, 0.1))[2], , drop = F]

ccle.cor <- ccle.cor.p <- c()

com <- intersect(rownames(ptypeOut), rownames(data_meta))

data_meta <- data_meta[com, , drop = F]

ptypeOut <- ptypeOut[com, ]

for (i in 1:ncol(ptypeOut)) {
    d <- colnames(ptypeOut)[i]
    a <- as.numeric(ptypeOut[, d])
    b <- as.numeric(data_meta[, 1])
    r <- cor.test(a, b, method = "spearman")$estimate
    names(r) <- d
    p <- cor.test(a, b, method = "spearman")$p.value
    names(p) <- d
    ccle.cor <- c(ccle.cor, r)
    ccle.cor.p <- c(ccle.cor.p, p)
}

p.values <- as.numeric(ccle.cor.p)

names <- names(ccle.cor.p)

drug_data_1 <- data.frame(Drug_Name = names, p.value = p.values)

correlation <- as.numeric(ccle.cor)

names <- names(ccle.cor)

drug_data_2 <- data.frame(Drug_Name = names, correlation = correlation)

drug_data <- merge(drug_data_2, drug_data_1, by.x = 1, by.y = 1)

drug_data$state <- ifelse(drug_data$cor > 0, "positive", "negative")

drug <- drug_data %>% dplyr::group_by(state) %>% dplyr::filter(p.value < 0.05) %>% dplyr::arrange(desc(abs(correlation))) %>% 
    dplyr::slice_head(n = number) %>% dplyr::pull(Drug_Name)

darkblue <- "#0772B9"

lightblue <- "#48C8EF"

cor.data <- data.frame(drug = drug, r = ccle.cor[drug], p = -log10(ccle.cor.p[drug]))

cor.data$p[cor.data$p == Inf] <- 15

ggplot(data = cor.data, aes(r, forcats::fct_reorder(drug, r, .desc = T))) + geom_segment(aes(xend = 0, yend = drug), 
    linetype = 2) + geom_point(aes(size = p), col = darkblue) + scale_size_continuous(range = c(2, 8)) + scale_x_reverse(breaks = c(0, 
    -0.3, -0.5), expand = expansion(mult = c(0.01, 0.1))) + theme_classic() + labs(x = "Correlation coefficient", 
    y = "", size = bquote("-log"[10] ~ "(" ~ italic(P) ~ "-value)")) + theme(legend.position = "bottom", axis.line.y = element_blank())

invisible(NULL)

ccle.boxdata <- NULL

for (d in drug) {
    a <- as.numeric(ptypeOut[rownames(top.pps), d])
    b <- as.numeric(ptypeOut[rownames(bot.pps), d])
    p <- wilcox.test(a, b)$p.value
    s <- as.character(cut(p, c(0, 0.001, 0.01, 0.05, 1), labels = c("***", "**", "*", "")))
    ccle.boxdata <- rbind.data.frame(ccle.boxdata, data.frame(drug = d, auc = c(a, b), p = p, s = s, group = rep(c("High_signature", 
        "Low_signature"), c(nrow(top.pps), nrow(bot.pps))), stringsAsFactors = F), stringsAsFactors = F)
}

p2 <- ggplot(ccle.boxdata, aes(drug, auc, fill = group)) + geom_boxplot(aes(col = group), outlier.shape = NA) + 
    geom_text(aes(drug, y = max(auc)), label = ccle.boxdata$s, data = ccle.boxdata, inherit.aes = F) + scale_fill_manual(values = c(darkblue, 
    lightblue)) + scale_color_manual(values = c(darkblue, lightblue)) + xlab(NULL) + ylab("Estimated AUC value") + 
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 10), legend.position = "bottom", 
    legend.title = element_blank())

dat <- ggplot_build(p2)$data[[1]]

p2 + geom_segment(data = dat, aes(x = xmin, xend = xmax, y = middle, yend = middle), color = "white", inherit.aes = F)

invisible(NULL)

dat_exp <- res_model$dat_risk$TCGA_PRAD %>% tibble::remove_rownames() %>% column_to_rownames("ID") %>% dplyr::select("RS")

dat_drug <- Drug_list$PRISM$TCGA_PRAD %>% as.data.frame() %>% log1p()

ptypeOut <- dat_drug

data_meta <- dat_exp

top.pps <- data_meta[data_meta[, 1] >= quantile(data_meta[, 1], probs = seq(0, 1, 0.1))[10], , drop = F]

bot.pps <- data_meta[data_meta[, 1] <= quantile(data_meta[, 1], probs = seq(0, 1, 0.1))[2], , drop = F]

ccle.cor <- ccle.cor.p <- c()

com <- intersect(rownames(ptypeOut), rownames(data_meta))

data_meta <- data_meta[com, , drop = F]

ptypeOut <- ptypeOut[com, ]

for (i in 1:ncol(ptypeOut)) {
    d <- colnames(ptypeOut)[i]
    a <- as.numeric(ptypeOut[, d])
    b <- as.numeric(data_meta[, 1])
    r <- cor.test(a, b, method = "spearman")$estimate
    names(r) <- d
    p <- cor.test(a, b, method = "spearman")$p.value
    names(p) <- d
    ccle.cor <- c(ccle.cor, r)
    ccle.cor.p <- c(ccle.cor.p, p)
}

p.values <- as.numeric(ccle.cor.p)

names <- names(ccle.cor.p)

drug_data_1 <- data.frame(Drug_Name = names, p.value = p.values)

correlation <- as.numeric(ccle.cor)

names <- names(ccle.cor)

drug_data_2 <- data.frame(Drug_Name = names, correlation = correlation)

drug_data <- merge(drug_data_2, drug_data_1, by.x = 1, by.y = 1)

drug_data$state <- ifelse(drug_data$cor > 0, "positive", "negative")

drug <- drug_data %>% dplyr::group_by(state) %>% dplyr::filter(p.value < 0.05) %>% dplyr::arrange(desc(abs(correlation))) %>% 
    dplyr::slice_head(n = number) %>% dplyr::pull(Drug_Name)

darkblue <- "#0772B9"

lightblue <- "#48C8EF"

cor.data <- data.frame(drug = drug, r = ccle.cor[drug], p = -log10(ccle.cor.p[drug]))

cor.data$p[cor.data$p == Inf] <- 15

ggplot(data = cor.data, aes(r, forcats::fct_reorder(drug, r, .desc = T))) + geom_segment(aes(xend = 0, yend = drug), 
    linetype = 2) + geom_point(aes(size = p), col = darkblue) + scale_size_continuous(range = c(2, 8)) + scale_x_reverse(breaks = c(0, 
    -0.3, -0.5), expand = expansion(mult = c(0.01, 0.1))) + theme_classic() + labs(x = "Correlation coefficient", 
    y = "", size = bquote("-log"[10] ~ "(" ~ italic(P) ~ "-value)")) + theme(legend.position = "bottom", axis.line.y = element_blank())

invisible(NULL)

ccle.boxdata <- NULL

for (d in drug) {
    a <- as.numeric(ptypeOut[rownames(top.pps), d])
    b <- as.numeric(ptypeOut[rownames(bot.pps), d])
    p <- wilcox.test(a, b)$p.value
    s <- as.character(cut(p, c(0, 0.001, 0.01, 0.05, 1), labels = c("***", "**", "*", "")))
    ccle.boxdata <- rbind.data.frame(ccle.boxdata, data.frame(drug = d, auc = c(a, b), p = p, s = s, group = rep(c("High_signature", 
        "Low_signature"), c(nrow(top.pps), nrow(bot.pps))), stringsAsFactors = F), stringsAsFactors = F)
}

p2 <- ggplot(ccle.boxdata, aes(drug, auc, fill = group)) + geom_boxplot(aes(col = group), outlier.shape = NA) + 
    geom_text(aes(drug, y = max(auc)), label = ccle.boxdata$s, data = ccle.boxdata, inherit.aes = F) + scale_fill_manual(values = c(darkblue, 
    lightblue)) + scale_color_manual(values = c(darkblue, lightblue)) + xlab(NULL) + ylab("Estimated AUC value") + 
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 10), legend.position = "bottom", 
    legend.title = element_blank())

dat <- ggplot_build(p2)$data[[1]]

p2 + geom_segment(data = dat, aes(x = xmin, xend = xmax, y = middle, yend = middle), color = "white", inherit.aes = F)

invisible(NULL)

invisible(NULL)

