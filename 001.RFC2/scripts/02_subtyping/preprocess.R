# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
if (TRUE) {
    library(rtracklayer)
    library(data.table)
    library(tidyverse)
    library(survival)
    library(impute)
}

prad_restore("input_019")

pdata <- total_clin_list$TCGA_PRAD %>% as.data.frame() %>% remove_rownames() %>% column_to_rownames("ID")

rownames(pdata) <- substr(rownames(pdata), 1, 15)

run_cox_filter <- function(data, pdata, p_threshold = 0.05) {
    cox_results <- apply(data, 1, function(gene_expr) {
        coxph(Surv(PFS.time, PFS) ~ gene_expr, data = pdata)
    })
    p_values <- sapply(cox_results, function(x) summary(x)$coefficients[5]) %>% na.omit()
    significant_features <- names(p_values[p_values < p_threshold])
    return(significant_features)
}

match_samples <- function(data, pdata) {
    com <- intersect(colnames(data), rownames(pdata))
    list(data = data[, com], pdata = pdata[com, ])
}

miRNA <- prad_input("input_020") %>% column_to_rownames("sample")

na_ratio <- rowMeans(is.na(miRNA))

mirna_filtered <- miRNA[na_ratio < 0.3, ]

dat_mirna <- impute.knn(as.matrix(mirna_filtered))$data %>% as.data.frame()

matched <- match_samples(dat_mirna, pdata)

res_mirna <- run_cox_filter(matched$data, matched$pdata, p_threshold = 0.1)

dat_mirna <- matched$data

methylation <- prad_input("input_021") %>% column_to_rownames("sample")

na_ratio <- rowMeans(is.na(methylation))

dat_methylation <- methylation[na_ratio == 0, ]

matched <- match_samples(dat_methylation, pdata)

res_methylation <- run_cox_filter(matched$data, matched$pdata)

dat_methylation <- matched$data

exp <- prad_input("input_022")$dat_exp %>% as.data.frame()

gtf <- prad_input("input_023") %>% as.data.frame() %>% select(gene_id, gene_name, gene_type) %>% distinct(gene_name, 
    .keep_all = TRUE)

lncRNA_genes <- gtf %>% filter(gene_type == "lncRNA") %>% pull(gene_name)

dat_lnc <- exp[, intersect(colnames(exp), lncRNA_genes)] %>% na.omit() %>% t() %>% as.data.frame() %>% log1p()

pdata_lnc <- pdata

rownames(pdata_lnc) <- gsub("-", ".", rownames(pdata_lnc))

matched <- match_samples(dat_lnc, pdata_lnc)

res_lnc <- run_cox_filter(matched$data, matched$pdata)

dat_lnc <- matched$data

dat_mutation <- prad_input("input_024") %>% filter(!is.na(sample)) %>% column_to_rownames("sample") %>% mutate(across(everything(), 
    as.numeric))

matched <- match_samples(dat_mutation, pdata)

res_mutation <- run_cox_filter(matched$data, matched$pdata)

dat_mutation <- matched$data

colnames(dat_mirna) <- gsub("-", ".", colnames(dat_mirna))

colnames(dat_methylation) <- gsub("-", ".", colnames(dat_methylation))

colnames(dat_mutation) <- gsub("-", ".", colnames(dat_mutation))

common_samples <- Reduce(intersect, list(colnames(dat_lnc), colnames(dat_mirna), colnames(dat_methylation), colnames(dat_mutation)))

dat_MOVICS <- list(dat_lnc = dat_lnc[, common_samples], dat_mirna = dat_mirna[, common_samples], dat_methylation = dat_methylation[, 
    common_samples], dat_mutation = dat_mutation[, common_samples], res_lnc = res_lnc, res_mirna = res_mirna, res_methylation = res_methylation, 
    res_mutation = res_mutation)

invisible(NULL)

