# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
library(tidyverse)

prad_restore("input_066")

res_model <- prad_input("input_067")

dat_risk <- res_model$dat_risk$TCGA_PRAD %>% as.data.frame() %>% tibble::remove_rownames() %>% column_to_rownames("ID") %>% 
    dplyr::select("RS") %>% dplyr::mutate(Group = ifelse(RS > median(RS), "High_RS", "Low_RS"))

dat_cli <- total_clin_list$TCGA_PRAD %>% as.data.frame() %>% dplyr::select(c(PFS, PFS.time)) %>% na.omit() %>% 
    dplyr::rename(event = "PFS", time = "PFS.time")

dat_exp <- total_expr_list$TCGA_PRAD %>% as.data.frame()

com_sam <- Reduce(intersect, list(rownames(dat_risk), rownames(dat_cli), colnames(dat_exp)))

dat_risk <- dat_risk[com_sam, , drop = FALSE]

dat_exp <- dat_exp[, com_sam]

dat_cli <- dat_cli[com_sam, ]

high_sam <- rownames(subset(dat_risk, dat_risk$Group == "High_RS"))

dat_cli_high <- dat_cli[high_sam, ]

dat_exp_high <- dat_exp[, high_sam]

res_TimiGP_high <- do_bulk_TimiGP(dat_exp = dat_exp_high, pdata = dat_cli_high, dataset = "Charoentong2017", log = FALSE, 
    GMNorm = FALSE, QV = FALSE)

TimiDotplot(resdata = res_TimiGP_high$res, select = c(1:25), condition = "P.Value", cutoff = 0.05)

TimiCellChord(resdata = res_TimiGP_high$res, dataset = "Charoentong2017")

TimiGeneChord(resdata = res_TimiGP_high$res, select = 1)

TimiFSBar(score = res_TimiGP_high$score, select = 1:25)

invisible(NULL)

low_sam <- rownames(subset(dat_risk, dat_risk$Group == "Low_RS"))

dat_cli_low <- dat_cli[low_sam, ]

dat_exp_low <- dat_exp[, low_sam]

res_TimiGP_low <- do_bulk_TimiGP(dat_exp = dat_exp_low, pdata = dat_cli_low, dataset = "Charoentong2017", log = FALSE, 
    GMNorm = FALSE, QV = FALSE)

TimiDotplot(resdata = res_TimiGP_low$res, select = c(1:25), condition = "P.Value", cutoff = 0.05)

TimiCellChord(resdata = res_TimiGP_low$res, dataset = "Charoentong2017", condition = "P.Value", cutoff = 0.05)

TimiGeneChord(resdata = res_TimiGP_low$res, select = 1)

TimiFSBar(score = res_TimiGP_low$score, select = 1:9)

invisible(NULL)

