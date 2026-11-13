# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
if (TRUE) {
    library(tidyverse)
    library(MOVICS)
    library(purrr)
    invisible(NULL)
    prad_restore("input_016")
    gene <- prad_input("input_017")
    dat_all <- prad_input("input_018")
}

do_MOVICS_valid <- function(dataset, event, time, unit = "y") {
    sur_valid <- total_clin_list[[dataset]] %>% as.data.frame() %>% dplyr::select(c(ID, event, time)) %>% na.omit() %>% 
        tibble::remove_rownames() %>% column_to_rownames("ID")
    dat_valid <- total_expr_list[[dataset]] %>% as.data.frame() %>% na.omit()
    dat_valid <- dat_valid[, rownames(sur_valid)]
    NTP.pred <- runNTP(expr = dat_valid, templates = marker.up$templates, scaleFlag = FALSE, centerFlag = FALSE, 
        doPlot = TRUE, fig.name = paste0("NTP HEATMAP FOR ", dataset))
    colnames(sur_valid)[colnames(sur_valid) == event] <- "fustat"
    colnames(sur_valid)[colnames(sur_valid) == time] <- "futime"
    if (unit == "y") {
        sur_valid$futime <- round(sur_valid$futime * 365)
    }
    else if (unit == "m") {
        sur_valid$futime <- round(sur_valid$futime * 30.4)
    }
    else if (unit == "d") {
        sur_valid$futime <- sur_valid$futime
    }
    else {
        stop(print("unit must choose y or m or d"))
    }
    surv.prad <- compSurv(moic.res = NTP.pred, surv.info = sur_valid, convt.time = "y", surv.median.line = "h", 
        clust.col = c("#00599F", "#D01910"), xyrs.est = c(1, 3, 5), fig.name = paste0("output/", "KAPLAN-MEIER CURVE OF CONSENSUSMOIC for ", 
            dataset))
}

optk.stad <- getClustNum(data = mo.data_cleaned, is.binary = c(F, F, F, F, T), try.N.clust = 2:8, fig.name = "CLUSTER NUMBER OF TCGA-STAD")

moic.res.list <- getMOIC(data = mo.data_cleaned, methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", 
    "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster", "iClusterBayes"), N.clust = 2, type = c("gaussian", 
    "gaussian", "gaussian", "gaussian", "binomial"))

cmoic.prad <- getConsensusMOIC(moic.res.list = moic.res.list, fig.name = "CONSENSUS HEATMAP", distance = "euclidean", 
    linkage = "average", clust.col = c("#00599F", "#D01910"))

getSilhouette(sil = cmoic.prad$sil, fig.path = invisible(NULL), fig.name = "SILHOUETTE", height = 5.5, width = 5, 
    clust.col = c("#00599F", "#D01910"))

indata <- mo.data

indata$dat_methylation <- log2(indata$dat_methylation/(1 - indata$dat_methylation))

plotdata <- getStdiz(data = indata, halfwidth = c(2, 2, 2, 2, NA), centerFlag = c(T, T, T, T, F), scaleFlag = c(T, 
    T, T, T, F))

feat <- moic.res.list$iClusterBayes$feat.res

feat1 <- feat[which(feat$dataset == "dat_mrna"), ][1:10, "feature"]

feat2 <- feat[which(feat$dataset == "dat_lnc"), ][1:10, "feature"]

feat3 <- feat[which(feat$dataset == "dat_mirna"), ][1:10, "feature"]

feat4 <- feat[which(feat$dataset == "dat_methylation"), ][1:10, "feature"]

feat5 <- feat[which(feat$dataset == "dat_mutation"), ][1:10, "feature"]

annRow <- list(feat1, feat2, feat3, feat4, feat5)

mRNA.col <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")

lncRNA.col <- c("#6699CC", "white", "#FF3C38")

miRNA.col <- c("#74ADD1", "#FFFFBF", "#FDAE61", "#D73027")

meth.col <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")

mut.col <- c("grey90", "black")

col.list <- list(mRNA.col, lncRNA.col, miRNA.col, meth.col, mut.col)

plotdata$dat_mrna = plotdata$dat_mrna[apply(plotdata$dat_mrna, 1, function(x) sd(x) != 0), ]

plotdata$dat_lnc = plotdata$dat_lnc[apply(plotdata$dat_lnc, 1, function(x) sd(x) != 0), ]

plotdata$dat_mirna = plotdata$dat_mirna[apply(plotdata$dat_mirna, 1, function(x) sd(x) != 0), ]

plotdata$dat_methylation = plotdata$dat_methylation[apply(plotdata$dat_methylation, 1, function(x) sd(x) != 0), 
    ]

plotdata$dat_mutation = plotdata$dat_mutation[apply(plotdata$dat_mutation, 1, function(x) sd(x) != 0), ]

getMoHeatmap(data = plotdata, row.title = c("mRNA", "lncRNA", "miRNA", "Methylation", "Mutation"), is.binary = c(F, 
    F, F, F, T), legend.name = c("mRNA.FPKM", "lncRNA.FPKM", "miRNA.FPKM", "M value", "Mutated"), clust.res = cmoic.prad$clust.res, 
    clust.dend = NULL, show.rownames = c(F, F, F, F, F), show.colnames = FALSE, annRow = annRow, color = col.list, 
    clust.col = c("#00599F", "#D01910"), annCol = NULL, annColors = NULL, width = 10, height = 5, fig.name = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")

colnames(dat_all$dat_cli)[5:6] <- c("fustat", "futime")

dat_all$dat_cli$futime <- round(dat_all$dat_cli$futime * 365)

surv.prad <- compSurv(moic.res = cmoic.prad, surv.info = dat_all$dat_cli, convt.time = "y", surv.median.line = "h", 
    clust.col = c("#00599F", "#D01910"), xyrs.est = c(1, 3, 5), fig.name = "output/KAPLAN-MEIER CURVE OF CONSENSUSMOIC for TCGA-PRAD")

dat_exp <- total_expr_list$TCGA_PRAD

colnames(dat_exp) <- substr(colnames(dat_exp), 1, 15)

colnames(dat_exp) <- gsub("-", ".", colnames(dat_exp))

com <- intersect(colnames(dat_exp), rownames(cmoic.prad$clust.res))

dat_exp <- dat_exp[, rownames(cmoic.prad$clust.res)]

runDEA(dea.method = "limma", expr = dat_exp, moic.res = cmoic.prad, prefix = "TCGA-PRAD")

marker.up <- runMarker(moic.res = cmoic.prad, dea.method = "limma", prefix = "TCGA-PRAD", dirct = "up", n.marker = 100, 
    doplot = TRUE, norm.expr = dat_exp, fig.name = "output/DOWNREGULATED BIOMARKER HEATMAP")

do_MOVICS_valid(dataset = "DKFZ2018", event = "RFS", time = "RFS.time", unit = "y")

do_MOVICS_valid(dataset = "GSE107299", event = "RFS", time = "RFS.time", unit = "y")

do_MOVICS_valid(dataset = "GSE116918", event = "RFS", time = "RFS.time", unit = "y")

do_MOVICS_valid(dataset = "GSE21034", event = "RFS", time = "RFS.time", unit = "y")

do_MOVICS_valid(dataset = "GSE21034", event = "OS", time = "OS.time", unit = "y")

do_MOVICS_valid(dataset = "GSE54460", event = "RFS", time = "RFS.time", unit = "y")

do_MOVICS_valid(dataset = "GSE70768", event = "RFS", time = "RFS.time", unit = "y")

do_MOVICS_valid(dataset = "GSE70769", event = "RFS", time = "RFS.time", unit = "y")

res_cluster <- cmoic.prad$clust.res %>% as.data.frame()

dat_exp <- total_expr_list$TCGA_PRAD %>% as.data.frame()

colnames(dat_exp) <- gsub("-", ".", colnames(dat_exp))

colnames(dat_exp) <- substr(colnames(dat_exp), 1, 15)

com_sam <- intersect(colnames(dat_exp), rownames(res_cluster))

dat_exp <- dat_exp[, com_sam]

res_cluster <- res_cluster[com_sam, ]

dat_use <- do_bulk_gsva_ssgsea(dat_exp = dat_exp, geneset = gene$candidata, id = "RGDPs")

dat_use <- dat_use$res_gsva %>% as.data.frame() %>% dplyr::mutate(Cluster = res_cluster$clust)

dat_use$Cluster[dat_use$Cluster == 1] <- "CS1"

dat_use$Cluster[dat_use$Cluster == 2] <- "CS2"

plot_expr_diff_Violin(data = dat_use, Group = "Cluster", variable = "RGDPs", colors = c(CS1 = "#00599F", CS2 = "#D01910"), 
    method = "t.test")

invisible(NULL)

mut.PRAD <- compMut(moic.res = cmoic.prad, mut.matrix = dat_all$dat_mutation, doWord = TRUE, doPlot = TRUE, freq.cutoff = 0.05, 
    p.adj.cutoff = 0.05, innerclust = TRUE, width = 6, height = 3, clust.col = c("#00599F", "#D01910"), fig.name = "output/ONCOPRINT FOR SIGNIFICANT MUTATIONS", 
    tab.name = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")

MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)

gsea.up <- runGSEA_modify(moic.res = cmoic.prad, dea.method = "limma", prefix = "TCGA-PRAD", dat.path = "./data/", 
    res.path = "./output/", msigdb.path = MSIGDB.FILE, norm.expr = dat_exp, dirct = "up", n.path = 50, p.cutoff = 0.05, 
    p.adj.cutoff = 0.25, gsva.method = "gsva", norm.method = "mean", fig.name = "output/UPREGULATED PATHWAY HEATMAP")

invisible(NULL)

Subtype <- cmoic.prad$clust.res

SNF <- moic.res.list[["SNF"]]$clust.res

Subtype$SNF <- SNF$clust

PINSPlus <- moic.res.list[["PINSPlus"]]$clust.res

Subtype$PINSPlus <- PINSPlus$clust

NEMO <- moic.res.list[["NEMO"]]$clust.res

Subtype$NEMO <- NEMO$clust

COCA <- moic.res.list[["COCA"]]$clust.res

Subtype$COCA <- COCA$clust

LRAcluster <- moic.res.list[["LRAcluster"]]$clust.res

Subtype$LRAcluster <- LRAcluster$clust

ConsensusClustering <- moic.res.list[["ConsensusClustering"]]$clust.res

Subtype$ConsensusClustering <- ConsensusClustering$clust

IntNMF <- moic.res.list[["IntNMF"]]$clust.res

Subtype$IntNMF <- IntNMF$clust

CIMLR <- moic.res.list[["CIMLR"]]$clust.res

Subtype$CIMLR <- CIMLR$clust

MoCluster <- moic.res.list[["MoCluster"]]$clust.res

Subtype$MoCluster <- MoCluster$clust

iClusterBayes <- moic.res.list[["iClusterBayes"]]$clust.res

Subtype$iClusterBayes <- iClusterBayes$clust

invisible(NULL)

invisible(NULL)

