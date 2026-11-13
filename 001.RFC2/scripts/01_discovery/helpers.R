# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
do_convert_h5ad <- function(seurat_obj, output_path) {
    anndata <- prad_input("input_001")
    np <- prad_input("input_002")
    pd <- prad_input("input_003")
    counts_matrix <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
    if (!inherits(counts_matrix, "matrix") && !inherits(counts_matrix, "dgCMatrix")) {
        counts_matrix <- as.matrix(counts_matrix)
    }
    counts_matrix <- Matrix::t(counts_matrix)
    counts_matrix <- as(counts_matrix, "dgCMatrix")
    meta_data <- seurat_obj@meta.data
    meta_data$barcode <- rownames(meta_data)
    adata <- anndata$AnnData(X = counts_matrix, obs = pd$DataFrame(meta_data))
    gene_names <- colnames(counts_matrix)
    adata$var_names <- np$array(gene_names)
    if (length(seurat_obj@reductions) > 0) {
        for (reduction_name in names(seurat_obj@reductions)) {
            embeddings <- Seurat::Embeddings(seurat_obj, reduction = reduction_name)
            obsm_key <- paste0("X_", reduction_name)
            adata$obsm[obsm_key] <- np$array(embeddings)
        }
    }
    invisible(NULL)
}

clinical_heatmap_plot <- function(data, score_var = "RMGs", cutoff_method = "median", output_file = "Clinical_Heatmap.pdf") {
    suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(tidyr)
        library(scales)
    })
    cutoff <- median(data[[score_var]], na.rm = TRUE)
    data$Group <- factor(ifelse(data[[score_var]] > cutoff, "High", "Low"), levels = c("Low", "High"))
    clinical_vars <- c("Age", "Race", "T", "N", "Gleason", "PSA", "Radiotherapy")
    clinical_vars <- clinical_vars[clinical_vars %in% colnames(data)]
    heatmap_data <- data.frame()
    chi_pvalues <- c()
    for (var in clinical_vars) {
        temp_data <- data[!is.na(data$Group) & !is.na(data[[var]]), ]
        if (nrow(temp_data) > 0 && length(unique(temp_data[[var]])) > 1) {
            perc_data <- temp_data %>% group_by(Group, !!sym(var)) %>% summarise(count = n(), .groups = "drop") %>% 
                group_by(Group) %>% mutate(percentage = count/sum(count) * 100) %>% ungroup() %>% mutate(Variable = var, 
                Category = as.character(!!sym(var))) %>% select(Variable, Category, Group, percentage)
            heatmap_data <- rbind(heatmap_data, perc_data)
            tryCatch({
                chi_test <- chisq.test(table(temp_data$Group, temp_data[[var]]))
                chi_pvalues <- c(chi_pvalues, chi_test$p.value)
            }, error = function(e) {
                chi_pvalues <<- c(chi_pvalues, NA)
            })
        }
    }
    pval_df <- data.frame(Variable = clinical_vars[1:length(chi_pvalues)], P_value = chi_pvalues, Sig = ifelse(chi_pvalues < 
        0.001, "***", ifelse(chi_pvalues < 0.01, "**", ifelse(chi_pvalues < 0.05, "*", "ns"))))
    heatmap_data <- heatmap_data %>% left_join(pval_df, by = "Variable") %>% mutate(Variable_label = paste0(Variable, 
        " (P", ifelse(P_value < 0.001, "<0.001", sprintf("=%.3f", P_value)), ")"))
    p <- ggplot(heatmap_data, aes(x = Group, y = Category, fill = percentage)) + geom_tile(color = "white", size = 1) + 
        geom_text(aes(label = sprintf("%.1f%%", percentage)), color = "black", size = 3.5, fontface = "bold") + 
        facet_grid(Variable ~ ., scales = "free_y", space = "free_y", switch = "y") + scale_fill_gradient2(low = "#3288BD", 
        mid = "white", high = "#D53E4F", midpoint = 50, limits = c(0, 100), name = "Percentage\n(%)") + labs(title = "", 
        x = paste0(score_var, " Group"), y = NULL) + theme_minimal() + theme(text = element_text(size = 12, family = "sans"), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14), axis.title.x = element_text(face = "bold", 
            size = 12), axis.text.x = element_text(color = "black", size = 11, face = "bold"), axis.text.y = element_text(color = "black", 
            size = 10, face = "bold"), strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 10), 
        strip.background = element_rect(fill = "grey90", color = "black", size = 0.5), legend.title = element_text(face = "bold", 
            size = 10), legend.text = element_text(size = 9), panel.grid = element_blank(), panel.spacing = unit(0.5, 
            "lines"))
    invisible(NULL)
    return(list(plot = p, chi_square_results = pval_df, heatmap_data = heatmap_data))
}

do_anndata_seurat <- function(filename) {
    library(sceasy)
    library(reticulate)
    use_condaenv("sceasy")
    sceasy::convertFormat(obj = filename, from = "anndata", to = "seurat", outFile = "seurat.rds")
}

do_sc_Read_TXT_CSV <- function(data_dir) {
    suppressPackageStartupMessages({
        library(data.table)
        library(tidyverse)
        library(Seurat)
    })
    file_paths <- prad_input("input_004")
    seurat_objects <- list()
    for (file_path in file_paths) {
        sample_id <- strsplit(basename(file_path), "_")[[1]][1]
        message("Processing sample: ", sample_id)
        expr_data <- prad_input("input_005") %>% as.data.frame() %>% distinct(Symbol, .keep_all = TRUE) %>% column_to_rownames("Symbol") %>% 
            .[, !duplicated(colnames(.))] %>% dplyr::select(-Gene_ID)
        if (ncol(expr_data) > nrow(expr_data)) {
            expr_data <- t(expr_data)
        }
        seurat_obj <- CreateSeuratObject(counts = expr_data, project = sample_id, assay = "RNA")
        seurat_objects[[sample_id]] <- seurat_obj
    }
    if (length(seurat_objects) == 1) {
        merged_seurat <- seurat_objects[[1]]
    }
    else {
        merged_seurat <- merge(seurat_objects[[1]], seurat_objects[-1])
    }
    if (packageVersion("Seurat") >= "5.0.0") {
        merged_seurat <- JoinLayers(merged_seurat)
    }
    message("Successfully processed", length(seurat_objects), "samples")
    return(merged_seurat)
}

runGSEA_modify <- function(moic.res = NULL, dea.method = c("deseq2", "edger", "limma"), norm.expr = NULL, prefix = NULL, 
    dat.path = NULL, res.path = NULL, dirct = "up", n.path = 10, msigdb.path = NULL, nPerm = 1000, minGSSize = 10, 
    maxGSSize = 500, p.cutoff = 0.05, p.adj.cutoff = 0.05, gsva.method = "gsva", norm.method = "mean", clust.col = c("#2EC4B6", 
        "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627", "#023E8A", "#9D4EDD"), color = NULL, fig.name = NULL, 
    fig.path = NULL, width = 15, height = 10) {
    comsam <- intersect(moic.res$clust.res$samID, colnames(norm.expr))
    if (length(comsam) == nrow(moic.res$clust.res)) {
        message("--all samples matched.")
    }
    else {
        message(paste0("--", (nrow(moic.res$clust.res) - length(comsam)), " samples mismatched from current subtypes."))
    }
    moic.res$clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
    norm.expr <- norm.expr[, comsam]
    n.moic <- length(unique(moic.res$clust.res$clust))
    mo.method <- moic.res$mo.method
    DEpattern <- paste(mo.method, "_", ifelse(is.null(prefix), "", paste0(prefix, "_")), dea.method, ".*._vs_Others.txt$", 
        sep = "")
    DEfiles <- dir(dat.path, pattern = DEpattern)
    if (length(DEfiles) == 0) {
        stop("no DEfiles!")
    }
    if (length(DEfiles) != n.moic) {
        stop("not all multi-omics clusters have DEfile!")
    }
    if (!is.element(dirct, c("up", "down"))) {
        stop("dirct type error! Allowed value contains c('up', 'down').")
    }
    if (!is.element(gsva.method, c("gsva", "ssgsea", "zscore", "plage"))) {
        stop("GSVA method error! Allowed value contains c('gsva', 'ssgsea', 'zscore', 'plage').")
    }
    if (dirct == "up") {
        outlabel <- "unique_upexpr_pathway.txt"
    }
    if (dirct == "down") {
        outlabel <- "unique_downexpr_pathway.txt"
    }
    gsea.list <- list()
    gseaidList <- c()
    for (filek in DEfiles) {
        DEres <- prad_input("input_006")
        DEres <- DEres[!duplicated(DEres[, 1]), ]
        DEres <- DEres[!is.na(DEres[, 1]), ]
        rownames(DEres) <- DEres[, 1]
        DEres <- DEres[, -1]
        geneList <- DEres$log2fc
        names(geneList) <- rownames(DEres)
        geneList <- sort(geneList, decreasing = TRUE)
        msigdb <- try(prad_input("input_007"), silent = TRUE)
        if (class(msigdb) == "try-error") {
            stop("please provide correct ABSOLUTE PATH for MSigDB file.")
        }
        moic.lab <- paste0("CS", gsub("_vs_Others.txt", "", sub(".*CS", "", filek)))
        gsea.list[[moic.lab]] <- suppressWarnings(clusterProfiler::GSEA(geneList = geneList, TERM2GENE = msigdb, 
            nPerm = nPerm, minGSSize = minGSSize, maxGSSize = maxGSSize, seed = TRUE, verbose = FALSE, pvalueCutoff = 1))
        gsea.dat <- as.data.frame(gsea.list[[moic.lab]])
        invisible(NULL)
        if (dirct == "up") {
            gseaidList <- c(gseaidList, rownames(gsea.dat[which(gsea.dat$NES > 0 & gsea.dat$pvalue < p.cutoff & 
                gsea.dat$p.adjust < p.adj.cutoff), ]))
        }
        if (dirct == "down") {
            gseaidList <- c(gseaidList, rownames(gsea.dat[which(gsea.dat$NES < 0 & gsea.dat$pvalue < p.cutoff & 
                gsea.dat$p.adjust < p.adj.cutoff), ]))
        }
        unqlist <- setdiff(gseaidList, gseaidList[duplicated(gseaidList)])
    }
    message("GSEA done...")
    GSEApattern <- paste(mo.method, "_", ifelse(is.null(prefix), "", paste0(prefix, "_")), dea.method, ".*._gsea_all_results.txt$", 
        sep = "")
    GSEAfiles <- dir(res.path, pattern = GSEApattern)
    pathway <- pathcore <- list()
    pathnum <- c()
    for (filek in GSEAfiles) {
        GSEAres <- prad_input("input_008")
        if (dirct == "up") {
            outk <- intersect(unqlist, rownames(GSEAres[which(GSEAres$NES > 0 & GSEAres$pvalue < p.cutoff & GSEAres$p.adjust < 
                p.adj.cutoff), ]))
            outk <- GSEAres[outk, ]
            outk <- outk[order(outk$NES, decreasing = TRUE), ]
            if (nrow(outk) > n.path) {
                pathway[[filek]] <- outk[1:n.path, ]
            }
            else {
                pathway[[filek]] <- outk
            }
            pathnum <- c(pathnum, nrow(pathway[[filek]]))
            pathway$dirct <- "up"
            for (i in rownames(pathway[[filek]])) {
                pathcore[[i]] <- msigdb[which(msigdb[, 1] %in% i), "gene"]
            }
        }
        if (dirct == "down") {
            outk <- intersect(unqlist, rownames(GSEAres[which(GSEAres$NES < 0 & GSEAres$pvalue < p.cutoff & GSEAres$p.adjust < 
                p.adj.cutoff), ]))
            outk <- GSEAres[outk, ]
            outk <- outk[order(outk$NES, decreasing = FALSE), ]
            if (nrow(outk) > n.path) {
                pathway[[filek]] <- outk[1:n.path, ]
            }
            else {
                pathway[[filek]] <- outk
            }
            pathnum <- c(pathnum, nrow(pathway[[filek]]))
            pathway$dirct <- "down"
            for (i in rownames(pathway[[filek]])) {
                pathcore[[i]] <- msigdb[which(msigdb[, 1] %in% i), "gene"]
            }
        }
        invisible(NULL)
    }
    standarize.fun <- function(indata = NULL, halfwidth = NULL, centerFlag = TRUE, scaleFlag = TRUE) {
        outdata = t(scale(t(indata), center = centerFlag, scale = scaleFlag))
        if (!is.null(halfwidth)) {
            outdata[outdata > halfwidth] = halfwidth
            outdata[outdata < (-halfwidth)] = -halfwidth
        }
        return(outdata)
    }
    rowmean <- function(x) {
        return(apply(x, 1, mean))
    }
    rowmedian <- function(x) {
        return(apply(x, 1, median))
    }
    if (max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 0)) {
        message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
        gset <- norm.expr
    }
    if (max(norm.expr) >= 25 & min(norm.expr) >= 0) {
        message("--log2 transformation done for expression data.")
        gset <- log2(norm.expr + 1)
    }
    sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, decreasing = FALSE), "samID"]
    colvec <- clust.col[1:n.moic]
    names(colvec) <- paste0("CS", 1:n.moic)
    annCol <- data.frame(Subtype = paste0("CS", moic.res$clust.res[sam.order, "clust"]), row.names = sam.order, 
        stringsAsFactors = FALSE)
    annColors <- list(Subtype = colvec)
    if (gsva.method == "gsva") {
        gsvaParamObj <- GSVA::gsvaParam(exprData = as.matrix(gset[, rownames(annCol), drop = FALSE]), geneSets = pathcore, 
            minSize = 1)
    }
    else if (gsva.method == "ssgsea") {
        gsvaParamObj <- GSVA::ssgseaParam(exprData = as.matrix(gset[, rownames(annCol), drop = FALSE]), geneSets = pathcore, 
            minSize = 1)
    }
    else if (gsva.method == "zscore") {
        gsvaParamObj <- GSVA::zscoreParam(exprData = as.matrix(gset[, rownames(annCol), drop = FALSE]), geneSets = pathcore, 
            minSize = 1)
    }
    else if (gsva.method == "plage") {
        gsvaParamObj <- GSVA::plageParam(exprData = as.matrix(gset[, rownames(annCol), drop = FALSE]), geneSets = pathcore, 
            minSize = 1)
    }
    es <- GSVA::gsva(gsvaParamObj)
    es.backup <- es
    es <- standarize.fun(es, halfwidth = 1, centerFlag = TRUE, scaleFlag = TRUE)
    message(gsva.method, " done...")
    esm <- data.frame(row.names = rownames(es))
    if (norm.method == "mean") {
        for (i in paste0("CS", 1:n.moic)) {
            esm <- cbind.data.frame(esm, data.frame(rowmean(es[, rownames(annCol[which(annCol$Subtype == i), , 
                drop = FALSE])])))
        }
    }
    if (norm.method == "median") {
        for (i in paste0("CS", 1:n.moic)) {
            esm <- cbind.data.frame(esm, data.frame(rowmedian(es[, rownames(annCol[which(annCol$Subtype == i), 
                , drop = FALSE])])))
        }
    }
    colnames(esm) <- paste0("CS", 1:n.moic)
    annRow <- data.frame(Subtype = rep(paste0("CS", 1:n.moic), pathnum), row.names = rownames(esm), stringsAsFactors = FALSE)
    if (is.null(color)) {
        mapcolor <- (grDevices::colorRampPalette(c("#0000FF", "#8080FF", "#FFFFFF", "#FF8080", "#FF0000")))(64)
    }
    else {
        mapcolor <- (grDevices::colorRampPalette(color))(64)
    }
    hm <- ComplexHeatmap::pheatmap(mat = as.matrix(esm), cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, 
        show_colnames = TRUE, annotation_row = annRow, annotation_colors = annColors, annotation_names_row = FALSE, 
        legend = TRUE, color = mapcolor, border_color = "black", legend_breaks = c(-1, -0.5, 0, 0.5, 1), legend_labels = c(-1, 
            -0.5, 0, 0.5, 1), cellwidth = 15, cellheight = 10)
    draw(hm, annotation_legend_side = "left", heatmap_legend_side = "left")
    if (is.null(fig.name)) {
        outFig <- paste0("gseaheatmap_using_", dirct, "regulated_pathways.pdf")
    }
    else {
        outFig <- paste0(fig.name, "_using_", dirct, "regulated_pathways.pdf")
    }
    invisible(NULL)
    draw(hm, annotation_legend_side = "left", heatmap_legend_side = "left")
    invisible(invisible(NULL))
    message("heatmap done...")
    return(list(gsea.list = gsea.list, raw.es = es.backup, scaled.es = es, grouped.es = esm, heatmap = hm))
}

