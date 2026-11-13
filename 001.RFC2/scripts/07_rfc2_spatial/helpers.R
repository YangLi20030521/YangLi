# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
plot_sc_markers_heatmap <- function(object = NULL, markerGene = NULL, group.by = "ident", assays = "RNA", slot = "data", 
    layer = NULL, scale.data = TRUE, htCol = c("#0099CC", "white", "#CC0033"), htRange = c(-2, 0, 2), annoCol = FALSE, 
    myanCol = NULL, colseed = 666, annoColType = "light", annoColTypeAlpha = 0, row_title = "Cluster top Marker genes", 
    clusterAnnoName = TRUE, showRowNames = TRUE, row_names_side = "left", fontsize = 10, column_names_rot = 45, 
    markGenes = NULL, mark_fontsize = NULL, border = FALSE, width = NULL, height = NULL, cluster.order = NULL, 
    gene.order = NULL, cluster_columns = FALSE, cluster_rows = FALSE, show_legend = TRUE, legend_title = "Z-score", 
    return_data = FALSE, ...) {
    if (is.null(object)) {
        stop("Please provide a Seurat object!")
    }
    if (is.null(markerGene)) {
        stop("Please provide marker genes!")
    }
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' is required but not installed.")
    }
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        stop("Package 'ComplexHeatmap' is required but not installed.")
    }
    message("Calculating average expression...")
    seurat_version <- utils::compareVersion(as.character(utils::packageVersion("Seurat")), "5")
    if (is.null(layer)) {
        layer <- slot
    }
    if (seurat_version == 1) {
        mean_gene_exp <- as.matrix(data.frame(Seurat::AverageExpression(object, features = markerGene, group.by = group.by, 
            assays = assays, layer = layer)))
    }
    else {
        mean_gene_exp <- as.matrix(data.frame(Seurat::AverageExpression(object, features = markerGene, group.by = group.by, 
            assays = assays, slot = slot)))
    }
    colnames(mean_gene_exp) <- levels(Seurat::Idents(object))
    if (scale.data) {
        message("Scaling data (Z-score normalization)...")
        htdf <- t(scale(t(mean_gene_exp), scale = TRUE, center = TRUE))
    }
    else {
        htdf <- mean_gene_exp
    }
    if (!is.null(cluster.order)) {
        message("Applying custom cluster order...")
        htdf <- htdf[, cluster.order, drop = FALSE]
    }
    if (!is.null(gene.order)) {
        message("Applying custom gene order...")
        htdf <- htdf[gene.order, , drop = FALSE]
    }
    col_fun <- circlize::colorRamp2(htRange, htCol)
    if (annoCol == FALSE) {
        set.seed(colseed)
        anno_col <- circlize::rand_color(ncol(htdf), luminosity = annoColType, transparency = annoColTypeAlpha)
        message("Generated annotation colors:")
        print(setNames(anno_col, colnames(htdf)))
    }
    else if (annoCol == TRUE) {
        if (is.null(myanCol)) {
            stop("When annoCol=TRUE, please provide myanCol parameter!")
        }
        anno_col <- myanCol
    }
    else {
        stop("annoCol parameter must be TRUE or FALSE!")
    }
    names(anno_col) <- colnames(htdf)
    column_ha <- ComplexHeatmap::HeatmapAnnotation(cluster = colnames(htdf), show_legend = FALSE, show_annotation_name = clusterAnnoName, 
        col = list(cluster = anno_col))
    right_annotation <- NULL
    if (!is.null(markGenes)) {
        message("Adding gene markers...")
        rowGene <- rownames(htdf)
        annoGene <- markGenes
        index <- match(annoGene, rowGene)
        index <- index[!is.na(index)]
        if (length(index) > 0) {
            if (is.null(mark_fontsize)) {
                mark_fontsize <- fontsize
            }
            right_annotation <- ComplexHeatmap::rowAnnotation(gene = ComplexHeatmap::anno_mark(at = index, labels = rowGene[index], 
                labels_gp = grid::gpar(fontface = "italic", fontsize = mark_fontsize)), ...)
        }
        else {
            warning("None of the specified markGenes found in the data!")
        }
    }
    message("Generating heatmap...")
    ht_params <- list(matrix = htdf, show_row_dend = cluster_rows, show_column_dend = cluster_columns, name = legend_title, 
        cluster_columns = cluster_columns, cluster_rows = cluster_rows, row_title = row_title, right_annotation = right_annotation, 
        show_row_names = showRowNames, row_names_gp = grid::gpar(fontface = "italic", fontsize = fontsize), row_names_side = row_names_side, 
        border = border, column_names_side = "top", column_names_rot = column_names_rot, top_annotation = column_ha, 
        col = col_fun, show_heatmap_legend = show_legend)
    if (!is.null(width) && !is.null(height)) {
        ht_params$width <- ggplot2::unit(width, "cm")
        ht_params$height <- ggplot2::unit(height, "cm")
    }
    ht_params <- c(ht_params, list(...))
    ht <- do.call(ComplexHeatmap::Heatmap, ht_params)
    if (return_data) {
        message("Returning heatmap and data...")
        return(list(heatmap = ht, scaled_data = htdf, raw_data = mean_gene_exp, anno_colors = anno_col))
    }
    else {
        return(ht)
    }
}

extract_top_markers <- function(markers, n = 5, order_by = "avg_log2FC") {
    if (!requireNamespace("dplyr", quietly = TRUE)) {
        clusters <- unique(markers$cluster)
        top_genes <- character()
        for (cl in clusters) {
            cluster_markers <- markers[markers$cluster == cl, ]
            cluster_markers <- cluster_markers[order(-cluster_markers[[order_by]]), ]
            top_genes <- c(top_genes, head(cluster_markers$gene, n))
        }
        return(unique(top_genes))
    }
    else {
        markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = !!rlang::sym(order_by), n = n) %>% 
            dplyr::pull(gene) %>% unique()
    }
}

