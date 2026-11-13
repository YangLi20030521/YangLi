# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
do_bulk_TimiGP <- function(dat_exp, pdata, dataset = "Charoentong2017", log = FALSE, GMNorm = FALSE, QV = FALSE) {
    library(TimiGP)
    if (dataset == "Bindea2013_Cancer") {
        prad_restore("input_055")
        geneset <- CellType_Bindea2013_cancer
    }
    else if (dataset == "Charoentong2017") {
        prad_restore("input_056")
        geneset <- CellType_Charoentong2017
    }
    else if (dataset == "Newman2015") {
        prad_restore("input_057")
        geneset <- CellType_Newman2015_LM22
    }
    else if (dataset == "Xu2018") {
        prad_restore("input_058")
        geneset <- CellType_Xu2018
    }
    else {
        stop("dataset must in one of Bindea2013_Cancer, Charoentong2017, Newman2015 or Xu2018")
    }
    marker <- unique(geneset$Gene)
    info <- TimiCheckEvent(pdata)
    com_sam <- intersect(rownames(info), colnames(dat_exp))
    info <- info[com_sam, ]
    dat_exp <- dat_exp[, com_sam]
    rna <- TimiPrePropress(marker = marker, rna = dat_exp, cohort = rownames(info), log = log, GMNorm = GMNorm)
    mps <- TimiGenePair(rna)
    res <- TimiCOX(mps = mps, info = info, p.adj = "BH", parallel = TRUE, core = 10)
    mps <- res[[1]]
    cox_res <- res[[2]]
    cell_pair <- TimiCellPair(geneset = geneset, core = 1, dataset = dataset)
    background <- TimiBG(marker.pair = row.names(cox_res))
    if (QV) {
        GP <- rownames(cox_res)[which(cox_res$QV < 0.05)]
    }
    else {
        GP <- rownames(cox_res)[which(cox_res$PV < 0.05)]
    }
    res <- TimiEnrich(gene = GP, background = background, geneset = cell_pair, p.adj = "BH", core = 1)
    score <- TimiFS(res)
    res_TimiGP <- list(res = res, score = score)
    return(res_TimiGP)
}

do_convert_h5ad <- function(seurat_obj, output_path) {
    anndata <- prad_input("input_059")
    np <- prad_input("input_060")
    pd <- prad_input("input_061")
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

do_anndata_seurat <- function(filename) {
    library(sceasy)
    library(reticulate)
    use_condaenv("sceasy")
    sceasy::convertFormat(obj = filename, from = "anndata", to = "seurat", outFile = "seurat.rds")
}

