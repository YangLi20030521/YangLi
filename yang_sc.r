do_h5ad_seurat <- function(input,output){
  library(sceasy)
  library(reticulate)
  library(Seurat)
  library(BiocParallel)
  register(MulticoreParam(workers = 4, progressbar = TRUE)) 
  
  sceasy::convertFormat(obj = input, 
                        from="anndata",
                        to="seurat",
                        outFile = output)
}


do_sc_read10X <- function(dir,version="v5"){
  library(Seurat)
  samples_name=as.vector(list.files(dir))
  samples_name
  sampes_name_path = paste0(dir,"/",samples_name)
  sampes_name_path
  scRNAlist <- list()
  for(i in sampes_name_path){
    project <- gsub(paste0(dir,"/"),"",i)
    counts <- Read10X(data.dir = i)
    scRNAlist[[i]] <- CreateSeuratObject(counts, 
                                         project = project,
                                         min.cells = 3, 
                                         min.features = 200)
    scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = project)
  }
  names(scRNAlist) <- samples_name
  scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
  scRNA <- JoinLayers(scRNA)
  if(version=="v5"){
    return(scRNA)
  }else{
    scRNA[["RNA"]] <- as(scRNA[["RNA"]], "Assay")
    return(scRNA)
  }
}


do_sc_qc <- function(scRNA,
                     pctMT=10,
                     minGene=500,
                     maxGene=4000,
                     maxUMI=15000){
  scRNA[["percent.MT"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
  scRNA <- subset(scRNA, subset = percent.MT < pctMT &
                    nCount_RNA < maxUMI &
                    nFeature_RNA > minGene &
                    nFeature_RNA < maxGene)
}


do_sc_pipeline <- function(scRNA,method = "LogNormalize",resolution = 0.1,integration = NULL){
  if(method=="LogNormalize"){
    scRNA <- NormalizeData(object = scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
    scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = 2000)
    scRNA <- ScaleData(scRNA)
    scRNA <- RunPCA(scRNA, npcs=20,pc.genes=VariableFeatures(object = scRNA), verbose=FALSE)
    pct <- scRNA[["pca"]]@stdev/sum(scRNA[["pca"]]@stdev) * 100
    cumu <- cumsum(pct)
    co1 <- which(cumu > 90 & pct < 5)[1]
    co1
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    co2
    pcSelect <- min(co1, co2)
    print(pcSelect)
    if(is.null(integration)){
      scRNA <- FindNeighbors(scRNA, reduction="pca",dims = 1:pcSelect)
      scRNA <- FindClusters(object = scRNA,reduction="pca",resolution = resolution)
      scRNA <- RunUMAP(object = scRNA, reduction="pca",dims = 1:pcSelect,check_duplicates = FALSE)
      return(scRNA)
    }else{
      library(harmony)
      scRNA <- RunHarmony(object=scRNA, group.by.vars="orig.ident",reduction.use = "pca",assay.use="RNA",max.iter.harmony = 20)
      scRNA <- FindNeighbors(scRNA, reduction="harmony",dims = 1:pcSelect)
      scRNA <- FindClusters(object = scRNA,reduction="harmony",resolution = resolution)
      scRNA <- RunUMAP(object = scRNA, reduction="harmony",dims = 1:pcSelect,check_duplicates = FALSE)
      return(scRNA)
    }
  }else if(method=="SCT"){
    scRNA <- SCTransform(scRNA)
    scRNA <- RunPCA(scRNA,assay = "SCT", npcs=20, pc.genes=VariableFeatures(object = scRNA), verbose=FALSE)
    pct <- scRNA[["pca"]]@stdev/sum(scRNA[["pca"]]@stdev) * 100
    cumu <- cumsum(pct)
    co1 <- which(cumu > 90 & pct < 5)[1]
    co1
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    co2
    pcSelect <- min(co1, co2)
    print(pcSelect)
    if(is.null(integration)){
      scRNA <- FindNeighbors(scRNA, reduction="pca",dims = 1:pcSelect)
      scRNA <- FindClusters(object = scRNA,reduction="pca",resolution = resolution)
      scRNA <- RunUMAP(object = scRNA, reduction="pca",dims = 1:pcSelect,check_duplicates = FALSE)
      return(scRNA)
    }else{
      library(harmony)
      scRNA <- RunHarmony(object=scRNA, group.by.vars="orig.ident",reduction.use = "pca",assay.use="RNA",max.iter.harmony = 20)
      scRNA <- FindNeighbors(scRNA, reduction="harmony",dims = 1:pcSelect)
      scRNA <- FindClusters(object = scRNA,reduction="harmony",resolution = resolution)
      scRNA <- RunUMAP(object = scRNA, reduction="harmony",dims = 1:pcSelect,check_duplicates = FALSE)
      return(scRNA)
    }
  }
}


get_sc_markers = function(label = "main_lable"){
  if(label == "main_lable"){
    message(paste0(
      "Epithelial_cells = c('EPCAM','KRT19','CDH1','KRT18') #PMID: 28474673/31067475\n",
      "T_lymphocytes = c('CD3D','CD3E','CD3G','TRAC') #PMID: 28475900/31209336\n",
      "B_lymphocytes = c('CD79A','IGHM','IGHG3','IGHA2') #PMID: 31712411/30523328\n",
      "Myeloid_cells = c('CD68',' MARCO','FCGR3A','LYZ') #PMID: 28475900/29967419\n",
      "NK_cells = c('NCAM1','NKG7','GNLY','KLRD1') #PMID: 28475900/31477722\n",
      "Mast_cells = c('KIT','MS4A2','GATA2') #PMID: 30979687\n",
      "Fibroblasts = c('DCN','COL1A1','COL1A2','THY1') #PMID: 31209336/29198524\n",
      "Endothelial_cells = c('PECAM1','CLDN5','FLT1','RAMP2') #PMID: 30674341/21460247/23355623\n",
      "Oligodendrocytes = c('OLIG1','OLIG2','MOG','CLDN11') #PMID: 29615592/26628089\n")
    )
  }else if(label == "fine_lable_t"){
    message(paste0(
      "NK = c('XCL1','FCGR3A','KLRD1','KLRF1') #PMID: 31477722\n",
      "T_lymphocytes = c('CD3D','CD3E','CD3G','TRAC') #PMID: 28475900/31209336\n",
      "CD4T = c('Il7R','CD4') #PMID: 14662907/28475900\n",
      "CD8T = c('CD8A','CD8B') #PMID: 28475900\n",
      "Naïve = c('TCF7','SELL','LEF1','CCR7') #PMID: 29942094\n",
      "Treg = c('IL2RA','FOXP3','IKZF2','TGFB1','TGFB3','TGFBI','TGFBR1') #PMID: 29942094/28474673\n",
      "Tex = c('LAG3','TIGIT','PDCD1','HAVCR2') #PMID: 29942094\n",
      "Teff = c('IL2','GZMA','GNLY','PRF1','GZMB','GZMK','IFNG','NKG7') #PMID: 29942094\n",
      "Th1 = c('STAT4','IFNG','IL12RB2') #PMID: 24987392/21685955\n",
      "Th2 = c('GATA3', 'STAT6','IL4') #PMID: 24987392\n",
      "Th17 = c('IRF4','CREM','NR4A2') #PMID: 21381156/27680869/23437182\n",
      "Tgd = c('TRDC','TRGC2','TRGC1') #PMID: 31118283\n",
      "Tfh = c('MAF','CXCR5','PDCD1','CXCL13') #PMID: 28265271/28570278\n")
    )
  }else if(label == "fine_lable_b"){
    message(paste0(
      "B_lymphocytes = c('CD79A','IGHM','IGHG3','IGHA2') ##PMID: 31712411/30523328\n",
      "B_GC_DZ = c('STMN1','AICDA','MKI67','BIRC5') #PMID: 30104629\n",
      "B_GC_LZ = c('LMO2','BCL2A1') #PMID: 30104629\n",
      "GrB_secreting = c('GZMB') #PMID: 21808264\n",
      "B_follicular = c('MS4A1','HLA-DRA') #PMID: 29988129\n",
      "B_MALT = c('JCHAIN','IGHA1') #PMID: 29988129\n",
      "Plasma = c('IGHG1') #PMID: 29988129\n"
    )
    )
  }else if(label == "fine_lable_dc"){
    message(paste0(
      "DCs = c('CLEC10A','CD1C','CLEC4C','PTCRA','CCR7','LAMP3') #PMID: 28475900\n",
      "DCs_CD1c = c('CD1C','ITGAX') #PMID: 24744755\n",
      "DCs_CD141 = c('CLEC9A','XCR1') #PMID: 24744755\n",
      "DCs_CD207_CD1a = c('CD207','CD1A') #PMID: 24744755\n",
      "DCs_Activated = c('CCR7','LAMP3') #PMID: 17312119\n",
      "pDCs = c('IL3RA','CLEC4C') #PMID: 28428369\n",
      "DCs_CD163_CD14 = c('CD14','CD163') #PMID: 31474513\n"))
  }else if(label == "fine_lable_myeloid"){
    message(paste0(
      "Monocyte = c('CTSS','FCN1','S100A8','S100A9','LYZ','VCAN') #PMID: 29967419\n",
      "Macrophage = c('LGMN','CTSB','CD14','FCGR3A') #PMID: 29967419\n",
      "mo_lineage = c('MAFB','MAF','CX3CR1','ITGAM','CSF1R') #PMID: 28257233\n",
      "Alveolar_Mac = c('MARCO','FABP4','MCEMP1') #PMID: 28257233\n",
      "Anti_inflammatory = c('CD163','APOE','SEPP1','C1QA','C1QB','C1QC') #PMID: 27381735/21350196/26053663/22523386\n",
      "Pro_inflammatory = c('CXCL8', 'IL1B') #PMID: 25339958\n",
      "Cycling = c('STMN1','MKI67','TOP2A','CDK1') #PMID: 29967419\n",
      "DC = c('CLEC10A','CD1C', 'CLEC4C','PTCRA','CCR7','LAMP3') #PMID: 28475900\n"
    ))
  }else if(label == "fine_lable_fibroblast"){
    message(paste0(
      "Fibroblasts = c('DCN','COL1A1','COL1A2','THY1') #PMID: 31209336/29198524\n",
      "FBs_COL13A1 = c('COL13A1','TCF21','ITGAB','CXCL14','NPNT') #PMID: 29590628\n",
      "FBs_COL14A1 = c('COL14A1','GSN', 'PI16','CYGB','PRRX1') #PMID: 29590628\n",
      "Myofibroblasts = c('ACTA2','MYH11','TAGLN', 'ACTG2','MYLK') #PMID: 29590628\n",
      "SMCs = c('CNN1','SYNPO2','CRYAB','DES') #PMID: 28564607\n",
      "Mesothelial_cells = c('UPK3B','MSLN','CALB2','WT1') #PMID: 29590628\n",
      "Pericytes = c('RGS5','CSPG4','ABCC9','KCNJ8') #PMID: 28564607\n",
      "Perivascular_FBs = c('CYP1B1','APOD') #PMID: 29443965\n",
      "Lipofibroblast = c('FABP4','FABP5','PPARG') #PMID: 29590628\n"))
  }else if(label == "fine_lable_endothelial"){
    message(paste0(
      "Endothelial_cells = c('PECAM1','CLDN5','FLT1','RAMP2') #PMID: 30674341/21460247/23355623\n",
      "Tip_like_ECs = c('RAMP3','RGCC','ADM') #PMID: 29449267\n",
      "Stalk_like_ECs = c('SELP', 'ACKR1') #PMID: 29449267\n",
      "Lymphatic_ECs = c('CCL21', 'LYVE1') #PMID: 29449267\n",
      "EPCs = c('TYROBP','C1QB') #PMID: 29449267\n",
      "Tumor_ECs = c('HSPG2','INSR','VWA1') #PMID: 29988129/30559346/10629090\n"
    ))
  }else if(label == "fine_lable_epithelial"){
    message(paste0(
      "Epithelial_cells = c('EPCAM','KRT19','CDH1','KRT18') #PMID: 28474673/31067475\n",
      "AT1 = c('AGER') #PMID: 30554520\n",
      "AT2 = c('SFTPC','LAMP3') #PMID: 30554520\n",
      "Club = c('SCGB1A1') #PMID: 30554520\n",
      "Ciliated = c('FOXJ1','RFX2') #PMID: 30554520\n"
    ))
  }
}


plot_ann_dotplot_easy_human <- function(scRNA,
                                        markers_list = list(Epith=c('EPCAM','CDH1','KRT5','AMACR','KRT8','KRT17'),
                                                            Meso = c("MSLN"),
                                                            Fibro=c("COL1A2","FGF7","DCN"),
                                                            Endo=c("TEK","PECAM1","FLT1","VWF"),
                                                            Macro= c("CD68","LYZ"),
                                                            Mono = c("CD14","FCGR3A"),
                                                            Mast=c("TPSAB1","TPSB2"),
                                                            B= c("CD79A","CD79B","CD19","MS4A1"),
                                                            Plasma=c("JCHAIN","XBP1"),
                                                            NK= c("KLRD1","NKG7","NCAM1",'XCL2','XCL1'),
                                                            T=c("CD3G","CD3D","CD3E"),
                                                            CD8=c("CD8A","CD8B"),
                                                            CD4=c("CD4"),
                                                            Granulocyte=c('LAMP5','CD147'),
                                                            cDC1=c('CD11c','HLA-DR','CD141'),
                                                            cDC2=c('CD1c','CD11b','FCER1A'),
                                                            pDC=c('CD303','CD123'),
                                                            moDC=c('CD11b','CD64','CD88')),
                                        Labels="seurat_clusters"){
  library(Seurat)
  library(ggplot2)
  
  plot <- DotPlot(scRNA,
                  features = markers_list,
                  group.by = Labels) +
    scale_colour_gradient2(low = c('lightgrey','#330066'),
                           mid = '#336699',
                           high =c('#66CC66','#FFCC33',"red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 设置角度为45度，并调整水平对齐方式
  return(plot)
}


plot_ann_dotplot_easy_mouse <- function(scRNA,
                                        markers_list = list(Cardiomyocytes = c('Myh6','Tnnt2','Actc1','Ckmt2','Nppa','Actn2','Tnni3'),
                                                            Endothelial_cell = c('Pecam1','Vwf','Cdh5','Nos3','Vegfa'),
                                                            Fibroblast_cell = c('Vim','Col1a1','Col3a1','Ckap4','Ddr2','Tcf21','Pdgfra'),
                                                            Macrophage = c('Cd68','Adgre1','Fcgr1','Lgals3','Itgam'),
                                                            Epicardial_cell = c('Upk3b','Upk1b','Muc16','Krt19','Krt18','Krt14','Wt1'),
                                                            T_cell = c('Cd3e','Cd3d','Il7r','Cd3g','Lat'),
                                                            Adipocytes = c('Pparg','Cebpa','Lepr','Adipoq'),
                                                            Pericytes = c('Pdgfrb','Nes','Cspg4','Gja1','Adamts1'),
                                                            B_cell = c('Cd79a','Cd79b','Ighm','Cd19','Igha2'),
                                                            # Myeloid_cell = c('Cd68','Csf1r','C3ar1','Adgre1','Itgam'),
                                                            Granulocytes = c('S100a9','Csf3r','Hcar2','Slpi'),
                                                            Smooth_muscle_cell = c('Cnn1','Myh11','Mylk','Mustn1','Tagln'),
                                                            Mast_cell = c('Fcer1g','Ptprc','Kit','Tpsab1')),
                                        Labels="seurat_clusters"){
  library(Seurat)
  library(ggplot2)
  
  plot <- DotPlot(scRNA,
                  features = markers_list,
                  group.by = Labels) +
    scale_colour_gradient2(low = c('lightgrey','#330066'),
                           mid = '#336699',
                           high =c('#66CC66','#FFCC33',"red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 设置角度为45度，并调整水平对齐方式
  return(plot)
}


plot_sc_bar <- function(obj,group,sample){
  library(reshape2)
  library(tidyr)
  library(Seurat)
  library(gplots)
  library(ggthemes)
  library(RColorBrewer)
  
  tb <- table(obj@meta.data[,sample], obj@meta.data[,group])
  bar_data <- as.data.frame(tb)
  bar_per <- bar_data %>% 
    group_by(Var1) %>%
    mutate(sum(Freq)) %>%
    mutate(percent = Freq / `sum(Freq)`)
  # col <- c("#CE2820","#F78000","#3FA116","#3176B7","#9265C1","#885649","#DD76C5","#BBBE00","#41BED1")
  col1 <- sample(RColorBrewer::brewer.pal(9,'Set1'))
  col2 <- sample(RColorBrewer::brewer.pal(8,'Set2'))
  col3 <- sample(RColorBrewer::brewer.pal(12,'Set3'))
  col4 <- c('skyblue','red')
  col <- c(col1,col2,col3,col4)
  ggplot(bar_per, aes(y = percent, x = Var1)) +
    geom_bar(aes(fill = Var2) , stat = "identity") + coord_flip() +
    theme(axis.ticks = element_line(linetype = "blank"),
          legend.position = "top",
          panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
          panel.background = element_rect(fill = NA),
          plot.background = element_rect(colour = NA)) +
    labs(y = " ", fill = NULL)+labs(x = 'Relative proportion(%)')+
    scale_fill_manual(values=col)+
    theme_few()+
    theme(plot.title = element_text(size=12,hjust=0.5)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename = 'barplot.pdf',width = 8,height = 5)
  
  pdf(file = 'tb.pdf',width = 10,height = 8)
  balloonplot(tb)
  dev.off()
}


plot_sc_bar_plus <- function(obj = sc,sample = 'sample',celltype = 'CellType'){
  library(reshape2)
  library(ggplot2)
  library(ggalluvial)
  library(ggh4x)
  library(Seurat)
  
  colors=c("#FBB463","#80B1D3","#F47F72","#BDBADB","#FBF8B4","#8DD1C6","#CE2820","#F78000","#3FA116","#3176B7","#9265C1","#885649","#DD76C5","#BBBE00","#41BED1")
  tb <- table(obj@meta.data[,sample], obj@meta.data[,celltype])
  bar_data <- as.data.frame(tb)
  bar_per <- bar_data %>% 
    group_by(Var1) %>%
    mutate(sum(Freq)) %>%
    mutate(percent = Freq / `sum(Freq)`)
  colnames(bar_per)[1] <- 'sample';colnames(bar_per)[2] <- 'celltype';colnames(bar_per)[5] <- 'Percent'
  
  p <- ggplot(bar_per, aes(x = sample, y = Percent, fill = celltype,
                           stratum = celltype, alluvium = celltype)) +
    geom_col(position = 'stack', width = 0.6) +
    geom_stratum(width = 0.6, color = 'white') +
    geom_alluvium(alpha = 0.4, width = 0.6, color = 'white', linewidth = 1, curve_type = "linear") +
    scale_fill_manual(values = colors) +
    xlab('') + 
    ylab('') +
    scale_y_continuous(expand = c(0, 0))+
    theme_bw(base_size = 12) + 
    theme(
      axis.text = element_text(color = "black"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12)
    )
  p
  ggsave(filename = 'bar_plus.pdf',width = 7.5,height = 6)
}


do_positive_negative <- function(sc,gene){
  library(Seurat)
  df <- FetchData(sc, vars = gene)
  df$status <- ifelse(df[,1] > 0, paste0(gene,'+'), paste0(gene,'-'))
  sc <- AddMetaData(object = sc,metadata = df)
  sc$status <- factor(sc$status,levels = c(paste0(gene,'+'), paste0(gene,'-')))
  Idents(sc) <- 'status'
  return(sc)
}


do_sc_sample_cell <- function(obj,group.by="seurat_clusters",sp.size=NULL,diet="true",sp.total=1000) {
  library(Seurat)
  
  # 1.使用DietSeurat对seurat对象瘦身，加快后续操作
  all <- obj
  if (diet=="true") {
    all <- DietSeurat(all)
  }
  
  # 2.计算每个分组需要选取多少细胞，即sp.total/nlen
  if (is.null(sp.size)) {
    nlen <- length(unique(all@meta.data[,group.by]))
    sp.size <- ceiling(sp.total/nlen)
  }
  
  # 3.执行挑取细胞，这里会报错，还不清楚怎样解决！！！
  seob_list <- list()
  i <- 1
  for (scobj in unique(all@meta.data[,group.by])){
    cellist <- colnames(all)[which(all@meta.data[,group.by] == scobj)]
    ob <- subset(all, cells=cellist)
    if (length(colnames(ob)) > sp.size) {
      ob <- subset(ob,cells=sample(colnames(ob), sp.size))
    }
    seob_list[[i]] <- ob
    i <- i+1
  }
  all <- Reduce(merge,seob_list)
  return(all)
}


do_sc_inferCNV <- function(dat,ref,group){
  library(Seurat)
  library(tidyverse)
  library(infercnv)
  
  sc <- dat
  # V5运行下面的代码！！！
  sc <- JoinLayers(sc)
  exprMatrix <- LayerData(object = sc,assay = 'RNA',layer = 'counts') %>% as.matrix()
  
  # 非V5直接运行下面代码！！！
  exprMatrix <- as.matrix(GetAssayData(sc, slot='counts'))
  
  cellAnnota <- subset(sc@meta.data, select=group)
  groupFiles <- 'groupFiles.txt'
  write.table(cellAnnota,file ="groupFiles.txt",sep = '\t',col.names = F)
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = exprMatrix,
                                       annotations_file = "groupFiles.txt",
                                       delim = "\t",
                                       gene_order_file = "A:/data/geneLocate.txt",
                                       ref_group_names = ref)
  infercnv_obj2 <- infercnv::run(infercnv_obj,
                                 cutoff = 0.1,            
                                 out_dir = 'cnv.REF/' ,   
                                 cluster_by_groups = T,  
                                 hclust_method = "ward.D2", 
                                 plot_steps = F,
                                 denoise = TRUE,
                                 HMM = F,
                                 # denoise = F,             
                                 num_threads = 6,
                                 write_expr_matrix = T)
  infer_CNV_obj<-readRDS('cnv.REF/run.final.infercnv_obj')
  expr<-infer_CNV_obj@expr.data
  data_cnv <- as.data.frame(expr)
  meta <- sc@meta.data
  if(T){
    tmp1=expr[,infer_CNV_obj@reference_grouped_cell_indices$T_cell]
    tmp=tmp1
    tmp2=expr[,infer_CNV_obj@reference_grouped_cell_indices$B_cell]
    tmp=cbind(tmp1,tmp2)
    down=mean(rowMeans(tmp))-2*mean(apply(tmp,1,sd))
    up=mean(rowMeans(tmp))+2*mean(apply(tmp,1,sd))
    oneCopy=up-down
    oneCopy
    a1=down-2*oneCopy
    a2=down-1*oneCopy
    down;up
    a3=up+1*oneCopy
    a4=up+2*oneCopy
    cnv_score_table<-infer_CNV_obj@expr.data
    cnv_score_table[1:4,1:4]
    cnv_score_mat<-as.matrix(cnv_score_table)
    cnv_score_table[cnv_score_mat>0&cnv_score_mat<a2]<-"A"#complete loss.2pts
    cnv_score_table[cnv_score_mat>=a2&cnv_score_mat<down]<-"B"#loss of one copy.1pts
    cnv_score_table[cnv_score_mat>=down&cnv_score_mat<up]<-"C"#Neutral.0pts
    cnv_score_table[cnv_score_mat>=up&cnv_score_mat<=a3]<-"D"#addition of one copy.1pts
    cnv_score_table[cnv_score_mat>a3&cnv_score_mat<=a4]<-"E"#addition of two copies.2pts
    cnv_score_table[cnv_score_mat>a4]<-"F"#addition of more than two copies.2pts
    table(cnv_score_table[,1])
    cnv_score_table_pts<-cnv_score_mat
    rm(cnv_score_mat)
    cnv_score_table_pts[cnv_score_table=="A"]<-2
    cnv_score_table_pts[cnv_score_table=="B"]<-1
    cnv_score_table_pts[cnv_score_table=="C"]<-0
    cnv_score_table_pts[cnv_score_table=="D"]<-1
    cnv_score_table_pts[cnv_score_table=="E"]<-2
    cnv_score_table_pts[cnv_score_table=="F"]<-2
    cnv_score_table_pts[1:4,1:4]
    str(as.data.frame(cnv_score_table_pts[1:4,1:4]))
    cell_scores_CNV<-as.data.frame(colSums(cnv_score_table_pts))
    colnames(cell_scores_CNV)<-"cnv_score"
  }
  score <- cell_scores_CNV
  meta$totalCNV <- score[match(colnames(sc),rownames(score)),1]
  identical(rownames(sc@meta.data),rownames(meta))
  sc@meta.data$TotalCNV <- meta$totalCNV
  return(sc)
}


do_sc_GSEA <- function(dat,group,org){
  library(Seurat)
  library(SCP)
  library(BiocParallel)
  library(ggplot2)
  
  if (org == 'mouse') {
    specie <- 'Mus_musculus'
  }else{
    specis <- 'Homo_sapiens'
  }
  sc <- dat
  register(MulticoreParam(workers = 8, progressbar = TRUE))
  sc <- RunDEtest(srt = sc, group_by = group, fc.threshold = 1, only.pos = FALSE)
  sc <- RunGSEA(
    srt = sc, 
    group_by = group,
    db = c("KEGG",'WikiPathway','Reactome','GO_BP'),
    species = specie,
    DE_threshold = "p_val_adj < 0.05",
    db_combine = T
  )
  result_all <- sc@tools[2]$GSEA_cell_type_wilcox$enrichment
  result_filter <- subset(result_all,result_all$p.adjust < 0.05)
  write.table(result_all,file = 'result_all.txt',sep = '\t',quote = F,col.names = NA)
  write.table(result_filter,file = 'result_filter.txt',sep = '\t',quote = F,col.names = NA)
  reture(sc)
}


plot_sc_feature_violin <- function(scRNA,GeneName,width=4.5,height=4){
  library(Seurat)
  library(ggplot2)
  pdf(paste0(GeneName,"_VlnPlot.pdf"), width = width, height = height)
  plot = VlnPlot(scRNA, features = GeneName, pt.size = 0)+ 
    scale_color_manual(values = alpha(c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF",
                                        "#1f77b4", "#ff7f0e", "#279e68", "#d62728","#aa40fc", "#8c564b", "#e377c2", "#b5bd61","#17becf","#aec7e8")
                                      ,0.5)) +
    scale_fill_manual(values = alpha(c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF",
                                       "#1f77b4", "#ff7f0e", "#279e68", "#d62728","#aa40fc", "#8c564b", "#e377c2", "#b5bd61","#17becf","#aec7e8")
                                     ,0.5))+
    scale_y_continuous(expand = c(0, 0))+
    xlab("Chemotherapy") +
    theme(
      panel.border = element_rect(colour = "black", fill=NA, size=0.4), 
      axis.ticks = element_line(size=0.2, color="black"),
      axis.ticks.length = unit(0.2,"cm"),
      legend.position = "none",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10))
  print(plot)
  dev.off()
}


get_sc_feature_exp <- function(scRNA,GeneName){
  expression_data <- GetAssayData(object = scRNA, assay = "RNA") %>% .[c(GeneName),] %>% as.matrix()
  if(length(GeneName)==1){
    colnames(expression_data) = GeneName
    expression_data = as.data.frame(expression_data)
  }else{
    expression_data = t(expression_data)
    expression_data = as.data.frame(expression_data)
  }
  return(expression_data)
}


get_gs_marker <- function(scRNA, idents) {
  library(genesorteR)
  library(Seurat)
  library(dplyr)
  
  # 一种能够代替findallmarker的算法！！！
  # 设置Seurat对象中的细胞身份
  Idents(scRNA) <- idents
  
  # 对基因进行排序并获取标记基因
  gs <- sortGenes(scRNA@assays$RNA@data, Idents(scRNA))
  gs_marker <- getMarkers(gs, quant = 0.99)
  
  # 提取特异性分数和基因名
  spec_scores <- gs$specScore
  genes <- rownames(spec_scores)
  
  # 初始化一个空的数据框来存储每个簇的前10个基因
  top10_genes_df <- data.frame()
  
  # 初始化一个列表来存储每个簇的cluster_scores
  cluster_scores_list <- list()
  
  # 遍历每个簇，获取前10个基因
  for (cluster in colnames(spec_scores)) {
    cluster_scores <- data.frame(
      gene = genes,
      specScore = spec_scores[, cluster],
      cluster = rep(cluster, length(genes))
    )
    top10_genes <- cluster_scores %>%
      arrange(desc(specScore)) %>%
      head(10)
    
    # 将每个簇的前10个基因添加到总的数据框中
    top10_genes_df <- rbind(top10_genes_df, top10_genes)
    
    # 将每个簇的cluster_scores添加到列表中
    cluster_scores_list[[cluster]] <- cluster_scores
  }
  
  # 返回一个包含top10_genes_df和cluster_scores_list的列表
  return(list(top10_genes_df = top10_genes_df, cluster_scores_list = cluster_scores_list))
}


get_sc_sub_cluster <- function(dat,resolution){
  # example
  # dat <- get_sc_sub_cluster(dat = sc,resolution = 1.5)
  
  library(Seurat)
  library(clustree)
  
  sc <- dat
  sc <- FindVariableFeatures(object = sc,nfeatures = 3000)
  sc <- Seurat::RunPCA(sc,assay = 'RNA',npcs = 50)
  sc <- FindNeighbors(sc, dims = 1:15, reduction = "pca")
  for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1,1.2,1.5,2,2.5,3)) {
    sc <- FindClusters(sc, graph.name = "RNA_snn", resolution = res, algorithm = 1)
  }
  apply(sc@meta.data[,grep("RNA_snn_res",colnames(sc@meta.data))],2,table)
  p2_tree <- clustree(sc@meta.data, prefix = "RNA_snn_res.")
  pdf(file = 'cutree.pdf',width = 10,height = 10)
  p2_tree
  dev.off()
  sc <- FindClusters(sc,resolution = resolution)
  print(table(sc$seurat_clusters))
  return(sc)
}


do_sc_AUCell = function(scRNA,markers,cores_use = 5){
  library(Seurat)
  library(AUCell)
  markers = intersect(markers,rownames(scRNA))
  geneSets = list(Custom = markers)
  cells_rankings = AUCell_buildRankings(as.matrix(GetAssayData(object = scRNA, assay = "RNA")),
                                        nCores = cores_use)
  
  merge_row <- function(data1,data2){
    samesample = intersect(rownames(data1),rownames(data2))
    data1 = data1[samesample,,drop=FALSE]
    data2 = data2[samesample,,drop=FALSE]
    data3 = cbind(data1,data2)
    return(data3)
  }
  
  cells_AUC = AUCell_calcAUC(geneSets,
                             cells_rankings,
                             nCores = cores_use,
                             aucMaxRank = nrow(cells_rankings)*0.05)
  #AUC_Exp = as.numeric(getAUC(cells_AUC)["Custom",])
  #scRNA$AUC_Exp = AUC_Exp
  AUC_Score = cells_AUC@assays@data@listData[["AUC"]]
  AUC_Score = t(AUC_Score)
  colnames(AUC_Score) = "AUC"
  scRNA@meta.data = merge_row(scRNA@meta.data,AUC_Score)
  return(scRNA)
}


do_sc_UCell = function(scRNA,markers,cores_use = 5, name_use = ""){
  library(Seurat)
  library(UCell)
  markers = intersect(markers,rownames(scRNA))
  geneSets = list(UCell = markers)
  scRNA = AddModuleScore_UCell(scRNA,features = geneSets,ncores = cores_use,name = name_use)
  return(scRNA)
}


do_sc_AddModuleScore = function(scRNA,markers,name_use = 'AddModuleScore'){
  library(Seurat)
  markers = intersect(markers,rownames(scRNA))
  geneSets = list(Custom = markers)
  scRNA = AddModuleScore(object = scRNA,features = geneSets, ctrl = 100, name = name_use)
  colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == paste0(name_use,"1")] = name_use
  return(scRNA)
}


do_sc_nmf = function(scRNA, nfeatures = 2000,rank_use = 12){
  library(Seurat)
  library(NMF)
  library(tidyverse)
  
  scRNA = FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = nfeatures)
  scRNA = ScaleData(scRNA,do.center = FALSE)
  vm = scRNA@assays$RNA$scale.data
  res = nmf(vm,rank_use,method = "snmf/r")
  
  fs = extractFeatures(res, 30L)
  fs = lapply(fs, function(x) rownames(res)[x])
  fs = do.call("rbind", fs)
  
  rownames(fs) = paste0("cluster",1:rank_use)
  write.csv(t(fs), "c_NMF_TopGenes.csv")
  DT::datatable(t(fs))
  scRNA = RunPCA(scRNA,verbose = F)
  scRNA@reductions$nmf = scRNA@reductions$pca
  scRNA@reductions$nmf@cell.embeddings = t(coef(res))
  scRNA@reductions$nmf@feature.loadings = basis(res)
  scRNA = RunUMAP(scRNA, reduction = 'nmf', dims = 1:rank_use)
  
  ## 基于NMF降维矩阵的聚类
  scRNA = FindNeighbors(scRNA, reduction='nmf', dims = 1:rank_use) %>% FindClusters()
  
  ## 基于因子最大载荷分类
  scRNA$nmf_cluster = apply(NMF::coefficients(res)[1:rank_use,], 2, which.max)
  
  return(scRNA)
}


do_sc_cellchat = function(scRNA,
                          gene_use = NULL,
                          cell_use,
                          label_use = "CellType",
                          is_v5 = FALSE,
                          height_1 = 8,
                          height_2 = 8,
                          org) {
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(Matrix)
  library(limma)
  library(NMF)
  library(ggplot2)
  library(ggalluvial)
  library(svglite)
  library(CellChat)
  
  if(is_v5 == TRUE){
    scRNA[["RNA"]] = as(scRNA[["RNA"]], "Assay")
  }
  colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == label_use] = "labels"
  cellchat = createCellChat(object = scRNA@assays$RNA@data, meta = scRNA@meta.data, group.by = "labels")
  
  #V5需要运行下面的两行代码
  #rownames(cellchat@data) = rownames(scRNA@assays$RNA)
  #colnames(cellchat@data) = colnames(scRNA@assays$RNA)
  
  cellchat = setIdent(cellchat, ident.use="labels")
  groupSize = as.numeric(table(cellchat@idents))
  
  if (org == 'human') {
    CellChatDB = CellChatDB.human
    CellChatDB.use = subsetDB(CellChatDB, search = "Secreted Signaling")
    cellchat@DB = CellChatDB.use
    
    cellchat = subsetData(cellchat)
    cellchat = identifyOverExpressedGenes(cellchat)
    cellchat = identifyOverExpressedInteractions(cellchat)
    cellchat = projectData(cellchat, PPI.human)
    cellchat = computeCommunProb(cellchat)
    cellchat = filterCommunication(cellchat, min.cells = 10)
  } else if (org == 'mouse') {
    CellChatDB = CellChatDB.mouse
    CellChatDB.use = subsetDB(CellChatDB, search = "Secreted Signaling")
    cellchat@DB = CellChatDB.use
    
    cellchat = subsetData(cellchat)
    cellchat = identifyOverExpressedGenes(cellchat)
    cellchat = identifyOverExpressedInteractions(cellchat)
    cellchat = projectData(cellchat, PPI.mouse)
    cellchat = computeCommunProb(cellchat)
    cellchat = filterCommunication(cellchat, min.cells = 10)
  }
  
  df.net = subsetCommunication(cellchat)
  write.table(df.net,file="1.CellChat.Comm.network.xls",  sep="\t", row.names=F, quote=F)
  
  cellchat = computeCommunProbPathway(cellchat)
  cellchat = aggregateNet(cellchat)
  pdf("2.CellChat_cellNetwork_Count.pdf", width = 8, height = 7)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge= F, title.name = "Number of interactions")
  dev.off()
  
  pdf("3.CellChat_cellNetwork_Weight.pdf", width = 8, height = 7)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge= F, title.name = "Interaction strength")
  dev.off()
  
  pdf("4.CellChat_singleCell_weight.pdf", width = 10, height = 6)
  weight_mat = cellchat@net$weight
  par(mfrow = c(2,4), mgp=c(0,0,0), xpd=TRUE)
  for (cel in unique(cellchat@idents)){
    cir_mat = matrix(0, nrow = nrow(weight_mat),
                     ncol = ncol(weight_mat),
                     dimnames = dimnames(weight_mat)
    )
    cir_mat[cel, ] = weight_mat[cel, ]
    netVisual_circle(cir_mat,
                     vertex.weight= groupSize,
                     weight.scale= T,
                     edge.weight.max = max(weight_mat),
                     vertex.label.cex = 0.8,
                     title.name=cel
    )
  }
  dev.off()
  
  pdf("5.CellChat_singleCell_count.pdf", width = 10, height = 6)
  weight_mat = cellchat@net$count
  par(mfrow = c(2,4), mgp = c(0,0,0), xpd=TRUE)
  for (cel in unique(cellchat@idents)){
    cir_mat = matrix(0,
                     nrow = nrow(weight_mat),
                     ncol = ncol(weight_mat),
                     dimnames = dimnames(weight_mat)
    )
    cir_mat[cel, ] = weight_mat[cel, ]
    netVisual_circle(cir_mat,
                     vertex.weight= groupSize,
                     weight.scale= T,
                     edge.weight.max = max(weight_mat),
                     vertex.label.cex=0.8,title.name=cel
    )
  }
  dev.off()
  
  # pdf("6.CellChat.bubble.pdf", width = 12, height = 12)
  # p = netVisual_bubble(cellchat, remove.isolate = FALSE, angle.x = 45)
  # print(p)
  # dev.off()#已损坏
  
  cellchat = netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  pdf("7.CellChat_outgoing-incoming-importance.pdf", width=6, height=5)
  gg1 = netAnalysis_signalingRole_scatter(cellchat)
  print(gg1)
  dev.off()
  
  pdf("8.CellChat.outgoing-incoming-pathway.pdf", width = 10, height = height_1)
  ht1 = netAnalysis_signalingRole_heatmap(cellchat,
                                          font.size = 6,
                                          height = height_1,
                                          color.heatmap = "Reds",
                                          pattern = "outgoing"
  )
  ht2 = netAnalysis_signalingRole_heatmap(cellchat,
                                          font.size = 6,
                                          height = height_1,
                                          color.heatmap = "BuGn",
                                          pattern = "incoming"
  )
  
  ht = ht1 + ht2
  print(ht)
  dev.off()
  if(!is.null(gene_use)){
    p1 = netVisual_bubble(cellchat,
                          sources.use = levels(cellchat@idents),
                          targets.use = paste0(gene_use,"+ ",cell_use),
                          remove.isolate = FALSE)
    
    p2 = netVisual_bubble(cellchat,
                          sources.use = levels(cellchat@idents),
                          targets.use = paste0(gene_use,"- ",cell_use), remove.isolate = FALSE)
    
    
    p3 = netVisual_bubble(cellchat,
                          sources.use = paste0(gene_use,"+ ",cell_use),
                          targets.use = levels(cellchat@idents),
                          remove.isolate = FALSE)
    
    p4 = netVisual_bubble(cellchat,
                          sources.use = paste0(gene_use,"- ",cell_use),
                          targets.use = levels(cellchat@idents),
                          remove.isolate = FALSE)
    
    p5 = p1 + p3
    pdf(paste0("9.CellChat_sources_targets_positive.pdf"),  width = 16, height = height_2)
    print(p5)
    dev.off()
    
    p6 = p2 + p4
    pdf(paste0("9.CellChat_sources_targets_negative.pdf"),  width = 16, height = height_2)
    print(p6)
    dev.off()
  }else if(is.null(gene_use)){
    p1 = netVisual_bubble(cellchat,
                          sources.use = levels(cellchat@idents),
                          targets.use = cell_use,
                          remove.isolate = FALSE)
    
    p3 = netVisual_bubble(cellchat,
                          sources.use = cell_use,
                          targets.use = levels(cellchat@idents),
                          remove.isolate = FALSE)
    
    p5 = p1 + p3
    pdf(paste0("9.CellChat_sources_targets.pdf"),  width = 16, height = height_2)
    print(p5)
    dev.off()
  }
  
  for(i in cellchat@netP$pathways){
    pathways.show = i
    pdf(paste0("10.CellChat_", pathways.show,"_circle.pdf"), width=8, height=6)
    circle = netVisual_aggregate(cellchat, 
                                 signaling=pathways.show, 
                                 layout="circle")
    print(circle)
    dev.off()
    
    # pdf(paste0("10.CellChat_", pathways.show, "_hierarchy.pdf"), width = 12, height = 6)
    # hierarchy = netVisual_aggregate(cellchat,
    #                                 signaling = pathways.show,
    #                                 layout = "hierarchy",
    #                                 vertex.receiver=c(1,2,3,4))
    # print(hierarchy)
    # dev.off()#会报错
    
    pdf(paste0("10.CellChat_", pathways.show, "_heatmap.pdf"), width=8, height=6)
    heatmap = netVisual_heatmap(cellchat, signaling=pathways.show, color.heatmap = "Reds", measure= 'weight')
    print(heatmap)
    dev.off()
    
    pdf(paste0("10.CellChat_", pathways.show,"_netAnalysis.pdf"), width=6, height=5)
    cellchat = netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    netAnalysis=netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 5, color.heatmap = "Reds", font.size = 12)
    print(netAnalysis)
    dev.off()
    
    pdf(paste0("10.CellChat.", pathways.show,".contribution.pdf"), width=8, height=2)
    contribution = netAnalysis_contribution(cellchat, signaling= pathways.show)
    print(contribution)
    dev.off()
    
    pdf(paste0("10.CellChat.", pathways.show, ".geneExp.pdf"), width=8, height=6)
    geneExp = plotGeneExpression(cellchat, signaling=pathways.show)
    print(geneExp)
    dev.off()
    
    pairLR = extractEnrichedLR(cellchat, signaling=pathways.show, geneLR.return=FALSE)
    pdf(paste0("10.CellChat.", pathways.show, ".pairLR.pdf"), width=9, height=8)
    pairCircos=netVisual_individual(cellchat, signaling=pathways.show, pairLR.use=pairLR[1] , layout="circle" )
    print(pairCircos)
    dev.off()
    
    for(i in 1:nrow(pairLR)){
      pdf(paste0("10.CellChat_", pathways.show,"_", pairLR[i,], "_pairLR.pdf"), width = 8, height = 6)
      pairChord = netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = pairLR[i,], layout = "chord")
      print(pairChord)
      dev.off()
    }
  }
  
  return(cellchat)
}


do_human_mouse <- function(human_gene){
  library(homologene)
  
  mouse_gene <- homologene(human_gene,inTax = 9606, outTax = 10090)
  mouse_gene <- mouse_gene[,'10090']
  return(mouse_gene)
}



do_mouse_human <- function(mouse_gene){
  library(homologene)
  
  human_gene <- homologene(mouse_gene,inTax = 10090,outTax = 9606)
  human_gene <- human_gene[,'9606']
  return(human_gene)
}