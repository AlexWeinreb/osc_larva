

# Run on cluster: load CellRanger output and run DoubltFinder


# Inits ----
library(Seurat)
library(DoubletFinder)

library(wbData)

gids <- wb_load_gene_coords(295)

# output
dir_save <- "intermediates/2502/250304_doubletFinder"

# inputs
dir_soupx <- "intermediates/2502/soupx"
samples <- list.files(dir_soupx) |>
  stringr::str_remove("_soupx_feature_bc_matrix$")




# predict nb of dblts (https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/76)

exp_rate <- read.table(text="cells	dblt
500	0.4
1000	0.8
2000	1.6
3000	2.3
4000	3.1
5000	3.9
6000	4.6
7000	5.4
8000	6.1
9000	6.9
10000	7.6
", header = TRUE)
db_rate_mod <- lm(dblt ~ cells, data = exp_rate)




# For each sample ----

for(sample in samples){
  
  message("-----   sample: ", sample)
  
  # load and do minimal preprocessing
  seu <- file.path(dir_soupx, paste0(sample, "_soupx_feature_bc_matrix")) |>
    Read10X(gene.column = 1) |>
    CreateSeuratObject(project = sample)
  
  # gene sets
  mt_genes <- gids |>
    dplyr::filter(chr == "MtDNA") |>
    dplyr::filter(gene_id %in% rownames(seu)) |>
    dplyr::pull(gene_id)
  
  
  
  seu[["percent.mt"]] <- PercentageFeatureSet(seu,
                                              features = mt_genes)
  
  
  
  
  
  subseu <- subset(seu, percent.mt < 15 & nFeature_RNA > 250)
  
  
  
  
  #~ preprocess ----
  
  subseu <- NormalizeData(subseu)
  subseu <- FindVariableFeatures(subseu)
  subseu <- ScaleData(subseu)
  
  subseu <- RunPCA(subseu, npcs = 20, verbose = FALSE)
  
  subseu <- FindNeighbors(subseu, dims = 1:20)
  subseu <- FindClusters(subseu,
                         resolution = .8)
  
  
  #~ Doublets ----
  
  sweep_dat1 <- paramSweep(subseu,
                           PCs = 1:20,
                           sct = FALSE)
  
  summ_sweep_dat1 <- summarizeSweep(sweep_dat1)
  
  bcmvn <- find.pK(summ_sweep_dat1)
  
  
  
  png(file.path(dir_save, paste0(sample, "_bcmvn.png")))
  plot(x = bcmvn$ParamID, y = bcmvn$BCmetric, pch = 16, 
       col = "#41b6c4", cex = 0.75)
  lines(x = bcmvn$ParamID, y = bcmvn$BCmetric, 
        col = "#41b6c4")
  abline(v = bcmvn$pK[which.max(bcmvn$BCmetric)],
         lty = 'dotted', lwd = 2, col = 'grey')
  dev.off()
  
  
  
  mpK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  db_rate_pred <- predict(db_rate_mod, newdata = data.frame(cells=ncol(subseu)))
  
  message("-----   Predicted pK: ", mpK)
  message("-----   Predicted dblt rate: ", round(db_rate_pred, 2))
  
  
  homotypic.prop <- modelHomotypic(Idents(subseu))
  nExp_poi <- round( db_rate_pred * ncol(subseu) / 100 )
  nExp_poi.adj <- round( nExp_poi * (1-homotypic.prop) )
  
  
  subseu <- doubletFinder(subseu,
                          PCs = 1:20,
                          pK = mpK,
                          nExp = nExp_poi,
                          reuse.pANN = FALSE,
                          sct = FALSE)
  
  message("-----   Doublet estimate done")
  
  # bring the doublet info into the unsubsetted object
  df_colname <- colnames(subseu[[]])[startsWith(colnames(subseu[[]]), "DF.classifications")]
  
  df_res <- subseu[[]][,df_colname, drop = FALSE]
  names(df_res) <- "doubletFinder"
  
  seu[["doubletFinder"]] <- df_res[rownames(seu[[]]), "doubletFinder"]
  
  table(seu$doubletFinder, useNA = 'always')
  
  
  message("-----   Saving ", sample)
  qs::qsave(seu, file.path(dir_save, paste0("250304_", sample,"_dblts.qs")) )
  
}



message("-----   sessioninfo ")

sessionInfo()





