
# Inits ----
library(tidyverse) |> suppressPackageStartupMessages()
library(Seurat) |> suppressPackageStartupMessages()


source("https://raw.githubusercontent.com/yanailab/knn-smoothing/refs/heads/master/knn_smooth.R")


library(wbData)

gids_230 <- wb_load_gene_ids(230)


dir_third_processed <- "intermediates/2502/250313_third_processed"
dir_out <- "intermediates/2502/250319_step1"



osc <- readxl::read_excel("../10x_grl18/data/oscillating/msb209498-sup-0003-datasetev1.xlsx",
                          sheet = "Dataset EV1 WBidToGeneNames_Osc",
                          na = "NA") |>
  mutate(gene_name = i2s(WB_ID, gids_230),
         PeakPhase_rad = PeakPhase*pi / 180 )






# Load ----

files_list <- list.files(dir_third_processed) |>
  enframe(value = "filename",
          name = NULL) |>
  mutate(filename2 = filename |>
           str_replace("seu_1","g1") |>
           str_replace("seu_2","g2")) |>
  filter(str_detect(filename2, "2[35]031[0-9]_(g1|g2)")) |>
  separate_wider_regex(filename2,
                       patterns = c(
                         "^2[35]031[0-9]_",
                         group = "(?:g1|g2)",
                         "_",
                         cell_type = "[a-zA-Z0-9_]+",
                         "\\.qs$"
                       ))

files_list |>
  count(cell_type) |>
  arrange(desc(n)) |>
  count(n)


files_list_by_type <- files_list |>
  summarize(both_names = list(filename),
            n = n(),
            .by = cell_type) |>
  filter(n == 2) |>
  select(-n)



#~ Next cell type ----

for(i in seq_along(files_list_by_type$cell_type)){
  
  
  .ct <- files_list_by_type$cell_type[[i]]
  .files <- files_list_by_type$both_names[[i]]
  
  
  message("---------  ", i,": ", .ct, "  ---------")
  .ct
  
  stopifnot(length(.files) == 2L)
  
  
  #~ load data ----
  seu_list <- .files |>
    map(~ qs::qread(file.path( dir_third_processed, .x ))) |>
    map(~{
      DefaultAssay(.x) <- "RNA"
      .x[["SCT"]] <- NULL
      
      .x <- JoinLayers(.x)
      .x
    })
  
  seu <- merge(
    seu_list[[1]],
    seu_list[[2]]
  ) |>
    JoinLayers()
  
  
  
  
  
  #~ add phase ----
  
  mat <- LayerData(seu)[osc$gene_name[osc$Class == "Osc"] |> intersect(rownames(seu)),]
  dim(mat)
  
  
  osc_ordered <- osc[match(rownames(mat), osc$gene_name),]
  
  
  stopifnot(all.equal( osc_ordered$gene_name, rownames(mat) ))
  
  mean_angles <- tibble(cell_bc = colnames(mat),
                        mean_angle = NA_real_,
                        mean_rho = NA_real_)
  
  
  for(rr in which(is.na(mean_angles$mean_angle))){
    
    expression <- mat[, rr ]
    
    mean_angles$mean_angle[[rr]] <- circhelp::weighted_circ_mean(osc_ordered$PeakPhase_rad,
                                                                 expression) %% (2*pi)
    mean_angles$mean_rho[[rr]] <- max(expression)*circhelp::weighted_circ_rho(osc_ordered$PeakPhase_rad,
                                                                              expression) %% (2*pi)
    
  }
  
  
  
  stopifnot(identical(seu[[]] |> rownames(),
                      mean_angles$cell_bc))
  
  # list(seu = seu[[]] |> rownames(),
  #      ang = mean_angles$cell_bc) |>
  #   eulerr::euler() |>
  #   plot(quantities = TRUE)
  
  
  seu$mean_angle <- mean_angles$mean_angle
  seu$mean_rho <- mean_angles$mean_rho
  
  
  
  
  
  #~ gene expression ----
  cnts_raw <- GetAssayData(seu, assay = "RNA", layer = "counts")
  
  gene_expressions <- data.frame(
    gene_name = rownames(cnts_raw),
    nb_cells = rowSums(cnts_raw > 0),
    prop_cells = rowMeans(cnts_raw > 0)
  )
  
  
  
  
  
  #~ smooth ----
  
  smoothed <- knn_smoothing(cnts_raw,
                            k = 5)
  
  seu <- CreateSeuratObject(smoothed |> as("dgCMatrix"),
                            meta.data = seu[[]])
  
  
  
  
  #~ PCA ----
  seu <- SCTransform(seu)
  
  seu <- RunPCA(seu, npcs = 2, verbose = FALSE)
  
  
  
  #~ save ----
  
  gg_batch <- DimPlot(seu,
                      reduction = "pca",
                      group.by = "orig.ident",
                      pt.size = 2,
                      alpha = .2) +
    NoLegend()
  
  
  ggsave(paste0(.ct,"_batch.png"), gg_batch,
         path = dir_out,
         width = 7, height = 5, units = "in")
  
  
  
  gg_phase <- ggplot(FetchData(seu,
                               vars = c("PC_1", "PC_2",
                                        "mean_angle", "mean_rho"))) +
    theme_classic() +
    scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                          limits = c(0, 2*pi)) +
    geom_point(aes(x = PC_1, y = PC_2,
                   color = mean_angle,
                   alpha = mean_rho))
  
  ggsave(paste0(.ct,"_phase.png"), gg_phase,
         path = dir_out,
         width = 7, height = 5, units = "in")
  
  
  
  
  
  qs::qsave(
    seu,
    file.path(dir_out,
              paste0(.ct, "_seu.qs"))
  )
  
  qs::qsave(
    gene_expressions,
    file.path(dir_out,
              paste0(.ct, "_gene_expressions.qs"))
  )
  
  
}


# End ----
message("Done  ", date())

sessionInfo()





# visual comparison for several values of k

# res <- map(c(3,5,7,10,20) |> set_names(),
#            ~{
#              smoothed <- knn_smoothing(cnts_raw,
#                                        k = .x)
#              
#              seu1 <- CreateSeuratObject(smoothed |> as("dgCMatrix"),
#                                         meta.data = seu[[]])
#              
#              
#              
#              # reprocess
#              seu1 <- SCTransform(seu1)
#              
#              nps_max <- pmin(200, ncol(seu1) - 2L)
#              seu1 <- RunPCA(seu1, npcs = nps_max, verbose = FALSE)
#              seu1
#            })
# 
# iwalk(res,
#       ~{
#         k <- .y
#         g1 <- DimPlot(.x,
#                       reduction = "pca",
#                       group.by = "orig.ident",
#                       pt.size = 2,
#                       alpha = .2) +
#           NoLegend()
#         ggsave(paste0("tmp/",.y,"_batch.png"), plot = g1)
#         
#         g2 <- ggplot(FetchData(.x,
#                                vars = c("PC_1", "PC_2",
#                                         "mean_angle", "mean_rho"))) +
#           theme_classic() +
#           scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
#                                 limits = c(0, 2*pi)) +
#           geom_point(aes(x = PC_1, y = PC_2,
#                          color = mean_angle,
#                          alpha = mean_rho))
#         ggsave(paste0("tmp/",.y,"_phase.png"), plot = g2)
#       })


