
# Inits ----
library(tidyverse) |> suppressPackageStartupMessages()
library(Seurat) |> suppressPackageStartupMessages()


# cloned for reproducibility
# source("https://raw.githubusercontent.com/yanailab/knn-smoothing/refs/heads/master/knn_smooth.R")
source("https://raw.githubusercontent.com/AlexWeinreb/knn-smoothing/refs/heads/master/knn_smooth.R")





dir_assembled <- "intermediates/2502/250605_assembled/"
dir_out <- "intermediates/2502/250609_step1"







# Load ----

seu <- qs::qread( file.path(dir_assembled, "250606_seu_all_herma.qs"))



cell_types <- seu[[]] |>
  count(cell_type, stage,
        name = "nb_cells") |>
  summarize(nb_stages = n(),
            nb_cells = sum(nb_cells),
            .by = cell_type) |>
  filter(nb_stages >= 1,
         nb_cells > 20) |>
  pull(cell_type)


message("Processing ", length(cell_types)," cell types")

for(.ct in cell_types ){
  
  message("---------  ", .ct, "  ---------")
  
  subseu <- subset(seu, cell_type == .ct)
  
  
  
  #~ gene expression ----
  cnts_raw <- GetAssayData(subseu, assay = "RNA", layer = "counts")
  
  gene_expressions <- data.frame(
    gene_name = rownames(cnts_raw),
    nb_cells = rowSums(cnts_raw > 0),
    prop_cells = rowMeans(cnts_raw > 0)
  )
  
  
  
  
  #~ unsmoothed ----
  
  message("    --- unsmoothed")
  sub_unsmoothed <- SCTransform(subseu)|>
    RunPCA(npcs = 2, verbose = FALSE)
  
  
  gg_batch_unsmoothed <- DimPlot(sub_unsmoothed,
                      reduction = "pca",
                      group.by = "orig.ident",
                      pt.size = 2,
                      alpha = .2) +
    NoLegend()
  
  
  ggsave(paste0(.ct,"_batch_unsmoothed.png"), gg_batch_unsmoothed,
         path = dir_out,
         width = 7, height = 5, units = "in")
  
  
  gg_phase_unsmoothed <- ggplot(FetchData(sub_unsmoothed,
                               vars = c("PC_1", "PC_2",
                                        "cell_phase_masked", "cell_rho"))) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                          limits = c(0, 360)) +
    geom_point(aes(x = PC_1, y = PC_2,
                   color = cell_phase_masked),
               alpha = .7, size = 3, shape = 16)
  
  ggsave(paste0(.ct,"_phase_unsmoothed.png"), gg_phase_unsmoothed,
         path = dir_out,
         width = 7, height = 5, units = "in")
  
  
  
  qs::qsave(
    sub_unsmoothed,
    file.path(dir_out,
              paste0(.ct, "_seu_unsmoothed.qs"))
  )
  
  
  
  
  
  
  #~ smooth ----
  
  message("    --- smoothed")
  cnts_smoothed <- knn_smoothing(cnts_raw,
                            k = 5)
  
  
  sub_smoothed <- CreateSeuratObject(cnts_smoothed |> as("dgCMatrix"),
                                     meta.data = subseu[[]])
  
  
  
  
  #~ PCA ----
  sub_smoothed <- SCTransform(sub_smoothed)|>
    RunPCA(npcs = 2, verbose = FALSE)
  
  
  
  #~ save ----
  
  gg_batch <- DimPlot(sub_smoothed,
                      reduction = "pca",
                      group.by = "orig.ident",
                      pt.size = 2,
                      alpha = .2) +
    NoLegend()
  
  
  ggsave(paste0(.ct,"_batch.png"), gg_batch,
         path = dir_out,
         width = 7, height = 5, units = "in")
  
  
  # gg_stage <-  ggplot(FetchData(sub,
  #                               vars = c("PC_1", "PC_2",
  #                                        "stage"))) +
  #   theme_classic() +
  #   geom_point(aes(x = PC_1, y = PC_2,
  #                  color = stage),
  #              size = 2,
  #              alpha = .2)
  # 
  # ggsave(paste0(.ct,"_stage.png"), gg_stage,
  #        path = "presentations/",
  #        width = 7, height = 5, units = "in")
  
  
  gg_phase <- ggplot(FetchData(sub_smoothed,
                               vars = c("PC_1", "PC_2",
                                        "cell_phase_masked", "cell_rho"))) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                          limits = c(0, 360)) +
    geom_point(aes(x = PC_1, y = PC_2,
                   color = cell_phase_masked),
               alpha = .7, size = 3, shape = 16)
  
  ggsave(paste0(.ct,"_phase.png"), gg_phase,
         path = dir_out,
         width = 7, height = 5, units = "in")
  
  
  
  
  # # same with a continuous alpha scale
  # gg_phase_with_legend <- ggplot(FetchData(sub,
  #                  vars = c("PC_1", "PC_2",
  #                           "cell_phase_masked", "cell_rho"))) +
  #   theme_classic() +
  #   scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
  #                         limits = c(0, 360)) +
  #   geom_point(aes(x = PC_1, y = PC_2,
  #                  color = cell_phase_masked,
  #                  alpha = cell_rho),
  #              shape = 16,
  #              size = 3,
  #              show.legend = c(color = TRUE, alpha = FALSE)) +
  #   scale_fill_gradient(high = "black", low = "grey90",
  #                       limits = c(0,1)) +
  #   geom_point(aes(x = PC_1, y = PC_2,
  #                  fill = cell_rho),
  #              alpha = 0)
  # ggsave(paste0(.ct,"_phase.pdf"), gg_phase_with_legend,
  #        path = dir_out,
  #        width = 7, height = 5, units = "in")

  qs::qsave(
    sub_smoothed,
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
message("======================================================================")
message("                             sessionInfo                              ")
message("======================================================================")

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


