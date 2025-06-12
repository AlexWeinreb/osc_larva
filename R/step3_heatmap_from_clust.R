# Inits ----
library(tidyverse)

library(wbData)


gids <- wb_load_gene_ids(295) |>
  add_row(X = "a",
          gene_id = "nsIs198",
          symbol = "GFP",
          sequence = "GFP",
          status = "Live",
          biotype = "protein_coding_gene",
          name = "GFP"
  )

source("R/utils_heatmap_processing.R")

dir_clust <- "intermediates/2502/250609_cluster"

dir_step2 <- "intermediates/2502/250606_step2/"
# dir_step2 <- "E:/backups/Projects_june2025/glia/osc_larva/intermediates/2502/250609_step2/"


dir_step3 <- "intermediates/2502/250606_step3_genes_by_celltype/"



# ct2tissue <- read_csv("data/cell_type2tissue.csv")



# Load ----
cluster_results <- read_csv(file.path(dir_clust, "250610_cluster_results.csv"))

smooth_noncentered <- list.files(dir_step2,
                                 pattern = "_preds\\.qs$") |>
  enframe(value = "filename",
          name = NULL) |>
  separate_wider_regex(filename,
                       patterns = c(
                         cell_type = "^.+",
                         "_preds\\.qs"
                       ),
                       cols_remove = FALSE) |>
  pmap(\(cell_type, filename){
    mat_preds <- qs::qread(file.path(dir_step2,
                                     filename))
    colnames(mat_preds) <- paste0(cell_type, "|", colnames(mat_preds))
    
    mat_preds
  }) |>
  do.call(cbind, args = _)


# preprocess cell types
by_cell_type <- cluster_results |>
  summarize(n_puls = sum(shape == "pulsatile"),
            n_tot = n(),
            .by = "cell_type") |>
  mutate(prop_puls = round( 100 * n_puls / n_tot )) |>
  arrange(desc(prop_puls)) |>
  mutate(cell_type = fct_inorder(cell_type)) |>
  filter(n_puls > 0)

cluster_results <- cluster_results |>
  mutate(cell_type = factor(cell_type,
                            levels = levels(by_cell_type$cell_type))) |>
  mutate(cellgene = paste0(cell_type, "|", gene_name))





stopifnot(all.equal(
  colnames(smooth_noncentered),
  cluster_results$cellgene
))




# heatmap ----

# printMat::matimage(log1p(smooth_noncentered)[,sample(ncol(smooth_noncentered), 200)])



# create a heatmap for each cell type

heatmaps_list <- cluster_results |>
  filter(shape == "pulsatile") |>
  summarize(cellgenes = list(cellgene),
            .by = cell_type) |>
  arrange(cell_type) |>
  deframe() |>
  map(\(.cellgenes){
    
    smooth_ct <- smooth_noncentered[ , .cellgenes ]
    
    peak_location <- apply(smooth_ct, 2, which.max)
    
    smooth_ct[, order(peak_location)]
  })
  


# sapply(heatmaps_list[1:4], dim)

heatmaps_scaled <- map(heatmaps_list,
                       ~{
                         apply(log10( 1 + .x), 2,
                               \(col){
                                 (col - min(col)) / (max(col) - min(col))
                               })
                       })




#~ uniformity ----

stopifnot(all.equal(
  as.character(by_cell_type$cell_type),
  names(heatmaps_list)
))


by_cell_type$uniformity_index <- map_dbl(heatmaps_list,
                                         \(hm) {
                                           
                                           if(ncol(hm) == 0L) return(0)
                                           
                                           hm_norm <- apply(hm, 2, \(x){
                                             (x - min(x)) / (max(x) - min(x))
                                           })
                                           
                                           diffrange <- hm_norm |>
                                             rowMeans() |>
                                             range() |>
                                             diff()
                                           
                                           1 - diffrange
                                         })



hist(by_cell_type$uniformity_index, breaks = 50)



#~ entropy ----

# 
# ct <- "ILso"
# mat <- heatmaps_scaled[[ct]]
# y <- 2*pi*apply(mat, 2, which.max) / nrow(mat)
# 
# circular::rao.spacing.test(circular::circular(y))
# 
# pvals <- heatmaps_list |>
#   map_dbl(\(.hm){
#     pos_peaks <- apply(.hm, 2, which.max) / nrow(.hm)
#     goftest::cvm.test(pos_peaks, null = "punif")$p.value
#   })
# 
# hist(pvals, breaks = 50)
# 
# 
# ct <- "hypodermis"
# mat <- heatmaps_scaled[[ct]]
# pos_peaks <- apply(mat, 2, which.max) / nrow(mat)
# 
# printMat::matimage(mat)
# points(seq_len(ncol(mat)) / ncol(mat) ,
#        1 - apply(mat, 2, which.max) / nrow(mat),
#        col = 'purple', cex = 2, pch = "-")
# 
# hist(apply(mat, 2, which.max)/nrow(mat), breaks = 50)
# hist(runif(ncol(mat)), breaks = 50, add = TRUE, col = alpha('lightgreen', .2))
# 
# 
# rose.diag(circular::circular(pos_peaks*2*pi), bins = 24)
# 
# 
# qqplot(qunif(ppoints(ncol(mat))), apply(mat, 2, which.max)/nrow(mat))
# abline(0, 1, col = "red")


stopifnot(identical(
  as.character(by_cell_type$cell_type),
  names(heatmaps_list)
))


circ_entropy <- function(angles, num_bins = 50) {
  
  bins <- seq(0, 2*pi, length.out = num_bins + 1)
  counts <- hist(angles, breaks = bins, plot = FALSE)$counts
  entropy::entropy(counts)
}

by_cell_type$perplexity <- heatmaps_list |>
  map_dbl(\(.hm, .ct){
    pos_peaks <- apply(.hm, 2, which.max) / nrow(.hm)
    
    angles <- circular::circular(pos_peaks * 2 * pi)
    
    exp( circ_entropy(angles) )
  })


hist(by_cell_type$perplexity, breaks = 50)



#~ distance ----


by_cell_type$similarity_diag <- map_dbl(heatmaps_list,
                                        \(hm){
                                          
                                          if(ncol(hm) == 0L) return(0)
                                          
                                          hm_norm <- apply(hm, 2, \(x){
                                            (x - min(x)) / (max(x) - min(x))
                                          })
                                          
                                          n_t <- nrow(hm_norm)
                                          n_g <- ncol(hm_norm)
                                          
                                          sig_cent <- make_ref_sig(n_t)
                                          
                                          null_mat <- sapply(seq_len(n_g),
                                                             \(i){
                                                               circ_perm(sig_cent, floor(n_t * i/n_g))
                                                             })
                                          
                                          
                                          dist <- norm(hm_norm - null_mat, "1")
                                          similarity <- dist/n_t
                                          1 - similarity
                                        })


hist(by_cell_type$similarity_diag, breaks = 30)


ggplot(by_cell_type) +
  theme_classic() +
  geom_text(aes(x = uniformity_index, y = similarity_diag, label = cell_type))

ggplot(by_cell_type) +
  theme_classic() +
  geom_text(aes(x = perplexity, y = similarity_diag, label = cell_type))



# # dists with pvals
# 
# by_cell_type$p_val <- map_dbl(
#   heatmaps_list,
#   \(hm){
#     
#     if(ncol(hm) == 0L) return(1)
#     
#     hm_norm <- apply(hm, 2, \(x){
#       (x - min(x)) / (max(x) - min(x))
#     })
#     
#     
#     n_t <- nrow(hm_norm)
#     n_g <- ncol(hm_norm)
#     
#     sig_cent <- make_ref_sig(n_t)
#     
#     null_mat <- sapply(seq_len(n_g),
#                        \(i){
#                          circ_perm(sig_cent, floor(n_t * i/n_g))
#                        })
#     
#     
#     
#     similarity_perms <- c(
#       {
#         dist <- norm(hm_norm - null_mat, "1")
#         similarity <- dist/n_t
#         1 - similarity
#       },
#       replicate(10000,{
#         hm_perm <- hm_norm[,sample(ncol(hm_norm))]
#         dist <- norm(hm_perm - null_mat, "1")
#         similarity <- dist/n_t
#         1 - similarity
#       })
#     )
#     
#     pval <- mean(similarity_perms >= similarity_perms[[1]])
#     pval
#   },
#   .progress = TRUE
# )
# 
# 
# by_cell_type$p_adj <- by_cell_type$p_val |>
#   p.adjust(method = "BH")
# hist(by_cell_type$p_val, breaks = 30)
# hist(by_cell_type$p_adj, breaks = 30)
# 
# 
# 
# ggplot(by_cell_type) +
#   theme_classic() +
#   aes(x = uniformity_index, y = similarity_diag, label = cell_type, shape = p_adj < .05, color = p_adj < .05) +
#   geom_point() +
#   ggrepel::geom_text_repel()
# 
# 
# 
# by_cell_type |>
#   filter(prop_puls >= 10) |>
#   ggplot() +
#   theme_classic() +
#   aes(x = uniformity_index, y = similarity_diag, label = cell_type, shape = p_adj < .05, color = p_adj < .05) +
#   geom_point() +
#   ggrepel::geom_text_repel()
# 
# 
# ggplot(by_cell_type) +
#   theme_classic() +
#   aes(x = similarity_diag, y = -log10(p_adj),
#       label = cell_type, shape = p_adj < .05, color = p_adj < .05) +
#   geom_point() +
#   ggrepel::geom_text_repel()
# 
# by_cell_type |>
#   filter(prop_puls >= 10) |>
#   ggplot() +
#   theme_classic() +
#   aes(x = similarity_diag, y = -log10(p_adj),
#       label = cell_type, shape = p_adj < .05, color = p_adj < .05) +
#   geom_point() +
#   ggrepel::geom_text_repel()











# # only label every n-th column
# ct <- "BWM"
# hmp_sparsified <- heatmaps_list[[ct]]
# 
# spar_index <- 7
# colnames_to_sparsify <- setdiff(seq_len(ncol(hmp_sparsified)),
#                                 spar_index * seq_len( ncol(hmp_sparsified) / spar_index ) )
# 
# colnames(hmp_sparsified)[colnames_to_sparsify] <- ""
# head(colnames(hmp_sparsified), 20)
# 
# 
# pheatmap::pheatmap(hmp_sparsified,
#                    cluster_rows = FALSE,
#                    cluster_cols = FALSE,
#                    show_rownames = FALSE,
#                    fontsize = 7,
#                    main = ct)






# Compare bulk ----

dir_assembled <- "intermediates/2502/250605_assembled/"
mean_dotprod_by_celltype_res_perm <- qs::qread(file.path(dir_assembled, "250606_coherence_unnorm_perm10000.qs"))
dotprod_by_cell <- qs::qread(file.path(dir_assembled, "250610_dotprod_by_cell.qs"))


p_vals <- mean_dotprod_by_celltype_res_perm |>
  group_by(tissue, cell_type) |>
  nest() |>
  summarize(p_val = map_dbl(data,
                            \(dat){
                              mean(dat$mean_coherence >= dat$mean_coherence[[1]])
                            }),
            .groups = 'drop') |>
  mutate(p_adj = p.adjust(p_val, method = "holm"))


cell_types_bulk <- dotprod_by_cell |>
  summarize(mean_coherence = mean(coherence),
            .by = "cell_type") |>
  left_join(p_vals,
            by = c("cell_type")) |>
  mutate(p_adj = if_else(is.na(p_adj), 1, p_adj))

cell_types_both <- inner_join(
  by_cell_type |>
    select(cell_type, n_puls, n_tot, prop_puls, perplexity)
  ,
  cell_types_bulk |>
    select(cell_type, tissue, mean_coherence, p_coherence_adj = p_adj)
)


cell_types_both |>
  ggplot() +
  theme_classic() +
  xlab("Mean local phase coherence") +
  ylab("Perplexity") +
  scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = 19)) +
  geom_point(aes(x = mean_coherence, y = perplexity, color = tissue,
                 shape = p_coherence_adj < .05),
             size = 3)

# ggsave("phasic_cell_types_unannot.png", path = "presentations/figures/250612_celltype_osc",
#        width = 80, height = 50, units = "mm",
#        scale = 2)
# ggsave("phasic_cell_types_unannot.pdf", path = "presentations/figures/250612_celltype_osc",
#        width = 80, height = 50, units = "mm",
#        scale = 2)


cell_types_both |>
  ggplot() +
  theme_classic() +
  xlab("Mean local phase coherence (bulk)") +
  ylab("Perplexity (sc)") +
  scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = 19)) +
  geom_point(aes(x = mean_coherence, y = perplexity, color = tissue,
                 shape = p_coherence_adj < .05),
             size = 3) +
  ggrepel::geom_text_repel(aes(x = mean_coherence, y = perplexity, label = cell_type))


# ggsave("phasic_cell_types_annot.png", path = "presentations/figures/250612_celltype_osc",
#        width = 80, height = 50, units = "mm",
#        scale = 2)
# ggsave("phasic_cell_types_annot.pdf", path = "presentations/figures/250612_celltype_osc",
#        width = 80, height = 50, units = "mm",
#        scale = 2)








# Save ----



# by_cell_type |>
#   qs::qsave(file.path(dir_step3, "cell_types_sc.qs"))


# cell_types_both |>
#   qs::qsave(file.path(dir_step3, "cell_types.qs"))






iwalk(heatmaps_list,
      \(.hm, .ct){
        
        if(ncol(.hm) == 0) return()
        
        hm_norm <- apply(.hm, 2, \(x){
          (x - min(x)) / (max(x) - min(x))
        })
          
        squash::savemat(t(hm_norm)[, nrow(hm_norm):1],
                        filename = file.path(dir_step3, "heatmaps_cts",
                                             paste0(.ct, "_heatmap.png")))
          
        })


iwalk(heatmaps_list,
      ~ qs::qsave(.x,
                  file.path(dir_step3, "heatmaps_cts",
                            paste0(.y, "_heatmap.qs")))
)

# heatmaps_list <- map(filtered_data$cell_type |> set_names(),
#                      ~ qs::qread(
#                   file.path(dir_step3,
#                             paste0(.x, "_heatmap.qs"))
#                   )
# )


# # example null matrix
# null_mat <- sapply(seq_len(800),
#                    \(i){
#                      circ_perm(sig_cent, floor(200 * i/800))
#                    })
# squash::savemat(t(null_mat)[, nrow(null_mat):1],
#                 filename = file.path(dir_step3,
#                                      paste0("null", "_heatmap.png")))




# Save prettier heatmaps



# hmp_sparsified <- heatmaps_list[[9]]
# ct <- names(heatmaps_list)[[9]]

iwalk(heatmaps_list[c("gonad_1", "ILso")],
      \(hmp_sparsified, ct){
        
        if(ncol(hmp_sparsified) <= 200){
          message("skipping ", ct)
          return()
        }
        
        if(is_null(colnames(hmp_sparsified))){
          message("saved not sparsified: ", ct)
          pheatmap::pheatmap(hmp_sparsified,
                             cluster_rows = FALSE,
                             cluster_cols = FALSE,
                             show_rownames = FALSE,
                             show_colnames = FALSE,
                             fontsize = 7,
                             filename = paste0("presentations/figures/250612_celltype_osc/heatmap_", ct, ".png"),
                             width = 4,
                             height = 2,
                             main = ct)
          
          return()
        }
        
        hmp_sparsified <- apply(hmp_sparsified, 2, \(x){
          (x - min(x)) / (max(x) - min(x))
        })
        
        
        colnames(hmp_sparsified) <- str_split_i(colnames(hmp_sparsified), fixed("|"), 2)
        
        spar_index <- round(ncol(hmp_sparsified)/40)
        colnames_to_sparsify <- setdiff(seq_len(ncol(hmp_sparsified)),
                                        spar_index * seq_len( ncol(hmp_sparsified) / spar_index ) )
        
        colnames(hmp_sparsified)[colnames_to_sparsify] <- ""
        head(colnames(hmp_sparsified), 20)
        
        
        pheatmap::pheatmap(hmp_sparsified,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           show_rownames = FALSE,
                           fontsize = 6,
                           filename = paste0("presentations/figures/250612_celltype_osc/heatmap_", ct, ".png"),
                           width = 5,
                           height = 2.5,
                           main = ct)
        
        pheatmap::pheatmap(hmp_sparsified,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           show_rownames = FALSE,
                           fontsize = 6,
                           filename = paste0("presentations/figures/250612_celltype_osc/heatmap_", ct, ".pdf"),
                           width = 5,
                           height = 2.5,
                           main = ct)
        
        message("saved: ", ct)
      })

# dev.off()




















