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

dir_clust <- "intermediates/2502/250624_cluster"

# dir_step2 <- "intermediates/2502/250624_step2/"
# dir_step2 <- "E:/2025-06-27/Projects/glia/osc_larva/intermediates/2502/250624_step2/"
dir_step2 <- "D:/2025-06-27/Projects/glia/osc_larva/intermediates/2502/250624_step2/"


dir_step3 <- "intermediates/2502/250624_step3_genes_by_celltype/"
# dir.create(dir_step3)

dir_figures3 <- "presentations/figures/250624_celltype_osc"
# dir.create(dir_figures3)


# Load ----
cluster_results <- read_csv(file.path(dir_clust, "250624_cluster_results.csv"))

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

## tests and explorations
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
# printMat::matimage(mat)
# points(seq_len(ncol(mat)) / ncol(mat) ,
#        1 - apply(mat, 2, which.max) / nrow(mat),
#        col = 'purple', cex = 2, pch = "-")
# 
# hist(apply(mat, 2, which.max)/nrow(mat), breaks = 50)
# hist(runif(ncol(mat)), breaks = 50, add = TRUE, col = alpha('lightgreen', .2))
# 
# qqplot(qunif(ppoints(ncol(mat))), apply(mat, 2, which.max)/nrow(mat))
# abline(0, 1, col = "red")
# 
# 
# circular::rose.diag(circular::circular(runif(10000)*2*pi),
#                     bins = 24,
#                     axes = FALSE,
#                     col = "lightblue", border = "lightblue",
#                     prop = 1.5)
# 
# circular::rose.diag(circular::circular(pos_peaks*2*pi),
#                     bins = 24,
#                     axes = FALSE,col = "grey",prop = 1.5,
#                     add = TRUE)



ct <- "ILso"
ct <- "BWM"
mat <- heatmaps_scaled[[ct]]
pos_peaks <- apply(mat, 2, which.max) / nrow(mat)



nbins <- 30

enframe(pos_peaks,
        name = "genecell",
        value = "peak_pt") |>
  mutate(peak_deg = peak_pt) |>
  ggplot() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 10),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  coord_polar() +
  xlab(NULL) +
  scale_y_sqrt() +
  scale_x_continuous(breaks = c(0, .25, .5, .75)) +
  geom_ribbon(aes(x = seq(from = 0, to = 1, length.out = length(peak_deg)),
                  ymin = 0,
                  ymax = length(peak_deg) / nbins),
              fill = "purple3",
              alpha = .35) +
  geom_hline(
    aes(yintercept = length(peak_deg) / nbins),
    linewidth = 1,
    color = "purple4",
    alpha = .6
  ) +
  geom_histogram(aes(x = peak_deg),
                 color = "black",
                 alpha = .9,
                 breaks = seq(0, 1, length.out = nbins + 1))

# ggsave(paste0("histogram_circ_",ct,".pdf"),
#        path = dir_figures3,
#        width = 35, height = 35, units = "mm")





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


# for annotation
circ_entropy(2*pi * (apply(heatmaps_list[["ILso"]], 2, which.max) / 128))
circ_entropy(2*pi * (apply(heatmaps_list[["BWM"]], 2, which.max) / 128))
by_cell_type |> filter(cell_type %in% c("ILso", "BWM"))





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

cell_types_bulk <- qs::qread(file.path(dir_assembled, "250610_coherence_by_ct.qs"))

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
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 10),
    plot.margin = unit(c(0,0,0,0), "mm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.margin = margin(),
    legend.box.margin = margin(),
    legend.key.size = unit(1, "mm")
  ) +
  xlab("Mean local phase coherence") +
  ylab("Perplexity") +
  scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = 19),
                     name = expression(p[adj] < .05)) +
  geom_point(aes(x = mean_coherence, y = perplexity, color = tissue,
                 shape = p_coherence_adj < .05),
             size = 1)


# ggsave("phasic_cell_types_unannot.pdf", path = dir_figures3,
#        width = 90, height = 60, units = "mm")


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



# ggsave("phasic_cell_types_annot.pdf", path = dir_figures3,
#        width = 80, height = 50, units = "mm",
#        scale = 2)




# Color by number of genes (power)

cell_types_both |>
  ggplot() +
  theme_classic() +
  # xlab("Mean local phase coherence (bulk)") +
  ylab("Perplexity (sc)") +
  scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = 19)) +
  geom_point(aes(x = n_tot, y = perplexity, color = n_tot,
                 shape = p_coherence_adj < .05),
             size = 3) +
  geom_label(aes(x = n_tot, y = perplexity, label = cell_type),
             data = cell_types_both |> filter(cell_type %in% c(
               "ILso", "seam",
               "glia_4", "glia_1", "glia_sheath_2", "early_gonad",
               "coelomocyte", "PHsh"
             )))


# by number of cells
seu_all <- qs::qread(file.path(dir_assembled, "250606_seu_all_herma.qs"))

by_cell_type |>
  left_join(tibble(cell_type = seu_all$cell_type) |> count(cell_type, name = "n_cells"),
            by = "cell_type") |>
  ggplot() +
  theme_classic() +
  # xlab("Mean local phase coherence (bulk)") +
  ylab("Perplexity (sc)") +
  scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = 19)) +
  geom_point(aes(x = n_cells, y = perplexity, color = n_tot),
             size = 3) +
  geom_text(aes(x = n_cells, y = perplexity, label = cell_type,
                 size = cell_type %in% c(
               "ILso", "seam",
               "glia_4", "glia_1", "glia_sheath_2", "early_gonad",
               "coelomocyte", "PHsh"
             )),
             show.legend = FALSE) +
  scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 0))

by_cell_type |>
  left_join(tibble(cell_type = seu_all$cell_type) |> count(cell_type, name = "n_cells"),
            by = "cell_type") |>
  filter(n_cells > 1000, perplexity > 30)



# Save ----



# by_cell_type |>
#   qs::qsave(file.path(dir_step3, "cell_types_sc.qs"))
# 
# cell_types_both |>
#   qs::qsave(file.path(dir_step3, "cell_types.qs"))


# For supp table
# cell_types_both |>
#   write_csv(file.path(dir_figures3, "cell_types_osc.csv"))



# dir.create(file.path(dir_step3, "heatmaps_cts"))

iwalk(heatmaps_list,
      \(.hm, .ct){
        
        if(ncol(.hm) == 0) return()
        
        hm_norm <- apply(.hm, 2, \(x){
          (x - min(x)) / (max(x) - min(x))
        })
          
        squash::savemat(
          t(hm_norm)[, nrow(hm_norm):1],
          filename = file.path(dir_step3, "heatmaps_cts",
                               paste0(.ct, "_heatmap.png"))
        )
          
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

mm_to_in <- 0.03937008

iwalk(heatmaps_list[c("BWM", "ILso")],
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
        
        spar_index <- round(ncol(hmp_sparsified)/35)
        colnames_to_sparsify <- setdiff(seq_len(ncol(hmp_sparsified)),
                                        spar_index * seq_len( ncol(hmp_sparsified) / spar_index ) )
        
        colnames(hmp_sparsified)[colnames_to_sparsify] <- ""
        
        
        
        png(paste0(dir_figures3, "/heatmap_", ct, ".png"),
            width = 80, height = 50, units = "mm", res = 500)
        pheatmap::pheatmap(hmp_sparsified,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           show_rownames = FALSE,
                           fontsize = 6,
                           # filename = paste0(dir_figures3, "/heatmap_", ct, ".png"),
                           # width = 4,
                           # height = 1.2,
                           main = ct)
        dev.off()
        
        
        pheatmap::pheatmap(hmp_sparsified,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           show_rownames = FALSE,
                           fontsize = 5,
                           filename = paste0(dir_figures3, "/heatmap_", ct, ".pdf"),
                           width = 80 * mm_to_in,
                           height = 50 * mm_to_in,
                           main = ct)
        
        
        
        message("saved: ", ct)
      })

# dev.off()




















