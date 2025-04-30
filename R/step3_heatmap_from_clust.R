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


dir_step1 <- "intermediates/2502/250330_step1/"
dir_step2_clust <- "intermediates/2502/250422_step2_nb_centered/"
dir_step2_noncentered <- "intermediates/2502/250422_step2_nb/"

dir_step3 <- "intermediates/2502/250425_step3_genes_by_celltype/"








cluster_results <- read_csv(file.path(dir_step2_clust, "250424_cluster_results.csv")) |>
  mutate(cellgene = paste0(cell_type, "|", gene_name))

by_cell_type <- cluster_results |>
  summarize(n_puls = sum(shape == "pulsatile"),
            n_tot = n(),
            .by = "cell_type") |>
  mutate(prop_puls = round( 100 * n_puls / n_tot )) |>
  arrange(desc(prop_puls))



all_preds <- list.files(dir_step2_noncentered,
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
    mat_preds <- qs::qread(file.path(dir_step2_noncentered,
                                     filename))
    colnames(mat_preds) <- paste0(cell_type, "|", colnames(mat_preds))
    
    mat_preds
  }) |>
  do.call(cbind, args = _)


stopifnot(all.equal(
  colnames(all_preds),
  cluster_results$cellgene
))




# heatmap ----

# printMat::matimage(log1p(all_preds))


# create a heatmap for each cell type

heatmaps_list <- by_cell_type$cell_type |>
  set_names() |>
  map(\(.ct){
    
    to_keep <- cluster_results$cell_type == .ct & cluster_results$shape == "pulsatile"
    
    preds <- all_preds[ , cluster_results$cellgene[to_keep] ]
    
    peak_location <- apply(preds, 2, which.max)
    
    preds[, order(peak_location)]
  })


# sapply(heatmaps_list, dim)





#~ uniformity ----


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



# dists with pvals

by_cell_type$p_val <- map_dbl(
  heatmaps_list,
  \(hm){
    
    if(ncol(hm) == 0L) return(1)
    
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
    
    
    
    similarity_perms <- c(
      {
        dist <- norm(hm_norm - null_mat, "1")
        similarity <- dist/n_t
        1 - similarity
      },
      replicate(10000,{
        hm_perm <- hm_norm[,sample(ncol(hm_norm))]
        dist <- norm(hm_perm - null_mat, "1")
        similarity <- dist/n_t
        1 - similarity
      })
    )
    
    pval <- mean(similarity_perms >= similarity_perms[[1]])
    pval
  },
  .progress = TRUE
)


by_cell_type$p_adj <- by_cell_type$p_val |>
  p.adjust(method = "BH")
hist(by_cell_type$p_val, breaks = 30)
hist(by_cell_type$p_adj, breaks = 30)



ggplot(by_cell_type) +
  theme_classic() +
  aes(x = uniformity_index, y = similarity_diag, label = cell_type, shape = p_adj < .05, color = p_adj < .05) +
  geom_point() +
  ggrepel::geom_text_repel()



by_cell_type |>
  filter(prop_puls >= 10) |>
  ggplot() +
  theme_classic() +
  aes(x = uniformity_index, y = similarity_diag, label = cell_type, shape = p_adj < .05, color = p_adj < .05) +
  geom_point() +
  ggrepel::geom_text_repel()


ggplot(by_cell_type) +
  theme_classic() +
  aes(x = similarity_diag, y = -log10(p_adj),
      label = cell_type, shape = p_adj < .05, color = p_adj < .05) +
  geom_point() +
  ggrepel::geom_text_repel()

by_cell_type |>
  filter(prop_puls >= 10) |>
  ggplot() +
  theme_classic() +
  aes(x = similarity_diag, y = -log10(p_adj),
      label = cell_type, shape = p_adj < .05, color = p_adj < .05) +
  geom_point() +
  ggrepel::geom_text_repel()


# # only label every n-th column
# ct <- "BWM"
# hmp_sparsified <- heatmaps_list[[ct]]
# hmp_sparsified <- null_mat
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

dir_assembled <- "intermediates/2502/250328_assembled"
mean_dotprod_by_celltype_res_perm <- qs::qread(file.path(dir_assembled, "250330_coherence_perm10000.qs"))
dotprod_by_cell <- qs::qread(file.path(dir_assembled, "250329_dotprod_by_cell.qs"))
# dotprod_by_cell$tissue[dotprod_by_cell$cell_type == "pharyngeal"] <- "neuron"

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

inner_join(
  by_cell_type |>
    select(cell_type, similarity_diag, p_diag = p_adj)
  ,
  cell_types_bulk |>
    select(cell_type, tissue, mean_coherence, p_coherence = p_adj)
) |>
  mutate(signif = case_when(
    p_diag < .05 & p_coherence < .05 ~ "*",
    p_diag < .05 & p_coherence >= .05 ~ "#",
    p_diag >= .05 & p_coherence < .05 ~ "$",
    p_diag >= .05 & p_coherence >= .05 ~ "o"
  )) |>
  ggplot() +
  theme_classic() +
  xlab("Mean local phase coherence") +
  ylab("Diagonal similarity") +
  scale_shape_manual(values = c("*" = 8, "#" = 7, "$" = 9, "o" = 16)) +
  geom_point(aes(x = mean_coherence, y = similarity_diag, color = tissue,
                 shape = signif),
             size = 3)


inner_join(
  by_cell_type |>
    select(cell_type, similarity_diag, p_diag = p_adj)
  ,
  cell_types_bulk |>
    select(cell_type, tissue, mean_coherence, p_coherence = p_adj)
) |>
  mutate(signif = case_when(
    p_diag < .05 & p_coherence < .05 ~ "both",
    p_diag < .05 & p_coherence >= .05 ~ "diag only",
    p_diag >= .05 & p_coherence < .05 ~ "bulk only",
    p_diag >= .05 & p_coherence >= .05 ~ "neither"
  )) |>
  ggplot() +
  theme_classic() +
  xlab("Mean local phase coherence (bulk)") +
  ylab("Diagonal similarity (sc)") +
  scale_shape_manual(values = c("both" = 8, "diag only" = 7, "bulk only" = 9, "neither" = 16)) +
  geom_point(aes(x = mean_coherence, y = similarity_diag, color = tissue,
                 shape = signif),
             size = 3) +
  ggrepel::geom_text_repel(aes(x = mean_coherence, y = similarity_diag, label = cell_type))


# ggsave("phasic_cell_types.pdf", path = "presentations/",
#        width = 60, height = 50, units = "mm",
#        scale = 2)








# Save ----



by_cell_type |>
  qs::qsave(file.path(dir_step3, "cell_types_sc.qs"))



inner_join(
  by_cell_type |>
    select(cell_type, n_puls, n_tot, prop_puls, similarity_diag, p_diag_adj = p_adj)
  ,
  cell_types_bulk |>
    select(cell_type, tissue, mean_coherence, p_coherence_adj = p_adj)
) |>
  qs::qsave(file.path(dir_step3, "cell_types.qs"))






iwalk(heatmaps_list,
      \(.hm, .ct){
        
        if(ncol(.hm) == 0)
        
        hm_norm <- apply(.hm, 2, \(x){
          (x - min(x)) / (max(x) - min(x))
        })
          
        squash::savemat(t(hm_norm)[, nrow(hm_norm):1],
                        filename = file.path(dir_step3,
                                             paste0(.ct, "_heatmap.png")))
          
        })


iwalk(heatmaps_list,
      ~ qs::qsave(.x,
                  file.path(dir_step3,
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




# Save pretty heatmaps
sig_cent <- make_ref_sig(200)
null_mat <- sapply(seq_len(800),
                   \(i){
                     circ_perm(sig_cent, floor(200 * i/800))
                   })

heatmaps_list <- append(heatmaps_list, list(null_mat = null_mat))

# hmp_sparsified <- heatmaps_list[[9]]
# ct <- names(heatmaps_list)[[9]]

iwalk(heatmaps_list,
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
                             filename = paste0("presentations/250411_heatmaps/", ct, ".png"),
                             width = 4,
                             height = 1.5,
                             main = ct)
          
          return()
        }
        
        
        spar_index <- round(ncol(hmp_sparsified)/40)
        colnames_to_sparsify <- setdiff(seq_len(ncol(hmp_sparsified)),
                                        spar_index * seq_len( ncol(hmp_sparsified) / spar_index ) )
        
        colnames(hmp_sparsified)[colnames_to_sparsify] <- ""
        head(colnames(hmp_sparsified), 20)
        
        
        pheatmap::pheatmap(hmp_sparsified,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           show_rownames = FALSE,
                           fontsize = 7,
                           filename = paste0("presentations/250411_heatmaps/", ct, ".png"),
                           width = 4,
                           height = 1.5,
                           main = ct)
        
        message("saved: ", ct)
      })





















