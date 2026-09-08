# Perplexity saturation analysis
# For each cell type, subsample cells at increasing sizes,
# recompute pseudotime and GAM fits for pulsatile genes only,
# then compute perplexity.

# Inits ----
library(tidyverse)
library(Seurat)
library(ElPiGraph.R)

source("R/utils_fit.R") # --> z03_utils_fit.R (for circ_perm_mat())

source("https://raw.githubusercontent.com/AlexWeinreb/knn-smoothing/refs/heads/master/knn_smooth.R") # -> knn_smoothing


set.seed(123)

# Params ----
dir_assembled <- "intermediates/2502/250605_assembled/"
dir_clust <- "intermediates/2502/250624_cluster"
dir_out <- "intermediates/2502/260527_saturation"
# dir.create(dir_out)

n_reps <- 20
len <- 128
n_rep_pt <- 10

# subsample_sizes <- c(50, 100, 200)
subsample_sizes <- c(25, 50, 75, 100, 150, 200, 300, 500, 750, 1000)


# cell_types_of_interest <- c("glia_4", "glia_1", "glia_sheath_2", "early_gonad")

  # cell_types_of_interest <- c(
  #   "ILso", "seam",
  #   "glia_4", "glia_1", "glia_sheath_2", "early_gonad", "PHsh",
  # "BWM",  "AM_PHso"
  # )

cell_types_of_interest <- c(
  "coelomocyte"
)


# Load data ----

seu_all <- qs::qread(file.path(dir_assembled, "250606_seu_all_herma.qs"))

seu_all$cell_type[seu_all$cell_type == "coelomyocyte"] <- "coelomocyte"



pulsatile_by_ct <- read_csv(file.path(dir_clust, "250624_cluster_results.csv"),
                            show_col_types = FALSE) |>
  filter(shape == "pulsatile") |>
  summarize(genes = list(gene_name), .by = cell_type) |>
  deframe() |>
  set_names(~ if_else(.x == "coelomyocyte", "coelomocyte", .x))



# Run pipeline on single subsample
compute_perplexity_from_cells <- function(seu_ct, puls_genes_ct, n_cells_sub, len, n_rep_pt) {
  
  
  # subsample
  if (n_cells_sub >= ncol(seu_ct)) {
    sub <- seu_ct
  } else {
    sub <- subset(seu_ct, cells = sample(colnames(seu_ct), n_cells_sub), seed = NULL)
  }
  
  
  cnts_raw <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  
  # KNN smoothing
  cnts_smoothed <- knn_smoothing(cnts_raw, k = 5, seed = NULL, verbose = FALSE)
  
  
  # Create smoothed Seurat, SCTransform, PCA
  sub_smooth <- cnts_smoothed |>
    as("dgCMatrix") |>
    CreateSeuratObject(meta.data = sub[[]]) |>
    SCTransform(verbose = FALSE, seed.use = NULL) |>
    RunPCA(npcs = 2, verbose = FALSE, seed.use = NULL)
  
  
  embding <- FetchData(sub_smooth, vars = c("PC_1", "PC_2")) |> as.matrix()
  
  
  # ElPiGraph pseudotime
  capture.output(
    CurveEPG <- computeElasticPrincipalCircle(
      X = embding, NumNodes = 20, Do_PCA = FALSE,
      nReps = n_rep_pt, drawPCAView = FALSE,
      drawAccuracyComplexity = FALSE, drawEnergy = FALSE,
      ProbPoint = 0.6
    ),
    file = nullfile()
  )
  
  
  targt <- CurveEPG[[n_rep_pt + 1L]]
  
  
  
  capture.output(
    subgraph <- GetSubGraph(
      Net = ConstructGraph(targt),
      Structure = "circle"
    )[["Circle_1"]],
    file = nullfile()
  )
  
  PartStruct <- PartitionData(X = embding, NodePositions = targt$NodePositions)
  ProjStruct <- project_point_onto_graph(
    X = embding, NodePositions = targt$NodePositions,
    Edges = targt$Edges$Edges, Partition = PartStruct$Partition
  )
  pseudotime <- getPseudotime(ProjStruct = ProjStruct, NodeSeq = names(subgraph))[["Pt"]]
  
  
  
  # GAM fits for pulsatile genes only
  # Use unsmoothed counts for GAM (as in original step 2)
  mat_cnt_unsmooth <- sub |>
    GetAssayData(assay = "RNA", layer = "counts")
  
  
  
  puls_genes <- rownames(mat_cnt_unsmooth)[rowSums(mat_cnt_unsmooth > 0) >= 5] |>
    intersect(puls_genes_ct)
  
  
  
  
  if (length(puls_genes) < 5) return(NA_real_)
  
  mat_cnt_puls <- mat_cnt_unsmooth[puls_genes, , drop = FALSE]
  
  nf <- edgeR::calcNormFactors(mat_cnt_unsmooth)
  size_factors <- colSums(mat_cnt_unsmooth) * nf
  mean_sf <- exp(mean(log(size_factors)))
  
  
  # Fit uncentered GAMs to get peak positions
  peak_positions <- vapply(puls_genes, \(.gene) {
    dat <- data.frame(
      expr = mat_cnt_puls[.gene, ],
      pseudotime = pseudotime / max(pseudotime)
    )
    
    mod <- tryCatch(
      mgcv::gam(
        expr ~ s(pseudotime, k = 6, bs = "cc") + offset(log(size_factors)),
        data = dat,
        family = mgcv::nb(link = "log")
      ),
      error = \(e) NULL
    )
    
    if (is.null(mod)) return(NA_real_)
    
    
    preds <- predict(mod, type = "response",
                     newdata = data.frame(
                       pseudotime = (0:(len - 1)) / len,
                       size_factors = rep(mean_sf, len)
                     ))
    which.max(preds) / len
  }, FUN.VALUE = double(1L))
  
  
  
  if(all(is.na(peak_positions))) return(NA_real_)
  
  stopifnot( !anyNA(peak_positions) )
  
  # peak_positions <- peak_positions[!is.na(peak_positions)]
  if (length(peak_positions) < 5) return(NA_real_)
  
  # Compute perplexity
  angles <- peak_positions * 2 * pi
  bins <- seq(0, 2 * pi, length.out = 51)
  counts <- hist(angles, breaks = bins, plot = FALSE)$counts
  exp(entropy::entropy(counts))
}



# Main loop ----
cell_types <- names(pulsatile_by_ct) |>
  intersect(cell_types_of_interest)



results_all <- list()

for (ct in cell_types) {
  message("===== ", ct, "(out of ", paste(cell_types, collapse = ", "), ")", " =====")
  
  seu_ct <- subset(seu_all, cell_type == ct)
  n_cells <- ncol(seu_ct)
  pulsatile_genes_ct <- pulsatile_by_ct[[ct]]
  
  message("  ", n_cells, " cells, ", length(pulsatile_genes_ct), " pulsatile genes")
  
  if (n_cells < 25) {
    message("  Too few cells, skipping")
    next
  }
  
  sizes_to_use <- subsample_sizes[subsample_sizes <= n_cells] |>
    c(n_cells) |>
    unique()
  
  
  ct_results <- list()
  
  for (sz in sizes_to_use) {
    n_reps_here <- if (sz == n_cells) 1L else n_reps
    
    message("  n=", sz, " (", n_reps_here, " reps)")
    
    perps <- vapply(
      seq_len(n_reps_here),
      \(rep_i) {
        message("    rep ", rep_i)
        # tryCatch(
        compute_perplexity_from_cells(
          seu_ct = seu_ct,
          puls_genes_ct = pulsatile_genes_ct,
          n_cells_sub = sz,
          len = len,
          n_rep_pt = n_rep_pt
        )
        # error = \(e) { message("    rep ", rep_i, " failed: ", e$message); NA_real_ })
      }, FUN.VALUE = double(1L))
    
    ct_results[[length(ct_results) + 1]] <- tibble(
      cell_type = ct,
      n_total_cells = n_cells,
      n_subsample = sz,
      replicate = seq_len(n_reps_here),
      perplexity = perps
    )
  }
  
  results_all[[ct]] <- bind_rows(ct_results)
  
  # Save incrementally
  qs::qsave(results_all[[ct]],
            file.path(dir_out, paste0(ct, "_saturation.qs")))
  
  message("  Done")
}


# Combine and save ----
results_all <- list.files("intermediates/2502/260527_saturation/",
                           pattern = "_saturation\\.qs",
                           full.names = TRUE) |>
  set_names(~ basename(.x) |> str_remove("_saturation\\.qs$")) |>
  map(qs::qread)

results_df <- bind_rows(results_all)
write_csv(results_df, file.path(dir_out, "saturation_results.csv"))



# Plot ----
results_summary <- results_df |>
  summarize(
    mean_perp = mean(perplexity, na.rm = TRUE),
    sd_perp = sd(perplexity, na.rm = TRUE),
    .by = c(cell_type, n_total_cells, n_subsample)
  ) |>
  mutate(sd_perp = replace_na(sd_perp, 0))

results_summary |>
  ggplot(aes(x = n_subsample, y = mean_perp)) +
  theme_classic() +
  geom_ribbon(aes(ymin = mean_perp - sd_perp,
                  ymax = mean_perp + sd_perp),
              alpha = 0.2) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~cell_type) +
  xlab("Number of cells subsampled") +
  ylab("Perplexity") +
  theme(strip.text = element_text(size = 7)) +
  scale_x_log10()


# > +++ Fig. EV5C +++ ----


results_summary |>
  filter(cell_type %in% c("ILso", "seam", "coelomocyte", "BWM", "glia_1","glia_4")) |>
  ggplot(aes(x = n_subsample, y = mean_perp)) +
  theme_classic() +
  geom_ribbon(aes(ymin = mean_perp - sd_perp,
                  ymax = mean_perp + sd_perp,
                  fill = cell_type),
              alpha = 0.2) +
  geom_line(aes(color = cell_type)) +
  geom_point(aes(color = cell_type),
             size = 1) +
  xlab("Number of cells subsampled") +
  ylab("Perplexity") +
  theme(strip.text = element_text(size = 7)) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0,NA)) +
  # theme(legend.position = "inside",
  #       legend.position.inside = c(.8,.3)) +
  theme(legend.position = "none") +
  geom_text(aes(label = cell_type),
            data = results_summary |>
              slice_tail(n = 1, by = "cell_type") |>
              filter(cell_type %in% c("ILso", "seam", "coelomocyte", "BWM", "glia_1","glia_4")),
            hjust = -.2)

# ggsave(filename = "260527_subsampled_perplexity.pdf",
#        path = "presentations/figures/260527_subsampled_perplexity",
#        width = 65, height = 55, units = "mm",
#        scale = 2)



results_df |>
  filter(cell_type == "glia_sheath_2") |>
  ggplot(aes(x = n_subsample, y = perplexity, color = as.factor(replicate), label = replicate)) +
  theme_classic() +
  ggbeeswarm::geom_quasirandom() +
  # ggrepel::geom_text_repel() +
  scale_x_log10()




