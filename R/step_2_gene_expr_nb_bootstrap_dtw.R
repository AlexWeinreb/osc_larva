# Second part ----

# Measure dtw distance to a narrow peak. Bootstraps to get CI


# Inits ----
library(ggplot2)
suppressPackageStartupMessages( library(Seurat) )
library("ElPiGraph.R")


library(getopt)


# center signal so the peak ends up in the middle
# note might be off by one
circ_perm <- function(preds){
  
  peak_loc <- apply(preds, 2, which.max)
  
  preds_permuted <- matrix(nrow = nrow(preds), ncol = ncol(preds))
  t_tot <- nrow(preds)
  t_cent <- floor(t_tot/2)
  
  for(i in seq_len(ncol(preds))){
    
    if(peak_loc[[i]] == t_cent){
      
      preds_permuted[, i] <- preds[,i]
      
    } else if(peak_loc[[i]] == t_tot){
      
      preds_permuted[, i] <- c(preds[ (t_cent+1):t_tot, i],
                               preds[ 1:t_cent ,i])
      
    } else if(peak_loc[[i]] < t_cent){
      
      #  1     a      t_c     b      t_t
      #  |_____|_______|______|______|
      
      
      a <- peak_loc[[i]]
      b <- t_cent + a
      
      preds_permuted[, i] <- c(preds[ (b+1):t_tot, i],
                               preds[ 1:b, i])
      
    } else{
      
      #  1     a      t_c     b      t_t
      #  |_____|_______|______|______|
      
      
      b <- peak_loc[[i]]
      a <- 1 + b - (t_cent+1)
      
      
      preds_permuted[, i] <- c(preds[ (a+1):b, i],
                               preds[ (b+1):t_tot, i],
                               preds[ 1:a, i])
    }
    
  }
  
  colnames(preds_permuted) <- colnames(preds)
  rownames(preds_permuted) <- rownames(preds)
  
  preds_permuted
}


if(! interactive()){
  spec <- matrix(c(
    'batch_rmed_dir', 'f', 1, 'character',
    'out_dir', 'o', 1, 'character',
    'model', 'm', 1, 'character',
    'i', 'i', 1, 'integer',
    'prop_thres', 'g', 1, 'double',
    'cnt_thres', 'c', 1, 'integer',
    'nb_subsamples_ptDE', 'n', 1, 'integer'
  ), byrow=TRUE, ncol=4)
  
  params <- getopt(spec)
  
} else{
  # Options for interactive
  params <- list(
    batch_rmed_dir = "intermediates/2502/250330_step1",
    out_dir = "intermediates/2502/250401_step2_boot",
    model = "auto",
    i = 3,
    prop_thres = 0.05,
    cnt_thres = 20,
    nb_subsamples_ptDE = 10
  )
}





n_rep_pt_global <- 50



set.seed(123)


cell_types <- list.files(params$batch_rmed_dir,
                         pattern = "_seu\\.qs$") |>
  stringr::str_remove("_seu\\.qs$")





# loop ----



cell_type <- cell_types[[params$i]]

cell_type


message(params$i,"/", length(cell_types), ": ", cell_type)

log_content <- paste0(cell_type, "\n")


message("---- load")

gene_expressions <- qs::qread(
  file.path(params$batch_rmed_dir,
            paste0(cell_type, "_gene_expressions.qs"))
)

seu_ref <- qs::qread(
  file.path(params$batch_rmed_dir,
            paste0(cell_type, "_seu.qs"))
)


# save the plot with phases
gg_phase <- FetchData(seu_ref,
                      vars = c("PC_1", "PC_2", "cell_phase_masked", "cell_rho")) |>
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 360)) +
  geom_point(aes(x = PC_1, y = PC_2,
                 color = cell_phase_masked, alpha = cell_rho)) +
  ggtitle(cell_type)

gg_phase

ggsave(paste0(cell_type,"_phase.png"), gg_phase,
       path = params$out_dir,
       width = 7, height = 5, units = "in")









#

# Pick genes ----
pos_high_genes <- which(
  (gene_expressions$prop_cells >= params$prop_thres) &
    (gene_expressions$nb_cells >= params$cnt_thres)
)

high_genes <- gene_expressions$gene_name[pos_high_genes]

length(high_genes)






run_once <- function(mat, i){
  
  message("Next bootstrap ##########")
  
  mat_resample <- t(mat)[,i]
  colnames(mat_resample) <- seq_along(colnames(mat_resample))
  
  
  message("   -- Normalize")
  
  seu <- CreateSeuratObject(counts = mat_resample) |>
    SCTransform(verbose = FALSE,
                seed.use = NULL) |>
    RunPCA(npcs = 2,
           verbose = FALSE,
           seed.use = NULL)
  
  
  embding <- seu |>
    FetchData(vars = c("PC_1", "PC_2")) |>
    as.matrix()
  
  stopifnot(all.equal( rownames(embding), rownames(seu[[]]) ))
  
  
  #~ pseudotime ----
  
  message("   -- pseudotime")
  capture.output(
    CurveEPG <- computeElasticPrincipalCircle(X = embding,
                                              NumNodes = 20,
                                              Do_PCA = FALSE,
                                              nReps = n_rep_pt_global,
                                              drawPCAView = FALSE,
                                              drawAccuracyComplexity = FALSE,
                                              drawEnergy = FALSE,
                                              ProbPoint = .6),
    file = nullfile()
  )
  
  stopifnot( length(CurveEPG) == n_rep_pt_global+1L )
  targt <- CurveEPG[[ n_rep_pt_global+1L ]]
  boots <- CurveEPG[ seq_len(n_rep_pt_global) ]
  
  
  gg_int_cca_circlepath <- PlotPG(embding,
                                  TargetPG = targt, BootPG = boots,
                                  Do_PCA = FALSE,VizMode = c("Target","Boot"),
                                  DimToPlot = 1:2)[[1]] +
    theme_classic() +
    scale_linewidth_manual(values = c(1,2,3))
  
  # gg_int_cca_circlepath 
  
  ggsave(paste0(cell_type,"_circlepath.png"), gg_int_cca_circlepath,
         path = params$out_dir,
         width = 6, height = 5, units = "in", dpi = 600)
  
  
  capture.output(
    subgraph <- GetSubGraph(Net = ConstructGraph(targt),
                            Structure = 'circle')[[ "Circle_1" ]],
    file = nullfile()
  )
  
  
  PartStruct <- PartitionData(X = embding,
                              NodePositions = targt$NodePositions)
  
  ProjStruct <- project_point_onto_graph(X = embding,
                                         NodePositions = targt$NodePositions,
                                         Edges = targt$Edges$Edges,
                                         Partition = PartStruct$Partition)
  
  
  pseudotime <- getPseudotime(ProjStruct = ProjStruct, NodeSeq = names(subgraph))[[ "Pt" ]]
  
  
  
  
  #~  GAM ----
  
  message("   -- fit GAM")
  mat_cnt <- GetAssayData(seu, assay = "RNA", layer = "count")[high_genes,]
  
  mods <- lapply(high_genes,
                 \(.gene){
                   dat <- data.frame(
                     cnt = mat_cnt[.gene,],
                     pseudotime = pseudotime/max(pseudotime)
                   )
                   
                   tryCatch({
                     mgcv::gam(cnt ~ s(pseudotime, k = 6, bs = 'cc'),
                               data = dat, family = mgcv::nb(link = "log"))
                   }, error = \(e) NA
                   )
                 }) |>
    setNames(high_genes)
  
  successful_fits <- which(!is.na(mods))
  mods <- mods[successful_fits]
  
  message("   -- get distance")
  
  len <- 100
  interp_pt <- (seq_len(len) - 1)/len
  
  all_preds <- vapply(mods,
                      \(.mod){
                        
                        predict(.mod,
                                type = "response",
                                newdata = data.frame(pseudotime = interp_pt ))
                        },
                      FUN.VALUE = double(len))
  
  all_preds_centered <- circ_perm(all_preds) |>
    apply(2, \(x) x/max(x))
  
  
  
  
  ref <- dnorm(interp_pt, mean = .5, sd = .01)
  ref <- (ref - min(ref))/max(ref - min(ref))
  ref <- circ_perm(matrix(ref, ncol = 1))
  
  
  # matplot(interp_pt, all_preds_centered, type = 'l',
  #         xlab = "Pseudotime",
  #         ylab = "log(count + 1)",
  #         lwd = .7)
  # lines(interp_pt, ref, lwd = 2.5)
  
  dist_dtw <- dtw::dtwDist(t(all_preds_centered), t(ref)) |> as.numeric()
  
  dist_with_na <- rep(NA_real_, length(high_genes))
  dist_with_na[successful_fits] <- dist_dtw
  
  message("Done ")
  
  dist_with_na |> setNames(high_genes)
}



#~ run bootstrap ----

# once <- run_once(seu_ref)
mat_ref <- LayerData(seu_ref, layer = "counts", assay = "RNA")

bootstraps <- boot::boot(t(mat_ref), statistic = run_once, R = 300)


qs::qsave(bootstraps,
          file.path(params$out_dir,
                    paste0("250424_bootstraps_",cell_type,".qs")))



boot_ci <- purrr::map_dfr(seq_along(bootstraps$t0),
                     ~ {
                       res <- boot::boot.ci(bootstraps,type = "perc", index = .x)
                       
                       if(!is.null(res)){
                         intervals <- tibble::as_tibble(res$percent)
                         colnames(intervals) <- c("level","index1","index2","lower","upper")
                       } else{
                         intervals <- data.frame(level = NA_real_,
                                                 index1 = NA_real_,
                                                 index2 = NA_real_,
                                                 lower = NA_real_,
                                                 upper = NA_real_)
                       }
                       
                       tibble::tibble(gene_name = names(bootstraps$t0)[[.x]],
                                      t0 = bootstraps$t0[[.x]]) |>
                         dplyr::bind_cols(intervals)
                       
                     }
)


qs::qsave(boot_ci,
          file.path(params$out_dir,
                    paste0("250424_bootstraps_ci_",cell_type,".qs")))


message("--------------------------------------")
sessionInfo()


message(date())




