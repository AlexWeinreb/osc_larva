# Second part ----

# Fit GAM (NB, other versions used binomial).
# We get a (uncentered) smooth prediction --> keep
#
# we can center and redo a fit --> centered smooth
#   * keep coefficients of centered GAM for clustering
#   * descriptors of centered smooth curve for clustering
#   * dtw of centered smooth curve for clustering (single number, other script used bootstrap for CI, doesn't appear useful)
# 


# Inits ----
library(ggplot2)
suppressPackageStartupMessages( library(Seurat) )
library("ElPiGraph.R")


source("R/utils_fit.R") # -> circ_perm_i()

library(getopt)



if(! interactive()){
  spec <- matrix(c(
    'dir_step1', 'f', 1, 'character',
    'out_dir', 'o', 1, 'character',
    'i', 'i', 1, 'integer',
    'prop_thres', 'g', 1, 'double',
    'cnt_thres', 'c', 1, 'integer',
    'nb_subsamples_ptDE', 'n', 1, 'integer'
  ), byrow=TRUE, ncol=4)
  
  params <- getopt(spec)
  
} else{
  # Options for interactive
  params <- list(
    dir_step1 = "intermediates/2502/250502_step1",
    out_dir = "intermediates/2502/250502_step2_nb",
    i = 8,
    prop_thres = 0.1,
    cnt_thres = 30,
    nb_subsamples_ptDE = 10
  )
}




# for ElPigraph
n_rep_pt_global <- 50

# predictions length
len <- 128




set.seed(123)


cell_types <- list.files(params$dir_step1,
                         pattern = "_seu\\.qs$") |>
  stringr::str_remove("_seu\\.qs$")





# loop ----



cell_type <- cell_types[[params$i]]

cell_type


message(params$i,"/", length(cell_types), ": ", cell_type)



message("---- load")

gene_expressions <- qs::qread(
  file.path(params$dir_step1,
            paste0(cell_type, "_gene_expressions.qs"))
)

seu <- qs::qread(
  file.path(params$dir_step1,
            paste0(cell_type, "_seu.qs"))
)

message("Number of cells: ", ncol(seu))

if(ncol(seu) < 25){
  message("####   Not enough cells, skipping   ####")
  quit()
}




# save the plot with phases
gg_phase <- FetchData(seu,
                      vars = c("PC_1", "PC_2", "cell_phase_masked", "cell_rho")) |>
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 360)) +
  geom_point(aes(x = PC_1, y = PC_2,
                 color = cell_phase_masked, alpha = cell_rho)) +
  ggtitle(cell_type)

# gg_phase

ggsave(paste0(cell_type,"_phase.png"), gg_phase,
       path = params$out_dir,
       width = 7, height = 5, units = "in")




# Pick genes ----
pos_high_genes <- which(
  (gene_expressions$prop_cells >= params$prop_thres) &
    (gene_expressions$nb_cells >= params$cnt_thres)
)

high_genes <- gene_expressions$gene_name[pos_high_genes]


message("Number of genes: ", length(high_genes))

if(length(high_genes) < 5){
  message("####   No gene to test, skip   ####")
  quit()
}



embding <- seu |>
  FetchData(vars = c("PC_1", "PC_2")) |>
  as.matrix()

all.equal(rownames(embding), rownames(seu[[]]))






#~ pseudotime ----


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




#~  GAM uncentered ----
message("---- fit GAM uncentered")

mat_cnt <- GetAssayData(seu, assay = "RNA", layer = "count")[high_genes,]


mods_uncentered <- lapply(
  high_genes,
  \(.gene){
    
    dat <- data.frame(
      prob = mat_cnt[.gene,],
      pseudotime = pseudotime/max(pseudotime)
    )
    
    
    mgcv::gam(prob ~ s(pseudotime, k = 6, bs = 'cc'),
              data = dat, family = mgcv::nb(link = "log"))
    
  }) |>
  setNames(high_genes)

# predictions
preds_uncentered <- vapply(mods_uncentered,
                    \(.mod) predict(.mod,
                                    type = "response",
                                    newdata = data.frame(pseudotime = (0:(len-1))/len )),
                    FUN.VALUE = double(len))



#~  GAM centered ----
message("---- fit GAM centered")

pos_peak <- apply(preds_uncentered, 2, which.max) / len

# ngenes <- ncol(preds_uncentered)
# printMat::matimage(log1p(preds_uncentered))
# points((1:ngenes - 1)/(ngenes - 1), (1 - pos_peak), pch = "-", cex = 3.5, col = 'purple')


mods_centered <- lapply(
  high_genes,
  \(.gene){
    
    dat <- data.frame(
      prob = mat_cnt[.gene,],
      pseudotime = pseudotime/max(pseudotime)
    )
    
    # shift pseudotime
    dat$pseudotime_centered <- ( .5 + dat$pseudotime - pos_peak[[.gene]] ) %% 1
    
    # plot(dat$pseudotime, log1p(dat$prob)); abline(v = pos_peak[[.gene]])
    # lines((0:(len-1))/len , log1p(preds_uncentered[,.gene]))
    # 
    # plot(dat$pseudotime_centered, dat$prob, lab = c(10,5,7)); abline(v = .5)
    
    
    mgcv::gam(prob ~ s(pseudotime_centered, k = 6, bs = 'cc'),
              data = dat, family = mgcv::nb(link = "log"))
    
  }) |>
  setNames(high_genes)



preds_centered <- vapply(mods_centered,
                           \(.mod) predict(.mod,
                                           type = "response",
                                           newdata = data.frame(pseudotime_centered = (0:(len-1))/len )),
                           FUN.VALUE = double(len))



# apply(preds_centered, 2, which.max)
# len / 2

# ngenes <- ncol(preds_centered)
# printMat::matimage(log1p(preds_centered))
# points((1:ngenes - 1)/(ngenes - 1),
#        y = (1 - apply(preds_centered, 2, which.max) / len), pch = "-", cex = 3.5, col = 'purple')






#~ metrics ----

message("---- compute metrics")


stopifnot(
  identical(
    high_genes,
    names(mods_centered)
  ) &&
    identical(
      high_genes,
      colnames(preds_centered)
    ) &&
    identical(
      high_genes,
      names(mods_uncentered)
    ) &&
    identical(
      high_genes,
      colnames(preds_uncentered)
    )
)






res <- data.frame(
  cell_type = cell_type,
  gene_name = names(mods_centered)
)


#~~ coefficients of centered GAM ----
message("  ---- coefs GAM")

coefs_gam <- vapply(
  mods_centered,
  coef,
  double(5L)
) |>
  t() |>
  as.data.frame() |>
  setNames(c("intercept", "s1", "s2", "s3", "s4"))


res <- cbind(res, coefs_gam)

res$mse <- vapply(
  mods_centered,
  \(.mod) mean( residuals(.mod)^2 ),
  double(1L)
)


#~~ curve descriptors ----
message("  ---- curve descriptors")

res$amplitude <- apply( preds_centered, 2,
                        \(.x) diff(range(.x)) )

res$dev_expl <- vapply(
  mods_centered,
  \(.mod) summary(.mod)[["dev.expl"]],
  FUN.VALUE = double(1L)
)

res$area_under_curve <- apply( preds_centered, 2,
                               \(.y) pracma::trapz(seq_along(.y), .y) )




#~~ dtw ----
message("  ---- dtw")

preds_scaled <- apply(preds_centered, 2, \(x) x/max(x) )

# printMat::matimage(log1p(preds_centered))
# printMat::matimage(preds_scaled)


# reference curve
ref <- dnorm( (seq_len(len) - 1)/len , mean = .5, sd = .01)
ref <- (ref - min(ref))/max(ref - min(ref))
ref <- circ_perm_mat(matrix(ref, ncol = 1))


# matplot((seq_len(len) - 1)/len, preds_scaled, type = 'l',
#         xlab = "Pseudotime",
#         ylab = "log(count + 1)",
#         lwd = .7)
# lines((seq_len(len) - 1)/len, ref, lwd = 2.5)

res$dist_dtw <- dtw::dtwDist(t(preds_scaled), t(ref)) |>
  as.numeric()


head(res)





# save results ----
message("---- save results")

# use for clustering
qs::qsave(res,
          file.path(params$out_dir, paste0(cell_type, "_descriptors.qs")) )

# use for representations
qs::qsave(preds_uncentered,
          file.path(params$out_dir, paste0(cell_type, "_preds.qs")))



# in case
qs::qsave(mods_uncentered,
          file.path(params$out_dir, paste0(cell_type, "_mods_uncentered.qs")))

qs::qsave(mods_centered,
          file.path(params$out_dir, paste0(cell_type, "_mods_centered.qs")))




message("--------------------------------------")
sessionInfo()



message("Finished ---- ", date())







