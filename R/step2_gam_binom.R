# Second part ----

# Resampling: take same cells, compute pseudotime,
# permute gene and compute GAM (not using PseudotimeDE),
# keep GAM parameters for perm test


# Inits ----
library(ggplot2)
suppressPackageStartupMessages( library(Seurat) )
library("ElPiGraph.R")


library(getopt)


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
    batch_rmed_dir = "intermediates/2502/250319_step1",
    out_dir = "intermediates/2502/250319_step2",
    model = "auto",
    i = 3,
    prop_thres = 0.1,
    cnt_thres = 30,
    nb_subsamples_ptDE = 10
  )
}





n_rep_pt_global <- 50



set.seed(123)


cell_types <- readLines(file.path(params$batch_rmed_dir,
                                  "cell_types.txt"))




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

seu <- qs::qread(
  file.path(params$batch_rmed_dir,
            paste0(cell_type, "_seu.qs"))
)


# save the plot with phases
gg_phase <- FetchData(seu,
                      vars = c("PC_1", "PC_2", "mean_angle", "mean_rho")) |>
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 2*pi)) +
  geom_point(aes(x = PC_1, y = PC_2,
                 color = mean_angle, alpha = mean_rho)) +
  ggtitle(cell_type)

ggsave(paste0(cell_type_path,"_harmony_phase.png"), gg_phase,
       path = params$out_dir,
       width = 7, height = 5, units = "in")




# Pick genes ----
pos_high_genes <- which(
  (gene_expressions$prop_cells >= params$prop_thres) &
    (gene_expressions$nb_cells >= params$cnt_thres)
)

high_genes <- gene_expressions$gene_name[pos_high_genes]

length(high_genes)




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

ggsave(paste0(cell_type_path,"_circlepath.png"), gg_int_cca_circlepath,
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
message("---- fit GAM")

mat_cnt <- GetAssayData(seu, assay = "RNA", layer = "count")[high_genes,]

mods <- lapply(high_genes,
               \(.gene){
                 dat <- data.frame(
                   prob = mat_cnt[.gene,] > 0,
                   pseudotime = pseudotime/max(pseudotime)
                 )
                 
                 mgcv::gam(prob ~ s(pseudotime, k = 6, bs = 'cc'),
                           data = dat,family = "binomial")
               }) |>
  setNames(high_genes)


all_preds <- vapply(mods,
                    \(.mod) predict(.mod,
                                    type = "response",
                                    newdata = data.frame(pseudotime = (0:199)/200)),
                    FUN.VALUE = double(200))



#~ metrics ----
message("---- compute metrics")

res <- tibble::tibble(
  cell_type = cell_type,
  gene_name = names(mods),
  amplitude = apply(all_preds, 2, \(.x) diff(range(.x))),
  dev_expl = vapply(mods,
                    \(.mod) summary(.mod)[["dev.expl"]],
                    FUN.VALUE = double(1L)),
  gam_fit = mods
)






# save results ----
message("---- save results")


qs::qsave(res,
          file.path(params$out_dir, paste0(cell_type_path, "_res_gam.qs")) )

qs::qsave(cells_pt_global,
          file.path(params$out_dir, paste0(cell_type_path, "_cells_pt_global.qs")))



log_content <- paste0(log_content,
                      "__\n\nSessionInfo:\n",
                      paste(toLatex(sessionInfo()), collapse = "\n"),
                      "\n")


writeLines(log_content,
           file.path(params$out_dir, paste0(cell_type_path, "_log.txt")) )


message("--------------------------------------")
sessionInfo()

