# Second part ----

# Inits ----
suppressPackageStartupMessages( library(tidyverse) )
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
    'cnt_thres', 'c', 1, 'integer'
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
    cnt_thres = 30
  )
}




n_rep_pt_global <- 200
nb_subsamples_ptDE <- 100
n_rep_pt_subsamples <- 50


# n_rep_pt_global <- 5
# nb_subsamples_ptDE <- 3
# n_rep_pt_subsamples <- 2


set.seed(123)


cell_types <- readLines(file.path(params$batch_rmed_dir,
                                  "cell_types.txt"))






# loop ----



cell_type <- cell_types[[params$i]]

cell_type_path <- fs::path_sanitize(cell_type)
cell_type



message(params$i,"/", length(cell_types), ": ", cell_type)

log_content <- paste0(cell_type, "\n")


message("---- load")

gene_expressions <- qs::qread(
  file.path(params$batch_rmed_dir,
            paste0(cell_type_path, "_gene_expressions.qs"))
)

subseu <- qs::qread(
  file.path(params$batch_rmed_dir,
            paste0(cell_type_path, "_seu.qs"))
)


# save the plot with phases
gg_phase <- FetchData(subseu,
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



cnts_raw <- GetAssayData(subseu, assay = "RNA", layer = "counts")


log_content <- paste0(log_content,
                      "Nb cells: ", ncol(cnts_raw), "\n")
log_content <- paste0(log_content,
                      "Nb genes raw: ", nrow(cnts_raw), "\n")



# filter genes worth testing
pos_high_genes <- which(
  (gene_expressions$prop_cells >= params$prop_thres) &
    (gene_expressions$nb_cells >= params$cnt_thres)
)

high_genes <- gene_expressions$gene_name[pos_high_genes]
cnts_filtered <- cnts_raw[high_genes, ]


message("---- filtering, keep ", length(high_genes), " genes")



log_content <- paste0(log_content,
                      "Nb tested genes: ", length(high_genes), "\n")



#~ Main pseudotime ----

message("---- pseudotime")

embding <- FetchData(subseu,
                     vars = c("PC_1", "PC_2")) |>
  as.matrix()



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

#gg_int_cca_circlepath

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









#~ Subsamples for PseudotimeDE ----

# note, pseudotimeDE expects between 0 and 1 --> divide by max
cells_pt_global <- tibble(cell = colnames(cnts_filtered),
                          pseudotime = pseudotime / max(pseudotime))



## Generate random subsamples
idx_subsamples <- lapply(seq_len(nb_subsamples_ptDE),
                         function(x) {
                           sample(x = nrow(cells_pt_global),
                                  size = 0.8*nrow(cells_pt_global),
                                  replace = FALSE)
                         })


sub_tbl <- map(idx_subsamples,
               \(idx_one_sub, embding) {
                 
                 
                 sub_embding <- embding[idx_one_sub,]
                 
                 capture.output(
                   CurveEPG <- computeElasticPrincipalCircle(X = sub_embding,
                                                             NumNodes = 20,
                                                             Do_PCA = FALSE,
                                                             nReps = n_rep_pt_subsamples,
                                                             drawPCAView = FALSE,
                                                             drawAccuracyComplexity = FALSE,
                                                             drawEnergy = FALSE,
                                                             ProbPoint = .6,
                                                             verbose = FALSE),
                   file = nullfile()
                 )
                 
                 
                 stopifnot( length(CurveEPG) == n_rep_pt_subsamples+1L )
                 targt <- CurveEPG[[ n_rep_pt_subsamples+1L ]]
                 boots <- CurveEPG[ seq_len(n_rep_pt_subsamples) ]
                 
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
                 
                 
                 pseudotime_sub <- getPseudotime(ProjStruct = ProjStruct,
                                                 NodeSeq = names(subgraph))[[ "Pt" ]]
                 
                 
                 cells_pt_sub <- tibble(cell = rownames(embding),
                                        pseudotime = pseudotime_sub / max(pseudotime_sub) )
                 
                 
                 
                 ## Make sure the subsample pseudotime is comparable to global pseudotime
                 
                 stopifnot(all.equal( cells_pt_sub$cell, cells_pt_global$cell ))
                 
                 
                 # synchronize starting point (pt = 0)
                 pt_offset <- cells_pt_sub$pseudotime[ which.min(cells_pt_global$pseudotime) ]
                 permuted <- (cells_pt_sub$pseudotime - pt_offset) %% 1
                 
                 # normalize to 0-1
                 cells_pt_sub$pseudotime <- permuted
                 
                 # general direction
                 if(cor(cells_pt_sub$pseudotime, cells_pt_global$pseudotime) < 0) {
                   cells_pt_sub <- dplyr::mutate(cells_pt_sub,
                                                 pseudotime = 1-pseudotime)
                 }
                 
                 cells_pt_sub
                 
               }, embding = embding)


stopifnot(all( high_genes %in% rownames(cnts_filtered) ))


message("---- checking pseudotime subsamples")


# check
png(file.path(params$out_dir, paste0(cell_type_path, "_pseudotime_subsamples.png")),
    width = 800, height = 400)
prepar <- par(no.readonly = TRUE)
par(mfrow = c(1,2))

plot(cells_pt_global$pseudotime, sub_tbl[[1]]$pseudotime)
plot(sub_tbl[[1]]$pseudotime, sub_tbl[[2]]$pseudotime)

par(prepar)
dev.off()


message("---- computing DE")


#~ DE ----
tictoc::tic()
res_pseudotimeDE <- PseudotimeDE::runPseudotimeDE(gene.vec = high_genes,
                                                  ori.tbl = cells_pt_global,
                                                  sub.tbl = sub_tbl,
                                                  mat = round(cnts_filtered),
                                                  model = params$model,
                                                  formula = expv ~ s(pseudotime, k = 6, bs = 'cc'),
                                                  mc.cores = 1)
tictoc::toc()


message("---- saving results")


qs::qsave(res_pseudotimeDE,
          file.path(params$out_dir, paste0(cell_type_path, "_res_pseudotimeDE.qs")) )

qs::qsave(cells_pt_global,
          file.path(params$out_dir, paste0(cell_type_path, "_cells_pt_global.qs")))


log_content <- paste0(log_content,
                      "__\n\nSessionInfo:\n",
                      paste(toLatex(sessionInfo()), collapse = "\n"),
                      "\n")


writeLines(log_content,
           file.path(params$out_dir, paste0(cell_type_path, "_log.txt")) )





