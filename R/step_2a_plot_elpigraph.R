


# Extracted code from `step_2a_gam_nb.R` for plotting figure




# Inits ----
library(ggplot2)
suppressPackageStartupMessages( library(Seurat) )
library("ElPiGraph.R")
source("R/utils_fit.R") # -> circ_perm_mat()

library(getopt)
if(! interactive()){} else{
  # Options for interactive
  params <- list(
    dir_step1 = "intermediates/2502/250609_step1",
    out_dir = "intermediates/2502/250624_step2",
    i = 8,
    prop_thres = 0.1,
    cnt_thres = 30
  )
}
# for ElPigraph
n_rep_pt_global <- 50
# predictions length
len <- 128
set.seed(123)


params$dir_step1 <- "E:/backups/Projects_june2025/glia/osc_larva/intermediates/2502/250609_step1/"

cell_types <- list.files(params$dir_step1,
                         pattern = "_seu\\.qs$") |>
  stringr::str_remove("_seu\\.qs$")


cell_type="ILso"
message(params$i,"/", length(cell_types), ": ", cell_type)
gene_expressions <- qs::qread(
  file.path(params$dir_step1,
            paste0(cell_type, "_gene_expressions.qs"))
)
seu <- qs::qread(
  file.path(params$dir_step1,
            paste0(cell_type, "_seu.qs"))
)

seu_unsmoothed <- qs::qread(
  file.path(params$dir_step1,
            paste0(cell_type, "_seu_unsmoothed.qs"))
)


# save the plot with phases
embding <- seu |>
  FetchData(vars = c("PC_1", "PC_2")) |>
  as.matrix()
all.equal(rownames(embding), rownames(seu[[]]))
head(embding)





# for ILso, invert axes for easier interpretation

embding[,1] <- -embding[,1]

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


PlotPG(embding,
       TargetPG = targt, BootPG = boots,
       Do_PCA = FALSE,VizMode = c("Target","Boot"),
       DimToPlot = 1:2)[[1]] +
  theme_classic() +
  scale_linewidth_manual(values = c(1,2,3))






## Modified PlotPG()
X = embding
TargetPG = targt
BootPG = boots
Do_PCA = FALSE
VizMode = c("Target","Boot")
DimToPlot = 1:2


PGCol = ""
PlotProjections = "none" 
GroupsLab = NULL
PointViz = "points"
Main = ""
p.alpha = 0.3
PointSize = NULL
NodeLabels = NULL
LabMult = 1

######

if (length(PGCol) == 1) {
  PGCol = rep(PGCol, nrow(TargetPG$NodePositions))
}
if (is.null(GroupsLab)) {
  GroupsLab = factor(rep("N/A", nrow(X)))
}
levels(GroupsLab) <- c(levels(GroupsLab), unique(PGCol))


BaseData <- X
DataVarPerc <- apply(X, 2, var)/sum(apply(X, 2, var))
DimToPlot <- intersect(DimToPlot, 1:ncol(X))


# only one combination (1 col of AllComb)
AllComb <- combn(DimToPlot, 2)
i <- 1

Idx1 <- AllComb[1, i]
Idx2 <- AllComb[2, i]


df1 <- data.frame(PCA = BaseData[, Idx1],
                  PCB = BaseData[, Idx2],
                  Group = GroupsLab)

p <- ggplot2::ggplot(data = df1,
                     mapping = ggplot2::aes(x = PCA,
                                            y = PCB),
                     environment = environment()) +
  ggplot2::geom_point(alpha = p.alpha,
                      mapping = ggplot2::aes(color = Group))
Inizialized <- TRUE



RotData <- TargetPG$NodePositions



tEdg <- t(sapply(
  1:nrow(TargetPG$Edges$Edges),
  function(i) {
    Node_1 <- TargetPG$Edges$Edges[i, 1]
    Node_2 <- TargetPG$Edges$Edges[i, 2]
    if (PGCol[Node_1] == PGCol[Node_2]) {
      tCol = paste("ElPiG", PGCol[Node_1])
    }
    if (PGCol[Node_1] != PGCol[Node_2]) {
      tCol = "ElPiG Multi"
    }
    if (any(PGCol[c(Node_1, Node_2)] == "None")) {
      tCol = "ElPiG None"
    }
    c(RotData[Node_1, c(Idx1, Idx2)], RotData[Node_2, 
                                              c(Idx1, Idx2)], tCol)
  }
))

AllEdg <- cbind(tEdg, 0)

TarPGVarPerc <- apply(TargetPG$NodePositions, 2, 
                      var)/sum(apply(TargetPG$NodePositions, 2, var))


df2 <- data.frame(x = as.numeric(AllEdg[, 1]),
                  y = as.numeric(AllEdg[, 2]),
                  xend = as.numeric(AllEdg[, 3]),
                  yend = as.numeric(AllEdg[, 4]),
                  Col = AllEdg[, 5],
                  Rep = as.numeric(AllEdg[, 6]),
                  stringsAsFactors = FALSE)



if (!is.null(BootPG) & ("Boot" %in% VizMode)) {
  AllEdg <- lapply(1:length(BootPG), function(i) {
    tTree <- BootPG[[i]]
    if (Do_PCA) {
      RotData <- t(t(tTree$NodePositions) - CombPCA$center) %*% 
        CombPCA$rotation
    }
    else {
      RotData <- tTree$NodePositions
    }
    tEdg <- t(sapply(1:nrow(tTree$Edges$Edges), 
                     function(i) {
                       c(RotData[tTree$Edges$Edges[i, 1], c(Idx1, 
                                                            Idx2)], RotData[tTree$Edges$Edges[i, 2], 
                                                                            c(Idx1, Idx2)])
                     }))
    cbind(tEdg, i)
  })
  AllEdg <- do.call(rbind, AllEdg)
  df3 <- data.frame(x = AllEdg[, 1], y = AllEdg[, 
                                                2], xend = AllEdg[, 3], yend = AllEdg[, 4], 
                    Rep = AllEdg[, 5])
  p <- p + ggplot2::geom_segment(data = df3, mapping = ggplot2::aes(x = x, 
                                                                    y = y, xend = xend, yend = yend), inherit.aes = FALSE, 
                                 alpha = 0.2, color = "black")
}







if ("Target" %in% VizMode) {
  if (is.factor(GroupsLab)) {
    p <- p + ggplot2::geom_segment(data = df2, mapping = ggplot2::aes(x = x, 
                                                                      y = y, xend = xend, yend = yend, col = Col), 
                                   inherit.aes = TRUE) + ggplot2::labs(linetype = "")
  }
  else {
    p <- p + ggplot2::geom_segment(data = df2, mapping = ggplot2::aes(x = x, 
                                                                      y = y, xend = xend, yend = yend), inherit.aes = FALSE)
  }
}




df4 <- data.frame(PCA = TargetPG$NodePositions[, 
                                               Idx1], PCB = TargetPG$NodePositions[, Idx2])




if ("Target" %in% VizMode) {
  if (!is.null(PointSize)) {
    if (!is.na(PointSize)) {
      p <- p + ggplot2::geom_point(mapping = ggplot2::aes(x = PCA, 
                                                          y = PCB, size = PointSize), data = df4, 
                                   inherit.aes = FALSE)
    }
    else {
      p <- p + ggplot2::geom_point(mapping = ggplot2::aes(x = PCA, 
                                                          y = PCB), data = df4, size = 0, inherit.aes = FALSE)
    }
  }
  else {
    p <- p 
  }
}





#~ replot final ----

ggplot(
  data = df1,
  mapping = aes(x = PCA, y = PCB)
) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  xlab("PC 1") + ylab("PC 2") +
  geom_point(
    alpha = .2,
    shape = 16
  ) +
  geom_segment(
    data = df3[1:(20*5),],
    mapping = aes(x = x, y = y,
                  xend = xend, yend = yend),
    inherit.aes = FALSE, 
    alpha = 0.4,
    linewidth = .7,
    color = "#cdb974"
  ) +
  geom_segment(
    data = df2,
    mapping = aes(x = x,y = y, xend = xend, yend = yend),
    color = "#4c72b1",
    linewidth = 1.5,
    inherit.aes = TRUE
  ) +
  geom_point(
    mapping = aes(x = PCA, 
                  y = PCB),
    data = df4, inherit.aes = FALSE,
    size = 5,
    shape = "+",
    color = "#c54e53"
  )

ggsave("elpigraph.pdf", path = "presentations/figures/250625_gam_illustrations",
       width = 85, height = 85, units = "mm")









