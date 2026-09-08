
library(SoupX)
library(Seurat)
library(wbData)

gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = NA, gene_id = "nsIs198",
                  symbol = "nsIs198", sequence = "nsIs198",
                  status = "Live", biotype = "protein_coding_gene",
                  name = "nsIs198")

dir_aligned <- "data/250225_aligned"
dir_filtered <- "intermediates/2502/empytdrops_filt"
dir_soupx <- "intermediates/2502/soupx"

samples <- list.files(dir_filtered, pattern = "filt_feature_bc_matrix$") |>
  stringr::str_remove("filt_feature_bc_matrix$")


# a list of genes that are good markers of specific tissues/cell types
tissue_marks <- readxl::read_excel("/gpfs/ycga/work/hammarlund/aw853/references/data/tissue_markers_250304.xlsx")
tissue_marks <- tissue_marks[tissue_marks$gene != "Y17D7C.4", ]

marks_by_tissue <- split(tissue_marks$gene, tissue_marks$tissue)

# samples <- samples[10:17]


for(sample in samples){
  # sample <- "L2_3"
  message("--- ", sample, " ---")
  
  tab_counts = Seurat::Read10X(file.path(dir_filtered, paste0(sample,"filt_feature_bc_matrix") ))
  rownames(tab_counts) <- i2s(rownames(tab_counts), gids, warn_missing = TRUE)
  
  tab_droplets = Seurat::Read10X(file.path(dir_aligned, sample, "outs", 'raw_feature_bc_matrix'))
  
  
  
  seu <- CreateSeuratObject(tab_counts)
  
  # seu
  
  # Generating initial clustering with Seurat
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, nfeatures = 1000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = 70, verbose = FALSE)
  
  # ElbowPlot(seu, ndims = 70)
  
  
  
  seu <- FindNeighbors(seu,
                       dims = 1:30,
                       verbose = FALSE)
  
  seu <- FindClusters(seu,
                      resolution = 2,
                      verbose = FALSE)
  
  seu <- RunUMAP(seu,
                 dims = 1:30,
                 verbose = FALSE)
  
  
  
  # SoupX
  sc <- SoupChannel(tab_droplets,tab_counts)
  sc <- setClusters(sc, Idents(seu))
  
  
  
  # use list of markers to estimate conta
  
  # unlist(marks_by_tissue)[!unlist(marks_by_tissue) %in% rownames(tab_counts)]
  useToEst = estimateNonExpressingCells(sc,
                                        nonExpressedGeneList = marks_by_tissue)
  
  sc = calculateContaminationFraction(sc,
                                      nonExpressedGeneList = marks_by_tissue,
                                      useToEst=useToEst,
                                      forceAccept = TRUE)
  
  out <- adjustCounts(sc)
  
  
  DropletUtils::write10xCounts(file.path(dir_soupx,
                                         paste0(sample,"_soupx_feature_bc_matrix")),
                               x = out,
                               version = "3")
  
}


sessionInfo()

