# save filtered matrix
library(Seurat)

dir_alignements <- "/vast/palmer/scratch/hammarlund/aw853/250331_align"

seu <- qs::qread("intermediates/2502/250328_assembled/250329_seu_all_herma.qs")
samples <- unique(seu$orig.ident)

stopifnot(all(
  samples %in% list.files(dir_alignements)
))


for(samp in samples){
  
  message(samp)
  
  file.rename(file.path(dir_alignements, samp, "outs", "filtered_feature_bc_matrix"),
              file.path(dir_alignements, samp, "outs", "cr_filtered_feature_bc_matrix"))
  
  sub <- subset(seu,
                orig.ident == samp)
  
  mat <- LayerData(sub, layer = "counts", assay = "RNA")
  
  match_bc <- stringr::str_match(colnames(mat),
                                 "^s[0-9]+_([ATCG]{16})\\-1_[1-9]$")
  
  stopifnot(anyDuplicated(match_bc[,2]) == 0)
  
  colnames(mat) <- match_bc[,2]
  
  DropletUtils::write10xCounts(
    file.path(dir_alignements, samp, "outs", "filtered_feature_bc_matrix"),
    x = mat,
    version = "3"
  )
}
