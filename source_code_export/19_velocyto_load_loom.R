
message("--- Starting ", date())

# Inits ----

library(Seurat)
library("velocyto.R")
library(tidyverse)



dir_loom <- "intermediates/2502/250409_loom/"
dir_velocyto <- "intermediates/2502/250825_velocyto"

seu <- qs::qread("intermediates/2502/250605_assembled/250606_seu_all_herma.qs")


samples <- unique(seu$orig.ident)

stopifnot(all(
  samples %in% (list.files(dir_loom) |> str_remove("\\.loom$"))
))



message("--- Load")
bc_table <- seu[[]] |>
  select(orig.ident, cell_type) |>
  rownames_to_column("bc_seu") |>
  separate_wider_regex(bc_seu,
                       patterns = c(sample = "^s[0-9]+",
                                    "_",
                                    cell_bc = "[ATCG]{16}",
                                    "\\-1_[1-9]$"),
                       cols_remove = FALSE) |>
  mutate(bc_velocyto = paste0(orig.ident,":",cell_bc,"x"))



mats_list <- map(samples,
                 \(samp){
                   path <- file.path(dir_loom,
                                     paste0(samp, ".loom"))
                   
                   read.loom.matrices(path)
                 }) |>
  transpose()

message("--- Extract matrices")

emat_tot <- do.call(cbind,
                    mats_list[["spliced"]])
nmat_tot <- do.call(cbind,
                    mats_list[["unspliced"]])

dim(emat_tot); dim(nmat_tot)

# stopifnot(all( colnames(emat_tot) %in% bc_table$bc_velocyto ))
stopifnot(all.equal( colnames(emat_tot), colnames(nmat_tot) ))

# colnames(emat_tot) <- column_to_rownames(bc_table, "bc_velocyto")[colnames(emat_tot), "bc_seu"] 
# colnames(nmat_tot) <- column_to_rownames(bc_table, "bc_velocyto")[colnames(nmat_tot), "bc_seu"] 

bc_seu_recoded <- enframe(colnames(emat_tot)) |> left_join(bc_table, by = c(value = "bc_velocyto"), multiple = "any") |> pull(bc_seu)

colnames(emat_tot) <- colnames(nmat_tot) <- bc_seu_recoded

message("--- Save")

qs::qsave(emat_tot,
          file.path(dir_velocyto, "emat_tot.qs"))
qs::qsave(nmat_tot,
          file.path(dir_velocyto, "nmat_tot.qs"))



message("--------------------------------------------------")
message("-------  Done ", date(), " -------")

sessionInfo()



