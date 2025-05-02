# Inits ----

library(Seurat)
library("velocyto.R")
library(tidyverse)

reticulate::use_virtualenv("anndata")


dir_loom <- "intermediates/2502/250409_loom/"
dir_out_anndata <- "intermediates/2502/250409_anndata/"


seu <- qs::qread("intermediates/2502/250328_assembled/250329_seu_all_herma.qs")
samples <- unique(seu$orig.ident)

stopifnot(all(
  samples %in% (list.files(dir_loom) |> str_remove("\\.loom$"))
))

# select "easy" sample for testing
# table(seu$orig.ident) |> sort()
# which(samples == "230421_AM")

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

emat_tot <- do.call(cbind,
                    mats_list[["spliced"]])
nmat_tot <- do.call(cbind,
                    mats_list[["unspliced"]])

dim(emat_tot); dim(nmat_tot)

stopifnot(all( colnames(emat_tot) %in% bc_table$bc_velocyto ))
stopifnot(all.equal( colnames(emat_tot), colnames(nmat_tot) ))

colnames(emat_tot) <- column_to_rownames(bc_table, "bc_velocyto")[colnames(emat_tot), "bc_seu"] 
colnames(nmat_tot) <- column_to_rownames(bc_table, "bc_velocyto")[colnames(nmat_tot), "bc_seu"] 




# by cell type
cell_types <- unique(seu$cell_type)

for (ct in cell_types){
  
  message("---- Cell type: ", ct)
  
  
  sub <- subset(seu,
                cell_type == ct)
  
  
  
  
  emat <- emat_tot[,colnames(sub)]
  nmat <- nmat_tot[,colnames(sub)]
  
  
  dim(emat)
  
  emat[1:3,1:4]
  
  anndata::AnnData(X = t(emat),
                   layers = list(spliced = t(emat),
                                 unspliced = t(nmat))) |>
    anndata::write_h5ad(filename = file.path(dir_out_anndata,
                                             paste0(ct,".h5ad")))
  
  message("done...")
  
}








# explore example cell type ----
# which(cell_types == "BWM")
ct <- cell_types[[81]]

ct

sub <- subset(seu,
              cell_type == ct)


emat <- emat_tot[,colnames(sub)]
nmat <- nmat_tot[,colnames(sub)]

#~ filter genes ----
dim(emat); dim(nmat)

emat |> rowMeans() |> log10() |> hist(breaks = 50); abline(v = log10(.02), col = 'red3')
nmat |> rowMeans() |> log10() |> hist(breaks = 50); abline(v = log10(.02), col = 'red3')

emat_f <- emat[rowMeans(emat) >= .02, ]
nmat_f <- nmat[rowMeans(nmat) >= .02, ]

dim(emat_f); dim(nmat_f)

length(intersect(rownames(nmat_f),rownames(nmat_f)))



#~ process PCA ----
sub <- SCTransform(sub)
sub <- RunPCA(sub, npcs = 2, verbose = FALSE)

dat <- FetchData(sub, vars = c("PC_1","PC_2","cell_phase_masked", "cell_rho"))

dat |>
  ggplot() +
  theme_classic() +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 360)) +
  geom_point(aes(x = PC_1, y = PC_2,
                 color = cell_phase_masked,
                 alpha = cell_rho))



#~ Process velocity ----

rvel.cd <- gene.relative.velocity.estimates(emat_f,
                                            nmat_f)


stopifnot(all.equal(rownames(dat), colnames(rvel.cd$current)))



# color-code cell phase per cell
phase_per_bc <- dat[["cell_phase_masked"]]
ref_cols <- pals::kovesi.cyclic_mrybm_35_75_c68(50) |> colorRamp()



phase_per_bc[!is.na(phase_per_bc)] <- (phase_per_bc / 360) |>
  na.exclude() |>
  ref_cols() |>
  (\(x) x / 255 )() |>
  rgb()

phase_per_bc[is.na(phase_per_bc)] <- 'grey'

names(phase_per_bc) <- rownames(dat)

show.velocity.on.embedding.cor(as.matrix(dat[,1:2]),
                               rvel.cd,
                               n = 300,
                               scale = 'sqrt',
                               cex = 0.8,
                               arrow.scale = 2,
                               show.grid.flow = TRUE,
                               cell.colors = phase_per_bc,
                               min.grid.cell.mass = 0.5,
                               grid.n = 40,
                               arrow.lwd = 1,
                               do.par = F,
                               cell.border.alpha = 0.1)












