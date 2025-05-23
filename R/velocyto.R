# Inits ----

library(Seurat)
library("velocyto.R")
library(tidyverse)



dir_loom <- "intermediates/2502/250409_loom/"
dir_velocyto <- "intermediates/2502/250522_velocyto"

seu <- qs::qread("intermediates/2502/250509_assembled/250509_seu_all_herma.qs")
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

# qs::qsave(emat_tot,
#           file.path(dir_velocyto, "emat_tot.qs"))
# qs::qsave(nmat_tot,
#           file.path(dir_velocyto, "nmat_tot.qs"))



# # by cell type, annData
# dir_out_anndata <- "intermediates/2502/250521_anndata/"
# 
# cell_types <- unique(seu$cell_type)
# 
# length(cell_types)
# for (ct in cell_types){
#   
#   message("---- Cell type: ", ct)
#   
#   
#   sub <- subset(seu,
#                 cell_type == ct)
#   
#   
#   
#   
#   emat <- emat_tot[,colnames(sub)]
#   nmat <- nmat_tot[,colnames(sub)]
#   
#   
#   dim(emat)
#   
#   emat[1:3,1:4]
#   
#   anndata::AnnData(X = t(emat),
#                    layers = list(spliced = t(emat),
#                                  unspliced = t(nmat))) |>
#     anndata::write_h5ad(filename = file.path(dir_out_anndata,
#                                              paste0(ct,".h5ad")))
#   
#   message("done...")
#   
# }










# explore example cell type ----

library(Seurat)
library("velocyto.R")
library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = NA, gene_id = "nsIs198",
                  symbol = "nsIs198", sequence = "nsIs198",
                  status = "Live", biotype = "protein_coding_gene",
                  name = "nsIs198") |>
  mutate(name = if_else(name == "E_BE45912.2",
                        "E-BE45912.2",
                        name))

# reticulate::use_virtualenv("anndata")


# dir_loom <- "intermediates/2502/250409_loom/"

dir_step1 <- "intermediates/2502/250522_step1"
dir_velocyto <- "intermediates/2502/250521_velocyto"


emat_tot <- qs::qread(file.path(dir_velocyto, "emat_tot.qs"))
nmat_tot <- qs::qread(file.path(dir_velocyto, "nmat_tot.qs"))

rownames(emat_tot) <- i2s(rownames(emat_tot), gids, warn_missing = TRUE)
rownames(nmat_tot) <- i2s(rownames(nmat_tot), gids, warn_missing = TRUE)


# cell_types <- unique(seu$cell_type)



# which(cell_types == "ILso")
# ct <- cell_types[[45]]

ct <- "ILso"
ct

sub <- qs::qread(file.path(dir_step1,
                           paste0(ct, "_seu_unsmoothed.qs")) )

stopifnot(all(
  rownames(sub) %in% rownames(emat_tot)
))
stopifnot(all(
  colnames(sub) %in% colnames(emat_tot)
))


emat <- emat_tot[rownames(sub), colnames(sub)]
nmat <- nmat_tot[rownames(sub), colnames(sub)]


# mat_cnt <- GetAssayData(sub, assay = "RNA", layer = "counts")
# 
# ex_rows <- sample(rownames(emat), 100)
# ex_cols <- sample(colnames(emat), 100)
# tibble(
#   emat = as.numeric(emat[ex_rows,ex_cols]),
#   seu = as.numeric(mat_cnt[ex_rows,ex_cols])
# ) |>
#   ggplot() +
#   theme_classic() +
#   scale_x_continuous(transform = "log1p") +
#   scale_y_continuous(transform = "log1p") +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point(aes(x = emat, y = seu),
#               alpha = .1)
#   # geom_jitter(aes(x = emat, y = seu),
#   #            alpha = .1,width = .01, height = 0)



# #~ filter genes ----
dim(emat); dim(nmat)

emat |> rowMeans() |> log10() |> hist(breaks = 50); abline(v = log10(.02), col = 'red3')
nmat |> rowMeans() |> log10() |> hist(breaks = 50); abline(v = log10(.02), col = 'red3')

emat_f <- emat[rowMeans(emat) >= .02, ]
nmat_f <- nmat[rowMeans(nmat) >= .02, ]

dim(emat_f); dim(nmat_f)

length(intersect(rownames(nmat_f),rownames(nmat_f)))



#~ process PCA ----

dat <- FetchData(sub, vars = c("PC_1","PC_2","cell_phase_masked", "cell_rho"))

# for ILso, invert axes for easier interpretation
if(ct == "ILso"){
  dat$PC_1 <- -dat$PC_1
}


dat |>
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 360)) +
  geom_point(aes(x = PC_1, y = PC_2,
                 color = cell_phase_masked,
                 alpha = cell_rho))



#~ Process velocity ----

rvel.cd <- gene.relative.velocity.estimates(emat,
                                            nmat)

# qs::qsave(rvel.cd, file.path(dir_velocyto, "rvel.cd.qs"))
rvel.cd <- qs::qread(file.path(dir_velocyto, "rvel.cd.qs"))

stopifnot(identical(rownames(dat), colnames(rvel.cd$current)))



# # color-code cell phase per cell
# phase_per_bc <- dat[["cell_phase_masked"]]
# ref_cols <- pals::kovesi.cyclic_mrybm_35_75_c68(50) |> colorRamp()
# 
# 
# 
# phase_per_bc[!is.na(phase_per_bc)] <- (phase_per_bc / 360) |>
#   na.exclude() |>
#   ref_cols() |>
#   (\(x) x / 255 )() |>
#   rgb()
# 
# phase_per_bc[is.na(phase_per_bc)] <- 'grey'
# 
# names(phase_per_bc) <- rownames(dat)




#~ plots ----

res <- show.velocity.on.embedding.cor(as.matrix(dat[,1:2]),
                               rvel.cd,
                               arrow.scale = 3,
                               show.grid.flow = TRUE,
                               grid.n = 40,
                               arrow.lwd = 1.5,
                               do.par = T,
                               cell.border.alpha = 0.1)

# opar <- par(no.readonly = T)
pdf(file.path(dir_velocyto, paste0(ct, "_velocity.pdf")),
    width = 6, height = 6)
show.velocity.on.embedding.cor(as.matrix(dat[,1:2]),
                               rvel.cd,
                               cc = res$cc,
                               arrow.scale = 3,
                               show.grid.flow = TRUE,
                               grid.n = 40,
                               arrow.lwd = 1.5,
                               do.par = T,
                               cell.border.alpha = 0.1)
dev.off()

png(file.path(dir_velocyto, paste0(ct, "_velocity.png")),
    width = 6, height = 6, units = "in",
    res = 300)
show.velocity.on.embedding.cor(as.matrix(dat[,1:2]),
                               rvel.cd,
                               cc = res$cc,
                               arrow.scale = 3,
                               show.grid.flow = TRUE,
                               grid.n = 40,
                               arrow.lwd = 1.5,
                               do.par = T,
                               cell.border.alpha = 0.1)
dev.off()






# Plot ILso ----

# adapted from celocyto_cell_type


message("Starting ", date())

# Inits ----

library(Seurat)
library("velocyto.R")
library(ggplot2)
library(wbData)

gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = NA, gene_id = "nsIs198",
                  symbol = "nsIs198", sequence = "nsIs198",
                  status = "Live", biotype = "protein_coding_gene",
                  name = "nsIs198") |>
  dplyr::mutate(name = dplyr::if_else(name == "E_BE45912.2",
                                      "E-BE45912.2",
                                      name))


library(getopt)



if(! interactive()){
  spec <- matrix(c(
    'dir_step1', 'f', 1, 'character',
    'dir_velocyto', 'o', 1, 'character',
    'i', 'i', 1, 'integer'
  ), byrow=TRUE, ncol=4)
  
  params <- getopt(spec)
  
} else{
  # Options for interactive
  params <- list(
    dir_step1 = "intermediates/2502/250522_step1",
    dir_velocyto = "intermediates/2502/250522_velocyto",
    i = 33
  )
}

# which(cell_types == "ILso")



emat_tot <- qs::qread(file.path(params$dir_velocyto, "emat_tot.qs"))
nmat_tot <- qs::qread(file.path(params$dir_velocyto, "nmat_tot.qs"))

rownames(emat_tot) <- i2s(rownames(emat_tot), gids, warn_missing = TRUE)
rownames(nmat_tot) <- i2s(rownames(nmat_tot), gids, warn_missing = TRUE)



cell_types <- list.files(params$dir_step1, "_seu_unsmoothed\\.qs$") |>
  stringr::str_remove("_seu_unsmoothed\\.qs$")




ct <- cell_types[[params$i]]


message("i: ", params$i,"; cell type: ", ct)

sub <- qs::qread(file.path(params$dir_step1,
                           paste0(ct, "_seu_unsmoothed.qs")) )

stopifnot(all(
  rownames(sub) %in% rownames(emat_tot)
))
stopifnot(all(
  colnames(sub) %in% colnames(emat_tot)
))


emat <- emat_tot[rownames(sub), colnames(sub)]
nmat <- nmat_tot[rownames(sub), colnames(sub)]




#~ process PCA ----

dat <- FetchData(sub, vars = c("PC_1","PC_2","cell_phase_masked", "cell_rho"))

x_axis <- "PC_1"
if(sum(!is.na(dat$cell_phase_masked)) > 20){
  
  
  # if counterclockwise, invert axis for easier interpretation
  cor_angles <- circular::cor.circular(circular::circular(dat$cell_phase_masked, units = "degrees"),
                                       circular::coord2rad(dat$PC_1,
                                                           dat$PC_2))
  
  
  if(cor_angles > 0){
    message("Counterclockwise correlation of angles: ", round(cor_angles, 2),"; inverting PC1")
    dat$PC_1 <- -dat$PC_1
    
    x_axis <- "PC_1 (inverted)"
  }
  
}


dat |>
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 360)) +
  geom_point(aes(x = PC_1, y = PC_2,
                 color = cell_phase_masked),
             alpha = .7, size = 3, shape = 16) +
  ggtitle(ct) +
  xlab(x_axis)

rvel.cd <- qs::qread(file.path(params$dir_velocyto, paste0( ct, "_rvel.cd.qs")))








stopifnot(identical(rownames(dat), colnames(rvel.cd$current)))





#~ plots ----
message("---  Process arrows")

# arrows_preproc <- show.velocity.on.embedding.cor(as.matrix(dat[,1:2]),
#                                                  rvel.cd,
#                                                  arrow.scale = 3,
#                                                  show.grid.flow = TRUE,
#                                                  grid.n = 40,
#                                                  arrow.lwd = 1.5,
#                                                  do.par = T,
#                                                  cell.border.alpha = 0.1)

arrows_preproc <- qs::qread(file.path(params$dir_velocyto, paste0( ct, "_arrows_preproc.qs")))


# # color-code cell phase per cell
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
                               cc = arrows_preproc$cc,
                               arrow.scale = 3,
                               show.grid.flow = TRUE,
                               grid.n = 40,
                               cell.colors = phase_per_bc,
                                 arrow.lwd = 1.5,
                               do.par = T,
                               cell.border.alpha = 0.1)

# color specific cells
cells_of_interest <- rownames(dat)[c(4,242)]

cell_colors <- rep("grey", nrow(dat)) |> setNames(rownames(dat))
cell_colors[c(4,242)] <- "red3"

show.velocity.on.embedding.cor(as.matrix(dat[,1:2]),
                               rvel.cd,
                               cc = arrows_preproc$cc,
                               cell.colors = cell_colors,
                               arrow.scale = 3,
                               show.grid.flow = TRUE,
                               grid.n = 40,
                               arrow.lwd = 1.5,
                               do.par = T,
                               cell.border.alpha = 0.1)


message("---  Save plots")

# pdf(file.path(params$dir_velocyto, paste0(ct, "_velocity.pdf")),
#     width = 6, height = 6)
# show.velocity.on.embedding.cor(as.matrix(dat[,1:2]),
#                                rvel.cd,
#                                cc = arrows_preproc$cc,
#                                arrow.scale = 3,
#                                show.grid.flow = TRUE,
#                                grid.n = 40,
#                                arrow.lwd = 1.5,
#                                do.par = T,
#                                cell.border.alpha = 0.1)
# dev.off()

png(file.path(params$dir_velocyto, paste0(ct, "_velocity.png")),
    width = 6, height = 6, units = "in",
    res = 300)
show.velocity.on.embedding.cor(as.matrix(dat[,1:2]),
                               rvel.cd,
                               cc = arrows_preproc$cc,
                               arrow.scale = 3,
                               show.grid.flow = TRUE,
                               grid.n = 40,
                               cell.colors = 
                                 arrow.lwd = 1.5,
                               do.par = T,
                               cell.border.alpha = 0.1)
dev.off()
















