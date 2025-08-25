
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
    dir_step1 = "intermediates/2502/250609_step1",
    dir_velocyto = "intermediates/2502/250825_velocyto",
    i = 42
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
                           paste0(ct, "_seu_unsmoothed.qs")) )[,colnames(emat_tot)]

stopifnot(all(
  rownames(sub) %in% rownames(emat_tot)
))
stopifnot(all(
  colnames(sub) %in% colnames(emat_tot)
))


emat <- emat_tot[rownames(sub), colnames(sub)]
nmat <- nmat_tot[rownames(sub), colnames(sub)]




#~ save PCA ----

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


gg_phase <- dat |>
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

ggsave(paste0(ct,"_phase.png"), gg_phase,
       path = params$dir_velocyto,
       width = 6, height = 6, units = "in")

ggsave(paste0(ct,"_phase.pdf"), gg_phase,
       path = params$dir_velocyto,
       width = 6, height = 6, units = "in")





# Process velocity ----
message("---  Process velocity")
rvel.cd <- gene.relative.velocity.estimates(emat,
                                            nmat)

qs::qsave(rvel.cd, file.path(params$dir_velocyto, paste0( ct, "_rvel.cd.qs")))


stopifnot(identical(rownames(dat), colnames(rvel.cd$current)))





#~ plots ----
message("---  Process arrows")

arrows_preproc <- show.velocity.on.embedding.cor(as.matrix(dat[,1:2]),
                                                 rvel.cd,
                                                 arrow.scale = 3,
                                                 show.grid.flow = TRUE,
                                                 grid.n = 40,
                                                 arrow.lwd = 1.5,
                                                 do.par = T,
                                                 cell.border.alpha = 0.1)

qs::qsave(arrows_preproc, file.path(params$dir_velocyto, paste0( ct, "_arrows_preproc.qs")))



message("---  Save plots")

#~~ all grey ----
pdf(file.path(params$dir_velocyto, paste0(ct, "_velocity.pdf")),
    width = 6, height = 6)
show.velocity.on.embedding.cor(as.matrix(dat[,1:2]),
                               rvel.cd,
                               cc = arrows_preproc$cc,
                               arrow.scale = 3,
                               show.grid.flow = TRUE,
                               grid.n = 40,
                               arrow.lwd = 1.5,
                               do.par = T,
                               cell.border.alpha = 0.1)
dev.off()

png(file.path(params$dir_velocyto, paste0(ct, "_velocity.png")),
    width = 6, height = 6, units = "in",
    res = 300)
show.velocity.on.embedding.cor(as.matrix(dat[,1:2]),
                               rvel.cd,
                               cc = arrows_preproc$cc,
                               arrow.scale = 3,
                               show.grid.flow = TRUE,
                               grid.n = 40,
                               arrow.lwd = 1.5,
                               do.par = T,
                               cell.border.alpha = 0.1)
dev.off()





#~~ phase-coded phase ----
message(" Color-coded phase")


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




png(file.path(params$dir_velocyto, paste0(ct, "_velocity_phase.png")),
    width = 6, height = 6, units = "in",
    res = 300)
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
dev.off()


pdf(file.path(params$dir_velocyto, paste0(ct, "_velocity_phase.pdf")),
    width = 6, height = 6)
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
dev.off()





message("Done ", date())
message("-------------------------------------")
sessionInfo()



