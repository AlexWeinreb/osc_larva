# replot velocyto

# Inits ----

library(Seurat)
library("velocyto.R")
library(ggplot2)

opar <- par(no.readonly = TRUE)

dir_fig_velocyto <- "presentations/figures/260521_velocyto"
in2mm <- 0.03937008

params <- list(dir_step1 = "intermediates/2502/250609_step1",
               dir_velocyto = "intermediates/2502/250825_velocyto")





cell_types <- c("ILso", "pharynx_epithelial", "seam", "BWM", "SMB")
for(ct in cell_types){
  
  # Load preprocessed data ----
  
  #~ load velocyto preproc ----
  
  rvel.cd <- qs::qread(file.path(params$dir_velocyto, paste0( ct, "_rvel.cd.qs")))
  
  arrows_preproc <- qs::qread(file.path(params$dir_velocyto, paste0( ct, "_arrows_preproc.qs")))
  
  
  
  
  #~ PCA ----
  sub <- qs::qread(file.path(params$dir_step1,
                             paste0(ct, "_seu_unsmoothed.qs")) )[,colnames(rvel.cd$cellKNN)]
  
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
  
  
  
  # replot PCA
  
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
  
  
  stopifnot(identical(rownames(dat), colnames(rvel.cd$current)))
  
  
  # # color-code cell phase per cell
  phase_per_bc <- dat[["cell_phase_masked"]]
  ref_cols <- pals::kovesi.cyclic_mrybm_35_75_c68(50) |> colorRamp()
  
  
  
  phase_per_bc[!is.na(phase_per_bc)] <- (phase_per_bc / 360) |>
    na.exclude() |>
    ref_cols() |>
    (\(x) x / 255 )() |>
    rgb() |>
    alpha(0.8)
  
  phase_per_bc[is.na(phase_per_bc)] <- 'grey'
  
  names(phase_per_bc) <- rownames(dat)
  
  phases_all_grey <- rep("grey", nrow(dat))
  names(phases_all_grey) <- rownames(dat)
  
  
  
  
  #~~ All grey
  pdf(file.path(dir_fig_velocyto, paste0(ct, "_velocity_grey.pdf")),
      width = 50 * in2mm, height = 47.7 * in2mm, pointsize = 10)
  
  par(
    mar = c(1.8, 1.8, 0.3, 0.3),
    mgp = c(0,.2,0),
    cex.axis = 0.7,
    cex.lab = 0.8,
    tck = -0.015
  )
  
  show.velocity.on.embedding.cor(as.matrix(dat[,1:2]),
                                 rvel.cd,
                                 cc = arrows_preproc$cc,
                                 arrow.scale = 5,
                                 show.grid.flow = TRUE,
                                 grid.n = 30,
                                 cell.colors = phases_all_grey,
                                 arrow.lwd = 1,
                                 do.par = F,
                                 cell.border.alpha = 0,
                                 cex = .5)
  dev.off()
  
  
  
  #~~ Color ----
  pdf(file.path(dir_fig_velocyto, paste0(ct, "_velocity_phase.pdf")),
      width = 50 * in2mm, height = 47.7 * in2mm, pointsize = 10)
  
  par(
    mar = c(1.8, 1.8, 0.3, 0.3),
    mgp = c(0,.2,0),
    cex.axis = 0.7,
    cex.lab = 0.8,
    tck = -0.015
  )
  
  show.velocity.on.embedding.cor(as.matrix(dat[,1:2]),
                                 rvel.cd,
                                 cc = arrows_preproc$cc,
                                 arrow.scale = 5,
                                 show.grid.flow = TRUE,
                                 grid.n = 30,
                                 cell.colors = phase_per_bc,
                                 arrow.lwd = 1,
                                 do.par = F,
                                 cell.border.alpha = 0,
                                 cex = .5)
  dev.off()
  
  
}






