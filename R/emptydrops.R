

message("Starting ", date())


library(DropletUtils) |> suppressPackageStartupMessages()
library(Seurat) |> suppressPackageStartupMessages()


dir_aligned <- "data/250225_aligned"
dir_emptydrops <- "intermediates/2502/emptydrops"
dir_filtered <- "intermediates/2502/empytdrops_filt"

gids <- wbData::wb_load_gene_coords(295)


samples <- list.files(dir_aligned)

# i <- 1

for(sample in samples){
  # sample <- "L4_1"
  
  #i <- i+1
  #sample <- samples[[i]]
  #sample
  
  message("### Processing sample ", sample, " ###")
  
  
  dat <- read10xCounts(file.path(dir_aligned, sample, "outs", "raw_feature_bc_matrix"))
  
  
  
  
  
  br.out <- barcodeRanks(dat)
  
  # Making a plot
  png(file.path(dir_emptydrops, paste0(sample, "_rank.png")))
   plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
   o <- order(br.out$rank)
   lines(br.out$rank[o], br.out$fitted[o], col="red")
   
   abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
   abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
   legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
          legend=c("knee", "inflection"))
  dev.off()
  

  set.seed(100)
  e.out <- qs::qread(file.path(dir_emptydrops, paste0(sample,"_e_out.qs")))
  # e.out <- emptyDrops(dat, lower = 100)
  # e.out
  
  is.cell <- e.out$FDR <= 0.01
  sum(is.cell, na.rm=TRUE)
  
  table(Limited=e.out$Limited, Significant=is.cell)
  
 # plot(log10(e.out$Total), log10(-e.out$LogProb), col=ifelse(is.cell, "red", "black"),
 #      xlab="Total UMI count", ylab="-Log Probability")

  # png(file.path(dir_emptydrops, paste0(sample, "_rank_total.png")))
  #   plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total",
  #        col=ifelse(is.cell, "red", "black"))
  # dev.off()
  
  
  mat <- counts(dat)
  colnames(mat) <- dat$Barcode
  
  

  write10xCounts(file.path(dir_filtered,
                           paste0(sample,"filt_feature_bc_matrix")),
                 x = mat[,which(is.cell)],
                 version = "3")
}




message("Done ", date())

message(toLatex(sessionInfo()))


