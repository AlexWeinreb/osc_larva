
# Inits ----

library(tidyverse)
library(Seurat)
library(wbData)

gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = NA, gene_id = "nsIs198",
                  symbol = "nsIs198", sequence = "nsIs198",
                  status = "Live", biotype = "protein_coding_gene",
                  name = "nsIs198")

source("R/utils.R")


dir_third <- "intermediates/2502/250311_third"
dir_processed <- "intermediates/2502/250529_third_processed"



rev_markers <- readxl::read_excel("/gpfs/ycga/work/hammarlund/aw853/references/data/tissue_markers_250305.xlsx",
                                  sheet = "selected") |>
  mutate(gene_id = s2i(gene, gids, warn_missing = TRUE) )


stopifnot(!any(duplicated(rev_markers$gene)))


# check later annotations
seu_annot <- qs::qread( file.path("intermediates/2502/250509_assembled", "250509_seu_all_herma.qs"))
table(seu_annot$cell_type[seu_annot$tissue == "reproductive"])

xx <- names(which(seu_annot$cell_type == "gonadal_sheath")) |> str_remove("^s[0-9]+_")


# transfer the files that are already atomic cell types
atomicfile <- paste0("250311_seu_",cond_here,"_", tissue_here,".qs")
atomicfile
file.exists(file.path(dir_third, atomicfile))
file.copy(file.path(dir_third, atomicfile),
          file.path(dir_processed, atomicfile))

# copy neurons
neuron_files <- c("230313_g2_ACh_motoneuron.qs", "230313_g2_ADL.qs", "230313_g2_AFD.qs", "230313_g2_AIA.qs",
                  "230313_g2_AIB.qs", "230313_g2_AIM_RIP.qs", "230313_g2_AIN.qs", "230313_g2_ALA.qs",
                  "230313_g2_ASEL.qs", "230313_g2_ASER.qs", "230313_g2_ASG.qs", "230313_g2_ASH.qs",
                  "230313_g2_ASJ.qs", "230313_g2_ASK.qs", "230313_g2_AVA.qs", "230313_g2_AVB.qs",
                  "230313_g2_AVD.qs", "230313_g2_AVE.qs", "230313_g2_AVH.qs", "230313_g2_AVJ.qs",
                  "230313_g2_AVK.qs", "230313_g2_AVL.qs", "230313_g2_AWA.qs", "230313_g2_AWC.qs",
                  "230313_g2_BAG.qs", "230313_g2_CEP.qs", "230313_g2_HSN.qs", "230313_g2_IL1.qs",
                  "230313_g2_IL2_LR.qs", "230313_g2_mechanosensory.qs", "230313_g2_OLQ.qs",
                  "230313_g2_pharyngeal.qs", "230313_g2_PVD.qs", "230313_g2_PVQ.qs", "230313_g2_RIA.qs",
                  "230313_g2_RIC.qs", "230313_g2_RID.qs", "230313_g2_RIH.qs", "230313_g2_RIM.qs",
                  "230313_g2_RIR.qs", "230313_g2_RIS.qs", "230313_g2_RIV.qs", "230313_g2_RMD_LR.qs",
                  "230313_g2_RMF.qs", "230313_g2_RMG.qs", "230313_g2_RMH.qs", "230313_g2_SAA.qs",
                  "230313_g2_SIA.qs", "230313_g2_SIB.qs", "230313_g2_SMB.qs", "230313_g2_SMD.qs",
                  "230313_g2_URX.qs", "230313_g2_URY.qs", "230313_g2_VD_DD.qs", "230318_g1_ACh_motoneuron.qs",
                  "230318_g1_AIN.qs", "230318_g1_ALN_PLN_SMB.qs", "230318_g1_ASK.qs", "230318_g1_CAN.qs",
                  "230318_g1_DVB.qs", "230318_g1_mechanosensory.qs", "230318_g1_PDA.qs", "230318_g1_PHC.qs",
                  "230318_g1_PVD.qs", "230318_g1_RMH.qs", "230318_g1_VD_DD.qs")

dir_third_old <- "intermediates/2502/250313_third_processed/"
dir_third_processed
walk(neuron_files,
     \(.file){
       stopifnot(file.exists(file.path(dir_third_old, .file)))
       stopifnot( !file.exists(file.path(dir_third_processed, .file)))
       file.copy(file.path(dir_third_old, .file),
                 file.path(dir_third_processed, .file))
     })


## Save ----

# if not atomic, use rest of script, save here individual objects

atomic_seu <- seu
atomic_seu <- subset(seu, idents = c(9,11,12))

stopifnot(
  ! file.exists(file.path(dir_processed,
                          "230529_g1_ILso.qs"))
)
qs::qsave(atomic_seu,
          file.path(dir_processed,
                    "230529_g1_ILso.qs"))
rm(atomic_seu)


## Load ----
list.files(dir_third, pattern = "250311_seu_1")

cond_here <- "1"
tissue_here <- "glia"

seu <- qs::qread(file.path( dir_third, paste0("250311_seu_",cond_here,"_", tissue_here,".qs") ))

# marks <- qs::qread(file.path( dir_third, paste0("250311_marks_",cond_here,"_", tissue_here,".qs") ))

DimPlot(seu,
        reduction = "umap",
        group.by = "orig.ident",
        label = TRUE,
        pt.size = 2,
        alpha = .2) +
  NoLegend()


DimPlot(seu,
        reduction = "umap",
        label = TRUE,
        pt.size = 2,
        alpha = .5) +
  NoLegend()

# ElbowPlot(seu, ndims = 200) + geom_vline(aes(xintercept = 30))
# seu <- FindNeighbors(seu,
#                      dims = 1:40)
# 
# seu <- FindClusters(seu,
#                     resolution = .6)

# check against previous annotation
xx <- colnames(seu)
xx <- colnames(seu)[seu$seurat_clusters == 0]
seu_annot$cell_type[which(str_remove(colnames(seu_annot), "^s[0-9]+_") %in% xx)] |> table()

xx <- names(which(seu_annot$cell_type == "BWM")) |> str_remove("^s[0-9]+_")
seu$seurat_clusters[which(colnames(seu) %in% xx)] |> table()

DimPlot(seu,
        reduction = "pca") +
  NoLegend()

seu$tmp <- Idents(seu) == "7"
seu$tmp <- colnames(seu) %in% xx
# seu$tmp <- Idents(seu) %in% c(0,2,4,7)
DimPlot(seu,
        reduction = "umap",
        group.by = "tmp",
        label = FALSE,
        pt.size = 1.5,
        alpha = .5) 
  # NoLegend()

# Specific markers ----


aggregate_tib(seu,
              gene_names = rev_markers$gene[rev_markers$tissue == "muscle"] ) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)


aggregate_tib(seu,
              gene_names = rev_markers$gene[(rev_markers$tissue == "muscle" &
                                               rev_markers$cell_type != "pharyngeal_muscle") |
                                              rev_markers$cell_type == "hmc"]) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)


aggregate_tib(seu,
              gene_names = c("nlp-35","nlp-52","nlp-49","nlp-70","F27B10.1","C12D5.3","B0432.14","nlp-55","nlp-18","nhr-236","R10H10.4","pnc-2","ceh-27") ) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)



# Gene sets ----

#~ germline ----
# # germline (see e.g. https://pmc.ncbi.nlm.nih.gov/articles/PMC6216576/)
genes <- c("rpl-30","Y6D1A.1","C07D8.6")

genes <- c("strl-1", "lag-2", "epi-1","aqp-2", "ckb-3")

genes <- rev_markers$gene[rev_markers$cell_type == "germline" ]

genes <- rev_markers$gene[startsWith(rev_markers$cell_type , "gonad") ]

genes <- rev_markers$gene[rev_markers$tissue == "reproductive" ]



#~ repro----
genes <- c("cey-3","pgl-1","fbf-1","fkh-6","lag-2")



#~ sperm ----
genes <- c("nspa-8","msp-59","msp-40","spch-3")


#~ neurons ----
genes <- rev_markers$gene[rev_markers$tissue == "neurons"]



#~ skin ----
genes <- rev_markers$gene[rev_markers$tissue == "skin" &
                            rev_markers$cell_type != "seam" ]

genes <- rev_markers$gene[rev_markers$cell_type == "seam" ]




#in hyp not seam (compiled on Wormbase)
genes <- c("semo-1","dpy-7","osm-8","ceh-14","nekl-3","nekl-2","lrp-1","unc-115","mlt-2",
           "mlt-4","mlt-3","rml-4","bus-5","gmap-1","sma-6","cut-6","nas-37","col-121",
           "vha-17","qua-1","nlp-29","rdy-2","rncs-1","mlt-10")

# the highest from above
genes <- c("mlt-10","nas-37","qua-1","lrp-1","gmap-1","nlp-29")


genes1 <- c("mlt-10","nas-37","qua-1","lrp-1","gmap-1","nlp-29")
genes2 <- rev_markers$gene[rev_markers$cell_type == "seam" ]

#~ XXX ----
genes <- c("eak-3","eak-4","sdf-9","daf-9")



#~ muscle ----
genes <- c("ttr-16","tni-3","tnt-3","egl-20","mls-1", "Y105E8B.9")

# BWM
genes <- c("hlh-1",
           "ttr-16",
           "D1086.12",
           "Y37D8A.2",
           "Y53H1B.2",
           "tnt-2",
           "tni-3"  # <- anterior BWM
)

# smooth muscles
genes <- c("egl-20",
           "ttr-10",
           "Y59C2A.1",
           "Y51A2D.8",
           "W03A5.1",
           "pxn-1",
           " mls-1",
           "C53B7.2",
           "Y17D7C.4"
)


#~ intestine ----
genes <- rev_markers$gene[rev_markers$tissue == "intestine" ] |>
  c("dct-16")


#~ rectal_gland ----
genes <- c("F56C3.8","Y75B12B.3", "ZC53.22","F13B6.3","ifb-1","nex-3")



unique(rev_markers$tissue)


genes <- rev_markers$gene_id[rev_markers$tissue == "pharynx"] |>
  set_names(~ i2s(.x, gids))


genes <- rev_markers$gene_id[rev_markers$tissue == "reproductive"] |>
  set_names(~ i2s(.x, gids))


#~ DTC ----
genes <- c("hlh-12","lag-2","T02E1.8","sex-1","nid-1")

#~ Anchor ----
genes <- c("cdh-3")

#~ hmc ----
genes <- c("glb-26","Y105E8B.9","W04G5.4","Y73B6A.3","twk-9","C28C12.11") |>
  s2i(gids) |>
  set_names(~ i2s(.x, gids))

#~ coelo ----
genes <- c("cup-4", "lgc-26", "unc-122")


#~ excr ----
genes <- c("klo-1", "Y55B1BL.1", "nspc-20","nspc-14")


#~ Pharynx ----

#~ marginal
genes <- c("marg-1", "nspb-11","nspb-12",
           "fipr-1","fipr-2","fipr-3","fipr-4","fipr-6","fipr-7",
           "sms-5", "pgp-14")

#~ arcade non-marg
genes <- c("F40E10.5", "F56D2.3","Y75B7AR.1","F30H5.3",
           "F14B4.1","cht-2","D1054.9","C23H3.9","dpf-6","R01E6.2")

#~ marg non-arcade
genes <- c("fipr-9", "T25B9.1","ugt-45","C54D2.1","fipr-7",
           "W01D2.6","nspb-6","fipr-1","T06E4.10","fipr-6","nspb-12")




#~ pha muscle ----
genes <- c("mlc-1",	"mlc-2", "myo-1", "myo-2", "tnc-2","clik-2","mlc-3","clik-3","hsp-12.2")


#~ pha gland ----
genes <- c("pqn-8","hlh-6","phat-1", "dod-6","phat-4","R06F6.14")

#~ pha epith ----
# note: here epith means `e epithelial` + marginal + arcade
genes <- c("hsp-43", "tat-3", "ajm-1", "marg-1", "pgp-14", "Y75B7AR.1","cht-2",
           "nspb-12","fipr-2")


genes <- c("hsp-43","tat-3","ajm-1","marg-1","let-23","tat-2","cdc-50b","ifa-1","inx-5","pgp-14","sms-5")






#~ glia ----
genes <- c("F16F9.3","F11C7.2") 

genes <- c("grl-18","col-53")

#AMso/PHso
genes <- c("grl-2", "unc-53","mab-31","lin-48")

#AMso only
genes <- c("unc-33","arx-2","grl-12","F44A2.3")

genes <- c("grl-2", "grl-12")

# CEPso/OLso
genes <- c("col-56","col-68")

## From wormglia
# ADEsh
genes <- c("T19B10.2","pugs-4","M03E7.2","lpr-4")
# ADE/PDEso
genes <- c("glam-3","F18E9.4","efhd-1")
# AM/PHsh
genes <- c("Y57G11C.38","K01A6.8","far-8")
# AMso/PHso1
genes <- c("Y11D7A.10","spp-24","pnc-1")
# PHso2
genes <- c("F52E1.5","spig-17","lin-44")
# CEPsh
genes <- c("C06E7.2","aqp-7")
# CEPso
genes <- c("col-174","col-113","mam-8")
# ILsh, OLsh
genes <- c("srr-1","spi-2","ttr-17","acp-7","nas-8")
# ILso
genes <- c("pugs-5","lips-17","srb-16","cyp-36A1")
# OLso
genes <- c("C04F12.5","col-53","tag-290","T26C5.4")
# ADEso
genes <- c("flp-24","F58F12.4","cpg-22")




# markers ILso in L4 herma
genes <- c("col-53","T05C3.6","pugs-3","grl-18","B0393.5","cutl-9","col-177","dpy-31",
           "R10E12.2","cutl-8","F58H1.5","ceh-31","cutl-3","cpg-23","dpy-18","C16B8.2",
           "cpg-22","dhs-5","col-155","col-174","F10D7.11","cdh-3","col-14","Y50E8A.1",
           "ram-5","oac-46","abch-1","igcm-1","col-97","col-33","msa-1","col-109",
           "col-48","cbn-1","cutl-26","dex-1","Y58A7A.2")


genes <- c("col-53","pugs-3","grl-18","B0393.5","cutl-8")

# markers OLso based on L4
genes <- c("C04F12.5","spig-13","cnc-10","lido-8","srp-6","grd-16","zipt-2.2","faxc-1","M60.4","gmap-2")


# markers hyp2 from L4
genes <- c("cutl-10","F49E11.2","nhr-67","F59B10.3","dpf-6","Y57A10A.24",
           "pugs-5","pde-4","F53F4.2","F45E4.5","cutl-16","lon-8","F55H12.4")

# markers tail hyp from L4
genes <- c("col-118","Y65B4BL.6","R12E2.7","col-166","grd-5","col-155")





aggregate_tib(seu2,
              gene_names = genes) |>
  summarize(prop_OLso = prop[cell_id == "11"],
            prop_ILso = median(prop[cell_id %in% 0:2]),
            .by = gene_name) |>
  filter(prop_OLso > .5,
         prop_ILso < .1)


genes <- c("amx-2")



#~~ plot----


genes <- genes[genes %in% rownames(seu)]

genes

FeaturePlot(seu,
            genes, alpha=.5, pt.size = 1.5,
            combine = FALSE) |># lapply(\(.x) .x + NoLegend()) |>
  patchwork::wrap_plots()


table(seu$doubletFinder, Idents(seu), useNA = 'ifany')


plotsum(genes, seu)
plotsum2(genes1, genes2, seu)
plotsum2(genes2, genes1, seu)



aggregate_tib(seu,
              gene_names = genes) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)



# glia
gids$name |> str_subset("^pugs-") |> s2i(gids) |> paste0(collapse = ",") |> message()
genes <- gids$name |> str_subset("^ceh-") |> paste0(collapse = ",") |> message()


wormglia_genes <- c("T19B10.2","pugs-4","M03E7.2","lpr-4", # ADEsh
                    "glam-3","F18E9.4","efhd-1",  # ADE/PDEso
                    "Y57G11C.38","K01A6.8","far-8",  # AM/PHsh
                    "Y11D7A.10","spp-24","pnc-1", # AMso/PHso1
                    "F52E1.5","spig-17","lin-44", # PHso2
                    "C06E7.2","aqp-7",  # CEPsh
                    "col-174","col-113","mam-8", # CEPso
                    "srr-1","spi-2","ttr-17","acp-7","nas-8", #ILsh, OLsh
                    "pugs-5","lips-17","srb-16","cyp-36A1",  # ILso
                    "C04F12.5","col-53","tag-290","T26C5.4",  # OLso
                    "flp-24","F58F12.4","cpg-22")

aggregate_tib(seu,
              gene_names = wormglia_genes) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)

aggregate_tib(seu,
              gene_names = rev_markers$gene[rev_markers$tissue == "glia"] |>
                union(c("gegf-1", "spig-2", "ttr-43", "txt-17", #sheath
                        "fig-1","nas-31","hmit-1.2","F52E1.2","vap-1",
                        "K02E11.4","R11D1.3","bgnt-1.1","toh-1",
                        "unc-53","mab-31","lite-1","alr-1","lin-48","arx-2","F44A2.3",
                        "hlh-17","inx-3","fmi-1","nhr-57","swip-10","twk-16"))) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)



aggregate_tib(seu,
              gene_names = genes) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)




# pharynx
aggregate_tib(seu,
              gene_names = c("mlc-1",	"mlc-2", "myo-1", "myo-2", "tnc-2","clik-2","mlc-3","clik-3","hsp-12.2",
                             "pqn-8","hlh-6","phat-1", "dod-6","phat-4","R06F6.14",
                             "hsp-43","tat-3","ajm-1","marg-1","let-23","tat-2","cdc-50b","ifa-1","inx-5","pgp-14","sms-5")) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)





# reverse approach ----


aggregate_tib(seu,
              gene_names = rev_markers$gene[1:57]) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)

aggregate_tib(seu,
              gene_names = rev_markers$gene[58:117]) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)




# plot annotation ----


annot <- read_csv(file.path(dir_second,
                            paste0("250124_",stage_here,"_tissues.csv"))) |>
  filter(first_pass == tissue_here) |>
  column_to_rownames("cluster")


head(annot)

seu$annot_second <- annot$cell_type[seu$seurat_clusters]


DimPlot(seu,
        reduction = "umap",
        group.by = "annot_second",
        label = TRUE,
        pt.size = 2,
        alpha = .2,
        repel = TRUE) +
  NoLegend() +
  ggtitle(NULL)





# pairwise comparison ----
mark_sing <- FindMarkers(seu,
                         ident.1 = c(0,5,13,16), ident.2 = c(17) )

mark_sing <- FindMarkers(seu,
                         ident.1 = c(20) )

mark_sing |>
  rownames_to_column("gene_name") |>
  as_tibble() |>
  filter(avg_log2FC > 1,
         pct.1 > .3,
         pct.2 < .2) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name) |> # head(6) -> genes
  paste0(collapse = ", ") |> message()

mark_sing |>
  rownames_to_column("gene_name") |>
  as_tibble() |>
  filter(avg_log2FC < -1,
         pct.1 < .15,
         pct.2 > .3) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name) |> paste0(collapse = " ") |> message()




# named genes
mark_sing |>
  rownames_to_column("gene_name") |>
  as_tibble() |>
  filter(avg_log2FC > 1,
         pct.1 > .3,
         pct.2 < .2) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name) |> str_subset("\\-") |> paste0(collapse = ", ") |> message()

mark_sing |>
  rownames_to_column("gene_name") |>
  as_tibble() |>
  filter(avg_log2FC > 1,
         pct.1 > .3,
         pct.2 < .2) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name) |> str_subset("\\-") |> head(6) -> genes


#~ TEA  ----

genelist <- mark_sing |>
  rownames_to_column("gene_name") |>
  as_tibble() |>
  filter(avg_log2FC > 1,
         pct.1 > .3,
         pct.2 < .2) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name)

length(genelist)

# wbe_dict <- wormbaseEnrich::fetch_dictionary("tissue")
# wbe_dict_go <- wormbaseEnrich::fetch_dictionary("GO")
tea_res <- wormbaseEnrich::enrichment_analysis(s2i(genelist, gids), wbe_dict,
                                               background_genes = rownames(seu) |> s2i(gids))



tea_res |>
  mutate(progenitor = startsWith(term_name, "AB")) |>
  # filter(! startsWith(term_name, "AB")) |>
  ggplot() +
  theme_classic() +
  scale_alpha_manual(values = c(`TRUE` = .2, `FALSE` = 1)) +
  aes(x = observed, y = -log10(FDR),
      alpha = progenitor,
      label = term_name) +
  geom_point(aes(size = enrichment_fc)) +
  ggrepel::geom_text_repel()

tea_res |> arrange(desc(observed))
tea_res |> arrange(FDR)


# Genes explanation
wbe_dict$wbid[as.logical(wbe_dict$`vulE WBbt:0006767`)] |>
  i2s(gids) |>
  intersect(genelist) |>
  paste(collapse = ", ")



#~ TEA from marker precomp ----
marks <- qs::qread(file.path( dir_third, paste0("250311_marks_",cond_here,"_", tissue_here,".qs") ))

genelist <- marks |>
  filter(cluster == 12) |>
  filter(avg_log2FC > 1,
         pct.1 > .3,
         pct.2 < .2) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name)

length(genelist)

# wbe_dict <- wormbaseEnrich::fetch_dictionary("tissue")
# wbe_dict_go <- wormbaseEnrich::fetch_dictionary("GO")
tea_res <- wormbaseEnrich::enrichment_analysis(s2i(genelist, gids), wbe_dict,
                                               background_genes = rownames(seu) |> s2i(gids))



tea_res |>
  mutate(progenitor = startsWith(term_name, "AB")) |>
  # filter(! startsWith(term_name, "AB")) |>
  ggplot() +
  theme_classic() +
  scale_alpha_manual(values = c(`TRUE` = .2, `FALSE` = 1)) +
  aes(x = observed, y = -log10(FDR),
      alpha = progenitor,
      label = term_name) +
  geom_point(aes(size = enrichment_fc)) +
  ggrepel::geom_text_repel()

tea_res |> arrange(desc(observed))
tea_res |> arrange(FDR)


# Genes explanation
wbe_dict$wbid[as.logical(wbe_dict$`excretory socket cell WBbt:0004534`)] |>
  i2s(gids) |>
  intersect(genelist) |>
  paste(collapse = ", ")






# test subset ----

sub <- subset(seu, idents = setdiff(0:8, c(0,6)))
sub <- subset(seu, idents = c(9,11,12))
# sub <- seu


sub <- SCTransform(sub)

maxnpcs <- pmin(200, ncol(sub) - 10)
sub <- RunPCA(sub, npcs = maxnpcs, verbose = FALSE)

npca <- 15

ElbowPlot(sub, ndims = maxnpcs) +
  geom_vline(aes(xintercept = npca))




sub <- RunUMAP(sub,
               dims = 1:npca,
               n.neighbors = 30)


DimPlot(
  sub,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 2,
  alpha = .3
) +
  NoLegend()

sub <- FindNeighbors(sub,
                     dims = 1:npca)

sub <- FindClusters(sub,
                    resolution = .5)

DimPlot(
  sub,
  reduction = "umap",
  label = TRUE,
  pt.size = 2,
  alpha = .2
) +
  NoLegend()

# qs::qsave(sub, file.path(dir_processed, "250318_subglia_g1.qs"))
# seu <- qs::qread(file.path(dir_processed, "250318_subglia_g1.qs"))
# seu <- qs::qread(file.path(dir_processed, "250317_subglia_g2.qs"))



#~ markers ----
genes <- genes[genes %in% rownames(sub)]

genes

FeaturePlot(sub,
            genes, alpha=.7, pt.size = 1.5,
            combine = FALSE) |># lapply(\(.x) .x + NoLegend()) |>
  patchwork::wrap_plots()



aggregate_tib(sub,
              gene_names = genes) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)


all_markers <- FindAllMarkers(sub,
                              logfc.threshold = 1,
                              min.pct = .3) |>
  as_tibble() |>
  rename(gene_name = gene)


all_markers |>
  filter(cluster == "2") |>
  filter(p_val_adj < 0.05,
         abs(pct.1 - pct.2) > 0.3,
         avg_log2FC > 0) |>
  pull(gene_name) |> paste0(collapse = " ")



table(sub$doubletFinder, Idents(sub), useNA = 'ifany')




#~ pairwise in sub ----


mark_sing <- FindMarkers(sub,
                         ident.1 = c(3) )

mark_sing |>
  rownames_to_column("gene_name") |>
  as_tibble() |>
  filter(avg_log2FC > 1,
         pct.1 > .3,
         pct.2 < .2) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name) |> paste0(collapse = ", ") |> message()


# named genes
mark_sing |>
  rownames_to_column("gene_name") |>
  as_tibble() |>
  filter(avg_log2FC > 1,
         pct.1 > .3,
         pct.2 < .2) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name) |> str_subset("\\-") |> head(15) |> paste0(collapse = ", ") |> message()

# as markers
genes <- mark_sing |>
  rownames_to_column("gene_name") |>
  as_tibble() |>
  filter(avg_log2FC > 1,
         pct.1 > .3,
         pct.2 < .2) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name) |> head()

FeaturePlot(sub,
            genes, alpha=.7, pt.size = 1.5,
            combine = FALSE) |># lapply(\(.x) .x + NoLegend()) |>
  patchwork::wrap_plots()

#~~ TEA  ----

genelist <- mark_sing |>
  rownames_to_column("gene_name") |>
  as_tibble() |>
  filter(avg_log2FC > 1,
         pct.1 > .3,
         pct.2 < .2) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name)

length(genelist)

# wbe_dict <- wormbaseEnrich::fetch_dictionary("tissue")
# wbe_dict_go <- wormbaseEnrich::fetch_dictionary("GO")
tea_res <- wormbaseEnrich::enrichment_analysis(s2i(genelist, gids), wbe_dict,
                                               background_genes = rownames(seu) |> s2i(gids))



tea_res |>
  mutate(progenitor = startsWith(term_name, "AB")) |>
  # filter(! startsWith(term_name, "AB")) |>
  ggplot() +
  theme_classic() +
  scale_alpha_manual(values = c(`TRUE` = .2, `FALSE` = 1)) +
  aes(x = observed, y = -log10(FDR),
      alpha = progenitor,
      label = term_name) +
  geom_point(aes(size = enrichment_fc)) +
  ggrepel::geom_text_repel()

tea_res |> arrange(desc(observed))
tea_res |> arrange(FDR)



# Genes explanation
wbe_dict$wbid[as.logical(wbe_dict$`ILshVL WBbt:0004525` | wbe_dict$`ILshDR WBbt:0004531`)] |>
  i2s(gids) |>
  intersect(genelist) |>
  paste(collapse = ", ")




# test integration ----

sub <- split(sub, f = seu$orig.ident)



sub <- NormalizeData(sub)
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub, verbose = FALSE)

ElbowPlot(sub, ndims = 50)



#~ Integrate ----

sub <- IntegrateLayers(object = sub,
                       method = HarmonyIntegration)



sub <- FindNeighbors(sub,
                     dims = 1:20,
                     reduction = "harmony")

sub <- FindClusters(sub,
                    cluster.name = "harmony_clusters",
                    resolution = .3)

sub <- RunUMAP(sub,
               dims = 1:20,
               reduction = "harmony",
               reduction.name = "umapharmony")


DimPlot(
  sub,
  reduction = "umapharmony",
  group.by = "orig.ident",
  pt.size = 2,
  alpha = .3
) +
  NoLegend()

DimPlot(
  sub,
  reduction = "umapharmony",
  group.by = "harmony_clusters",
  label = TRUE,
  pt.size = 2,
  alpha = .2
) +
  NoLegend()


genes <- rev_markers$gene_id[rev_markers$cell_type == "anal_muscle"] |>
  set_names(~ i2s(.x, gids))

genes <- c("T20G5.8", "F23H11.7", "unc-87","clik-1") |>
  s2i(gids) |>
  set_names(~ i2s(.x, gids))

#~~ plot----


genes <- genes[genes %in% rownames(sub)]

genes

FeaturePlot(sub,
            reduction = "umapharmony",
            genes, alpha=.2,pt.size = 1.5,
            combine = FALSE) |>
  map2(names(genes),
       ~ {.x + ggtitle(.y) + xlab(NULL) + ylab(NULL) + NoLegend()
       }) |>
  patchwork::wrap_plots()




aggregate_tib(sub,
              gene_names = rev_markers$gene[rev_markers$tissue == "muscle"]) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)



#DEGs----


sub <- JoinLayers(sub)



all_markers <- FindAllMarkers(sub,
                              logfc.threshold = 1,
                              min.pct = .3) |>
  as_tibble() |>
  rename(gene_id = gene) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE))


all_markers |>
  filter(cluster == "3") |>
  filter(p_val_adj < 0.05,
         abs(pct.1 - pct.2) > 0.3,
         avg_log2FC > 0) |>
  pull(gene_name) |> paste0(collapse = " ")


FeaturePlot(sub,
            features = c("egl-20", "col-118", "C39E9.8","W03F9.11") |> s2i(gids))




#~ pairwise comparison ----
mark_sing <- FindMarkers(sub,
                         ident.1 = c(5), ident.2 = c(0,2,10,11,6) )



mark_sing |>
  rownames_to_column("gene_id") |>
  mutate(gene_name = i2s(gene_id, gids)) |>
  as_tibble() |>
  filter(avg_log2FC > 1,
         pct.1 > .3,
         pct.2 < .15) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name) |> paste0(collapse = " ")

mark_sing |>
  rownames_to_column("gene_id") |>
  mutate(gene_name = i2s(gene_id, gids)) |>
  as_tibble() |>
  filter(avg_log2FC < -1,
         pct.1 < .15,
         pct.2 > .3) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name) |> paste0(collapse = " ")


# hclust ----


mat <- AggregateExpression(seu, assays = "RNA")[["RNA"]] |> as.matrix() |> log1p()

head(mat)
hist(mat)
dim(mat)
mat2 <- mat[VariableFeatures(seu),]
dim(mat2)



hc_genes <- hclust(as.dist((1 - cor(t(mat2)))/2))
hc_clusts <- hclust(as.dist((1 - cor(mat2))/2))

pheatmap::pheatmap(mat2,
                   cluster_rows = hc_genes,
                   cluster_cols = hc_clusts)

clg <- cutree(hc_genes, h=0.5)

table(clg)

mat3_list <- lapply(unique(clg),
                    \(.cl){
                      colSums(mat2[clg == .cl,])
                    })

mat3 <- do.call(rbind, mat3_list)
rownames(mat3) <- paste0("clg", unique(clg))



ann <- enframe(clg, value = "clg") |>
  column_to_rownames("name") |>
  mutate(clg = paste0("clg", clg)) |>
  mutate(is_8 = 1*(clg == "clg8"),
         is_49 = 1*(clg == "clg49"),
         is_47 = 1*(clg == "clg47"),
         is_37 = 1*(clg == "clg37"),
         is_35 = 1*(clg == "clg35"))


pheatmap::pheatmap(mat2,
                   cluster_rows = hc_genes,
                   cluster_cols = hc_clusts,
                   annotation_row = ann)


pheatmap::pheatmap(mat3,
                   cluster_cols = hc_clusts)


names(clg)[clg == 8]
names(clg)[clg == 35] |> paste0(collapse = ",") |> message()




tea_res <- wormbaseEnrich::enrichment_analysis(s2i(names(clg)[clg == 35], gids), wbe_dict,
                                               background_genes = rownames(seu) |> s2i(gids))

tea_res #|> arrange(desc(observed))








# add phase ----
osc <- readxl::read_excel("../10x_grl18/data/oscillating/msb209498-sup-0003-datasetev1.xlsx",
                          sheet = "Dataset EV1 WBidToGeneNames_Osc",
                          na = "NA")

mat <- LayerData(seu)[osc$WB_ID[osc$Class == "Osc"] |> i2s(gids) |> intersect(rownames(seu)),]
dim(mat)






osc$PeakPhase_rad <- osc$PeakPhase*pi/180
osc2 <- osc[match(s2i(rownames(mat), gids), osc$WB_ID),]
osc2$gene_name <- i2s(osc2$WB_ID, gids)

stopifnot(all.equal( osc2$gene_name, rownames(mat) ))

mean_angles <- tibble(cell_bc = colnames(mat),
                      mean_angle = NA_real_,
                      mean_rho = NA_real_)


for(i in which(is.na(mean_angles$mean_angle))){
  
  expression <- mat[, i ]
  
  mean_angles$mean_angle[[i]] <- circhelp::weighted_circ_mean(osc2$PeakPhase_rad,
                                                              expression) %% (2*pi)
  mean_angles$mean_rho[[i]] <- max(expression)*circhelp::weighted_circ_rho(osc2$PeakPhase_rad,
                                                                           expression) %% (2*pi)
  
}







stopifnot(identical(seu[[]] |> rownames(),
                    mean_angles$cell_bc))

# list(seu = seu[[]] |> rownames(),
#      ang = mean_angles$cell_bc) |>
#   eulerr::euler() |>
#   plot(quantities = TRUE)


seu$mean_angle <- mean_angles$mean_angle
seu$mean_rho <- mean_angles$mean_rho


ggplot(FetchData(seu,
                 vars = c("umap_1", "umap_2",
                          "mean_angle", "mean_rho"))) +
  theme_classic() +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 2*pi)) +
  geom_point(aes(x = umap_1, y = umap_2,
                 color = mean_angle,
                 alpha = mean_rho))







# Add identities ----


#~ group 1 glia ----

annot_g1_glia <- tibble(
  cell_bc = colnames(seu), glia_type = NA_character_
)
table(annot_g1_glia$glia_type, useNA = 'ifany')



# not subset
annot_g1_glia$glia_type[annot_g1_glia$cell_bc %in% colnames(seu)[Idents(seu) %in% c("0")]] <- "PHsh"

annot_g1_glia$glia_type[annot_g1_glia$cell_bc %in% colnames(seu)[Idents(seu) %in% c("6")]] <- "AMsh"



annot_g1_glia$glia_type[annot_g1_glia$cell_bc %in% colnames(seu)[Idents(seu) %in% c("11")]] <- "skin_2"





annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(seu)[Idents(seu) %in% c("9")]] <- "glia_socket_1"

# subset sockets setdiff(0:13, c(4,7,8,10,11))
annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(13)]] <- "glia_1"

annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(12)]] <- "glia_sheath_1"
annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(3,5,7)]] <- "glia_sheath_2"
annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(16)]] <- "glia_2"

annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(4)]] <- "unclear"
annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(6,9,10,15)]] <- "CEPso"

annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(0,1,2,8,11,14)]] <- "ILso"

# qs::qsave(annot_g2_glia, file.path(dir_processed, "250529_annot_g2_glia.qs"))


# #~ group 2 glia ----
# 
# annot_g2_glia <- tibble(
#   cell_bc = colnames(seu), glia_type = NA_character_
# )
# table(annot_g2_glia$glia_type, useNA = 'ifany')
# 
# 
# # subset AM/PHsh (seu cl 7)
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(1:4)]] <- "AMsh"
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c("0")]] <- "PHsh"
# 
# # not subset
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(seu)[Idents(seu) %in% c("8")]] <- "AM_PHso"
# 
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(seu)[Idents(seu) %in% c("10")]] <- "skin_3"
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(seu)[Idents(seu) %in% c("11")]] <- "skin_2"
# 
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(seu)[Idents(seu) %in% c("4")]] <- "CEPsh"
# 
# 
# 
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(seu)[Idents(seu) %in% c("9")]] <- "glia_socket_1"
# 
# # subset sockets setdiff(0:13, c(4,7,8,10,11))
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(13)]] <- "glia_1"
# 
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(12)]] <- "glia_sheath_1"
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(3,5,7)]] <- "glia_sheath_2"
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(16)]] <- "glia_2"
# 
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(4)]] <- "unclear"
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(6,9,10,15)]] <- "CEPso"
# 
# annot_g2_glia$glia_type[annot_g2_glia$cell_bc %in% colnames(sub)[Idents(sub) %in% c(0,1,2,8,11,14)]] <- "ILso"
# 
# # qs::qsave(annot_g2_glia, file.path(dir_processed, "250529_annot_g2_glia.qs"))







# save each glial type
annot_g2_glia <- annot_g2_glia |> column_to_rownames("cell_bc")


# seu$annot_third_glia <- annot_g2_glia[colnames(seu), "glia_type"]
Idents(seu) <- "annot_third_glia"



#~ save glia ----

levels(Idents(seu)) |>
  setdiff("unclear") |>
  walk(~{
    filename <- file.path(dir_processed,
                          paste0("230529_g2_", .x, ".qs"))

    sub <- subset(seu, idents = .x)


    if(!file.exists(filename)){
      qs::qsave(sub,
                filename)
      message("saved")
    } else{
      stop("exists!")
    }
  })




full_join(
seu_annot$cell_type[seu_annot$tissue == "glia"] |>
  enframe(value = "cell_type",
          name = "cell_bc") |>
  mutate(cell_bc = str_remove(cell_bc, "^s[0-9]+_"))
,
seu$annot_third_glia |>
  enframe(value = "cell_type",
          name = "cell_bc")
 ,
by = "cell_bc"
) |> filter(is.na(cell_type.y)) |> pull(cell_type.x) |> table()



seu_annot[[]] |>
  rownames_to_column("cell_bc") |>
  mutate(cell_bc = str_remove(cell_bc, "^s[0-9]+_")) |>
  as_tibble() |>
  select(cell_bc, cell_type, tissue, orig.ident) |>
  filter(tissue == "glia") |>
  anti_join(seu$annot_third_glia |>
              enframe(value = "cell_type",
                      name = "cell_bc"),
            by = "cell_bc"
  ) |>
  count(orig.ident)








