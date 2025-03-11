
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


dir_second <- "intermediates/2502/250306_per_condition_tissue"




rev_markers <- readxl::read_excel("/gpfs/ycga/work/hammarlund/aw853/references/data/tissue_markers_250305.xlsx",
                                  sheet = "selected") |>
  mutate(gene_id = s2i(gene, gids, warn_missing = TRUE) )


stopifnot(!any(duplicated(rev_markers$gene)))



## Load ----
cond_here <- "4"
tissue_here <- "skin"

seu <- qs::qread(file.path( dir_second, paste0("250306_seu_",cond_here,"_", tissue_here,".qs") ))

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




seu$tmp <- Idents(seu) == "56"
# seu$tmp <- Idents(seu) %in% c(0,2,4,7)
DimPlot(seu,
        reduction = "umap",
        group.by = "tmp",
        label = FALSE,
        pt.size = 1.5,
        alpha = .5) +
  NoLegend()

# Specific markers ----


aggregate_tib(seu,
              gene_names = rev_markers$gene[rev_markers$tissue == "reproductive"] ) |>
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


#~ sperm ----
genes <- c("nspa-8","msp-59","msp-40","spch-3")


#~ skin ----
genes <- rev_markers$gene[rev_markers$tissue == "skin" &
                            rev_markers$cell_type != "seam" ]

genes <- rev_markers$gene[rev_markers$cell_type == "seam" ]


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


#~ repro----
genes <- c("cey-3","pgl-1","fbf-1","fkh-6","lag-2") |>
  s2i(gids) |>
  set_names(~ i2s(.x, gids))

#~ DTC ----
genes <- c("hlh-12","lag-2","T02E1.8","sex-1","nid-1") |>
  s2i(gids) |>
  set_names(~ i2s(.x, gids))


#~ hmc ----
genes <- c("glb-26","Y105E8B.9","W04G5.4","Y73B6A.3","twk-9","C28C12.11") |>
  s2i(gids) |>
  set_names(~ i2s(.x, gids))

#~ coelo ----
genes <- c("cup-4", "lgc-26", "unc-122") |>
  s2i(gids) |>
  set_names(~ i2s(.x, gids))


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





genes <- c("rbm-3.1")



#~~ plot----


genes <- genes[genes %in% rownames(seu)]

genes

FeaturePlot(seu,
            genes, alpha=.7, pt.size = 1.5,
            combine = FALSE) |># lapply(\(.x) .x + NoLegend()) |>
  patchwork::wrap_plots()


table(seu$doubletFinder, Idents(seu))




# glia
gids$name |> str_subset("^pugs-") |> s2i(gids) |> paste0(collapse = ",") |> message()

aggregate_tib(seu,
              gene_names = rev_markers$gene[rev_markers$tissue == "glia"] |>
                union(c("gegf-1", "spig-2", "ttr-43",
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
              gene_names = c("nsy-7","ntc-1","K02A6.2","srt-42","srt-28","dac-1","srsx-5","srx-2")) |>
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
                         ident.1 = c(8,4,12), ident.2 = c(0) )

mark_sing <- FindMarkers(seu,
                         ident.1 = c(17) )

mark_sing |>
  rownames_to_column("gene_name") |>
  as_tibble() |>
  filter(avg_log2FC > 1,
         pct.1 > .3,
         pct.2 < .2) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name) |> paste0(collapse = " ") |> message()

mark_sing |>
  rownames_to_column("gene_name") |>
  as_tibble() |>
  filter(avg_log2FC < -1,
         pct.1 < .15,
         pct.2 > .3) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name) |> paste0(collapse = " ") |> message()



#~ TEA  ----

genelist <- mark_sing |>
  rownames_to_column("gene_name") |>
  as_tibble() |>
  filter(avg_log2FC > 1,
         pct.1 > .2,
         pct.2 < .1) |>
  arrange(desc(avg_log2FC)) |> pull(gene_name)

length(genelist)

# wbe_dict <- wormbaseEnrich::fetch_dictionary("tissue")
tea_res <- wormbaseEnrich::enrichment_analysis(s2i(genelist, gids), wbe_dict,
                                               background_genes = rownames(seu) |> s2i(gids))

tea_res |> arrange(desc(observed))
tea_res |> arrange(FDR)

tea_res |>
  filter(! startsWith(term_name, "AB")) |>
  ggplot() +
  theme_classic() +
  aes(x = observed, y = -log10(FDR), label = term_name) +
  geom_point() +
  ggrepel::geom_text_repel()



# test subset ----

sub <- subset(seu, idents = setdiff(0:6, c("4")))
sub <- subset(seu, idents = c("20","6","15","5","21"))



sub <- SCTransform(sub)

sub <- RunPCA(sub, npcs = 200, verbose = FALSE)

ElbowPlot(sub, ndims = 200)



sub <- FindNeighbors(sub,
                     dims = 1:40)

sub <- FindClusters(sub,
                    resolution = .8)

sub <- RunUMAP(sub,
               dims = 1:40)


DimPlot(
  sub,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 2,
  alpha = .3
) +
  NoLegend()

DimPlot(
  sub,
  reduction = "umap",
  label = TRUE,
  pt.size = 2,
  alpha = .2
) +
  NoLegend()


genes <- rev_markers$gene_id[rev_markers$tissue == "pharynx" | rev_markers$cell_type == "pharyngeal_muscle"] |>
  set_names(~ i2s(.x, gids))

genes <- genes[genes %in% rownames(sub)]



FeaturePlot(sub,
            reduction = "umap",
            genes,
            combine = FALSE) |>
  map2(names(genes),
       ~ {.x + ggtitle(.y) + xlab(NULL) + ylab(NULL) + NoLegend()}) |>
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






sub <- PrepSCTFindMarkers(sub)



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







