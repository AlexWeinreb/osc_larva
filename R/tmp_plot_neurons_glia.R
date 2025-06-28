library(Seurat)

dir_step1 <- "E:/backups/Projects_june2025/glia/osc_larva/intermediates/2502/250609_step1/"




ilso_subseu <- qs::qread( file.path(dir_step1,
                                    paste0("ILso", "_seu.qs")) )
samples_table <- read_tsv("data/samples_table.tsv")


sample2prom <- samples_table |>
  select(sample_name, promoter) |>
  column_to_rownames("sample_name")

all( unique(ilso_subseu$orig.ident) %in% rownames(sample2prom) )

ilso_subseu$promoter <- sample2prom[ilso_subseu$orig.ident, "promoter"]

# all cells
FeaturePlot(ilso_subseu,
            features = c("grl-18", "nsIs198"),
            pt.size = 2,
            alpha = .8)

# only from grl-18 sorts
ilso_subseu$promoter |> table()

ilso_subseu2 <- ilso_subseu |> subset(promoter == "grl-18")


FeaturePlot(ilso_subseu2,
            features = c("grl-18", "nsIs198"),
            pt.size = 2,
            alpha = .8)

GetAssayData(ilso_subseu2)["grl-18",] |> hist(breaks = 50)
GetAssayData(ilso_subseu2,assay = "RNA", layer = "counts")["grl-18",] |> hist(breaks = 50)


plot(GetAssayData(ilso_subseu2)["grl-18",],
     GetAssayData(ilso_subseu2,assay = "RNA", layer = "counts")["grl-18",] |> log1p())



plot(GetAssayData(ilso_subseu2)["nsIs198",],
     GetAssayData(ilso_subseu2,assay = "RNA", layer = "counts")["nsIs198",] )

plot(GetAssayData(ilso_subseu,assay = "RNA", layer = "counts")["grl-18",] |> log10(),
     GetAssayData(ilso_subseu,assay = "RNA", layer = "counts")["nsIs198",] |> log10() )


plot(GetAssayData(ilso_subseu,assay = "RNA", layer = "counts")["grl-18",] ,
     GetAssayData(ilso_subseu,assay = "RNA", layer = "counts")["nsIs198",]  )





