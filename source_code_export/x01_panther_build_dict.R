
# Using PANTHER.db (bioc package), build a dictionary for enrichment analysis


# First time: force download of annotation hub to cache
# library(AnnotationHub)
# ann_hub <- AnnotationHub()
# ann_hub_panther <- query(ann_hub, "PANTHER.db")[[1]]


# Inits ----
library(PANTHER.db)
library(tidyverse)





PANTHER.db
availablePthOrganisms(PANTHER.db)[7,]

pthOrganisms(PANTHER.db) <- "WORM"
columns(PANTHER.db)
keytypes(PANTHER.db)
keys(PANTHER.db, keytype = "FAMILY_ID") |> head()




family_to_uniprot <- mapIds(PANTHER.db,
                            keys=keys(PANTHER.db, keytype = "FAMILY_ID"),
                            column="UNIPROT",
                            keytype="FAMILY_ID",
                            multiVals="list") |>
  stack() |>
  rename(uniprot_id = values,
         panther_family = ind) |>
  as_tibble()


family_with_terms <- PANTHER.db::select(PANTHER.db,
                   keys = family_to_uniprot$panther_family,
                   column = c("CLASS_TERM","SUBFAMILY_TERM"),
                   keytype = "FAMILY_ID",
                   join = "left") |>
  mutate(CLASS_TERM = replace_na(CLASS_TERM, ""),
         SUBFAMILY_TERM = replace_na(SUBFAMILY_TERM, ""),
         family_terms = paste(CLASS_TERM, SUBFAMILY_TERM) |> str_trim()) |>
  summarize(family_terms = paste(family_terms, collapse = "; "),
            .by = FAMILY_ID) |>
  mutate(family_terms = paste(family_terms, FAMILY_ID))

family_to_uniprot <- left_join(family_to_uniprot,
          family_with_terms,
          by = c(panther_family = "FAMILY_ID"))




# Downloaded from Wormbase simplemine
id_table <- read_tsv("data/simplemine_uniprot_ids.txt",
                     show_col_types = FALSE) |>
  select(gene_id = `WormBase Gene ID`,
         uniprot_id = UniProt) |>
  mutate(uniprot_id = str_split(uniprot_id, fixed(", ")) ) |>
  unnest(uniprot_id) |>
  filter(uniprot_id != "N.A.")


# there are a few missing
table(family_to_uniprot$uniprot_id %in% id_table$uniprot_id)
# ignore them explicitly:
family_to_uniprot <- family_to_uniprot |>
  filter(uniprot_id %in% id_table$uniprot_id)



# note the same UniProt ID can correspond to several genes,
# e.g. 16 his-xx genes encode the same HIS-3/P09588 protein
table(duplicated(id_table$uniprot_id))
# in this case, ensure we have all rows in the dictionary (i.e. each his-xx gene has an identical row)


# conversely, many prots in Panther are in more than one family/subfamily
# hence many-to-many relationship


dict <- left_join(family_to_uniprot,
                id_table,
                by = "uniprot_id",
                relationship = "many-to-many") |>
  select(-uniprot_id) |>
  distinct() |>
  add_column(tmp = 1L) |>
  pivot_wider(id_cols = gene_id,
              names_from = "family_terms",
              values_from = tmp,
              values_fill = 0L) |>
  rename(wbid = gene_id) |>
  rename_with(~paste(.x,"WBbt:0000000"), -wbid)

# qs::qsave(dict, "intermediates/2502/250330_step3_genes_by_celltype/panther_dict.qs")


# wormbaseEnrich::enrichment_analysis(sample(id_table$gene_id, 5), dict)






























