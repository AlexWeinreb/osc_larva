# Is our sc-based method biased about peak width?

# Look at Meeuse 2020 peak width, is there a bias in that we only classify as pulsatile the genes with wide/narrow peaks?



# Inits ----

opar <- par(no.readonly = TRUE)

library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(295)

dir_out <- "presentations/figures/260622_peak_widths"


transform_log10p <- scales::new_transform(
  name      = "log10p",
  transform = \(x) log10(x + 1),
  inverse   = \(x) 10^x - 1,
  breaks    = \(range){
    log_range <- log10(range + 1)
    10^pretty(log_range, n = 6)
  },
  domain    = c(0, Inf)
)





# for plot annotations
stage_df <- tribble(
  ~stage, ~start, ~end,
  "L1",      0,   13,
  "L2",     13,   20,
  "L3",     20,   27,
  "L4",     27,   34,
  "adult",  34,   50
)

scale_fill_stages <- scale_fill_manual(
  values = c("L1" = "#e8dcc8",
             "L2" = "#c8d8e0",
             "L3" = "#d8e0c8",
             "L4" = "#e0d0c8",
             "adult" = "#d4c5a9")
)


#~ load ----

meeuse_raw <- read_tsv("data/GSE130811_expr_mRNA_CE10_coding.tab.gz",
                col_names = c("gene_id",
                              read_lines("data/GSE130811_expr_mRNA_CE10_coding.tab.gz", n_max = 1L) |> str_split_1("\\t")),
                skip = 1) |>
  select(-ends_with(".2")) |>
  pivot_longer(cols = -c(1,2),
               names_to = "time",
               names_pattern = c("^([0-9]{1,2})hr$")) |>
  mutate(time = time |> str_remove("hr$") |> as.numeric(),
         gene_name = i2s(gene_id, gids))

# note, not exact match, here we focus on the common set
list(in_raw = unique(meeuse_raw$gene_name),
     in_table = unique(wormOsc::table_osc_genes$gene_name)) |>
  eulerr::euler() |>
  plot(quantities = TRUE)

genes_in_meeuse <- intersect(
  unique(meeuse_raw$gene_name),
  unique(wormOsc::table_osc_genes$gene_name)
) |>
  setdiff(NA_character_)

meeuse_raw_filt <- meeuse_raw |>
  filter(gene_name %in% genes_in_meeuse)

tab_osc <- wormOsc::table_osc_genes |>
  filter(gene_name %in% genes_in_meeuse)


genes_osc <- tab_osc |> filter(osc_amplitude > 1.5) |> pull(gene_name)







# Normalization TPM/TMM ----

# TPM

normalize_to_tpm <- function(dat){
  
  norm_factors <- tibble(
    time = unique(dat$time),
    tmm_factor = dat |>
      select(gene_id, time, value) |>
      pivot_wider(names_from = time, values_from = value) |>
      column_to_rownames("gene_id") |>
      as.matrix() |>
      edgeR::DGEList() |>
      edgeR::calcNormFactors(method = "TMM") |>
      chuck("samples", "norm.factors")
  )
  
  dat |>
    left_join(norm_factors,
              by = "time") |>
    mutate(nf = sum(value / width), .by = time) |>
    mutate(
      TPM = 1e6 * value / (width * nf),
      TPM_TMM = TPM / tmm_factor
    ) |>
    select(gene_name, gene_id, time, TPM_TMM)
}

meeuse_norm <- normalize_to_tpm(meeuse_raw)




#~ Plot normalized expression ----

# dpy-6 should match coordinates of molt
meeuse_norm |>
  filter(gene_name == "dpy-6") |>
  ggplot(aes(x = time, y = log10(1+TPM_TMM), linetype = gene_name, shape = gene_name, color = gene_name)) +
  theme_bw() +
  scale_x_continuous(n.breaks = 10, minor_breaks = unique(meeuse_norm$time)) +
  scale_fill_stages +
  geom_rect(data = stage_df,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = stage),
            inherit.aes = FALSE, alpha = 0.6) +
  geom_point(size = 3, alpha = .5) +
  geom_line(linewidth = 1)


# random osc genes
genes_sel <- genes_osc |> sample(2)


meeuse_norm |>
  filter(gene_name %in% genes_sel) |>
  ggplot(aes(x = time, y = log10(1+TPM_TMM), linetype = gene_name, shape = gene_name, color = gene_name)) +
  theme_bw() +
  scale_x_continuous(n.breaks = 10, minor_breaks = unique(meeuse_norm$time)) +
  geom_rect(data = stage_df,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = stage),
            inherit.aes = FALSE, alpha = 0.15) +
  geom_point(size = 3, alpha = .5) +
  geom_line(linewidth = .8)



# # other sets of genes to plot
# 
# genes_sel <- c("lin-14", "lin-29", "lin-28", "lin-41", "dpy-6")
# genes_sel <- c("lin-14", "lin-29", "lin-28", "lin-46", "hbl-1")
# genes_sel <- c("him-3", "glp-1", "fog-1", "fog-3")
# 
# genes_sel <- "vit-" |> paste0(1:6)
# 
# 
# genes_sel <- sample(meeuse_norm$gene_name, 3)





# Peak width ----

## conversion phase <-> time
# see bottom of the script for conversion numbers
tab_peak_times <- tab_osc |>
  filter(gene_name %in% genes_osc) |>
  mutate(peak_time_L3_h = ((23.5 + (peak_phase_deg / 360) * 8 - 19.5) %% 8) + 19.5,
         peak_time_L2_h = ((16.5 + (peak_phase_deg / 360) * 8 - 12.5) %% 8) + 12.5)



fit_gaussian <- function(gene, mat_expr, peak_times, half_window = 4) {
  
  peak_time <- peak_times[[gene]] |>
    round()
  
  window <- seq(peak_time - half_window, peak_time + half_window)
  
  stopifnot(all(window %in% as.numeric(rownames(mat_expr))))
  
  
  window_data <- mat_expr[as.character(window), gene]
  
  # plot(window, window_data, type = "b"); abline(v = peak_time)
  
  fit <- tryCatch(
    nls(
      y ~ baseline + A * exp(-(t - mu)^2 / (2 * sigma^2)),
      data    = data.frame(y = window_data, t = window),
      start   = list(baseline = min(window_data),
                     A = max(window_data) - min(window_data),
                     mu = peak_time,
                     sigma = 1.5),
      lower = list(baseline = 0,
                   A = 0,
                   mu = peak_time - half_window,
                   sigma = 0.5),
      upper = list(baseline = Inf,
                   A = Inf,
                   mu = peak_time + half_window,
                   sigma = half_window),
      algorithm = "port",  # required for bounds
      control = nls.control(maxiter = 100, warnOnly = TRUE)
    ),
    error = \(e) NULL
  )
  
  
  # Return NA row on failed/non-converged fit
  if (is.null(fit) || !fit$convInfo$isConv) {
    return(tibble(gene_name = gene, baseline = NA_real_, amplitude = NA_real_,
                  mu = NA_real_, sigma = NA_real_, width_gaussian_h = NA_real_))
  }
  
  params <- coef(fit)
  
  tibble(
    gene_name      = gene,
    baseline       = params[["baseline"]],
    amplitude      = params[["A"]],
    mu             = params[["mu"]],
    sigma          = abs(params[["sigma"]]),
    width_gaussian_h = 2 * sigma * sqrt(2 * log(2))  # FWHM
  )
}



mat_expr <- meeuse_norm |>
  filter(gene_name %in% genes_osc) |>
  pivot_wider(id_cols = gene_name, names_from = time, values_from = TPM_TMM) |>
  column_to_rownames("gene_name") |>
  as.matrix() |>
  t() |>
  log1p()

width_fits_l2 <- map_dfr(genes_osc,
                         fit_gaussian,
                         mat_expr = mat_expr,
                         peak_times = tab_peak_times$peak_time_L2_h |> set_names(tab_peak_times$gene_name),
                         .progress = TRUE)



#~ Plot with width ----
genes_sel <- sample(genes_osc, 3)


meeuse_norm |>
  inner_join(width_fits_l2 |> select(gene_name, width_gaussian_h),
            by = "gene_name") |>
  filter(gene_name %in% genes_sel) |>
  ggplot(aes(x = time, y = log10(1+TPM_TMM), linetype = gene_name, shape = gene_name, color = gene_name)) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = margin(t = 5, r = 60, b = 5, l = 5)) +
  coord_cartesian(clip = "off",
                  xlim = c(0, 50)) +
  scale_x_continuous(n.breaks = 6,
                     minor_breaks = unique(meeuse_norm$time)) +
  scale_fill_stages +
  geom_rect(data = stage_df,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = stage),
            inherit.aes = FALSE, alpha = 0.25) +
  geom_rect(data = stage_df,
            aes(xmin = start, xmax = end, ymin = -0.1, ymax = 0, fill = stage),
            inherit.aes = FALSE, alpha = 1, show.legend = FALSE) +
  annotate("text",
           x = stage_df$start + (stage_df$end - stage_df$start) / 2,
           y = -.05,
           label = stage_df$stage) +
  geom_point(size = 3, alpha = .5) +
  geom_line(linewidth = .8) +
  ggrepel::geom_label_repel(data = meeuse_norm |>
               inner_join(width_fits_l2 |> select(gene_name, width_gaussian_h), by = "gene_name") |>
               filter(gene_name %in% genes_sel) |>
               slice_max(time,
                         by = c(gene_name, width_gaussian_h)) |>  # place label at the end of the line
               mutate(label = sprintf("%s: %.1fh", gene_name, width_gaussian_h)),
             aes(x = time, y = log10(1 + TPM_TMM), label = label, color = gene_name),
             hjust = 0,
             nudge_x = 10,
             xlim = c(0, 60),
             direction = "y") +
  ggtitle("Expression during larval development (Meeuse 2020 data)") +
  xlab("Developmental time (h)") +
  ylab(expression(Expression:~log[10](1 + TPM[TMM])))





# Comparison sc clustering ----

dir_sc_clust <- "intermediates/2502/250624_cluster"

sc_clust <- read_csv(file.path(dir_sc_clust, "250624_cluster_results.csv"),
                  show_col_types = FALSE) |>
  mutate(cell_type = if_else(cell_type == "coelomyocyte", "coelomocyte", cell_type))


# we have computed widths for all osc genes, and only for osc genes
stopifnot(identical(
  genes_osc,
  width_fits_l2 |>
    pull(gene_name)
))


list(detected_in_sc = unique(sc_clust$gene_name),
     width_from_meeuse = genes_osc) |>
  eulerr::euler() |>
  plot(quantities = TRUE)




#~ Any cell type mixed ----

# check we are not overcounting too much by taking "any" cell type pulsatile
sc_clust |>
  summarize(is_pulsatile = any(shape == "pulsatile"),
            .by = gene_name) |>
  inner_join(tab_osc,
             by = "gene_name") |>
  mutate(class = case_when(
    is.na(osc_amplitude) ~ "non-puls",
    osc_amplitude <= 1.5 ~ "low",
    osc_amplitude > 1.5 ~ "puls",
    .default = "BUG"
  )) |>
  summarize(n_puls = sum(is_pulsatile),
            prop_puls = mean(is_pulsatile),
            n_tot = n(),
            .by = class)



# sc_vs_width <- clust |>
#   summarize(is_pulsatile = sum(shape == "pulsatile") >= 3L,
#             .by = gene_name) |>
#   inner_join(width_fits_l2 |> select(gene_name, width_gaussian_h),
#              by = "gene_name")


sc_clust |>
  summarize(is_pulsatile = any(shape == "pulsatile"),
            .by = gene_name) |>
  inner_join(width_fits_l2 |> select(gene_name, width_gaussian_h),
             by = "gene_name") |>
  ggplot() +
  theme_classic() +
  geom_histogram(
    aes(x = width_gaussian_h, fill = is_pulsatile),
    color = "white"
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#BC7858", "FALSE" = "#C0ADD7"),
    labels = c("TRUE" = "Pulsatile", "FALSE" = "Non-pulsatile"),
    name = "Gene pulsatile in any cell type"
  ) +
  xlab(expression("Peak width bulk RNA-Seq (Meeuse " * italic("et al.") * ")")) +
  ylab("Number of genes") +
  theme(legend.position = "inside",
        legend.position.inside = c(.8, .6))



# sc_vs_width |>
#   ggplot() +
#   theme_classic() +
#   geom_density(
#     aes(x = width_gaussian_h, fill = is_pulsatile, y = after_stat(count)),
#     alpha = .5
#   ) +
#   scale_fill_manual(
#     values = c("TRUE" = "#BC7858", "FALSE" = "#C0ADD7"),
#     labels = c("TRUE" = "Pulsatile", "FALSE" = "Non-pulsatile"),
#     name = "Gene pulsatile in any cell type"
#   ) +
#   xlab(expression("Peak width bulk RNA-Seq (Meeuse " * italic("et al.") * ")")) +
#   ylab("Number of genes") +
#   theme(legend.position = "inside",
#         legend.position.inside = c(.8, .6))





#~ ILso only ----

list(
  in_sc_ILso = sc_clust |>
    filter(cell_type == "ILso") |>
    pull(gene_name)
  ,
  width_from_meeuse = genes_osc
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)

# still largely overlaps
list(
  puls_ILso = sc_clust |>
    filter(cell_type == "ILso",
           shape == "pulsatile") |>
    pull(gene_name)
  ,
  width_from_meeuse = genes_osc
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)




sc_vs_width_ILso <- sc_clust |>
  filter(cell_type == "ILso") |>
  mutate(shape = factor(shape, levels = c("nonpulsatile", "low", "pulsatile"))) |>
  left_join(width_fits_l2 |> select(gene_name, width_gaussian_h),
             by = "gene_name")


table(is.na(sc_vs_width_ILso$width_gaussian_h),
      sc_vs_width_ILso$shape)


# > +++ Fig. EV4A top +++ ----

sc_vs_width_ILso |>
  filter(!is.na(width_gaussian_h)) |>
  ggplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "mm"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.spacing.y = unit(0, "mm"),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0, vjust = -5),
    strip.clip = "off"
  ) +
  scale_fill_manual(
    values = c("nonpulsatile" = "#C0ADD7", "low" = "#D4B483", "pulsatile" = "#BC7858"),
    labels = c("nonpulsatile" = "Non-pulsatile", "low" = "Low-amplitude", "pulsatile" = "Pulsatile"),
    name = "Oscillating gene expressed in seam"
  ) +
  xlab(expression("Peak width from bulk RNA-Seq (Meeuse " * italic("et al.") * ", hours)")) +
  ylab("Number of genes") +
  ggtitle("ILso") +
  geom_histogram(
    aes(x = width_gaussian_h, fill = shape),
    color = "white",
    bins = 30
  )

# ggsave("peak_width_ILso.png", path = dir_out,
#        width = 10, height = 7, units = "cm",
#        scale = 1.5)
# 
# ggsave("peak_width_ILso.pdf", path = dir_out,
#        width = 6, height = 8, units = "cm",
#        scale = 1.5)



#~ Seam only ----


list(
  in_sc_seam = sc_clust |>
    filter(cell_type == "seam") |>
    pull(gene_name)
  ,
  width_from_meeuse = genes_osc
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)

sc_vs_width_seam <- sc_clust |>
  filter(cell_type == "seam") |>
  mutate(shape = factor(shape, levels = c("nonpulsatile", "low", "pulsatile"))) |>
  left_join(width_fits_l2 |> select(gene_name, width_gaussian_h),
            by = "gene_name")

table(is.na(sc_vs_width_seam$width_gaussian_h),
      sc_vs_width_seam$shape)


# > +++ Fig. EV4A bottom +++ ----

sc_vs_width_seam |>
  ggplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "mm"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.spacing.y = unit(0, "mm"),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0, vjust = -5),
    strip.clip = "off"
  ) +
  scale_fill_manual(
    values = c("nonpulsatile" = "#C0ADD7", "low" = "#D4B483", "pulsatile" = "#BC7858"),
    labels = c("nonpulsatile" = "Non-pulsatile", "low" = "Low-amplitude", "pulsatile" = "Pulsatile"),
    name = "Oscillating gene expressed in seam"
  ) +
  xlab(expression("Peak width from bulk RNA-Seq (Meeuse " * italic("et al.") * ", hours)")) +
  ylab("Number of genes") +
  ggtitle("Seam cells") +
  geom_histogram(
    aes(x = width_gaussian_h, fill = shape),
    color = "white",
    bins = 30
  )

# ggsave("peak_width_seam.png", path = dir_out,
#        width = 10, height = 7, units = "cm",
#        scale = 1.5)
# 
# ggsave("peak_width_seam.pdf", path = dir_out,
#        width = 6, height = 8, units = "cm",
#        scale = 1.5)


















##  _______________________   ----

# __Peak width optimization__ ----

# Code used to test and optimize the peak width computation


mat_expr <- meeuse_norm |>
  filter(! is.na(gene_name),
         gene_name %in% genes_osc) |>
  pivot_wider(id_cols = gene_name,
              names_from = time,
              values_from = TPM_TMM) |>
  column_to_rownames("gene_name") |>
  as.matrix() |>
  t() |>
  (\(mat) log10(1 + mat) )()





## conversion phase <-> time

# #> L2 (13:20), visually, peak_time_h = ((16.5 + (peak_phase_deg / 360) * 8 - 12.5) %% 8) + 12.5)

apply(mat_expr[as.character(13:20), genes_osc], 2, \(x) which.max(x)) |>
  (\(x) (13:20)[x])() |>
  set_names(genes_osc) |>
  enframe("gene_name", "peak_time_h") |>
  left_join(wormOsc::table_osc_genes |> select(gene_name, peak_phase_deg),
            by = "gene_name") |>
  ggplot() +
  geom_point(aes(peak_time_h, peak_phase_deg)) +
  geom_line(aes(time_h, phase_deg),
             data = tibble(
               phase_deg = 0:360,
               time_h = ((16.5 + (phase_deg / 360) * 8 - 12.5) %% 8) + 12.5
             ),
             color = "green4",
            linewidth = 1)


# #> L3 (20:27), visually, peak_time_h = ((23.5 + (peak_phase_deg / 360) * 8 - 20) %% 8) + 20

apply(mat_expr[as.character(20:27), genes_osc], 2, \(x) which.max(x)) |>
  (\(x) (20:27)[x])() |>
  set_names(genes_osc) |>
  enframe("gene_name", "peak_time_h") |>
  left_join(wormOsc::table_osc_genes |> select(gene_name, peak_phase_deg),
            by = "gene_name") |>
  ggplot() +
  geom_point(aes(peak_time_h, peak_phase_deg)) +
  geom_line(aes(time_h, phase_deg),
            data = tibble(
              phase_deg = 0:360,
              time_h = ((23.5 + (phase_deg / 360) * 8 - 20) %% 8) + 20
            ),
            color = "green4",
            linewidth = 1)







## check that from phase we predict correct peak time

genes_sel <- genes_osc |> sample(1)

meeuse_norm |>
  filter(gene_name %in% genes_sel) |>
  ggplot(aes(x = time, y = log10(1+TPM_TMM), linetype = gene_name, shape = gene_name, color = gene_name)) +
  theme_classic() +
  scale_fill_stages +
  geom_rect(data = stage_df,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = stage),
            inherit.aes = FALSE, alpha = 0.6) +
  geom_point(size = 3, alpha = .5) +
  geom_line() +
  geom_vline(
    data = tab_peak_times[genes_sel, ],
    mapping = aes(xintercept = peak_time_L3_h),
    color = "green3"
    ) +
  geom_vline(
    data = tab_peak_times[genes_sel, ],
    mapping = aes(xintercept = peak_time_L2_h),
    color = "blue"
  )



# For each gene:

# We fit a gaussian to either L2 or L3 peak, and recover the width of the peak in hours.
# This width is directly reflecting sigma of the gaussian (FWHM = 2sqrt(2ln2)*sigma) https://en.wikipedia.org/wiki/Gaussian_function#Properties
# We can also look at mu (the center, i.e. time of max expression) for QC
#
# Another approach is to look what fraction of the stage is above a threshold (rel_width), after QC it doesn't work as well
threshold_log10 <- 2

fit_gaussian <- function(gene, mat_expr, tab_peak_times, col = "peak_time_L2_h", half_window = 4) {
  
  peak_time <- tab_peak_times[gene, col] |> round()
  window <- (peak_time - half_window):(peak_time + half_window)
  stopifnot(all( window %in% as.numeric(rownames(mat_expr)) ))
  
  window_data <- mat_expr[as.character(window), gene]
  
  # Fit 4-parameter Gaussian
  fit <- nls(
    y ~ baseline + A * exp(-(t - mu)^2 / (2*sigma^2)),
    data = data.frame(y = window_data, t = window),
    start = list(baseline = min(window_data),
                 A = max(window_data) - min(window_data),
                 mu = peak_time,
                 sigma = 1.5),
    control = nls.control(maxiter = 100, warnOnly = TRUE)
  )
  
  if (is.null(fit) || !fit$convInfo$isConv) return(
    tibble(
      gene_name = gene,
      baseline = NA_real_,
      amplitude = NA_real_,
      mu = NA_real_,
      sigma = NA_real_,
      width_relative = NA_real_,
      width_gaussian_h = NA_real_
    )
  )
  
  params <- coef(fit)
  
  # plot(window, window_data, type = "p")
  # lines(window, predict(fit))
  # abline(v = peak_time, col = 'red3')
  # abline(v = params[["mu"]], col = 'blue3')
  
  
  # Get width
  threshold <- params[["baseline"]] + params[["A"]] / 2
  
  rel_width <- sum(window_data > threshold_log10) / length(window_data)
  
  
  width_gaussian_h <- 2 * abs(params[["sigma"]]) * sqrt(2 * log(2))
  
  # export results
  tibble(
    gene_name = gene,
    baseline = params[["baseline"]],
    amplitude = params[["A"]],
    mu = params[["mu"]],
    sigma = abs(params[["sigma"]]),
    width_relative = rel_width,
    width_gaussian_h = width_gaussian_h
  )
}



width_fits_l2 <- map_dfr(osc_genes,
                         fit_gaussian,
                         col = "peak_time_L2_h",
                         mat_expr = mat_expr,
                         tab_peak_times = tab_peak_times,
                         .progress = TRUE)

width_fits_l3 <- map_dfr(osc_genes,
                         fit_gaussian,
                         col = "peak_time_L3_h",
                         mat_expr = mat_expr,
                         tab_peak_times = tab_peak_times,
                         .progress = TRUE)

# some failed to fit
table(is.na(width_fits_l2$baseline))
table(is.na(width_fits_l3$baseline))

# we recover timing of peak
width_fits_l2 |>
  left_join(tab_peak_times |> rownames_to_column("gene_name"),
            by = "gene_name") |>
  ggplot() +
  geom_point(aes(x = peak_time_L2_h, y = mu)) +
  geom_abline(color = "red3")

width_fits_l3 |>
  left_join(tab_peak_times |> rownames_to_column("gene_name"),
            by = "gene_name") |>
  ggplot() +
  geom_point(aes(x = peak_time_L3_h, y = mu)) +
  geom_abline(color = "red3")




# the relative width above threshold is very sensitive to baseline, the gaussian width is more robust
width_fits_l2 |> ggplot() +
  geom_jitter(aes(x = baseline, y = width_relative))

pairs(width_fits_l2 |> select(baseline, amplitude, mu, sigma, width_relative, width_gaussian_h),
      pch = ".", cex = 2)
cor(width_fits_l2 |> select(baseline, amplitude, mu, sigma, width_relative, width_gaussian_h),
    use = "complete.obs") |>
  round(2)




# mu-sigma correlation?
width_fits_l2 |>
  left_join(tab_peak_times |> rownames_to_column("gene_name"),
            by = "gene_name") |>
  ggplot() +
  geom_point(aes(x = peak_phase_deg, y = sigma))

width_fits_l3 |>
  left_join(tab_peak_times |> rownames_to_column("gene_name"),
            by = "gene_name") |>
  ggplot() +
  geom_point(aes(x = peak_phase_deg, y = sigma)) +
  scale_y_continuous(limits = c(NA, 5))


# Merge L2 and L3 fits
compare <- inner_join(
  width_fits_l2 |> select(gene_name, sigma_l2 = sigma, mu_l2 = mu),
  width_fits_l3 |> select(gene_name, sigma_l3 = sigma, mu_l3 = mu),
  by = "gene_name"
) |> filter(!(sigma_l3 > 5))

plot(compare$sigma_l2, compare$sigma_l3,
     xlab = "sigma (L2 fit)", ylab = "sigma (L3 fit)")
abline(0, 1, col = "red")
cor(compare$sigma_l2, compare$sigma_l3, use = "complete.obs")

# L2 is a bit more robust

hist(width_fits_l2$width_gaussian_h, breaks = 30)






# Examples and extreme cases

genes_sel <- width_fits_l2 |> filter(width_gaussian_h > 4) |> pull(gene_name) |> sample(3)
genes_sel <- width_fits_l2 |> filter(width_gaussian_h < 1.5) |> pull(gene_name) |> sample(3)

i <- i+3
genes_sel <- genes_thin_cand[i:(i+2)]

dat_norm |>
  filter(gene_name %in% genes_sel) |>
  ggplot(aes(x = time, y = log10(1+TPM_TMM), linetype = gene_name, shape = gene_name, color = gene_name)) +
  theme_classic() +
  geom_rect(data = stage_df,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = stage),
            inherit.aes = FALSE, alpha = 0.15) +
  geom_point(size = 3, alpha = .5) +
  geom_line()

width_fits_l2 |>
  filter(gene_name %in% genes_sel)

# bcc-1 plateaus, but not in sc data (not expressed?? weird one paper reports it in hyp but sc data doesn't)

# ham-2 and Y11D7A.3 are nice examples of broad peaks that can be found in the sc data

# mam-2, C26B9.7 are thin, F45E4.5 is thin and very late

# genes thin (here) and expressed in my sc data
genes_thin_cand <- c('F45E4.5','T14A8.2','mam-2','R148.5','R74.2','fmo-3','Y43F4A.1','C02F12.5','ztf-16','C26B9.7','C45B2.8','T13C5.7','gpx-5','nhr-91','F58H10.1','ptr-6','bli-2','C30G12.4','col-38','E01G4.6','R05A10.6','Y43D4A.5','R07B1.5','K02H11.4','grl-27','nhr-181','C35A11.2','C54F6.3','F28G4.2','col-79')


genes_sel <- c('C26B9.7', 'mam-2', 'F45E4.5','ham-2', 'Y11D7A.3')
genes_sel <- c("ham-2", "mam-2")
genes_sel <- c('C26B9.7','Y11D7A.3')



dat_norm |>
  filter(gene_name %in% genes_sel) |>
  ggplot(aes(x = time, y = log10(1+TPM_TMM), linetype = gene_name, shape = gene_name, color = gene_name)) +
  theme_classic() +
  scale_fill_manual(values = c("L1" = "#e8dcc8",
                               "L2" = "#c8d8e0",
                               "L3" = "#d8e0c8",
                               "L4" = "#e0d0c8",
                               "adult" = "#d4c5a9")) +
  geom_rect(data = stage_df2,
            aes(xmin = start, xmax = end, ymin = 2.15, ymax = 2.25, fill = stage),
            inherit.aes = FALSE, alpha = 1, show.legend = FALSE) +
  annotate("text",
           x = stage_df2$start + (stage_df2$end - stage_df2$start) / 2,
           y = 2.2,
           label = stage_df2$stage) +
  geom_point(size = 3, alpha = .5) +
  geom_line() +
  ggtitle("Expression during larval development (Meeuse 2020 data)") +
  xlab("Developmental time (h)") +
  ylab(expression(Expression:~log[10](1 + TPM[TMM])))

