# velocycle results


# Inits ----
library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(295) |>
  add_row(X = "a",
          gene_id = "nsIs198",
          symbol = "GFP",
          sequence = "GFP",
          status = "Live",
          biotype = "protein_coding_gene",
          name = "GFP"
  )

osc_raw <- readxl::read_excel("../grl18/data/oscillating/msb209498-sup-0003-datasetev1.xlsx",
                              sheet = "Dataset EV1 WBidToGeneNames_Osc",
                              na = "NA") |>
  mutate(gene_id = wb_clean_gene_names(WB_ID),
         gene_name = i2s(gene_id, gids) )



# Load ----
dir_velocycle <- "ipynb/"
list.files(dir_velocycle)


res_vcy <- read_csv(file.path(dir_velocycle,
                              "velocycle_BWM.csv")) |>
  rename(gene_id = `...1`) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE),
         .before = 1)


osc_table <- osc_raw |>
  filter(gene_name %in% res_vcy$gene_name) |>
  mutate(OscAmplitude = if_else(Class == "nonOsc",
                                0,
                                OscAmplitude)) |>
  select(gene_name, gene_id,
         bulk_amplitude = OscAmplitude,
         bulk_peak_phase_deg = PeakPhase) |>
  group_by(gene_name) |>
  slice_sample(n = 1) |>
  ungroup() |>
  full_join(res_vcy |> select(1),
            by = "gene_name")




# reproduce angle/radius
res_vcy |>
  mutate(ang = Arg(`nu1cos mean` + 1i * `nu1sin mean`) ) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = peak_phase, y = ang))


res_vcy |>
  mutate(xs = `nu1cos mean`,
         ys = `nu1sin mean`,
         r = log10( sqrt(xs^2 + ys^2) /(`nu1sin std` + `nu1cos std`) ),
         angle = ( atan2(ys, xs))  ) |>
  select(gene_name, peak_phase, amplitude, r, angle) |>
  ggplot() +
  theme_classic() +
  # geom_abline() +
  geom_point(aes(x = peak_phase, y = angle))

res_vcy |>
  mutate(xs = `nu1cos mean`,
         ys = `nu1sin mean`,
         r = log10( sqrt(xs^2 + ys^2) /(`nu1sin std` + `nu1cos std`) ),
         angle = ( atan2(ys, xs))  ) |>
  select(gene_name, peak_phase, amplitude, r, angle) |>
  ggplot() +
  theme_classic() +
  geom_abline() +
  geom_point(aes(x = amplitude, y = r))






# Compare bulk ----

res_vcy |>
  select(gene_name, peak_phase, amplitude) |>
  left_join(osc_table,
            by = "gene_name") |>
  ggplot() +
  theme_classic() +
  geom_abline() +
  geom_point(aes(x = bulk_peak_phase_deg, y = peak_phase*180/pi,
                 alpha = log10(bulk_amplitude)))

res_vcy |>
  select(gene_name, peak_phase, amplitude) |>
  left_join(osc_table,
            by = "gene_name") |>
  filter(bulk_amplitude != 0) |>
  ggplot() +
  theme_classic() +
  # geom_abline() +
  aes(x = log10(bulk_amplitude), y = amplitude) +
  geom_point() +
  geom_smooth(method = lm)


# with uncertainty
res_vcy |>
  mutate(uncertainty = `nu0 std` + `nu1sin std` + `nu1cos std`) |>
  select(gene_name, peak_phase, amplitude, uncertainty) |>
  arrange(desc(amplitude)) |>
  left_join(osc_table,
            by = "gene_name") |>
  filter(bulk_amplitude != 0) |>
  ggplot() +
  theme_classic() +
  # geom_abline() +
  aes(x = log10(bulk_amplitude), y = amplitude) +
  geom_point(aes(alpha = uncertainty)) +
  geom_smooth(method = lm)


res_vcy |>
  mutate(uncertainty = `nu0 std` + `nu1sin std` + `nu1cos std`) |>
  select(gene_name, peak_phase, amplitude, uncertainty) |>
  left_join(osc_table,
            by = "gene_name") |>
  ggplot() +
  theme_classic() +
  geom_abline() +
  geom_point(aes(x = bulk_peak_phase_deg, y = peak_phase*180/pi,
                 color = uncertainty))


res_vcy |>
  mutate(uncertainty = `nu0 std` + `nu1sin std` + `nu1cos std`) |>
  select(gene_name, peak_phase, amplitude, uncertainty) |>
  left_join(osc_table,
            by = "gene_name") |>
  ggplot() +
  theme_classic() +
  aes(x = bulk_amplitude, y = uncertainty) +
  geom_point() + geom_smooth(method = "lm")

res_vcy |>
  mutate(uncertainty = `nu0 std` + `nu1sin std` + `nu1cos std`) |>
  select(gene_name, peak_phase, amplitude, uncertainty) |>
  left_join(osc_table,
            by = "gene_name") |>
  ggplot() +
  theme_classic() +
  aes(x = bulk_amplitude > 1.5, y = uncertainty) +
  geom_boxplot()



# look at them


res_vcy |>
  ggplot() +
  theme_minimal() +
  scale_y_continuous(limits = c(0,NA)) +
  scale_x_continuous(limits = c(0, 2*pi)) +
  coord_polar(start = -pi/2, direction = -1) +
  geom_point(aes(x = peak_phase, y = amplitude)) +
  ggrepel::geom_text_repel(aes(x = peak_phase, y = amplitude,
                               label = gene_name),
                           data = res_vcy |> filter(amplitude > 1.3))



res_vcy |>
  mutate(uncertainty = `nu0 std` + `nu1sin std` + `nu1cos std`) |>
  select(gene_name, peak_phase, amplitude, uncertainty) |>
  ggplot() +
  theme_classic() +
  aes(x = uncertainty, y = amplitude) +
  geom_point()


# AM/PHso
pred_manual <- readxl::read_excel("intermediates/2502/250330_step2/manual_AMPHso.xlsx")

res_vcy |>
  mutate(uncertainty = `nu0 std` + `nu1sin std` + `nu1cos std`) |>
  select(gene_name, peak_phase, amplitude, uncertainty) |>
  left_join(pred_manual,
            by = "gene_name") |>
  ggplot() +
  theme_classic() +
  aes(x = uncertainty, y = amplitude) +
  geom_point(aes(color = manual),
             alpha = .2) +
  geom_point(aes(color = manual),
             alpha = .8,
             size = 2,
             data = res_vcy |>
               mutate(uncertainty = `nu0 std` + `nu1sin std` + `nu1cos std`) |>
               select(gene_name, peak_phase, amplitude, uncertainty) |>
               left_join(pred_manual,
                         by = "gene_name") |> filter(!is.na(manual))) +
  ggrepel::geom_text_repel(aes(label = gene_name),
                           data = res_vcy |>
                             mutate(uncertainty = `nu0 std` + `nu1sin std` + `nu1cos std`) |>
                             select(gene_name, peak_phase, amplitude, uncertainty) |>
                             left_join(pred_manual,
                                       by = "gene_name") |> filter(!is.na(manual)))






# Load velocity fit ----
dir_velocycle <- "intermediates/2502/250418_velocycle/"
ct <- "AM_PHso"


# check files

list.files(file.path(dir_velocycle, ct), "fit1")

xx <- read_csv(file.path(dir_velocycle, ct,
                              paste0("fit0_ζ.csv")))
dim(xx)
matrixStats::colVars(as.matrix(xx)) |> table()

colMeans(xx) |> plot(type = "l")
colMeans(xx) |> hist()
plot(res_vcy$amplitude, log10(colMeans(xx)))

plot(phi_cells$phi, colMeans(xx))


tibble(gene_name = i2s(res_vcy$`...1`, gids),
       rho = colMeans(xx)) |>
  left_join(pred_manual,
            by = "gene_name") |> filter(!is.na(manual)) |>
  ggplot() + geom_density(aes(x = rho, fill = manual), alpha = .2)

# fit0_shape_inv.csv
tibble(gene_name = i2s(res_vcy$`...1`, gids),
       amplitude = res_vcy$amplitude,
       noise = colMeans(xx)) |>
  left_join(pred_manual,
            by = "gene_name") |>
  ggplot() +
  scale_y_log10() +
  geom_point(aes(x = amplitude, y = noise),alpha = .5) +
  geom_point(aes(x = amplitude, y = noise, color = manual),
             data = tibble(gene_name = i2s(res_vcy$`...1`, gids),
                           amplitude = res_vcy$amplitude,
                           noise = colMeans(xx)) |>
               left_join(pred_manual,
                         by = "gene_name") |> filter(!is.na(manual)),
             size = 2)

# fit0_shape_inv.csv
res_vcy |>
  add_column(noise = colMeans(xx)) |>
  mutate(uncertainty = nu0_std + nu1sin_std + nu1cos_std) |>
  ggplot() +
  scale_y_log10() + scale_x_log10() +
  geom_point(aes(x = uncertainty, y = noise),alpha = .5) 


#fit1_ζ.csv
n <- ncol(xx)/3; stopifnot(n == round(n))
zeta1 <- xx[(0:(n-1))*3 + 1] |> colMeans()
zeta2 <- xx[(0:(n-1))*3 + 2] |> colMeans()
zeta3 <- xx[(0:(n-1))*3 + 3] |> colMeans()
table(zeta1)

plot(zeta2, zeta3)




####
res_vcy <- read_csv(file.path(dir_velocycle, ct,
                              paste0("velocycle_",ct,".csv")))


#manifold fit (1 harmonic velocity)

phi_cells <- read_csv(file.path(dir_velocycle, ct,
                                paste0(ct,"_cell_phi1.csv")))

log_gamma <- read_csv(file.path(dir_velocycle, ct,"fit1_logγg.csv")) |>
  colMeans()


omega_raw <- read_csv(file.path(dir_velocycle, ct,"fit1_ω.csv"))

mean_gamma <- exp(mean(log_gamma))




tibble(phi = phi_cells$phi,
       omega = colMeans(omega_raw) / mean_gamma) |>
  mutate(perc_5 = apply(omega_raw, 2, \(x) quantile(x, .05)) / mean_gamma,
         perc_95 = apply(omega_raw, 2, \(x) quantile(x, .95)) / mean_gamma ) |>
ggplot() +
  theme_classic() +
  ggtitle(ct) +
  geom_ribbon(aes(x = phi,
                  ymin = perc_5, ymax = perc_95),
              alpha = .5) +
  geom_line(aes(x = phi, y = omega),
            color = 'black',
            linetype = "dashed",
            linewidth = 2)









#manifold fit (constant velocity)

phi_cells <- read_csv(file.path(dir_velocycle, ct,
                                paste0(ct,"_cell_phi0.csv")))

log_gamma <- read_csv(file.path(dir_velocycle, ct,"fit0_logγg.csv")) |>
  colMeans()


omega_raw <- read_csv(file.path(dir_velocycle, ct,"fit0_ω.csv"))

mean_gamma <- exp(mean(log_gamma))




tibble(phi = phi_cells$phi,
       omega = colMeans(omega_raw) / mean_gamma) |>
  mutate(perc_5 = apply(omega_raw, 2, \(x) quantile(x, .05)) / mean_gamma,
         perc_95 = apply(omega_raw, 2, \(x) quantile(x, .95)) / mean_gamma ) |>
  ggplot() +
  theme_classic() +
  ggtitle(ct) +
  geom_ribbon(aes(x = phi,
                  ymin = perc_5, ymax = perc_95),
              alpha = .5) +
  geom_line(aes(x = phi, y = omega),
            color = 'black',
            linetype = "dashed",
            linewidth = 2)



#~ several cell types ----
velocities <- lapply(c("BWM", "ILso", "AM_PHso"),
                     \(ct){
                       # res_vcy <- read_csv(file.path(dir_velocycle, ct,
                       #                               paste0("velocycle_",ct,".csv")),
                       #                     col_types = "cdddddddd")
                       
                       
                       #manifold fit (1 harmonic velocity)
                       
                       phi_cells <- read_csv(file.path(dir_velocycle, ct,
                                                       paste0(ct,"_cell_phi1.csv")),
                                             col_types = "id")
                       
                       log_gamma <- read_csv(file.path(dir_velocycle, ct,"fit1_logγg.csv"),
                                             col_types = "d") |>
                         colMeans()
                       
                       
                       omega_raw <- read_csv(file.path(dir_velocycle, ct,"fit1_ω.csv"),
                                             col_types = "d")
                       
                       mean_gamma <- exp(mean(log_gamma))
                       
                       
                       
                       
                       tibble(cell_type = ct,
                              phi = phi_cells$phi,
                              omega = colMeans(omega_raw) / mean_gamma) |>
                         mutate(perc_5 = apply(omega_raw, 2, \(x) quantile(x, .05)) / mean_gamma,
                                perc_95 = apply(omega_raw, 2, \(x) quantile(x, .95)) / mean_gamma )
                     }) |>
  bind_rows()





velocities |>
  ggplot() +
  theme_classic() +
  geom_ribbon(aes(x = phi,
                  ymin = perc_5, ymax = perc_95,
                  fill = cell_type),
              alpha = .5) +
  geom_line(aes(x = phi, y = omega,
                color = cell_type),
            linetype = "dashed",
            linewidth = 2)









#~ compare with batch ----

ct <- "AM_PHso"

dir_velocycle1 <- "intermediates/2502/250418_velocycle/"
dir_velocycle2 <- "intermediates/2502/250418_velocycle_with_batch/"

res_vcy1 <- read_csv(file.path(dir_velocycle1, ct,
                              paste0("velocycle_",ct,".csv"))) |>
  rename(gene_id = `...1`) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE),
         .before = 1)

res_vcy2 <- read_csv(file.path(dir_velocycle2,
                               paste0("velocycle_",ct,".csv"))) |>
  rename(gene_id = `...1`) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE),
         .before = 1)




# AM/PHso
pred_manual <- readxl::read_excel("intermediates/2502/250330_step2/manual_AMPHso.xlsx")

gd_trth1 <- res_vcy1 |>
  mutate(uncertainty = nu0_std + nu1sin_std + nu1cos_std) |>
  select(gene_name, peak_phase, amplitude, uncertainty) |>
  left_join(pred_manual,
            by = "gene_name") |> filter(!is.na(manual))


res_vcy1 |>
  mutate(uncertainty = nu0_std + nu1sin_std + nu1cos_std) |>
  select(gene_name, peak_phase, amplitude, uncertainty) |>
  left_join(pred_manual,
            by = "gene_name") |>
  ggplot() +
  theme_classic() +
  scale_x_continuous(limits = c(0, .5)) +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  aes(x = uncertainty, y = amplitude) +
  geom_point(aes(color = manual),
             alpha = .2) +
  geom_point(aes(color = manual),
             alpha = .8,
             size = 2,
             data = gd_trth1) +
  ggrepel::geom_text_repel(aes(label = gene_name),
                           data = gd_trth1)



gd_trth2 <- res_vcy2 |>
  mutate(uncertainty = nu0_std + nu1sin_std + nu1cos_std) |>
  # select(gene_name, peak_phase, amplitude, uncertainty) |>
  left_join(pred_manual,
            by = "gene_name") |> filter(!is.na(manual))


res_vcy2 |>
  mutate(uncertainty = nu0_std + nu1sin_std + nu1cos_std) |>
  select(gene_name, peak_phase, amplitude, uncertainty) |>
  left_join(pred_manual,
            by = "gene_name") |>
  ggplot() +
  theme_classic() +
  scale_x_continuous(limits = c(0, .5)) +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  aes(x = uncertainty, y = amplitude) +
  geom_point(aes(color = manual),
             alpha = .2) +
  geom_point(aes(color = manual),
             alpha = .8,
             size = 2,
             data = gd_trth2) +
  ggrepel::geom_text_repel(aes(label = gene_name),
                           data = gd_trth2)





gd_trth2 <- res_vcy2 |>
  mutate(uncertainty = nu0_std + nu1sin_std + nu1cos_std,
         amplitude_abs = log10(sqrt(nu1sin_mean^2 + nu1cos_mean^2))) |>
  # select(gene_name, peak_phase, amplitude, uncertainty) |>
  left_join(pred_manual,
            by = "gene_name") |> filter(!is.na(manual))


res_vcy2 |>
  mutate(uncertainty = nu0_std + nu1sin_std + nu1cos_std,
         amplitude_abs = log10(sqrt(nu1sin_mean^2 + nu1cos_mean^2) )) |>
  # select(gene_name, peak_phase, amplitude, uncertainty) |>
  left_join(pred_manual,
            by = "gene_name") |>
  ggplot() +
  theme_classic() +
  # scale_x_continuous(limits = c(0, .5)) +
  # scale_y_continuous(limits = c(-1.5, 1.5)) +
  aes(x = nu0_mean, y = amplitude_abs) +
  geom_point(aes(color = manual),
             alpha = .2) +
  geom_point(aes(color = manual),
             alpha = .8,
             size = 2,
             data = gd_trth2) +
  ggrepel::geom_text_repel(aes(label = gene_name),
                           data = gd_trth2)


res_vcy2 |>
  mutate(uncertainty = nu0_std + nu1sin_std + nu1cos_std,
         amplitude_abs = log10(sqrt(nu1sin_mean^2 + nu1cos_mean^2) )) |>
  # select(gene_name, peak_phase, amplitude, uncertainty) |>
  left_join(pred_manual,
            by = "gene_name") |>
  ggplot() +
  theme_classic() +
  # scale_x_continuous(limits = c(0, .5)) +
  scale_y_log10() +
  aes(x = 1, y = amplitude_abs / nu0_mean) +
  ggbeeswarm::geom_quasirandom(aes(color = manual),
             alpha = .2) +
  ggbeeswarm::geom_quasirandom(aes(color = manual),
             alpha = .8,
             size = 2,
             data = gd_trth2) +
  ggrepel::geom_text_repel(aes(label = gene_name),
                           data = gd_trth2)

