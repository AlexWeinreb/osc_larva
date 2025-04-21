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

osc_raw <- readxl::read_excel("../10x_grl18/data/oscillating/msb209498-sup-0003-datasetev1.xlsx",
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

library(tidyverse)


omega <- read_csv("ipynb/BWM/fit0_ω.csv")
log_gamma <- read_csv("ipynb/BWM/fit0_logγg.csv")
om <- ( omega/mean(colMeans(log_gamma)) )

hist(colMeans(om))
plot(res_vcy$peak_phase, as.numeric(omega[1,]))

plot(sort(as.numeric(omega[10,])))
plot(sort(res_vcy$peak_phase))

hist(as.matrix(read_csv("ipynb/BWM/fit0_logγg.csv")))
hist(as.matrix(read_csv("ipynb/AM_PHso/fit0_logγg.csv")))
hist(as.matrix(read_csv("ipynb/ILso/fit0_logγg.csv")))



# python
omega = full_pps_velo["ω"].squeeze().numpy() / torch.exp(torch.mean(full_pps_velo["logγg"].squeeze().mean(0).detach())).numpy()
phi = phase_pyro.phis

plot(res_vcy$peak_phase,
     colMeans(om)[1:664])

ids = np.array([n2n[i] for i in np.array(data_to_fit.obs["batch"])])
for i in range(len(data_to_fit.obs["batch"].unique())):
  omega1 = omega[:,np.where(ids == i)]
phi1 = phi[np.where(ids == i)]
omegas.append(omega1)
phis.append(phi1)

labels = np.array(data_to_fit.obs["batch"].unique())

colors = ["tab:blue"]
for i in range(len(omegas)):
  plt.plot(phis[i][np.argsort(phis[i])], omegas[i].mean(0)[0][np.argsort(phis[i])], c="black", linestyle='dashed')












