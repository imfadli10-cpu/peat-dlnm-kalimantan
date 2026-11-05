# =====================================================================
# 10b_gambar_inggris_cepat.R
# Membuat gambar paper berbahasa Inggris TANPA menjalankan INLA lagi.
# Semua gambar dibuat dari hasil yang sudah tersimpan (script 01-09).
# Waktu: ~1 menit.
#
#   fig1a_timeseries, fig1b_peat_map        (data)
#   fig2a_lag_peat, fig2b_cumulative_peat   (model M4, dari hasil_05)
#   fig3a_seasonal, fig3b_spatial_risk      (M3a/M3b, dari hasil_03)
#   fig5a_sim_ratio, fig5b_sim_coverage     (simulasi, dari hasil_09)
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(sf)
})
sf_use_s2(FALSE)
if (basename(getwd()) != "paper_karhutla" && dir.exists("~/paper_karhutla")) {
  setwd("~/paper_karhutla")
}
DIR <- "data_kalimantan"; FIG <- "figures_en"
dir.create(FIG, showWarnings = FALSE)
butuh <- c(file.path(DIR, "panel_kalimantan.rds"), file.path(DIR, "kab_kalimantan.gpkg"),
           file.path("hasil_02", "model_02.rds"), file.path("hasil_03", "m3_ringkas.rds"),
           file.path("hasil_05", "m4_ringkas.rds"), file.path("hasil_09", "simulasi_per_ulangan.csv"))
for (f in butuh) if (!file.exists(f)) stop("File ", f, " tidak ada.")
message("Folder kerja: ", normalizePath(getwd()))

tema <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"))
tema_peta <- theme_minimal(base_size = 11)

GAGAL <- character()
simpan <- function(p, nama, w, h) {
  ok <- tryCatch({
    ggsave(file.path(FIG, paste0(nama, ".png")), p, width = w, height = h, dpi = 300)
    ggsave(file.path(FIG, paste0(nama, ".pdf")), p, width = w, height = h)
    TRUE
  }, error = function(e) { message("  ! GAGAL menyimpan ", nama, ": ", conditionMessage(e)); FALSE })
  if (ok) message("  tersimpan: ", nama) else GAGAL <<- c(GAGAL, nama)
}

panel <- readRDS(file.path(DIR, "panel_kalimantan.rds"))
kab   <- st_read(file.path(DIR, "kab_kalimantan.gpkg"), quiet = TRUE)
kab   <- kab[, setdiff(names(kab), c("prop_gambut", "luas_gambut_km2"))] |>
  left_join(distinct(panel, id_kab, prop_gambut), by = "id_kab")
m02 <- readRDS(file.path("hasil_02", "model_02.rds"))
m3  <- readRDS(file.path("hasil_03", "m3_ringkas.rds"))
m4  <- readRDS(file.path("hasil_05", "m4_ringkas.rds"))

# ---- Figure 1 -----------------------------------------------------------
tot <- panel |> group_by(bulan) |> summarise(n = sum(hotspot), .groups = "drop")
simpan(ggplot(tot, aes(bulan, n / 1000)) + geom_col(fill = "firebrick", width = 25) +
         labs(x = NULL, y = "Active-fire detections (thousands)") + tema,
       "fig1a_timeseries", 7, 3.2)
simpan(ggplot(kab) + geom_sf(aes(fill = 100 * prop_gambut), colour = "grey40", linewidth = 0.1) +
         scale_fill_distiller(palette = "YlOrBr", direction = 1, name = "Peatland (%)") + tema_peta,
       "fig1b_peat_map", 5.5, 5)

# ---- Figure 2 (M4, hasil script 05) -------------------------------------
lab0 <- "No peat (0%)"; labt <- sprintf("Peat-rich (%.0f%%)", 100 * m4$P_TINGGI)
warna <- setNames(c("steelblue", "firebrick"), c(lab0, labt))
ke_en <- function(x) ifelse(grepl("^Tanpa", x), lab0, labt)

lag_df <- m4$lag_df |> mutate(group = factor(ke_en(kelompok), levels = c(lab0, labt)))
simpan(ggplot(lag_df, aes(lag, rr, colour = group, fill = group)) +
         geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
         geom_line(linewidth = 0.9) + geom_hline(yintercept = 1, linetype = 2) +
         scale_colour_manual(values = warna) + scale_fill_manual(values = warna) +
         scale_x_continuous(breaks = 0:6) +
         labs(x = "Lag (months)", y = sprintf("RR (%.0f mm vs %.0f mm)", m02$kering, m02$cen),
              colour = NULL, fill = NULL) + tema,
       "fig2a_lag_peat", 5.5, 4)

kum_df <- m4$kum_df |> mutate(group = factor(ke_en(kelompok), levels = c(lab0, labt)))
simpan(ggplot(kum_df, aes(hujan, rr, colour = group, fill = group)) +
         geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
         geom_line(linewidth = 0.9) + geom_hline(yintercept = 1, linetype = 2) +
         geom_vline(xintercept = m02$kering, linetype = 3) + scale_y_log10() +
         scale_colour_manual(values = warna) + scale_fill_manual(values = warna) +
         labs(x = "Monthly rainfall (mm)", y = "Cumulative RR, lag 0-6 (log scale)",
              colour = NULL, fill = NULL) + tema,
       "fig2b_cumulative_peat", 5.5, 4)

# ---- Figure 3 -----------------------------------------------------------
simpan(ggplot(m3$musim, aes(bl, rr, colour = model, fill = model)) +
         geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
         geom_line(linewidth = 0.9) + geom_hline(yintercept = 1, linetype = 2) +
         scale_x_continuous(breaks = 1:12, labels = month.abb) + scale_y_log10() +
         scale_colour_manual(values = c(M3a = "grey40", M3b = "darkorange")) +
         scale_fill_manual(values = c(M3a = "grey40", M3b = "darkorange")) +
         labs(x = NULL, y = "Seasonal RR (log scale)", colour = "Model", fill = "Model") + tema,
       "fig3a_seasonal", 6, 3.8)
peta <- kab |> left_join(m3$spasial[, c("id_kab", "rr_spasial")], by = "id_kab")
simpan(ggplot(peta) + geom_sf(aes(fill = rr_spasial), colour = "grey40", linewidth = 0.1) +
         scale_fill_distiller(palette = "RdYlBu", trans = "log10", name = "Spatial RR") + tema_peta,
       "fig3b_spatial_risk", 5.5, 5)

# ---- Figure 5 (simulasi) ------------------------------------------------
sim <- read.csv(file.path("hasil_09", "simulasi_per_ulangan.csv"))
sim <- sim[!sim$gagal & !is.na(sim$lr_est), ]
if (nrow(sim)) {
  sim$scen <- ifelse(sim$skenario == "A", "A: peat modification", "B: no peat modification")
  benar <- sim |> group_by(scen) |> summarise(true = exp(median(lr_est - lr_bias)), .groups = "drop")
  simpan(ggplot(sim, aes(model, exp(lr_est), fill = model)) +
           geom_boxplot(alpha = 0.6, outlier.size = 0.8, width = 0.6) +
           geom_hline(data = benar, aes(yintercept = true), linetype = 2, colour = "firebrick") +
           facet_wrap(~ scen, scales = "free_y") + scale_y_log10() +
           scale_fill_manual(values = c(M4 = "grey60", M5b = "steelblue")) +
           labs(x = NULL, y = "Estimated ratio R (log scale)") + tema + theme(legend.position = "none"),
         "fig5a_sim_ratio", 6, 3.8)
  cov <- sim |> group_by(scen, model) |> summarise(c = mean(lr_cov), n = n(), .groups = "drop") |>
    mutate(lo = pmax(0, c - 1.96 * sqrt(c * (1 - c) / n)), hi = pmin(1, c + 1.96 * sqrt(c * (1 - c) / n)))
  simpan(ggplot(cov, aes(model, c, colour = model)) +
           geom_hline(yintercept = 0.95, linetype = 2) + geom_pointrange(aes(ymin = lo, ymax = hi)) +
           facet_wrap(~ scen) + coord_cartesian(ylim = c(0, 1)) +
           scale_colour_manual(values = c(M4 = "grey35", M5b = "steelblue")) +
           labs(x = NULL, y = "Coverage of 95% CrI for R") + tema + theme(legend.position = "none"),
         "fig5b_sim_coverage", 6, 3.5)
}

if (length(GAGAL)) message("\nPERHATIAN: gambar yang gagal: ", paste(GAGAL, collapse = ", "))
message("\nSelesai. Gambar di folder: ", normalizePath(FIG))
