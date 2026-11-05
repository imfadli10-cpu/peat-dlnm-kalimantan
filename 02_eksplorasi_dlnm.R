# =====================================================================
# 02_eksplorasi_dlnm.R
# Ide 1: DLNM untuk titik panas Kalimantan - eksplorasi & model pertama
#
# Isi:
#   1. Eksplorasi: pola nol, musiman, hubungan hujan-titik panas per lag
#   2. Cross-basis DLNM: curah hujan bulanan, lag 0-6 bulan
#   3. Model pembanding (tanpa INLA):
#        M0 = Negative Binomial tanpa lag (hujan bulan yang sama saja)
#        M1 = Negative Binomial + DLNM (efek non-linear & tertunda)
#      keduanya dengan efek tetap kabupaten
#   4. Grafik exposure-lag-response + ringkasan angka
#
# Butuh: hasil script 01 (folder data_kalimantan/)
# Cara pakai: buka di RStudio, klik "Source". Waktu: ~1-3 menit.
# =====================================================================

# ---- 0. Paket & folder kerja -----------------------------------------
pkgs <- c("dplyr", "tidyr", "ggplot2", "dlnm", "MASS", "splines", "sf")
baru <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(baru)) install.packages(baru, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(dlnm); library(splines); library(sf)
})
sf_use_s2(FALSE)

if (basename(getwd()) != "paper_karhutla" && dir.exists("~/paper_karhutla")) {
  setwd("~/paper_karhutla")
}
DIR  <- "data_kalimantan"
OUT  <- "hasil_02"
dir.create(OUT, showWarnings = FALSE)
f_panel <- file.path(DIR, "panel_kalimantan.rds")
if (!file.exists(f_panel)) stop("File ", f_panel, " tidak ada. Jalankan script 01 dulu.")
message("Folder kerja: ", normalizePath(getwd()))

panel <- readRDS(f_panel) |> arrange(id_kab, bulan)   # urutan WAJIB: kabupaten lalu waktu
kab   <- st_read(file.path(DIR, "kab_kalimantan.gpkg"), quiet = TRUE)

LAG <- 6                                               # lag maksimum (bulan)

# ---- 1. Eksplorasi -----------------------------------------------------
sink_file <- file.path(OUT, "ringkasan_02.txt")
catat <- function(...) { cat(..., "\n"); cat(..., "\n", file = sink_file, append = TRUE) }
cat("", file = sink_file)

catat("================ EKSPLORASI ================")
catat("Observasi:", nrow(panel), "| Kabupaten:", n_distinct(panel$id_kab),
      "| Bulan:", n_distinct(panel$bulan))
catat("Proporsi nol keseluruhan:", round(mean(panel$hotspot == 0), 3))

nol_bulan <- panel |> group_by(bl) |>
  summarise(prop_nol = mean(hotspot == 0), rata2 = mean(hotspot), .groups = "drop")
catat("\nProporsi nol & rata-rata titik panas per bulan kalender:")
for (i in seq_len(nrow(nol_bulan))) {
  catat(sprintf("  %s : nol = %.2f | rata-rata = %8.1f",
                month.abb[nol_bulan$bl[i]], nol_bulan$prop_nol[i], nol_bulan$rata2[i]))
}

nol_kab <- panel |> group_by(id_kab, kabupaten) |>
  summarise(prop_nol = mean(hotspot == 0), total = sum(hotspot), .groups = "drop") |>
  arrange(desc(prop_nol))
catat("\n10 kabupaten/kota dengan proporsi nol tertinggi:")
for (i in 1:min(10, nrow(nol_kab))) {
  catat(sprintf("  %-28s nol = %.2f | total = %d",
                nol_kab$kabupaten[i], nol_kab$prop_nol[i], nol_kab$total[i]))
}

# Korelasi Spearman: hujan pada lag k vs titik panas bulan ini
kor <- sapply(0:LAG, function(k) {
  d <- panel |> group_by(id_kab) |>
    mutate(h_lag = dplyr::lag(hujan_mm, k)) |> ungroup()
  suppressWarnings(cor(d$h_lag, d$hotspot, method = "spearman", use = "complete.obs"))
})
catat("\nKorelasi Spearman hujan (lag k) vs titik panas:")
for (k in 0:LAG) catat(sprintf("  lag %d bulan : %+.3f", k, kor[k + 1]))

p_musim <- ggplot(nol_bulan, aes(factor(month.abb[bl], levels = month.abb), rata2)) +
  geom_col(fill = "firebrick") +
  labs(x = NULL, y = "Rata-rata titik panas per kabupaten",
       title = "Pola musiman titik panas (2012-2025)") +
  theme_minimal()
ggsave(file.path(OUT, "g1_musiman.png"), p_musim, width = 8, height = 4, dpi = 200)

p_sebar <- panel |> filter(hotspot > 0) |>
  ggplot(aes(hujan_mm, hotspot)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "firebrick") +
  scale_y_log10() +
  labs(x = "Curah hujan bulan yang sama (mm/bulan)", y = "Titik panas (skala log)",
       title = "Hujan vs titik panas (lag 0)") +
  theme_minimal()
ggsave(file.path(OUT, "g2_sebaran_lag0.png"), p_sebar, width = 7, height = 5, dpi = 200)

# ---- 2. Cross-basis DLNM -----------------------------------------------
# Dimensi paparan : natural spline, knot di persentil 10/50/90 curah hujan
# Dimensi lag     : natural spline, knot log-spaced (lebih rapat di lag awal)
# group = id_kab  : lag dihitung di dalam tiap kabupaten (tidak "bocor" antarkabupaten)
q_hujan <- quantile(panel$hujan_mm, c(.10, .50, .90), na.rm = TRUE)
cb <- crossbasis(panel$hujan_mm, lag = LAG,
                 argvar = list(fun = "ns", knots = q_hujan),
                 arglag = list(fun = "ns", knots = logknots(LAG, 2)),
                 group  = panel$id_kab)
catat("\nCross-basis: ", ncol(cb), "kolom (paparan df =", attr(cb, "df")[1],
      ", lag df =", attr(cb, "df")[2], ")")

# ---- 3. Model pembanding (Negative Binomial) --------------------------
dat <- panel |> mutate(kab_f = factor(id_kab))
ok  <- complete.cases(cb) & !is.na(dat$hujan_mm)
catat("Baris dipakai model:", sum(ok), "dari", nrow(dat),
      "(", LAG, "bulan pertama tiap kabupaten terpakai untuk lag )")

kontrol <- glm.control(maxit = 100)
m0 <- MASS::glm.nb(hotspot ~ ns(hujan_mm, knots = q_hujan) + kab_f,
                   data = dat[ok, ], control = kontrol)
m1 <- MASS::glm.nb(hotspot ~ cb[ok, ] + kab_f, data = dat[ok, ], control = kontrol)

# cb[ok, ] menghilangkan atribut crossbasis -> simpan koefisien & vcov secara eksplisit
idx_cb <- grep("^cb", names(coef(m1)))
coef_cb <- coef(m1)[idx_cb]
vcov_cb <- vcov(m1)[idx_cb, idx_cb]

catat("\n================ PERBANDINGAN MODEL ================")
catat(sprintf("M0 (NB, tanpa lag)   : AIC = %.1f | theta = %.3f", AIC(m0), m0$theta))
catat(sprintf("M1 (NB + DLNM lag %d) : AIC = %.1f | theta = %.3f", LAG, AIC(m1), m1$theta))
catat(sprintf("Selisih AIC (M0 - M1): %.1f  (positif = DLNM lebih baik)", AIC(m0) - AIC(m1)))

# ---- 4. Prediksi DLNM & grafik ----------------------------------------
cen <- unname(q_hujan[2])                        # pembanding: curah hujan median
at  <- seq(floor(min(panel$hujan_mm, na.rm = TRUE)),
           ceiling(quantile(panel$hujan_mm, .99, na.rm = TRUE)), length.out = 80)
pred <- crosspred(cb, coef = coef_cb, vcov = vcov_cb, model.link = "log",
                  at = at, cen = cen, bylag = 0.2)

kering  <- unname(quantile(panel$hujan_mm, .10, na.rm = TRUE))
pred_k  <- crosspred(cb, coef = coef_cb, vcov = vcov_cb, model.link = "log",
                     at = kering, cen = cen, bylag = 0.2)
pred_k1 <- crosspred(cb, coef = coef_cb, vcov = vcov_cb, model.link = "log",
                     at = kering, cen = cen)            # lag bulat untuk tabel

catat("\n================ HASIL DLNM (M1) ================")
catat(sprintf("Pembanding (cen) = hujan median = %.0f mm/bulan", cen))
catat(sprintf("Kondisi kering   = persentil 10 = %.0f mm/bulan", kering))
catat(sprintf("RR kumulatif lag 0-%d (kering vs median): %.2f (95%% CI %.2f - %.2f)",
              LAG, pred_k$allRRfit[1], pred_k$allRRlow[1], pred_k$allRRhigh[1]))
catat("RR per lag (kering vs median):")
for (l in 0:LAG) {
  catat(sprintf("  lag %d : %.2f (%.2f - %.2f)", l,
                pred_k1$matRRfit[1, l + 1], pred_k1$matRRlow[1, l + 1], pred_k1$matRRhigh[1, l + 1]))
}

png(file.path(OUT, "g3_permukaan_3d.png"), width = 1600, height = 1300, res = 200)
plot(pred, xlab = "Curah hujan (mm/bulan)", ylab = "Lag (bulan)", zlab = "RR",
     theta = 210, phi = 30, lphi = 30, main = "Exposure-lag-response: hujan -> titik panas")
dev.off()

png(file.path(OUT, "g4_kumulatif.png"), width = 1600, height = 1100, res = 200)
plot(pred, "overall", xlab = "Curah hujan (mm/bulan)", ylab = "RR kumulatif (lag 0-6)",
     main = "Efek kumulatif curah hujan terhadap titik panas", log = "y", col = "firebrick", lwd = 2)
abline(v = kering, lty = 2)
dev.off()

png(file.path(OUT, "g5_lag_kering.png"), width = 1600, height = 1100, res = 200)
plot(pred_k, "slices", var = kering, xlab = "Lag (bulan)", ylab = "RR",
     main = sprintf("Efek bulan kering (%.0f mm) menurut lag", kering), col = "firebrick", lwd = 2)
dev.off()

saveRDS(list(cb = cb, q_hujan = q_hujan, cen = cen, kering = kering, LAG = LAG,
             m0 = m0, m1 = m1, coef_cb = coef_cb, vcov_cb = vcov_cb),
        file.path(OUT, "model_02.rds"))

catat("\nSelesai. Semua hasil (teks & gambar) ada di folder:", normalizePath(OUT))
