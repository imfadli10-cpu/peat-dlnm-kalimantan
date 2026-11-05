# =====================================================================
# 06_diagnostik_validasi.R
# Ide 1: diagnostik model & validasi prakiraan (out-of-sample)
#
# A. Diagnostik M4 (model utama):
#    - histogram PIT (kalibrasi prediksi)
#    - autokorelasi residual Pearson di dalam kabupaten (lag 1-6 bulan)
#    - total titik panas bulanan: observasi vs model
# B. Validasi "hindcast" tahun kebakaran besar: 2015 (El Nino), 2019 (IOD+), 2023 (El Nino)
#    Data titik panas tahun target DISEMBUNYIKAN, model dilatih tanpa tahun itu,
#    lalu memprediksi tahun itu (curah hujan dianggap diketahui).
#    Model dibandingkan: M2, M3b, M4. Ukuran: log score, MAE log, korelasi Spearman,
#    rasio total prediksi/observasi.
#
# Butuh: script 01-05 sudah dijalankan
# Cara pakai: buka di RStudio, klik "Source". Waktu: ~5-15 menit (10 model INLA).
# Hasil: folder hasil_06/
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(dlnm); library(INLA); library(ggplot2)
})
if (basename(getwd()) != "paper_karhutla" && dir.exists("~/paper_karhutla")) {
  setwd("~/paper_karhutla")
}
DIR <- "data_kalimantan"; OUT <- "hasil_06"
dir.create(OUT, showWarnings = FALSE)
for (f in c(file.path(DIR, "panel_kalimantan.rds"), file.path(DIR, "kalimantan.graph"),
            file.path("hasil_02", "model_02.rds"))) {
  if (!file.exists(f)) stop("File ", f, " tidak ada. Jalankan script 01-05 dulu.")
}
message("Folder kerja: ", normalizePath(getwd()))

ringkas <- file.path(OUT, "ringkasan_06.txt")
catat <- function(...) { cat(..., "\n"); cat(..., "\n", file = ringkas, append = TRUE) }
cat("", file = ringkas)

# ---- 1. Data (sama dengan script 05) -----------------------------------
panel <- readRDS(file.path(DIR, "panel_kalimantan.rds")) |> arrange(id_kab, bulan)
if (!"prop_gambut" %in% names(panel)) stop("Kolom prop_gambut belum ada. Jalankan script 04 dulu.")
m02 <- readRDS(file.path("hasil_02", "model_02.rds"))
cb  <- m02$cb
if (nrow(cb) != nrow(panel)) stop("Jumlah baris cross-basis tidak cocok dengan panel. Jalankan ulang script 02.")
K <- ncol(cb)
cb_mat <- unclass(cb)[, seq_len(K), drop = FALSE]
cb_df  <- as.data.frame(cb_mat);                     names(cb_df)  <- paste0("cb",  seq_len(K))
cbg_df <- as.data.frame(cb_mat * panel$prop_gambut); names(cbg_df) <- paste0("cbg", seq_len(K))

ok <- complete.cases(cb_df) & !is.na(panel$hujan_mm)
info <- panel[ok, c("id_kab", "bulan", "tahun", "bl")]
dat0 <- data.frame(hotspot  = panel$hotspot,
                   id_kab   = as.integer(panel$id_kab),
                   bl       = as.integer(panel$bl),
                   tahun_id = as.integer(panel$tahun - min(panel$tahun) + 1L),
                   log_luas = log(panel$luas_km2))
dat0 <- cbind(dat0, cb_df, cbg_df)[ok, ]
g <- inla.read.graph(file.path(DIR, "kalimantan.graph"))

pc_prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
pc_bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                phi  = list(prior = "pc",      param = c(0.5, 2/3)))
f_bym2  <- "+ f(id_kab, model = 'bym2', graph = g, scale.model = TRUE, constr = TRUE, hyper = pc_bym2)"
f_musim <- "+ f(bl, model = 'rw2', cyclic = TRUE, scale.model = TRUE, constr = TRUE, hyper = pc_prec)"
f_tahun <- "+ f(tahun_id, model = 'iid', constr = TRUE, hyper = pc_prec)"
x_cb    <- paste(names(cb_df),  collapse = " + ")
x_cbg   <- paste(names(cbg_df), collapse = " + ")
rumus <- list(
  M2  = as.formula(paste("hotspot ~ 1 +", x_cb, "+ offset(log_luas)", f_bym2)),
  M3b = as.formula(paste("hotspot ~ 1 +", x_cb, "+ offset(log_luas)", f_bym2, f_musim, f_tahun)),
  M4  = as.formula(paste("hotspot ~ 1 +", x_cb, "+", x_cbg, "+ offset(log_luas)", f_bym2, f_musim, f_tahun)))

fit <- function(nama, data) {
  m <- tryCatch(
    inla(rumus[[nama]], family = "nbinomial", data = data,
         control.compute   = list(cpo = TRUE, config = FALSE),
         control.predictor = list(compute = TRUE, link = 1),
         verbose = FALSE),
    error = function(e) { message("ERROR INLA (", nama, "): ", conditionMessage(e)); NULL })
  if (is.null(m)) stop("INLA gagal untuk ", nama, ". Kirim pesan ERROR di atas ke Claude.")
  m
}
ukuran_nb <- function(m) {
  hp <- m$summary.hyperpar
  hp[grep("size", rownames(hp)), "mean"][1]
}

# ---- A. Diagnostik M4 (semua data) -------------------------------------
message("A. Diagnostik: menjalankan M4 dengan semua data ...")
m4 <- fit("M4", dat0)
mu   <- m4$summary.fitted.values$mean[seq_len(nrow(dat0))]
size <- ukuran_nb(m4)
res  <- (dat0$hotspot - mu) / sqrt(mu + mu^2 / size)          # residual Pearson NB

d_res <- info |> mutate(res = res) |> arrange(id_kab, bulan)
acf_kab <- sapply(1:6, function(k) {
  d <- d_res |> group_by(id_kab) |> mutate(r_lag = dplyr::lag(res, k)) |> ungroup()
  suppressWarnings(cor(d$res, d$r_lag, use = "complete.obs"))
})
catat("================ A. DIAGNOSTIK M4 ================")
catat(sprintf("Ukuran NB (size) = %.3f", size))
catat("Autokorelasi residual Pearson di dalam kabupaten:")
for (k in 1:6) catat(sprintf("  lag %d bulan : %+.3f", k, acf_kab[k]))
catat("  (|r| < 0.1 = baik; > 0.3 = masih ada autokorelasi kuat yang perlu ditangani)")

pit <- m4$cpo$pit
catat(sprintf("PIT: rata-rata %.3f (ideal 0.5) | proporsi < 0.05: %.3f | > 0.95: %.3f (ideal 0.05)",
              mean(pit, na.rm = TRUE), mean(pit < 0.05, na.rm = TRUE), mean(pit > 0.95, na.rm = TRUE)))

png(file.path(OUT, "g15_pit_m4.png"), width = 1400, height = 1000, res = 200)
hist(pit, breaks = 20, col = "grey80", border = "white", main = "Histogram PIT - M4 (ideal: rata)",
     xlab = "PIT", freq = FALSE)
abline(h = 1, lty = 2, col = "firebrick")
dev.off()

tot <- info |> mutate(obs = dat0$hotspot, fit = mu) |>
  group_by(bulan) |> summarise(obs = sum(obs), fit = sum(fit), .groups = "drop")
catat(sprintf("Korelasi total bulanan observasi vs model (log): %.3f",
              cor(log1p(tot$obs), log1p(tot$fit))))
g16 <- ggplot(tot, aes(bulan)) +
  geom_col(aes(y = obs), fill = "grey70") +
  geom_line(aes(y = fit), colour = "firebrick", linewidth = 0.7) +
  labs(x = NULL, y = "Titik panas per bulan", title = "Observasi (batang) vs model M4 (garis)") +
  theme_minimal()
ggsave(file.path(OUT, "g16_obs_vs_fit_m4.png"), g16, width = 9, height = 4, dpi = 200)
rm(m4); invisible(gc())

# ---- B. Validasi hindcast ----------------------------------------------
TARGET <- intersect(c(2015, 2019, 2023), unique(info$tahun))
hasil <- list(); pred_bulanan <- list()
for (th in TARGET) {
  uji <- info$tahun == th
  dat_l <- dat0; dat_l$hotspot[uji] <- NA
  for (nm in names(rumus)) {
    message("B. Hindcast ", th, " - ", nm, " ...")
    m  <- fit(nm, dat_l)
    mu <- m$summary.fitted.values$mean[seq_len(nrow(dat_l))][uji]
    sz <- ukuran_nb(m)
    y  <- dat0$hotspot[uji]
    hasil[[paste(th, nm)]] <- data.frame(
      tahun = th, model = nm,
      log_score = -mean(dnbinom(y, size = sz, mu = mu, log = TRUE)),
      mae_log   = mean(abs(log1p(y) - log1p(mu))),
      spearman  = suppressWarnings(cor(y, mu, method = "spearman")),
      rasio_total = sum(mu) / sum(y))
    pred_bulanan[[paste(th, nm)]] <- info[uji, ] |>
      mutate(obs = y, pred = mu) |>
      group_by(bulan) |> summarise(obs = sum(obs), pred = sum(pred), .groups = "drop") |>
      mutate(tahun = th, model = nm)
    rm(m); invisible(gc())
  }
}
hasil <- bind_rows(hasil); pred_bulanan <- bind_rows(pred_bulanan)
write.csv(hasil, file.path(OUT, "validasi_hindcast.csv"), row.names = FALSE)

catat("\n================ B. VALIDASI HINDCAST ================")
catat("log score & MAE: makin KECIL makin baik | Spearman: makin BESAR | rasio total: ideal 1")
catat(sprintf("%-6s %-5s %10s %9s %9s %12s", "Tahun", "Model", "LogScore", "MAE_log", "Spearman", "Pred/Obs"))
for (i in seq_len(nrow(hasil))) {
  catat(sprintf("%-6d %-5s %10.3f %9.3f %9.3f %12.2f", hasil$tahun[i], hasil$model[i],
                hasil$log_score[i], hasil$mae_log[i], hasil$spearman[i], hasil$rasio_total[i]))
}
rata <- hasil |> group_by(model) |>
  summarise(log_score = mean(log_score), mae_log = mean(mae_log),
            spearman = mean(spearman), .groups = "drop")
catat("\nRata-rata ketiga tahun:")
for (i in seq_len(nrow(rata))) {
  catat(sprintf("  %-4s log score = %.3f | MAE log = %.3f | Spearman = %.3f",
                rata$model[i], rata$log_score[i], rata$mae_log[i], rata$spearman[i]))
}
catat("Model terbaik (log score terkecil):", rata$model[which.min(rata$log_score)])

g17 <- ggplot(pred_bulanan, aes(bulan)) +
  geom_col(data = distinct(pred_bulanan, tahun, bulan, obs), aes(y = obs), fill = "grey75") +
  geom_line(aes(y = pred, colour = model), linewidth = 0.9) +
  facet_wrap(~ tahun, scales = "free", ncol = 1) +
  scale_colour_manual(values = c(M2 = "grey40", M3b = "steelblue", M4 = "firebrick")) +
  labs(x = NULL, y = "Titik panas per bulan", colour = "Model",
       title = "Hindcast: observasi (batang) vs prediksi tanpa data tahun tsb (garis)") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(OUT, "g17_hindcast.png"), g17, width = 8, height = 9, dpi = 200)

catat("\nSelesai. Hasil di folder:", normalizePath(OUT))
