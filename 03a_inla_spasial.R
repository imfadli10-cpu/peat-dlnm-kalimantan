# =====================================================================
# 03a_inla_spasial.R
# Ide 1: DLNM titik panas Kalimantan - model INLA pertama
#
#   M2 = Negative Binomial + DLNM (cross-basis yang SAMA dengan script 02)
#        + offset log(luas) + efek spasial BYM2 antar-kabupaten
#
# Butuh: hasil script 01 (data_kalimantan/) dan script 02 (hasil_02/model_02.rds)
# Cara pakai: buka di RStudio, klik "Source". Waktu: ~1-5 menit.
# Hasil: folder hasil_03/ (ringkasan_03a.txt, m2_ringkas.rds, gambar)
# =====================================================================

# ---- 0. Paket & folder kerja -----------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(dlnm); library(INLA)
})
if (basename(getwd()) != "paper_karhutla" && dir.exists("~/paper_karhutla")) {
  setwd("~/paper_karhutla")
}
DIR <- "data_kalimantan"; OUT <- "hasil_03"
dir.create(OUT, showWarnings = FALSE)
for (f in c(file.path(DIR, "panel_kalimantan.rds"), file.path(DIR, "kalimantan.graph"),
            file.path("hasil_02", "model_02.rds"))) {
  if (!file.exists(f)) stop("File ", f, " tidak ada. Jalankan script 01 dan 02 dulu.")
}
message("Folder kerja: ", normalizePath(getwd()))

ringkas <- file.path(OUT, "ringkasan_03a.txt")
catat <- function(...) { cat(..., "\n"); cat(..., "\n", file = ringkas, append = TRUE) }
cat("", file = ringkas)

# ---- 1. Data & cross-basis (sama persis dengan script 02) --------------
panel <- readRDS(file.path(DIR, "panel_kalimantan.rds")) |> arrange(id_kab, bulan)
m02   <- readRDS(file.path("hasil_02", "model_02.rds"))
cb    <- m02$cb
if (nrow(cb) != nrow(panel)) stop("Jumlah baris cross-basis tidak cocok dengan panel. Jalankan ulang script 02.")

cb_df <- as.data.frame(unclass(cb)[, , drop = FALSE])
names(cb_df) <- paste0("cb", seq_len(ncol(cb_df)))

dat <- bind_cols(panel, cb_df) |>
  mutate(log_luas = log(luas_km2))
ok  <- complete.cases(cb_df) & !is.na(dat$hujan_mm)
dat <- dat[ok, ]

g <- inla.read.graph(file.path(DIR, "kalimantan.graph"))
catat("================ CEK DATA ================")
catat("Baris model      :", nrow(dat), "| kabupaten:", n_distinct(dat$id_kab))
n_cc <- if (!is.null(g$cc$n)) g$cc$n else NA
catat("Graph BYM2       :", g$n, "simpul | komponen terpisah:", n_cc)
if (g$n != n_distinct(panel$id_kab)) stop("Jumlah simpul graph != jumlah kabupaten.")
if (!all(sort(unique(dat$id_kab)) == seq_len(g$n))) stop("id_kab harus 1..n sesuai graph.")

# ---- 2. Model M2: NB + DLNM + BYM2 ------------------------------------
# Prior PC (penalised complexity), lihat Simpson dkk. 2017
pc_prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
pc_bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                phi  = list(prior = "pc",      param = c(0.5, 2/3)))

rumus_m2 <- as.formula(paste(
  "hotspot ~ 1 +", paste(names(cb_df), collapse = " + "),
  "+ offset(log_luas)",
  "+ f(id_kab, model = 'bym2', graph = g, scale.model = TRUE, constr = TRUE, hyper = pc_bym2)"))

# Hanya kolom yang dipakai model (menghindari masalah kolom tanggal/teks di INLA)
dat_m <- as.data.frame(dat[, c("hotspot", "id_kab", "log_luas", names(cb_df))])
dat_m$id_kab <- as.integer(dat_m$id_kab)

message("Menjalankan INLA untuk M2 (bisa beberapa menit) ...")
t0 <- Sys.time()
m2 <- tryCatch(
  inla(rumus_m2, family = "nbinomial", data = dat_m,
       control.fixed     = list(correlation.matrix = TRUE),
       control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
       control.predictor = list(compute = TRUE, link = 1),
       verbose = FALSE),
  error = function(e) { message("ERROR INLA: ", conditionMessage(e)); NULL })
if (is.null(m2)) stop("INLA gagal. Kirim pesan ERROR di atas ke Claude.")
durasi <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

# ---- 3. Ringkasan ------------------------------------------------------
catat("\n================ MODEL M2 (NB + DLNM + BYM2) ================")
catat("Waktu hitung     :", durasi, "menit")
catat(sprintf("WAIC = %.1f | DIC = %.1f", m2$waic$waic, m2$dic$dic))
catat(sprintf("Gagal CPO (harus 0 atau sangat sedikit): %d", sum(m2$cpo$failure > 0, na.rm = TRUE)))
catat("\nHiperparameter:")
hp <- m2$summary.hyperpar
for (i in seq_len(nrow(hp))) {
  catat(sprintf("  %-45s mean = %8.3f (%.3f - %.3f)", rownames(hp)[i],
                hp[i, "mean"], hp[i, "0.025quant"], hp[i, "0.975quant"]))
}

# Koefisien & kovarians cross-basis dari INLA -> crosspred (seperti Lowe dkk. 2021)
idx <- match(names(cb_df), m2$names.fixed)
if (anyNA(idx)) stop("Nama koefisien cross-basis tidak ditemukan di hasil INLA.")
coef_cb <- m2$summary.fixed$mean[idx]
V <- m2$misc$lincomb.derived.covariance.matrix
if (!is.null(V) && all(dim(V) == length(m2$names.fixed))) {
  V <- as.matrix(V)
  vcov_cb <- if (!is.null(rownames(V)) && all(names(cb_df) %in% rownames(V))) {
    V[names(cb_df), names(cb_df)]
  } else V[idx, idx]
  catat("Kovarians koefisien: dari INLA (correlation.matrix)")
} else {
  # Cadangan: hitung kovarians dari 2000 sampel posterior
  sel <- setNames(as.list(rep(1, length(idx))), names(cb_df))
  smp <- inla.posterior.sample(2000, m2, selection = sel)
  mat <- t(sapply(smp, function(s) s$latent[, 1]))
  vcov_cb <- cov(mat)
  catat("Kovarians koefisien: dari 2000 sampel posterior")
}
vcov_cb <- unname(vcov_cb)

pred_k <- crosspred(cb, coef = coef_cb, vcov = vcov_cb, model.link = "log",
                    at = m02$kering, cen = m02$cen)
pred_k_halus <- crosspred(cb, coef = coef_cb, vcov = vcov_cb, model.link = "log",
                          at = m02$kering, cen = m02$cen, bylag = 0.2)
pred_m1 <- crosspred(cb, coef = m02$coef_cb, vcov = m02$vcov_cb, model.link = "log",
                     at = m02$kering, cen = m02$cen)

catat(sprintf("\nBulan kering (%.0f mm) vs median (%.0f mm):", m02$kering, m02$cen))
catat(sprintf("  RR kumulatif lag 0-%d  M2 (INLA): %.2f (%.2f - %.2f)",
              m02$LAG, pred_k$allRRfit[1], pred_k$allRRlow[1], pred_k$allRRhigh[1]))
catat(sprintf("  RR kumulatif lag 0-%d  M1 (GLM) : %.2f (%.2f - %.2f)",
              m02$LAG, pred_m1$allRRfit[1], pred_m1$allRRlow[1], pred_m1$allRRhigh[1]))
catat("\n  RR per lag       M2 (INLA)              M1 (GLM, script 02)")
for (l in 0:m02$LAG) {
  catat(sprintf("  lag %d :   %.2f (%.2f - %.2f)      %.2f (%.2f - %.2f)", l,
                pred_k$matRRfit[1, l + 1],  pred_k$matRRlow[1, l + 1],  pred_k$matRRhigh[1, l + 1],
                pred_m1$matRRfit[1, l + 1], pred_m1$matRRlow[1, l + 1], pred_m1$matRRhigh[1, l + 1]))
}

png(file.path(OUT, "g6_lag_kering_M2.png"), width = 1600, height = 1100, res = 200)
plot(pred_k_halus, "slices", var = m02$kering, xlab = "Lag (bulan)", ylab = "RR",
     main = sprintf("M2 (INLA + BYM2): efek bulan kering (%.0f mm)", m02$kering),
     col = "firebrick", lwd = 2)
dev.off()

# Efek spasial (risiko relatif kabupaten setelah memperhitungkan hujan)
re <- m2$summary.random$id_kab[seq_len(g$n), ]
spasial <- data.frame(id_kab = seq_len(g$n), re_mean = re$mean,
                      rr_spasial = exp(re$mean))
nm <- panel |> distinct(id_kab, kabupaten)
spasial <- left_join(spasial, nm, by = "id_kab") |> arrange(desc(rr_spasial))
catat("\n5 kabupaten dengan efek spasial tertinggi (risiko lebih dari yang dijelaskan hujan):")
for (i in 1:5) catat(sprintf("  %-25s RR spasial = %.2f", spasial$kabupaten[i], spasial$rr_spasial[i]))

saveRDS(list(coef_cb = coef_cb, vcov_cb = vcov_cb, waic = m2$waic$waic, dic = m2$dic$dic,
             hyper = m2$summary.hyperpar, spasial = spasial, durasi = durasi),
        file.path(OUT, "m2_ringkas.rds"))

catat("\nSelesai. Hasil di folder:", normalizePath(OUT))
