# =====================================================================
# 09_simulasi.R
# Ide 1: STUDI SIMULASI
#
# Data dibangkitkan dari model "sebenarnya" yang parameternya diambil dari
# hasil analisis data asli (geometri 55 kabupaten, curah hujan CHIRPS asli,
# proporsi gambut asli, efek DLNM & gambut dari M4, hiperparameter dari M5b).
#
#   Skenario A (ADA modifikasi gambut + AR1 per kabupaten):
#     dipasang M3b (abaikan gambut & AR1), M4 (abaikan AR1), M5b (benar)
#     -> bias, RMSE, cakupan (coverage) selang 95%, kuasa uji (power)
#   Skenario B (TIDAK ADA modifikasi gambut, delta = 0, tetap ada AR1):
#     dipasang M4 dan M5b
#     -> galat tipe I: seberapa sering model "menemukan" efek gambut palsu
#
# Butuh: script 01-07 sudah dijalankan
# Cara pakai: buka di RStudio, klik "Source".
# WAKTU: LAMA. Sekitar 12-15 menit per ulangan (5 model INLA), jadi
#        30 ulangan ~ 6-8 jam. Bisa dihentikan kapan saja (tombol STOP):
#        ulangan yang sudah selesai TERSIMPAN, jalankan ulang untuk melanjutkan.
# Hasil: folder hasil_09/
# =====================================================================

N_REP <- 30          # jumlah ulangan per skenario (boleh dinaikkan)
SEED  <- 20261104

if (!requireNamespace("sn", quietly = TRUE)) {      # dipakai INLA saat mengulang model yang sulit
  install.packages("sn", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(dlnm); library(INLA); library(ggplot2)
})
if (basename(getwd()) != "paper_karhutla" && dir.exists("~/paper_karhutla")) {
  setwd("~/paper_karhutla")
}
DIR <- "data_kalimantan"; OUT <- "hasil_09"; REPDIR <- file.path(OUT, "ulangan")
dir.create(REPDIR, recursive = TRUE, showWarnings = FALSE)
butuh <- c(file.path(DIR, "panel_kalimantan.rds"), file.path(DIR, "kalimantan.graph"),
           file.path("hasil_02", "model_02.rds"), file.path("hasil_03", "m3_ringkas.rds"),
           file.path("hasil_05", "m4_ringkas.rds"), file.path("hasil_07", "m5_ringkas.rds"))
for (f in butuh) if (!file.exists(f)) stop("File ", f, " tidak ada. Jalankan script 01-07 dulu.")
message("Folder kerja: ", normalizePath(getwd()))

# ---- 1. Bahan dari analisis data asli ----------------------------------
panel <- readRDS(file.path(DIR, "panel_kalimantan.rds")) |> arrange(id_kab, bulan)
m02 <- readRDS(file.path("hasil_02", "model_02.rds"))
m3  <- readRDS(file.path("hasil_03", "m3_ringkas.rds"))
m4  <- readRDS(file.path("hasil_05", "m4_ringkas.rds"))
m5  <- readRDS(file.path("hasil_07", "m5_ringkas.rds"))
g   <- inla.read.graph(file.path(DIR, "kalimantan.graph"))
cb  <- m02$cb; K <- ncol(cb)
if (nrow(cb) != nrow(panel)) stop("Jumlah baris cross-basis tidak cocok dengan panel. Jalankan ulang script 02.")
if (length(m4$theta) != K) stop("Koefisien M4 tidak cocok dengan cross-basis. Jalankan ulang script 05.")
P_TINGGI <- m4$P_TINGGI

ambil_hp <- function(hp, pola, cadangan) {
  r <- grep(pola, rownames(hp))
  if (length(r)) hp[r[1], "mean"] else { warning("Hiperparameter '", pola, "' tidak ditemukan, pakai ", cadangan); cadangan }
}
hp5 <- m5$M5b$hyper
SIZE     <- ambil_hp(hp5, "^size",                  0.8)
PREC_AR1 <- ambil_hp(hp5, "^Precision for t$",      1.0)
RHO_AR1  <- ambil_hp(hp5, "^Rho for t$",            0.47)
PREC_TH  <- ambil_hp(hp5, "^Precision for tahun_id$", 9.0)

# Efek spasial & musiman "sebenarnya" = estimasi dari M3b (data asli)
sp_true <- log(m3$spasial$rr_spasial[order(m3$spasial$id_kab)])
ms <- m3$musim[m3$musim$model == "M3b", ]
ms_true <- log(ms$rr[order(ms$bl)])
if (length(sp_true) != g$n || length(ms_true) != 12) stop("Efek spasial/musiman dari script 03b tidak lengkap.")

theta0 <- m4$theta; delta0 <- m4$delta
cbm <- unclass(cb)[, seq_len(K), drop = FALSE]
ok  <- complete.cases(cbm) & !is.na(panel$hujan_mm)
cb_df  <- as.data.frame(cbm);                     names(cb_df)  <- paste0("cb",  seq_len(K))
cbg_df <- as.data.frame(cbm * panel$prop_gambut); names(cbg_df) <- paste0("cbg", seq_len(K))
X  <- cbind(cb_df, cbg_df)[ok, ]
id_kab <- as.integer(panel$id_kab[ok]); bl <- as.integer(panel$bl[ok])
tahun_id <- as.integer(panel$tahun[ok] - min(panel$tahun[ok]) + 1L)
tt <- as.integer(factor(panel$bulan[ok])); tt <- tt - min(tt) + 1L
log_luas <- log(panel$luas_km2[ok])
n_t <- max(tt); n_th <- max(tahun_id)

# Nilai sebenarnya (target estimasi)
benar <- function(delta) {
  p0 <- crosspred(cb, coef = theta0, vcov = diag(1e-12, K), model.link = "log",
                  at = m02$kering, cen = m02$cen)
  pt <- crosspred(cb, coef = theta0 + P_TINGGI * delta, vcov = diag(1e-12, K), model.link = "log",
                  at = m02$kering, cen = m02$cen)
  list(cum0 = log(p0$allRRfit[1]), cumt = log(pt$allRRfit[1]),
       logratio = log(pt$allRRfit[1]) - log(p0$allRRfit[1]),
       lag0 = log(p0$matRRfit[1, ]), lagt = log(pt$matRRfit[1, ]))
}
TRUTH <- list(A = benar(delta0), B = benar(rep(0, K)))

# Intersep dikalibrasi supaya rata-rata titik panas simulasi ~ data asli
eta_tetap <- function(delta) {
  as.vector(cbm[ok, ] %*% theta0 + (cbm[ok, ] * panel$prop_gambut[ok]) %*% delta) +
    log_luas + sp_true[id_kab] + ms_true[bl]
}
# koreksi lognormal: efek acak tahun & AR1 menaikkan rata-rata sebesar exp(var/2)
koreksi  <- (1 / PREC_TH + 1 / PREC_AR1) / 2
INTERSEP <- list(A = log(sum(panel$hotspot[ok]) / sum(exp(eta_tetap(delta0)))) - koreksi,
                 B = log(sum(panel$hotspot[ok]) / sum(exp(eta_tetap(rep(0, K))))) - koreksi)

bangkitkan <- function(skenario, seed) {
  set.seed(seed)
  delta <- if (skenario == "A") delta0 else rep(0, K)
  v_th  <- rnorm(n_th, 0, 1 / sqrt(PREC_TH))
  sd_ar <- 1 / sqrt(PREC_AR1)                     # simpangan baku marginal AR1
  w <- matrix(0, g$n, n_t)
  w[, 1] <- rnorm(g$n, 0, sd_ar)
  for (j in 2:n_t) w[, j] <- RHO_AR1 * w[, j - 1] + rnorm(g$n, 0, sd_ar * sqrt(1 - RHO_AR1^2))
  eta <- INTERSEP[[skenario]] + eta_tetap(delta) + v_th[tahun_id] + w[cbind(id_kab, tt)]
  mu  <- exp(pmin(eta, 20))
  rnbinom(length(mu), size = SIZE, mu = mu)
}

# ---- 2. Model yang dipasang ---------------------------------------------
pc_prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
pc_bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                phi  = list(prior = "pc",      param = c(0.5, 2/3)))
pc_ar1  <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                rho  = list(prior = "pc.cor0", param = c(0.5, 0.5)))
x_cb  <- paste(names(cb_df),  collapse = " + ")
x_cbg <- paste(names(cbg_df), collapse = " + ")
acak  <- paste("+ offset(log_luas)",
               "+ f(id_kab, model = 'bym2', graph = g, scale.model = TRUE, constr = TRUE, hyper = pc_bym2)",
               "+ f(bl, model = 'rw2', cyclic = TRUE, scale.model = TRUE, constr = TRUE, hyper = pc_prec)",
               "+ f(tahun_id, model = 'iid', constr = TRUE, hyper = pc_prec)")
ar1   <- "+ f(t, model = 'ar1', replicate = id_rep, constr = TRUE, hyper = pc_ar1)"
rumus <- list(
  M3b = as.formula(paste("hotspot ~ 1 +", x_cb, acak)),
  M4  = as.formula(paste("hotspot ~ 1 +", x_cb, "+", x_cbg, acak)),
  M5b = as.formula(paste("hotspot ~ 1 +", x_cb, "+", x_cbg, acak, ar1)))
MODEL <- list(A = c("M3b", "M4", "M5b"), B = c("M4", "M5b"))

pasang <- function(nama, y) {
  d <- data.frame(hotspot = y, id_kab = id_kab, id_rep = id_kab, bl = bl,
                  tahun_id = tahun_id, t = tt, log_luas = log_luas, X)
  m <- tryCatch(
    inla(rumus[[nama]], family = "nbinomial", data = d,
         control.fixed   = list(correlation.matrix = TRUE),
         control.compute = list(config = TRUE),
         verbose = FALSE),
    error = function(e) { message("  ERROR INLA (", nama, "): ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)
  pakai_g <- nama != "M3b"
  nm_all <- if (pakai_g) c(names(cb_df), names(cbg_df)) else names(cb_df)
  idx <- match(nm_all, m$names.fixed)
  if (anyNA(idx)) return(NULL)
  V <- m$misc$lincomb.derived.covariance.matrix
  if (!is.null(V) && all(dim(V) == length(m$names.fixed))) {
    V <- as.matrix(V)
    Vall <- if (!is.null(rownames(V)) && all(nm_all %in% rownames(V))) V[nm_all, nm_all] else V[idx, idx]
  } else {
    sel  <- setNames(as.list(rep(1, length(nm_all))), nm_all)
    smp  <- inla.posterior.sample(1000, m, selection = sel)
    Vall <- cov(t(sapply(smp, function(s) s$latent[, 1])))
  }
  Vall <- unname(as.matrix(Vall)); b <- m$summary.fixed$mean[idx]
  it <- seq_len(K); th <- b[it]
  cp <- function(co, vc) crosspred(cb, coef = co, vcov = vc, model.link = "log",
                                   at = m02$kering, cen = m02$cen)
  p0 <- cp(th, Vall[it, it])
  out <- list(cum0 = log(p0$allRRfit[1]), cum0_se = p0$allse[1],
              lag0 = p0$matfit[1, ], lag0_se = p0$matse[1, ])
  if (pakai_g) {
    id <- K + seq_len(K); de <- b[id]
    pt <- cp(th + P_TINGGI * de,
             Vall[it, it] + P_TINGGI^2 * Vall[id, id] + P_TINGGI * (Vall[it, id] + Vall[id, it]))
    pd <- cp(de, Vall[id, id])
    out <- c(out, list(cumt = log(pt$allRRfit[1]), cumt_se = pt$allse[1],
                       lagt = pt$matfit[1, ], lagt_se = pt$matse[1, ],
                       logratio = unname(P_TINGGI * pd$allfit[1]),
                       logratio_se = unname(P_TINGGI * pd$allse[1])))
  }
  rm(m); invisible(gc())
  out
}

# ---- 3. Loop simulasi (dengan checkpoint) ------------------------------
message(sprintf("Simulasi: %d ulangan x skenario A (3 model) + B (2 model). Ulangan yang sudah ada dilewati.", N_REP))
for (r in seq_len(N_REP)) {
  for (sk in c("A", "B")) {
    f_rep <- file.path(REPDIR, sprintf("rep_%s_%03d.rds", sk, r))
    if (file.exists(f_rep)) next
    t0 <- Sys.time()
    y <- bangkitkan(sk, SEED + r * 10 + (sk == "B"))
    hasil_r <- lapply(setNames(MODEL[[sk]], MODEL[[sk]]), function(nm) {
      tryCatch(pasang(nm, y),
               error = function(e) { message("  ERROR (", nm, "): ", conditionMessage(e)); NULL })
    })
    saveRDS(list(skenario = sk, rep = r, hasil = hasil_r,
                 rata_y = mean(y), prop_nol = mean(y == 0)), f_rep)
    message(sprintf("  Ulangan %d/%d skenario %s selesai (%.1f menit)", r, N_REP, sk,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

# ---- 4. Ringkasan ------------------------------------------------------
ringkas <- file.path(OUT, "ringkasan_09.txt")
catat <- function(...) { cat(..., "\n"); cat(..., "\n", file = ringkas, append = TRUE) }
cat("", file = ringkas)

f_all <- list.files(REPDIR, pattern = "^rep_[AB]_\\d+\\.rds$", full.names = TRUE)
semua <- lapply(f_all, readRDS)
baris <- list()
for (s in semua) for (nm in names(s$hasil)) {
  h <- s$hasil[[nm]]; tr <- TRUTH[[s$skenario]]
  if (is.null(h)) { baris[[length(baris) + 1]] <- data.frame(skenario = s$skenario, rep = s$rep, model = nm, gagal = TRUE); next }
  cakup <- function(est, se, tru) abs(est - tru) <= 1.96 * se
  b <- data.frame(skenario = s$skenario, rep = s$rep, model = nm, gagal = FALSE,
                  cum0_bias = h$cum0 - tr$cum0, cum0_cov = cakup(h$cum0, h$cum0_se, tr$cum0),
                  lag0_cov = mean(cakup(h$lag0, h$lag0_se, tr$lag0)))
  if (!is.null(h$logratio)) {
    b$cumt_bias <- h$cumt - tr$cumt
    b$cumt_cov  <- cakup(h$cumt, h$cumt_se, tr$cumt)
    b$lagt_cov  <- mean(cakup(h$lagt, h$lagt_se, tr$lagt))
    b$lr_est    <- h$logratio
    b$lr_bias   <- h$logratio - tr$logratio
    b$lr_cov    <- cakup(h$logratio, h$logratio_se, tr$logratio)
    b$lr_tolak  <- abs(h$logratio / h$logratio_se) > 1.96
    b$lr_se     <- h$logratio_se
  } else {
    # M3b: satu kurva untuk semua -> dibandingkan dengan nilai sebenarnya gambut tinggi juga
    b$cumt_bias <- h$cum0 - tr$cumt
  }
  baris[[length(baris) + 1]] <- b
}
tab <- bind_rows(baris)
# pastikan semua kolom ada walaupun banyak model gagal
for (kol in c("cum0_bias", "cum0_cov", "lag0_cov", "cumt_bias", "cumt_cov", "lagt_cov",
              "lr_est", "lr_bias", "lr_cov", "lr_tolak", "lr_se")) {
  if (!kol %in% names(tab)) tab[[kol]] <- NA
}
write.csv(tab, file.path(OUT, "simulasi_per_ulangan.csv"), row.names = FALSE)

mcse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
catat("================ STUDI SIMULASI ================")
catat(sprintf("Ulangan selesai: A = %d, B = %d (target %d)",
              n_distinct(tab$rep[tab$skenario == "A"]), n_distinct(tab$rep[tab$skenario == "B"]), N_REP))
catat(sprintf("Parameter sebenarnya: size NB = %.3f | AR1 rho = %.3f, sd = %.3f | sd tahun = %.3f",
              SIZE, RHO_AR1, 1 / sqrt(PREC_AR1), 1 / sqrt(PREC_TH)))
catat(sprintf("Gambut tinggi = %.0f%% | bulan kering %.0f mm vs median %.0f mm",
              100 * P_TINGGI, m02$kering, m02$cen))
catat(sprintf("Nilai sebenarnya skenario A: RR kum non-gambut = %.2f | gambut = %.2f | rasio = %.2f",
              exp(TRUTH$A$cum0), exp(TRUTH$A$cumt), exp(TRUTH$A$logratio)))
catat("Model gagal:", sum(tab$gagal), "dari", nrow(tab))

ok_tab <- tab |> filter(!gagal)
catat("\n--- SKENARIO A: ada modifikasi gambut + AR1 ---")
catat("(bias & RMSE pada skala log; cakupan ideal 0.95; kuasa = proporsi rasio gambut signifikan)")
catat(sprintf("%-5s %10s %10s %8s %10s %10s %8s %8s %8s",
              "Model", "bias_rasio", "RMSE_rasio", "cak_ras", "bias_non", "bias_gamb",
              "cak_non", "cak_lag", "kuasa"))
for (nm in MODEL$A) {
  d <- ok_tab |> filter(skenario == "A", model == nm)
  if (!nrow(d)) next
  has_g <- "lr_bias" %in% names(d) && any(!is.na(d$lr_bias))
  catat(sprintf("%-5s %10s %10s %8s %+10.3f %+10.3f %8.2f %8.2f %8s", nm,
                if (has_g) sprintf("%+.3f", mean(d$lr_bias, na.rm = TRUE)) else "-",
                if (has_g) sprintf("%.3f", sqrt(mean(d$lr_bias^2, na.rm = TRUE))) else "-",
                if (has_g) sprintf("%.2f", mean(d$lr_cov, na.rm = TRUE)) else "-",
                mean(d$cum0_bias, na.rm = TRUE), mean(d$cumt_bias, na.rm = TRUE),
                mean(d$cum0_cov, na.rm = TRUE), mean(d$lag0_cov, na.rm = TRUE),
                if (has_g) sprintf("%.2f", mean(d$lr_tolak, na.rm = TRUE)) else "-"))
}
catat("Keterangan: cak_ras = cakupan selang 95% rasio gambut/non; cak_non = cakupan RR kum non-gambut;")
catat("            cak_lag = rata-rata cakupan kurva lag non-gambut (lag 0-6).")
catat("            M3b tidak punya efek gambut: bias_gamb = selisih kurva tunggalnya dari nilai gambut tinggi.")

catat("\n--- SKENARIO B: TIDAK ada modifikasi gambut (+AR1) ---")
catat("(galat tipe I ideal = 0.05; makin besar = makin sering 'menemukan' efek gambut palsu)")
for (nm in MODEL$B) {
  d <- ok_tab |> filter(skenario == "B", model == nm)
  if (!nrow(d)) next
  p <- mean(d$lr_tolak, na.rm = TRUE)
  catat(sprintf("  %-4s galat tipe I = %.2f (MCSE %.2f) | cakupan rasio = %.2f | bias rasio (log) = %+.3f",
                nm, p, sqrt(p * (1 - p) / nrow(d)), mean(d$lr_cov, na.rm = TRUE), mean(d$lr_bias, na.rm = TRUE)))
}

# ---- 5. Grafik ---------------------------------------------------------
if (!any(!is.na(ok_tab$lr_est))) {
  catat("\nBelum ada hasil model yang berhasil untuk digambar. Jalankan ulang script ini.")
} else {
gd <- ok_tab |> filter(!is.na(lr_est)) |>
  mutate(rasio = exp(lr_est),
         skenario = ifelse(skenario == "A", "A: ada efek gambut", "B: tanpa efek gambut"))
garis <- data.frame(skenario = c("A: ada efek gambut", "B: tanpa efek gambut"),
                    benar = c(exp(TRUTH$A$logratio), 1))
g20 <- ggplot(gd, aes(model, rasio, fill = model)) +
  geom_boxplot(alpha = 0.6, outlier.size = 0.8) +
  geom_hline(data = garis, aes(yintercept = benar), linetype = 2, colour = "firebrick") +
  facet_wrap(~ skenario, scales = "free_y") + scale_y_log10() +
  scale_fill_manual(values = c(M4 = "grey60", M5b = "steelblue")) +
  labs(x = NULL, y = "Estimasi rasio RR (gambut/non), skala log",
       title = "Simulasi: estimasi rasio efek gambut (garis merah = nilai sebenarnya)") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(OUT, "g20_simulasi_rasio.png"), g20, width = 9, height = 5, dpi = 200)

cov_df <- ok_tab |> filter(!is.na(lr_cov)) |> group_by(skenario, model) |>
  summarise(cakupan = mean(lr_cov), n = n(), .groups = "drop") |>
  mutate(lo = pmax(0, cakupan - 1.96 * sqrt(cakupan * (1 - cakupan) / n)),
         hi = pmin(1, cakupan + 1.96 * sqrt(cakupan * (1 - cakupan) / n)),
         skenario = ifelse(skenario == "A", "A: ada efek gambut", "B: tanpa efek gambut"))
g21 <- ggplot(cov_df, aes(model, cakupan, colour = model)) +
  geom_hline(yintercept = 0.95, linetype = 2) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.6) +
  facet_wrap(~ skenario) + coord_cartesian(ylim = c(0, 1)) +
  scale_colour_manual(values = c(M4 = "grey40", M5b = "steelblue")) +
  labs(x = NULL, y = "Cakupan selang 95% rasio gambut",
       title = "Simulasi: cakupan selang kredibel 95% (ideal = 0.95)") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(OUT, "g21_simulasi_cakupan.png"), g21, width = 9, height = 4.5, dpi = 200)
}

catat("\nSelesai. Hasil di folder:", normalizePath(OUT))
