# =====================================================================
# 05_inla_gambut.R
# Ide 1: MODEL UTAMA - DLNM dengan modifikasi efek oleh lahan gambut
#
#   M4 = M3b + interaksi cross-basis x proporsi gambut
#        log(mu) = ... + cb %*% theta + prop_gambut * (cb %*% delta) + ...
#   Efek hujan di kabupaten dengan proporsi gambut p:  theta + p * delta
#
# Butuh: script 01, 02, 03b, 04 sudah dijalankan
# Cara pakai: buka di RStudio, klik "Source". Waktu: ~1-5 menit.
# Hasil: folder hasil_05/
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(dlnm); library(INLA); library(ggplot2)
})
if (basename(getwd()) != "paper_karhutla" && dir.exists("~/paper_karhutla")) {
  setwd("~/paper_karhutla")
}
DIR <- "data_kalimantan"; OUT <- "hasil_05"
dir.create(OUT, showWarnings = FALSE)
for (f in c(file.path(DIR, "panel_kalimantan.rds"), file.path(DIR, "kalimantan.graph"),
            file.path("hasil_02", "model_02.rds"), file.path("hasil_03", "m3_ringkas.rds"))) {
  if (!file.exists(f)) stop("File ", f, " tidak ada. Jalankan script 01-03b dulu.")
}
message("Folder kerja: ", normalizePath(getwd()))

ringkas <- file.path(OUT, "ringkasan_05.txt")
catat <- function(...) { cat(..., "\n"); cat(..., "\n", file = ringkas, append = TRUE) }
cat("", file = ringkas)

# ---- 1. Data -----------------------------------------------------------
panel <- readRDS(file.path(DIR, "panel_kalimantan.rds")) |> arrange(id_kab, bulan)
if (!"prop_gambut" %in% names(panel)) stop("Kolom prop_gambut belum ada. Jalankan script 04 dulu.")
m02 <- readRDS(file.path("hasil_02", "model_02.rds"))
m3  <- readRDS(file.path("hasil_03", "m3_ringkas.rds"))
cb  <- m02$cb
if (nrow(cb) != nrow(panel)) stop("Jumlah baris cross-basis tidak cocok dengan panel. Jalankan ulang script 02.")

K <- ncol(cb)
cb_mat  <- unclass(cb)[, seq_len(K), drop = FALSE]
cb_df   <- as.data.frame(cb_mat);                  names(cb_df)  <- paste0("cb",  seq_len(K))
cbg_df  <- as.data.frame(cb_mat * panel$prop_gambut); names(cbg_df) <- paste0("cbg", seq_len(K))

ok <- complete.cases(cb_df) & !is.na(panel$hujan_mm)
dat_m <- data.frame(hotspot  = panel$hotspot,
                    id_kab   = as.integer(panel$id_kab),
                    bl       = as.integer(panel$bl),
                    tahun_id = as.integer(panel$tahun - min(panel$tahun) + 1L),
                    log_luas = log(panel$luas_km2))
dat_m <- cbind(dat_m, cb_df, cbg_df)[ok, ]
g <- inla.read.graph(file.path(DIR, "kalimantan.graph"))

catat("================ DATA ================")
catat("Baris model:", nrow(dat_m), "| kabupaten:", n_distinct(dat_m$id_kab))
catat(sprintf("Proporsi gambut: min %.2f | median %.2f | maks %.2f",
              min(panel$prop_gambut), median(unique(panel$prop_gambut)), max(panel$prop_gambut)))

# ---- 2. Model M4 -------------------------------------------------------
pc_prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
pc_bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                phi  = list(prior = "pc",      param = c(0.5, 2/3)))
rumus_m4 <- as.formula(paste(
  "hotspot ~ 1 +", paste(names(cb_df), collapse = " + "), "+",
  paste(names(cbg_df), collapse = " + "),
  "+ offset(log_luas)",
  "+ f(id_kab, model = 'bym2', graph = g, scale.model = TRUE, constr = TRUE, hyper = pc_bym2)",
  "+ f(bl, model = 'rw2', cyclic = TRUE, scale.model = TRUE, constr = TRUE, hyper = pc_prec)",
  "+ f(tahun_id, model = 'iid', constr = TRUE, hyper = pc_prec)"))

message("Menjalankan INLA untuk M4 ...")
t0 <- Sys.time()
m4 <- tryCatch(
  inla(rumus_m4, family = "nbinomial", data = dat_m,
       control.fixed     = list(correlation.matrix = TRUE),
       control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
       control.predictor = list(compute = TRUE, link = 1),
       verbose = FALSE),
  error = function(e) { message("ERROR INLA: ", conditionMessage(e)); NULL })
if (is.null(m4)) stop("INLA gagal. Kirim pesan ERROR di atas ke Claude.")
durasi <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

# ---- 3. Koefisien & kovarians gabungan (theta, delta) ------------------
nm_all <- c(names(cb_df), names(cbg_df))
idx <- match(nm_all, m4$names.fixed)
if (anyNA(idx)) stop("Nama koefisien tidak ditemukan di hasil INLA.")
V <- m4$misc$lincomb.derived.covariance.matrix
if (!is.null(V) && all(dim(V) == length(m4$names.fixed))) {
  V <- as.matrix(V)
  Vall <- if (!is.null(rownames(V)) && all(nm_all %in% rownames(V))) V[nm_all, nm_all] else V[idx, idx]
  sumber_v <- "INLA (correlation.matrix)"
} else {
  sel  <- setNames(as.list(rep(1, length(nm_all))), nm_all)
  smp  <- inla.posterior.sample(2000, m4, selection = sel)
  Vall <- cov(t(sapply(smp, function(s) s$latent[, 1])))
  sumber_v <- "2000 sampel posterior"
}
Vall  <- unname(as.matrix(Vall))
b_all <- m4$summary.fixed$mean[idx]
it <- seq_len(K); id <- K + seq_len(K)
theta <- b_all[it]; delta <- b_all[id]

# Efek pada tingkat gambut p: coef = theta + p*delta ; vcov = Vtt + p^2 Vdd + p(Vtd + Vdt)
efek_p <- function(p, at, bylag = 1) {
  co <- theta + p * delta
  vc <- Vall[it, it] + p^2 * Vall[id, id] + p * (Vall[it, id] + Vall[id, it])
  crosspred(cb, coef = co, vcov = vc, model.link = "log", at = at, cen = m02$cen, bylag = bylag)
}
P_TINGGI <- round(unname(quantile(unique(panel$prop_gambut), 0.9)), 2)   # kabupaten bergambut
P_LEVEL  <- c(tanpa_gambut = 0, gambut_tinggi = P_TINGGI)

# ---- 4. Ringkasan ------------------------------------------------------
catat("\n================ MODEL M4 (M3b + modifikasi gambut) ================")
catat("Waktu hitung:", durasi, "menit | kovarians dari:", sumber_v)
catat(sprintf("WAIC M4  = %.1f | DIC = %.1f | CPO gagal = %d",
              m4$waic$waic, m4$dic$dic, sum(m4$cpo$failure > 0, na.rm = TRUE)))
catat(sprintf("WAIC M3b = %.1f  (dari script 03b)", m3$waic[["M3b"]]))
catat(sprintf("Selisih WAIC (M3b - M4) = %.1f  (positif = modifikasi gambut memperbaiki model)",
              m3$waic[["M3b"]] - m4$waic$waic))

catat(sprintf("\nBulan kering (%.0f mm) vs median (%.0f mm):", m02$kering, m02$cen))
catat(sprintf("Tingkat gambut dibandingkan: 0%% vs %.0f%% (persentil 90 kabupaten)", 100 * P_TINGGI))
pk <- lapply(P_LEVEL, efek_p, at = m02$kering)
for (nm in names(pk)) {
  catat(sprintf("  RR kumulatif lag 0-%d, gambut %3.0f%% : %.2f (%.2f - %.2f)", m02$LAG,
                100 * P_LEVEL[[nm]], pk[[nm]]$allRRfit[1], pk[[nm]]$allRRlow[1], pk[[nm]]$allRRhigh[1]))
}

# Rasio RR kumulatif (gambut tinggi / tanpa gambut) = exp(p * b'delta)
pd <- crosspred(cb, coef = delta, vcov = Vall[id, id], model.link = "log",
                at = m02$kering, cen = m02$cen)
lf <- P_TINGGI * pd$allfit[1]; lse <- P_TINGGI * pd$allse[1]
z  <- lf / lse
catat(sprintf("  Rasio RR kumulatif (gambut %.0f%% / 0%%): %.2f (%.2f - %.2f) | z = %.2f, p = %.4f",
              100 * P_TINGGI, exp(lf), exp(lf - 1.96 * lse), exp(lf + 1.96 * lse), z, 2 * pnorm(-abs(z))))

catat("\n  RR per lag        gambut 0%                gambut tinggi")
for (l in 0:m02$LAG) {
  a <- pk$tanpa_gambut; b <- pk$gambut_tinggi
  catat(sprintf("  lag %d :   %.2f (%.2f - %.2f)      %.2f (%.2f - %.2f)", l,
                a$matRRfit[1, l + 1], a$matRRlow[1, l + 1], a$matRRhigh[1, l + 1],
                b$matRRfit[1, l + 1], b$matRRlow[1, l + 1], b$matRRhigh[1, l + 1]))
}

catat("\nHiperparameter M4:")
hp <- m4$summary.hyperpar
for (i in seq_len(nrow(hp))) {
  catat(sprintf("  %-55s %8.3f (%.3f - %.3f)", rownames(hp)[i],
                hp[i, "mean"], hp[i, "0.025quant"], hp[i, "0.975quant"]))
}

# ---- 5. Grafik ---------------------------------------------------------
lab <- c(tanpa_gambut = "Tanpa gambut (0%)",
         gambut_tinggi = sprintf("Gambut tinggi (%.0f%%)", 100 * P_TINGGI))
warna <- setNames(c("steelblue", "firebrick"), unname(lab))   # non-gambut biru, gambut merah

# Kurva lag untuk bulan kering
lag_df <- bind_rows(lapply(names(P_LEVEL), function(nm) {
  p <- efek_p(P_LEVEL[[nm]], at = m02$kering, bylag = 0.2)
  data.frame(kelompok = lab[[nm]], lag = as.numeric(sub("lag", "", colnames(p$matRRfit))),
             rr = p$matRRfit[1, ], lo = p$matRRlow[1, ], hi = p$matRRhigh[1, ])
}))
g13 <- ggplot(lag_df, aes(lag, rr, colour = kelompok, fill = kelompok)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_hline(yintercept = 1, linetype = 2) +
  scale_colour_manual(values = warna, breaks = unname(lab)) +
  scale_fill_manual(values = warna, breaks = unname(lab)) +
  labs(x = "Lag (bulan)", y = "RR", colour = NULL, fill = NULL,
       title = sprintf("Efek bulan kering (%.0f mm) menurut lag: gambut vs non-gambut", m02$kering)) +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(OUT, "g13_lag_gambut_vs_non.png"), g13, width = 8, height = 5, dpi = 200)

# Efek kumulatif menurut curah hujan
at <- seq(floor(min(panel$hujan_mm, na.rm = TRUE)),
          ceiling(quantile(panel$hujan_mm, .99, na.rm = TRUE)), length.out = 80)
kum_df <- bind_rows(lapply(names(P_LEVEL), function(nm) {
  p <- efek_p(P_LEVEL[[nm]], at = at)
  data.frame(kelompok = lab[[nm]], hujan = p$predvar,
             rr = p$allRRfit, lo = p$allRRlow, hi = p$allRRhigh)
}))
g14 <- ggplot(kum_df, aes(hujan, rr, colour = kelompok, fill = kelompok)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_hline(yintercept = 1, linetype = 2) +
  geom_vline(xintercept = m02$kering, linetype = 3) +
  scale_y_log10() +
  scale_colour_manual(values = warna, breaks = unname(lab)) +
  scale_fill_manual(values = warna, breaks = unname(lab)) +
  labs(x = "Curah hujan (mm/bulan)", y = "RR kumulatif lag 0-6 (skala log)",
       colour = NULL, fill = NULL,
       title = "Efek kumulatif curah hujan: gambut vs non-gambut") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(OUT, "g14_kumulatif_gambut_vs_non.png"), g14, width = 8, height = 5, dpi = 200)

saveRDS(list(theta = theta, delta = delta, Vall = Vall, K = K, P_TINGGI = P_TINGGI,
             waic = m4$waic$waic, dic = m4$dic$dic, hyper = m4$summary.hyperpar,
             lag_df = lag_df, kum_df = kum_df),
        file.path(OUT, "m4_ringkas.rds"))

catat("\nSelesai. Hasil di folder:", normalizePath(OUT))
