# =====================================================================
# 08_sensitivitas.R
# Ide 1: uji sensitivitas model final (M5b) terhadap pilihan DLNM
#
#   S1 = lag maks 6  (seperti model utama)
#   S2 = lag maks 9
#   S3 = lag maks 12
#   S4 = lag maks 6, knot paparan di persentil 25/50/75 (bukan 10/50/90)
#   S5 = lag maks 6, 3 knot lag (lebih lentur)
# Semua: NB + DLNM x gambut + BYM2 + musim + tahun + AR1 per kabupaten (= M5b)
# Semua model memakai BARIS DATA YANG SAMA (12 bulan pertama tiap kabupaten
# dibuang) supaya WAIC bisa dibandingkan.
#
# Butuh: script 01-05 sudah dijalankan
# Cara pakai: buka di RStudio, klik "Source". Waktu: ~15-40 menit (5 model).
# Hasil: folder hasil_08/
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(dlnm); library(INLA); library(ggplot2)
})
if (basename(getwd()) != "paper_karhutla" && dir.exists("~/paper_karhutla")) {
  setwd("~/paper_karhutla")
}
DIR <- "data_kalimantan"; OUT <- "hasil_08"
dir.create(OUT, showWarnings = FALSE)
for (f in c(file.path(DIR, "panel_kalimantan.rds"), file.path(DIR, "kalimantan.graph"),
            file.path("hasil_02", "model_02.rds"))) {
  if (!file.exists(f)) stop("File ", f, " tidak ada. Jalankan script 01-05 dulu.")
}
message("Folder kerja: ", normalizePath(getwd()))

ringkas <- file.path(OUT, "ringkasan_08.txt")
catat <- function(...) { cat(..., "\n"); cat(..., "\n", file = ringkas, append = TRUE) }
cat("", file = ringkas)

panel <- readRDS(file.path(DIR, "panel_kalimantan.rds")) |> arrange(id_kab, bulan)
if (!"prop_gambut" %in% names(panel)) stop("Kolom prop_gambut belum ada. Jalankan script 04 dulu.")
m02 <- readRDS(file.path("hasil_02", "model_02.rds"))
g   <- inla.read.graph(file.path(DIR, "kalimantan.graph"))
P_TINGGI <- round(unname(quantile(unique(panel$prop_gambut), 0.9)), 2)
LAG_MAKS <- 12

# ---- 1. Variasi cross-basis -------------------------------------------
q_10_90 <- quantile(panel$hujan_mm, c(.10, .50, .90), na.rm = TRUE)
q_25_75 <- quantile(panel$hujan_mm, c(.25, .50, .75), na.rm = TRUE)
varian <- list(
  S1 = list(lag = 6,  kv = q_10_90, nk_lag = 2, ket = "lag 6, knot paparan p10/50/90, 2 knot lag (utama)"),
  S2 = list(lag = 9,  kv = q_10_90, nk_lag = 2, ket = "lag 9"),
  S3 = list(lag = 12, kv = q_10_90, nk_lag = 2, ket = "lag 12"),
  S4 = list(lag = 6,  kv = q_25_75, nk_lag = 2, ket = "lag 6, knot paparan p25/50/75"),
  S5 = list(lag = 6,  kv = q_10_90, nk_lag = 3, ket = "lag 6, 3 knot lag"))

buat_cb <- function(v) {
  crossbasis(panel$hujan_mm, lag = v$lag,
             argvar = list(fun = "ns", knots = v$kv),
             arglag = list(fun = "ns", knots = logknots(v$lag, v$nk_lag)),
             group  = panel$id_kab)
}

# Baris yang sama untuk semua varian: buang LAG_MAKS bulan pertama tiap kabupaten
urutan <- ave(seq_len(nrow(panel)), panel$id_kab, FUN = seq_along)
pakai  <- urutan > LAG_MAKS & !is.na(panel$hujan_mm)
catat("================ SENSITIVITAS M5b ================")
catat("Baris dipakai (sama untuk semua varian):", sum(pakai), "| gambut tinggi =", 100 * P_TINGGI, "%")

pc_prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
pc_bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                phi  = list(prior = "pc",      param = c(0.5, 2/3)))
pc_ar1  <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                rho  = list(prior = "pc.cor0", param = c(0.5, 0.5)))

jalan <- function(nama) {
  v  <- varian[[nama]]
  cb <- buat_cb(v); K <- ncol(cb)
  cbm <- unclass(cb)[, seq_len(K), drop = FALSE]
  cb_df  <- as.data.frame(cbm);                     names(cb_df)  <- paste0("cb",  seq_len(K))
  cbg_df <- as.data.frame(cbm * panel$prop_gambut); names(cbg_df) <- paste0("cbg", seq_len(K))
  if (anyNA(cb_df[pakai, ])) stop("Ada NA di cross-basis untuk ", nama, " (tidak seharusnya).")
  d <- data.frame(hotspot  = panel$hotspot,
                  id_kab   = as.integer(panel$id_kab),
                  id_rep   = as.integer(panel$id_kab),
                  bl       = as.integer(panel$bl),
                  tahun_id = as.integer(panel$tahun - min(panel$tahun) + 1L),
                  t        = as.integer(factor(panel$bulan)),
                  log_luas = log(panel$luas_km2))
  d <- cbind(d, cb_df, cbg_df)[pakai, ]
  d$t <- d$t - min(d$t) + 1L
  rumus <- as.formula(paste(
    "hotspot ~ 1 +", paste(names(cb_df), collapse = " + "), "+", paste(names(cbg_df), collapse = " + "),
    "+ offset(log_luas)",
    "+ f(id_kab, model = 'bym2', graph = g, scale.model = TRUE, constr = TRUE, hyper = pc_bym2)",
    "+ f(bl, model = 'rw2', cyclic = TRUE, scale.model = TRUE, constr = TRUE, hyper = pc_prec)",
    "+ f(tahun_id, model = 'iid', constr = TRUE, hyper = pc_prec)",
    "+ f(t, model = 'ar1', replicate = id_rep, constr = TRUE, hyper = pc_ar1)"))
  message("Menjalankan INLA untuk ", nama, " (", v$ket, ") ...")
  t0 <- Sys.time()
  m <- tryCatch(
    inla(rumus, family = "nbinomial", data = d,
         control.fixed   = list(correlation.matrix = TRUE),
         control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
         verbose = FALSE),
    error = function(e) { message("ERROR INLA (", nama, "): ", conditionMessage(e)); NULL })
  if (is.null(m)) { message("  -> varian ", nama, " DILEWATI karena INLA gagal."); return(NULL) }
  menit <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

  nm_all <- c(names(cb_df), names(cbg_df)); it <- seq_len(K); id <- K + seq_len(K)
  idx <- match(nm_all, m$names.fixed)
  if (anyNA(idx)) stop("Nama koefisien tidak ditemukan di hasil INLA.")
  V <- m$misc$lincomb.derived.covariance.matrix
  if (!is.null(V) && all(dim(V) == length(m$names.fixed))) {
    V <- as.matrix(V)
    Vall <- if (!is.null(rownames(V)) && all(nm_all %in% rownames(V))) V[nm_all, nm_all] else V[idx, idx]
  } else {
    sel  <- setNames(as.list(rep(1, length(nm_all))), nm_all)
    smp  <- inla.posterior.sample(2000, m, selection = sel)
    Vall <- cov(t(sapply(smp, function(s) s$latent[, 1])))
  }
  Vall <- unname(as.matrix(Vall)); b <- m$summary.fixed$mean[idx]
  th <- b[it]; de <- b[id]
  efek <- function(p, bylag = 1) crosspred(cb, coef = th + p * de,
    vcov = Vall[it, it] + p^2 * Vall[id, id] + p * (Vall[it, id] + Vall[id, it]),
    model.link = "log", at = m02$kering, cen = m02$cen, bylag = bylag)
  pd  <- crosspred(cb, coef = de, vcov = Vall[id, id], model.link = "log", at = m02$kering, cen = m02$cen)
  lf  <- unname(P_TINGGI * pd$allfit[1]); lse <- unname(P_TINGGI * pd$allse[1])
  out <- list(nama = nama, ket = v$ket, lag = v$lag, waic = m$waic$waic, menit = menit,
              p0 = efek(0), pt = efek(P_TINGGI), p0h = efek(0, 0.2), pth = efek(P_TINGGI, 0.2),
              rasio = c(fit = exp(lf), lo = exp(lf - 1.96 * lse), hi = exp(lf + 1.96 * lse),
                        p = 2 * pnorm(-abs(lf / lse))))
  rm(m); invisible(gc())
  out
}

hasil <- lapply(setNames(names(varian), names(varian)), jalan)
gagal <- names(hasil)[sapply(hasil, is.null)]
hasil <- Filter(Negate(is.null), hasil)
if (!length(hasil)) stop("Semua varian gagal. Kirim pesan ERROR di atas ke Claude.")
if (length(gagal)) catat("PERHATIAN: varian gagal & dilewati:", paste(gagal, collapse = ", "))

# ---- 2. Tabel ----------------------------------------------------------
wmin <- min(sapply(hasil, `[[`, "waic"))
catat(sprintf("\nBulan kering (%.0f mm) vs median (%.0f mm); gambut 0%% vs %.0f%%",
              m02$kering, m02$cen, 100 * P_TINGGI))
catat(sprintf("%-3s %-48s %9s %20s %20s %26s", "", "Varian", "dWAIC",
              "RRkum non-gambut", "RRkum gambut", "Rasio gambut/non"))
for (h in hasil) {
  catat(sprintf("%-3s %-48s %9.1f %20s %20s %26s", h$nama, h$ket, h$waic - wmin,
                sprintf("%.2f (%.2f-%.2f)", h$p0$allRRfit, h$p0$allRRlow, h$p0$allRRhigh),
                sprintf("%.2f (%.2f-%.2f)", h$pt$allRRfit, h$pt$allRRlow, h$pt$allRRhigh),
                sprintf("%.2f (%.2f-%.2f)", h$rasio["fit"], h$rasio["lo"], h$rasio["hi"])))
}
catat("(dWAIC = selisih dari varian terbaik; selisih < ~10 dianggap setara)")

catat("\nRR per lag di kabupaten bergambut tinggi (lag 0-12):")
catat(sprintf("  %-4s %s", "lag", paste(sprintf("%-18s", names(hasil)), collapse = "")))
for (l in 0:LAG_MAKS) {
  sel <- sapply(hasil, function(h) {
    if (l > h$lag) return(sprintf("%-18s", "-"))
    sprintf("%-18s", sprintf("%.2f (%.2f-%.2f)", h$pt$matRRfit[1, l + 1],
                             h$pt$matRRlow[1, l + 1], h$pt$matRRhigh[1, l + 1]))
  })
  catat(sprintf("  %-4d %s", l, paste(sel, collapse = "")))
}
catat("\nRR per lag di kabupaten TANPA gambut (lag 0-12):")
catat(sprintf("  %-4s %s", "lag", paste(sprintf("%-18s", names(hasil)), collapse = "")))
for (l in 0:LAG_MAKS) {
  sel <- sapply(hasil, function(h) {
    if (l > h$lag) return(sprintf("%-18s", "-"))
    sprintf("%-18s", sprintf("%.2f (%.2f-%.2f)", h$p0$matRRfit[1, l + 1],
                             h$p0$matRRlow[1, l + 1], h$p0$matRRhigh[1, l + 1]))
  })
  catat(sprintf("  %-4d %s", l, paste(sel, collapse = "")))
}

# ---- 3. Grafik ---------------------------------------------------------
lv <- c("Tanpa gambut (0%)", sprintf("Gambut tinggi (%.0f%%)", 100 * P_TINGGI))
lag_df <- bind_rows(lapply(hasil, function(h) {
  bind_rows(lapply(c("p0h", "pth"), function(k) {
    p <- h[[k]]
    data.frame(varian = paste0(h$nama, ": ", h$ket),
               kelompok = if (k == "p0h") lv[1] else lv[2],
               lag = as.numeric(sub("lag", "", colnames(p$matRRfit))),
               rr = p$matRRfit[1, ], lo = p$matRRlow[1, ], hi = p$matRRhigh[1, ])
  }))
}))
lag_df$kelompok <- factor(lag_df$kelompok, levels = lv)
g19 <- ggplot(lag_df, aes(lag, rr, colour = kelompok, fill = kelompok)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) + geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ varian, ncol = 3, labeller = label_wrap_gen(30)) +
  scale_colour_manual(values = setNames(c("steelblue", "firebrick"), lv)) +
  scale_fill_manual(values = setNames(c("steelblue", "firebrick"), lv)) +
  labs(x = "Lag (bulan)", y = "RR", colour = NULL, fill = NULL,
       title = sprintf("Sensitivitas: efek bulan kering (%.0f mm) menurut lag (model M5b)", m02$kering)) +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(OUT, "g19_sensitivitas_lag.png"), g19, width = 11, height = 7, dpi = 200)

saveRDS(lapply(hasil, function(h) h[c("nama", "ket", "lag", "waic", "rasio")]),
        file.path(OUT, "sensitivitas_ringkas.rds"))
catat("\nSelesai. Hasil di folder:", normalizePath(OUT))
