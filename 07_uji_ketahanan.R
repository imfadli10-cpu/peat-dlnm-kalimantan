# =====================================================================
# 07_uji_ketahanan.R
# Ide 1: uji ketahanan (robustness) temuan gambut terhadap autokorelasi waktu
#
#   M4  = model utama (script 05)
#   M5a = M4 + AR1 bersama antar-bulan (guncangan regional)
#   M5b = M4 + AR1 per kabupaten (memori bulan sebelumnya di tiap kabupaten)
#
# Dicek: WAIC, autokorelasi residual, PIT, dan apakah efek gambut tetap kuat.
# Butuh: script 01-05 sudah dijalankan
# Cara pakai: buka di RStudio, klik "Source". Waktu: ~3-15 menit (3 model).
# Hasil: folder hasil_07/
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(dlnm); library(INLA); library(ggplot2)
})
if (basename(getwd()) != "paper_karhutla" && dir.exists("~/paper_karhutla")) {
  setwd("~/paper_karhutla")
}
DIR <- "data_kalimantan"; OUT <- "hasil_07"
dir.create(OUT, showWarnings = FALSE)
for (f in c(file.path(DIR, "panel_kalimantan.rds"), file.path(DIR, "kalimantan.graph"),
            file.path("hasil_02", "model_02.rds"))) {
  if (!file.exists(f)) stop("File ", f, " tidak ada. Jalankan script 01-05 dulu.")
}
message("Folder kerja: ", normalizePath(getwd()))

ringkas <- file.path(OUT, "ringkasan_07.txt")
catat <- function(...) { cat(..., "\n"); cat(..., "\n", file = ringkas, append = TRUE) }
cat("", file = ringkas)

# ---- 1. Data (sama dengan script 05/06) --------------------------------
panel <- readRDS(file.path(DIR, "panel_kalimantan.rds")) |> arrange(id_kab, bulan)
if (!"prop_gambut" %in% names(panel)) stop("Kolom prop_gambut belum ada. Jalankan script 04 dulu.")
m02 <- readRDS(file.path("hasil_02", "model_02.rds"))
cb  <- m02$cb
if (nrow(cb) != nrow(panel)) stop("Jumlah baris cross-basis tidak cocok dengan panel. Jalankan ulang script 02.")
K <- ncol(cb)
cb_mat <- unclass(cb)[, seq_len(K), drop = FALSE]
cb_df  <- as.data.frame(cb_mat);                     names(cb_df)  <- paste0("cb",  seq_len(K))
cbg_df <- as.data.frame(cb_mat * panel$prop_gambut); names(cbg_df) <- paste0("cbg", seq_len(K))

ok   <- complete.cases(cb_df) & !is.na(panel$hujan_mm)
info <- panel[ok, c("id_kab", "bulan")]
dat0 <- data.frame(hotspot  = panel$hotspot,
                   id_kab   = as.integer(panel$id_kab),
                   id_rep   = as.integer(panel$id_kab),          # salinan untuk replicate AR1
                   bl       = as.integer(panel$bl),
                   tahun_id = as.integer(panel$tahun - min(panel$tahun) + 1L),
                   t        = as.integer(factor(panel$bulan)),
                   log_luas = log(panel$luas_km2))
dat0 <- cbind(dat0, cb_df, cbg_df)[ok, ]
dat0$t <- dat0$t - min(dat0$t) + 1L
g <- inla.read.graph(file.path(DIR, "kalimantan.graph"))

pc_prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
pc_bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                phi  = list(prior = "pc",      param = c(0.5, 2/3)))
pc_ar1  <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                rho  = list(prior = "pc.cor0", param = c(0.5, 0.5)))
dasar <- paste("hotspot ~ 1 +", paste(names(cb_df), collapse = " + "), "+",
               paste(names(cbg_df), collapse = " + "),
               "+ offset(log_luas)",
               "+ f(id_kab, model = 'bym2', graph = g, scale.model = TRUE, constr = TRUE, hyper = pc_bym2)",
               "+ f(bl, model = 'rw2', cyclic = TRUE, scale.model = TRUE, constr = TRUE, hyper = pc_prec)",
               "+ f(tahun_id, model = 'iid', constr = TRUE, hyper = pc_prec)")
rumus <- list(
  M4  = as.formula(dasar),
  M5a = as.formula(paste(dasar, "+ f(t, model = 'ar1', constr = TRUE, hyper = pc_ar1)")),
  M5b = as.formula(paste(dasar, "+ f(t, model = 'ar1', replicate = id_rep, constr = TRUE, hyper = pc_ar1)")))

fit <- function(nama) {
  message("Menjalankan INLA untuk ", nama, " ...")
  t0 <- Sys.time()
  m <- tryCatch(
    inla(rumus[[nama]], family = "nbinomial", data = dat0,
         control.fixed     = list(correlation.matrix = TRUE),
         control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
         control.predictor = list(compute = TRUE, link = 1),
         verbose = FALSE),
    error = function(e) { message("ERROR INLA (", nama, "): ", conditionMessage(e)); NULL })
  if (is.null(m)) stop("INLA gagal untuk ", nama, ". Lihat pesan ERROR di atas.")
  m$durasi_menit <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  m
}

# ---- 2. Fungsi ringkasan -----------------------------------------------
nm_all <- c(names(cb_df), names(cbg_df)); it <- seq_len(K); id <- K + seq_len(K)
P_TINGGI <- round(unname(quantile(unique(panel$prop_gambut), 0.9)), 2)

ringkas_model <- function(m) {
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
  efek <- function(p, bylag = 1) {
    crosspred(cb, coef = th + p * de,
              vcov = Vall[it, it] + p^2 * Vall[id, id] + p * (Vall[it, id] + Vall[id, it]),
              model.link = "log", at = m02$kering, cen = m02$cen, bylag = bylag)
  }
  pd  <- crosspred(cb, coef = de, vcov = Vall[id, id], model.link = "log",
                   at = m02$kering, cen = m02$cen)
  lf  <- unname(P_TINGGI * pd$allfit[1]); lse <- unname(P_TINGGI * pd$allse[1])

  mu   <- m$summary.fitted.values$mean[seq_len(nrow(dat0))]
  hp   <- m$summary.hyperpar
  size <- hp[grep("size", rownames(hp)), "mean"][1]
  res  <- (dat0$hotspot - mu) / sqrt(mu + mu^2 / size)
  d    <- info |> mutate(res = res) |> arrange(id_kab, bulan)
  acf1 <- sapply(1:3, function(k) {
    dd <- d |> group_by(id_kab) |> mutate(r = dplyr::lag(res, k)) |> ungroup()
    suppressWarnings(cor(dd$res, dd$r, use = "complete.obs"))
  })
  pit <- m$cpo$pit
  list(waic = m$waic$waic, dic = m$dic$dic, menit = m$durasi_menit,
       cpo_gagal = sum(m$cpo$failure > 0, na.rm = TRUE),
       p0 = efek(0), pt = efek(P_TINGGI),
       p0h = efek(0, 0.2), pth = efek(P_TINGGI, 0.2),
       rasio = c(fit = exp(lf), lo = exp(lf - 1.96 * lse), hi = exp(lf + 1.96 * lse),
                 p = 2 * pnorm(-abs(lf / lse))),
       acf = acf1, pit_hi = mean(pit > 0.95, na.rm = TRUE), pit_lo = mean(pit < 0.05, na.rm = TRUE),
       hyper = hp)
}

hasil <- list()
for (nm in names(rumus)) {
  m <- fit(nm)
  hasil[[nm]] <- ringkas_model(m)
  rm(m); invisible(gc())
}

# ---- 3. Tabel ----------------------------------------------------------
catat("================ UJI KETAHANAN: AUTOKORELASI WAKTU ================")
catat(sprintf("%-5s %10s %8s %8s %8s %8s %8s %7s", "Model", "WAIC", "dWAIC",
              "ACF1", "ACF2", "ACF3", "PIT>.95", "menit"))
wmin <- min(sapply(hasil, `[[`, "waic"))
for (nm in names(hasil)) {
  h <- hasil[[nm]]
  catat(sprintf("%-5s %10.1f %8.1f %+8.3f %+8.3f %+8.3f %8.3f %7.1f", nm, h$waic, h$waic - wmin,
                h$acf[1], h$acf[2], h$acf[3], h$pit_hi, h$menit))
}
catat("(dWAIC = selisih dari model terbaik; ACF = autokorelasi residual lag 1-3; ideal |ACF| < 0.1)")

catat(sprintf("\nTEMUAN GAMBUT - bulan kering (%.0f mm) vs median (%.0f mm), gambut 0%% vs %.0f%%:",
              m02$kering, m02$cen, 100 * P_TINGGI))
catat(sprintf("%-5s %-26s %-26s %-30s", "Model", "RR kumulatif gambut 0%", "RR kumulatif gambut tinggi",
              "Rasio (gambut/non) & p"))
for (nm in names(hasil)) {
  h <- hasil[[nm]]
  catat(sprintf("%-5s %-26s %-26s %-30s", nm,
                sprintf("%.2f (%.2f-%.2f)", h$p0$allRRfit, h$p0$allRRlow, h$p0$allRRhigh),
                sprintf("%.2f (%.2f-%.2f)", h$pt$allRRfit, h$pt$allRRlow, h$pt$allRRhigh),
                sprintf("%.2f (%.2f-%.2f), p=%.4f", h$rasio["fit"], h$rasio["lo"], h$rasio["hi"], h$rasio["p"])))
}

catat("\nRR per lag di kabupaten bergambut tinggi (apakah efek tetap bertahan lama?):")
catat(sprintf("  %-5s %s", "lag", paste(sprintf("%-20s", names(hasil)), collapse = "")))
for (l in 0:m02$LAG) {
  sel <- sapply(hasil, function(h) sprintf("%-20s", sprintf("%.2f (%.2f-%.2f)",
                h$pt$matRRfit[1, l + 1], h$pt$matRRlow[1, l + 1], h$pt$matRRhigh[1, l + 1])))
  catat(sprintf("  %-5d %s", l, paste(sel, collapse = "")))
}

for (nm in c("M5a", "M5b")) {
  hp <- hasil[[nm]]$hyper; r <- grep("Rho|rho", rownames(hp))
  if (length(r)) catat(sprintf("\n%s: korelasi AR1 (rho) = %.3f (%.3f - %.3f)", nm,
                               hp[r[1], "mean"], hp[r[1], "0.025quant"], hp[r[1], "0.975quant"]))
}

# ---- 4. Grafik ---------------------------------------------------------
lag_df <- bind_rows(lapply(names(hasil), function(nm) {
  bind_rows(lapply(c("p0h", "pth"), function(k) {
    p <- hasil[[nm]][[k]]
    data.frame(model = nm,
               kelompok = if (k == "p0h") "Tanpa gambut (0%)" else sprintf("Gambut tinggi (%.0f%%)", 100 * P_TINGGI),
               lag = as.numeric(sub("lag", "", colnames(p$matRRfit))),
               rr = p$matRRfit[1, ], lo = p$matRRlow[1, ], hi = p$matRRhigh[1, ])
  }))
}))
lv <- c("Tanpa gambut (0%)", sprintf("Gambut tinggi (%.0f%%)", 100 * P_TINGGI))
lag_df$kelompok <- factor(lag_df$kelompok, levels = lv)
g18 <- ggplot(lag_df, aes(lag, rr, colour = kelompok, fill = kelompok)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ model, nrow = 1) +
  scale_colour_manual(values = setNames(c("steelblue", "firebrick"), lv)) +
  scale_fill_manual(values = setNames(c("steelblue", "firebrick"), lv)) +
  labs(x = "Lag (bulan)", y = "RR", colour = NULL, fill = NULL,
       title = sprintf("Uji ketahanan: efek bulan kering (%.0f mm) menurut lag", m02$kering)) +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(OUT, "g18_ketahanan_lag.png"), g18, width = 11, height = 4.5, dpi = 200)

saveRDS(lapply(hasil, function(h) h[c("waic", "dic", "rasio", "acf", "pit_hi", "pit_lo", "hyper")]),
        file.path(OUT, "m5_ringkas.rds"))
catat("\nSelesai. Hasil di folder:", normalizePath(OUT))
