# =====================================================================
# 03b_inla_musiman.R
# Ide 1: DLNM titik panas Kalimantan - model INLA dengan efek musiman & tahunan
#
#   M2  = NB + DLNM + offset log(luas) + BYM2                (sama dengan 03a)
#   M3a = M2 + efek musiman (RW2 siklik, 12 bulan)
#   M3b = M3a + efek tahunan (iid)
#
# Butuh: hasil script 01, 02, 03a
# Cara pakai: buka di RStudio, klik "Source". Waktu: ~3-10 menit (3 model).
# Hasil: folder hasil_03/ (ringkasan_03b.txt, gambar g7-g11, m3_ringkas.rds)
# =====================================================================

# ---- 0. Paket & folder kerja -----------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(dlnm); library(INLA); library(ggplot2); library(sf)
})
sf_use_s2(FALSE)
if (basename(getwd()) != "paper_karhutla" && dir.exists("~/paper_karhutla")) {
  setwd("~/paper_karhutla")
}
DIR <- "data_kalimantan"; OUT <- "hasil_03"
dir.create(OUT, showWarnings = FALSE)
for (f in c(file.path(DIR, "panel_kalimantan.rds"), file.path(DIR, "kalimantan.graph"),
            file.path(DIR, "kab_kalimantan.gpkg"), file.path("hasil_02", "model_02.rds"))) {
  if (!file.exists(f)) stop("File ", f, " tidak ada. Jalankan script 01 dan 02 dulu.")
}
message("Folder kerja: ", normalizePath(getwd()))

ringkas <- file.path(OUT, "ringkasan_03b.txt")
catat <- function(...) { cat(..., "\n"); cat(..., "\n", file = ringkas, append = TRUE) }
cat("", file = ringkas)

# ---- 1. Data & cross-basis (sama persis dengan script 02/03a) ---------
panel <- readRDS(file.path(DIR, "panel_kalimantan.rds")) |> arrange(id_kab, bulan)
kab   <- st_read(file.path(DIR, "kab_kalimantan.gpkg"), quiet = TRUE)
m02   <- readRDS(file.path("hasil_02", "model_02.rds"))
cb    <- m02$cb
if (nrow(cb) != nrow(panel)) stop("Jumlah baris cross-basis tidak cocok dengan panel. Jalankan ulang script 02.")

cb_df <- as.data.frame(unclass(cb)[, , drop = FALSE])
names(cb_df) <- paste0("cb", seq_len(ncol(cb_df)))
dat <- bind_cols(panel, cb_df) |> mutate(log_luas = log(luas_km2))
dat <- dat[complete.cases(cb_df) & !is.na(dat$hujan_mm), ]

dat_m <- data.frame(hotspot  = dat$hotspot,
                    id_kab   = as.integer(dat$id_kab),
                    bl       = as.integer(dat$bl),
                    tahun_id = as.integer(dat$tahun - min(dat$tahun) + 1L),
                    log_luas = dat$log_luas)
dat_m <- cbind(dat_m, as.data.frame(dat[, names(cb_df)]))
g <- inla.read.graph(file.path(DIR, "kalimantan.graph"))

# ---- 2. Model ----------------------------------------------------------
pc_prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
pc_bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                phi  = list(prior = "pc",      param = c(0.5, 2/3)))

dasar <- paste("hotspot ~ 1 +", paste(names(cb_df), collapse = " + "),
               "+ offset(log_luas)",
               "+ f(id_kab, model = 'bym2', graph = g, scale.model = TRUE, constr = TRUE, hyper = pc_bym2)")
musim <- "+ f(bl, model = 'rw2', cyclic = TRUE, scale.model = TRUE, constr = TRUE, hyper = pc_prec)"
tahun <- "+ f(tahun_id, model = 'iid', constr = TRUE, hyper = pc_prec)"

rumus <- list(M2  = as.formula(dasar),
              M3a = as.formula(paste(dasar, musim)),
              M3b = as.formula(paste(dasar, musim, tahun)))

jalan_inla <- function(nama) {
  message("Menjalankan INLA untuk ", nama, " ...")
  t0 <- Sys.time()
  m <- tryCatch(
    inla(rumus[[nama]], family = "nbinomial", data = dat_m,
         control.fixed     = list(correlation.matrix = TRUE),
         control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
         control.predictor = list(compute = TRUE, link = 1),
         verbose = FALSE),
    error = function(e) { message("ERROR INLA (", nama, "): ", conditionMessage(e)); NULL })
  if (is.null(m)) stop("INLA gagal untuk ", nama, ". Lihat pesan ERROR di atas.")
  m$durasi_menit <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  m
}
mods <- lapply(setNames(names(rumus), names(rumus)), jalan_inla)

# Koefisien & kovarians cross-basis dari hasil INLA
ambil_cb <- function(m) {
  idx <- match(names(cb_df), m$names.fixed)
  if (anyNA(idx)) stop("Nama koefisien cross-basis tidak ditemukan di hasil INLA.")
  V <- m$misc$lincomb.derived.covariance.matrix
  if (!is.null(V) && all(dim(V) == length(m$names.fixed))) {
    V <- as.matrix(V)
    vc <- if (!is.null(rownames(V)) && all(names(cb_df) %in% rownames(V))) {
      V[names(cb_df), names(cb_df)]
    } else V[idx, idx]
  } else {
    sel <- setNames(as.list(rep(1, length(idx))), names(cb_df))
    smp <- inla.posterior.sample(2000, m, selection = sel)
    vc  <- cov(t(sapply(smp, function(s) s$latent[, 1])))
  }
  list(coef = m$summary.fixed$mean[idx], vcov = unname(vc))
}
koef <- lapply(mods, ambil_cb)
pred_kering <- lapply(koef, function(k)
  crosspred(cb, coef = k$coef, vcov = k$vcov, model.link = "log", at = m02$kering, cen = m02$cen))

# ---- 3. Tabel perbandingan --------------------------------------------
catat("================ PERBANDINGAN MODEL ================")
catat(sprintf("%-4s %10s %10s %8s %7s  %s", "Model", "WAIC", "DIC", "CPOgagal", "menit",
              "RR kumulatif lag 0-6 (kering vs median)"))
waic_min <- min(sapply(mods, function(m) m$waic$waic))
for (nm in names(mods)) {
  m <- mods[[nm]]; p <- pred_kering[[nm]]
  catat(sprintf("%-5s %10.1f %10.1f %8d %7.1f  %.2f (%.2f - %.2f)", nm, m$waic$waic, m$dic$dic,
                sum(m$cpo$failure > 0, na.rm = TRUE), m$durasi_menit,
                p$allRRfit[1], p$allRRlow[1], p$allRRhigh[1]))
}
catat("Selisih WAIC terhadap model terbaik (0 = terbaik):")
for (nm in names(mods)) catat(sprintf("  %-4s : %.1f", nm, mods[[nm]]$waic$waic - waic_min))

catat(sprintf("\nRR per lag, bulan kering (%.0f mm) vs median (%.0f mm):", m02$kering, m02$cen))
catat(sprintf("  %-5s %-22s %-22s %-22s", "lag", "M2", "M3a (+musim)", "M3b (+musim+tahun)"))
for (l in 0:m02$LAG) {
  sel <- sapply(names(mods), function(nm) {
    p <- pred_kering[[nm]]
    sprintf("%.2f (%.2f-%.2f)", p$matRRfit[1, l + 1], p$matRRlow[1, l + 1], p$matRRhigh[1, l + 1])
  })
  catat(sprintf("  %-5d %-22s %-22s %-22s", l, sel[1], sel[2], sel[3]))
}

catat("\nHiperparameter M3b:")
hp <- mods$M3b$summary.hyperpar
for (i in seq_len(nrow(hp))) {
  catat(sprintf("  %-55s %8.3f (%.3f - %.3f)", rownames(hp)[i],
                hp[i, "mean"], hp[i, "0.025quant"], hp[i, "0.975quant"]))
}

# ---- 4. Grafik ---------------------------------------------------------
terbaik <- names(which.min(sapply(mods, function(m) m$waic$waic)))
catat("\nModel terbaik menurut WAIC:", terbaik, "-> dipakai untuk grafik g7-g9")
kb <- koef[[terbaik]]
at <- seq(floor(min(panel$hujan_mm, na.rm = TRUE)),
          ceiling(quantile(panel$hujan_mm, .99, na.rm = TRUE)), length.out = 80)
pred_all <- crosspred(cb, coef = kb$coef, vcov = kb$vcov, model.link = "log",
                      at = at, cen = m02$cen, bylag = 0.2)
pred_kh  <- crosspred(cb, coef = kb$coef, vcov = kb$vcov, model.link = "log",
                      at = m02$kering, cen = m02$cen, bylag = 0.2)

png(file.path(OUT, "g7_permukaan_3d_terbaik.png"), width = 1600, height = 1300, res = 200)
plot(pred_all, xlab = "Curah hujan (mm/bulan)", ylab = "Lag (bulan)", zlab = "RR",
     theta = 210, phi = 30, lphi = 30, main = paste0("Exposure-lag-response (", terbaik, ")"))
dev.off()

png(file.path(OUT, "g8_kumulatif_terbaik.png"), width = 1600, height = 1100, res = 200)
plot(pred_all, "overall", xlab = "Curah hujan (mm/bulan)", ylab = "RR kumulatif (lag 0-6)",
     main = paste0("Efek kumulatif curah hujan (", terbaik, ")"), log = "y", col = "firebrick", lwd = 2)
abline(v = m02$kering, lty = 2)
dev.off()

png(file.path(OUT, "g9_lag_kering_terbaik.png"), width = 1600, height = 1100, res = 200)
plot(pred_kh, "slices", var = m02$kering, xlab = "Lag (bulan)", ylab = "RR",
     main = sprintf("Efek bulan kering (%.0f mm) menurut lag (%s)", m02$kering, terbaik),
     col = "firebrick", lwd = 2)
dev.off()

# Efek musiman (dari M3a & M3b)
musim_df <- bind_rows(lapply(c("M3a", "M3b"), function(nm) {
  s <- mods[[nm]]$summary.random$bl
  data.frame(model = nm, bl = s$ID, rr = exp(s$mean),
             lo = exp(s$`0.025quant`), hi = exp(s$`0.975quant`))
}))
p_musim <- ggplot(musim_df, aes(bl, rr, colour = model, fill = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_hline(yintercept = 1, linetype = 2) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) + scale_y_log10() +
  labs(x = NULL, y = "RR musiman (skala log)",
       title = "Efek musiman setelah memperhitungkan curah hujan") +
  theme_minimal()
ggsave(file.path(OUT, "g10_efek_musiman.png"), p_musim, width = 8, height = 4.5, dpi = 200)

# Peta risiko spasial (model terbaik)
re <- mods[[terbaik]]$summary.random$id_kab[seq_len(g$n), ]
peta <- kab |> mutate(rr_spasial = exp(re$mean[id_kab]))
p_peta <- ggplot(peta) +
  geom_sf(aes(fill = rr_spasial), colour = "grey40", linewidth = 0.1) +
  scale_fill_distiller(palette = "RdYlBu", trans = "log10", name = "RR spasial") +
  labs(title = paste0("Risiko spasial titik panas (", terbaik, ")"),
       subtitle = c(M2  = "Setelah memperhitungkan curah hujan & luas wilayah",
                    M3a = "Setelah memperhitungkan curah hujan, luas wilayah & musim",
                    M3b = "Setelah memperhitungkan curah hujan, luas wilayah, musim & tahun")[[terbaik]]) +
  theme_minimal()
ggsave(file.path(OUT, "g11_peta_risiko_spasial.png"), p_peta, width = 7, height = 7, dpi = 200)

saveRDS(list(waic = sapply(mods, function(m) m$waic$waic),
             dic  = sapply(mods, function(m) m$dic$dic),
             koef = koef, terbaik = terbaik,
             hyper = lapply(mods, function(m) m$summary.hyperpar),
             musim = musim_df,
             spasial = st_drop_geometry(peta)[, c("id_kab", "kabupaten", "rr_spasial")]),
        file.path(OUT, "m3_ringkas.rds"))

catat("\nSelesai. Hasil di folder:", normalizePath(OUT))
