# =====================================================================
# 04_gambut.R
# Ide 1: hitung proporsi lahan gambut per kabupaten Kalimantan
#
# Sumber peta: PEATMAP (Xu dkk. 2018, CATENA 160:134-140), file SEA_Peatland.shp
#   https://archive.researchdata.leeds.ac.uk/251/  (lisensi CC BY 4.0)
# Script mencari SEA_Peatland.shp otomatis di folder Downloads (termasuk subfolder).
#
# Output:
#   - kolom prop_gambut ditambahkan ke data_kalimantan/panel_kalimantan.rds & .csv
#   - data_kalimantan/kab_kalimantan.gpkg (ikut diperbarui)
#   - hasil_04/gambut_per_kabupaten.csv, hasil_04/g12_peta_gambut.png
# Cara pakai: buka di RStudio, klik "Source". Waktu: ~1-5 menit.
# =====================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(ggplot2)
})
sf_use_s2(FALSE)
if (basename(getwd()) != "paper_karhutla" && dir.exists("~/paper_karhutla")) {
  setwd("~/paper_karhutla")
}
DIR <- "data_kalimantan"; OUT <- "hasil_04"
dir.create(OUT, showWarnings = FALSE)
message("Folder kerja: ", normalizePath(getwd()))

# ---- 1. Cari file PEATMAP ----------------------------------------------
lokasi <- c(getwd(), DIR, file.path(Sys.getenv("USERPROFILE"), "Downloads"))
lokasi <- unique(lokasi[dir.exists(lokasi)])
f_peat <- unlist(lapply(lokasi, list.files, pattern = "^SEA_Peatland\\.shp$",
                        recursive = TRUE, full.names = TRUE, ignore.case = TRUE))
f_peat <- f_peat[!grepl("\\.xml$", f_peat, ignore.case = TRUE)]
if (!length(f_peat)) {
  stop("File SEA_Peatland.shp tidak ditemukan di Downloads.\n",
       "Unduh 'Asia - PEATMAP shapefile' dari https://archive.researchdata.leeds.ac.uk/251/ ",
       "lalu ekstrak (Extract All) di folder Downloads.")
}
f_peat <- f_peat[1]
message("Peta gambut: ", normalizePath(f_peat))

# ---- 2. Baca kabupaten & potong peta gambut ke Kalimantan --------------
kab <- st_read(file.path(DIR, "kab_kalimantan.gpkg"), quiet = TRUE)
kab <- kab[, setdiff(names(kab), c("prop_gambut", "luas_gambut_km2"))]   # buang hasil lama (jika ada)
kab4326 <- st_transform(kab, 4326)
bb <- st_bbox(kab4326)

peat <- st_read(f_peat, quiet = TRUE)
message("Jumlah poligon PEATMAP (Asia Tenggara): ", nrow(peat))
if (is.na(st_crs(peat))) {
  warning("CRS peta gambut kosong -> diasumsikan WGS84 (EPSG:4326)")
  st_crs(peat) <- 4326
}
peat <- st_transform(peat, 4326)
# Saring dulu poligon yang bersinggungan dengan kotak Kalimantan (lebih cepat)
kotak <- st_as_sfc(bb)
dekat <- tryCatch(lengths(st_intersects(peat, kotak)) > 0,
                  error = function(e) rep(TRUE, nrow(peat)))
peat <- peat[dekat, ]
peat <- suppressWarnings(st_crop(st_make_valid(peat), bb))
peat <- peat[!st_is_empty(peat), ]
message("Poligon gambut di sekitar Kalimantan: ", nrow(peat))
if (!nrow(peat)) stop("Tidak ada poligon gambut di wilayah Kalimantan. Cek file peta.")

# ---- 3. Luas gambut per kabupaten (proyeksi equal-area EPSG:6933) ------
kab_ea  <- st_transform(kab, 6933)
peat_ea <- st_make_valid(st_transform(st_geometry(peat), 6933))
peat_ea <- peat_ea[st_is(peat_ea, c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION"))]

potong <- suppressWarnings(st_intersection(kab_ea[, "id_kab"], peat_ea))
geo <- st_geometry(potong)
# union per kabupaten dulu (hindari hitung ganda jika poligon gambut tumpang tindih)
luas_vec <- tapply(seq_along(geo), potong$id_kab, function(ix) {
  as.numeric(st_area(st_union(st_make_valid(geo[ix]))))
}) / 1e6
luas_g <- data.frame(id_kab = as.integer(names(luas_vec)),
                     luas_gambut_km2 = as.numeric(luas_vec))

kab <- kab |>
  left_join(luas_g, by = "id_kab") |>
  mutate(luas_gambut_km2 = coalesce(luas_gambut_km2, 0),
         prop_gambut     = pmin(1, luas_gambut_km2 / luas_km2))

tab <- st_drop_geometry(kab) |>
  select(id_kab, provinsi, kabupaten, luas_km2, luas_gambut_km2, prop_gambut) |>
  arrange(desc(prop_gambut))
write.csv(tab, file.path(OUT, "gambut_per_kabupaten.csv"), row.names = FALSE)

cat("\n================ GAMBUT PER KABUPATEN ================\n")
cat(sprintf("Total luas gambut Kalimantan (PEATMAP): %.0f km2\n", sum(tab$luas_gambut_km2)))
cat(sprintf("Kabupaten dengan gambut > 0        : %d dari %d\n", sum(tab$prop_gambut > 0), nrow(tab)))
cat(sprintf("Kabupaten dengan gambut >= 20%%     : %d\n", sum(tab$prop_gambut >= 0.2)))
cat("\n10 kabupaten dengan proporsi gambut tertinggi:\n")
for (i in 1:min(10, nrow(tab))) {
  cat(sprintf("  %-25s %5.1f%%  (%6.0f km2)\n", tab$kabupaten[i],
              100 * tab$prop_gambut[i], tab$luas_gambut_km2[i]))
}

# ---- 4. Simpan ke panel & peta kabupaten -------------------------------
panel <- readRDS(file.path(DIR, "panel_kalimantan.rds"))
panel <- panel[, setdiff(names(panel), c("prop_gambut", "luas_gambut_km2"))]
panel <- left_join(panel, tab[, c("id_kab", "prop_gambut")], by = "id_kab")
if (anyNA(panel$prop_gambut)) stop("Ada kabupaten tanpa nilai prop_gambut (tidak seharusnya).")
saveRDS(panel, file.path(DIR, "panel_kalimantan.rds"))
write.csv(panel, file.path(DIR, "panel_kalimantan.csv"), row.names = FALSE)
st_write(kab, file.path(DIR, "kab_kalimantan.gpkg"), delete_dsn = TRUE, quiet = TRUE)

p <- ggplot() +
  geom_sf(data = kab, aes(fill = 100 * prop_gambut), colour = "grey40", linewidth = 0.1) +
  geom_sf(data = peat, fill = NA, colour = "black", linewidth = 0.05, alpha = 0.3) +
  scale_fill_distiller(palette = "YlOrBr", direction = 1, name = "% gambut") +
  labs(title = "Proporsi lahan gambut per kabupaten (PEATMAP, Xu dkk. 2018)") +
  theme_minimal()
ggsave(file.path(OUT, "g12_peta_gambut.png"), p, width = 7, height = 7, dpi = 200)

cat("\nSelesai. prop_gambut sudah ditambahkan ke panel. Hasil di folder:", normalizePath(OUT), "\n")
