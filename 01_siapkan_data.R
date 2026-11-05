# =====================================================================
# 01_siapkan_data.R
# Ide 1: DLNM spasio-temporal untuk titik panas Kalimantan
#
# Output: panel kabupaten x bulan (2012-2025) berisi
#   - hotspot   : jumlah titik panas VIIRS S-NPP (NASA FIRMS)
#   - hujan_mm  : curah hujan bulanan rata-rata kabupaten (CHIRPS v3.0)
#   - hujan_anom: anomali hujan terstandar per kabupaten-per bulan kalender
#   - nino34, dmi : indeks ENSO dan IOD (NOAA PSL)
#   - luas_km2  : luas kabupaten (untuk offset)
#   - prop_gambut (opsional, jika file peta gambut tersedia)
# + file graph ketetanggaan untuk BYM2 di R-INLA
#
# Cara pakai: buka di RStudio, lalu klik "Source".
# Sebelum menjalankan: unduh file CHIRPS "data.nc" lewat browser (link di bagian 3).
# Kalau terputus, jalankan ulang saja: file yang sudah ada tidak diunduh lagi.
# =====================================================================

# ---- 0. Paket ---------------------------------------------------------
pkgs <- c("sf", "terra", "geodata", "dplyr", "tidyr", "readr",
          "lubridate", "spdep", "ggplot2", "ncdf4")
baru <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(baru)) install.packages(baru)
# geodata & ncdf4 dipanggil dengan "::" saja (tidak di-load)
invisible(lapply(setdiff(pkgs, c("geodata", "ncdf4")), library, character.only = TRUE))

# Pastikan folder kerja = paper_karhutla (setelah RStudio restart, folder kerja kembali ke Documents)
if (basename(getwd()) != "paper_karhutla" && dir.exists("~/paper_karhutla")) {
  setwd("~/paper_karhutla")
}
message("Folder kerja: ", normalizePath(getwd()))

sf_use_s2(FALSE)
options(timeout = 3600)

TAHUN <- 2012:2025            # VIIRS S-NPP tersedia sejak 2012
DIR   <- "data_kalimantan"
dir.create(file.path(DIR, "firms"),  recursive = TRUE, showWarnings = FALSE)

# ---- 1. Batas kabupaten/kota Kalimantan (GADM 4.1) -------------------
adm <- geodata::gadm("IDN", level = 2, path = DIR)
if (is.null(adm)) stop("Gagal mengunduh batas wilayah GADM. Cek internet lalu jalankan ulang.")
kab <- st_as_sf(adm) |>
  filter(grepl("Kalimantan", NAME_1)) |>
  st_make_valid() |>
  transmute(provinsi = NAME_1, kabupaten = NAME_2, GID_2) |>
  mutate(id_kab = row_number(), .before = 1)

kab$luas_km2 <- as.numeric(st_area(st_transform(kab, 6933))) / 1e6
message("Jumlah kabupaten/kota Kalimantan: ", nrow(kab))
print(table(kab$provinsi))

bb <- st_bbox(kab)

# ---- 2. Titik panas VIIRS S-NPP (NASA FIRMS) ---------------------------
# Script mencoba mengunduh file tahunan per negara secara otomatis.
# JIKA GAGAL, unduh manual dari https://firms.modaps.eosdis.nasa.gov/download/ :
#   (a) klik "Download yearly summary by country [.csv]", pilih VIIRS S-NPP,
#       ambil file Indonesia untuk tiap tahun (ekstrak dulu kalau berupa .zip); ATAU
#   (b) "Create New Request": pilih Country = Indonesia, sumber VIIRS S-NPP,
#       rentang tanggal 2012-01-01 s.d. 2025-12-31, format CSV (link dikirim via email).
# Simpan semua file .csv ke folder data_kalimantan/firms/ (nama file bebas).
unduh_firms <- function(th) {
  f <- file.path(DIR, "firms", sprintf("viirs-snpp_%d_Indonesia.csv", th))
  if (file.exists(f)) return(invisible(TRUE))
  url <- sprintf(paste0("https://firms.modaps.eosdis.nasa.gov/data/country/",
                        "viirs-snpp/%d/viirs-snpp_%d_Indonesia.csv"), th, th)
  ok <- tryCatch({
    download.file(url, f, mode = "wb", quiet = TRUE); TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  # Server kadang mengembalikan halaman HTML/error, bukan CSV -> anggap gagal
  if (ok) ok <- grepl("latitude", readLines(f, n = 1, warn = FALSE)[1])
  if (!isTRUE(ok)) {
    if (file.exists(f)) file.remove(f)
    message("  ! Gagal unduh otomatis FIRMS ", th,
            " (abaikan jika file tahun ini sudah Anda unduh manual)")
  }
  invisible(ok)
}
message("Mengunduh FIRMS ...")
invisible(lapply(TAHUN, unduh_firms))

f_firms <- list.files(file.path(DIR, "firms"), pattern = "\\.csv$", full.names = TRUE)
csv_valid <- vapply(f_firms, function(f) {
  h <- readLines(f, n = 1, warn = FALSE)
  length(h) == 1 && all(c("latitude", "longitude", "acq_date") %in% strsplit(gsub('"', "", h), ",")[[1]])
}, logical(1))
if (any(!csv_valid)) message("  ! File diabaikan (bukan CSV FIRMS): ",
                             paste(basename(f_firms[!csv_valid]), collapse = ", "))
f_firms <- f_firms[csv_valid]
if (!length(f_firms)) stop("Tidak ada file CSV FIRMS yang valid di ", file.path(DIR, "firms"),
                           ". Lihat petunjuk unduh manual di bagian 2 script.")

hs <- lapply(f_firms, readr::read_csv, show_col_types = FALSE,
             col_types = readr::cols(.default = readr::col_character())) |>
  bind_rows() |>
  mutate(latitude  = as.numeric(latitude),
         longitude = as.numeric(longitude),
         acq_date  = as.Date(acq_date)) |>
  filter(between(longitude, bb[["xmin"]], bb[["xmax"]]),
         between(latitude,  bb[["ymin"]], bb[["ymax"]]),
         year(acq_date) %in% TAHUN)

# Pastikan SEMUA tahun ada datanya. Tanpa cek ini, tahun yang gagal diunduh
# akan terisi nol semua dan merusak analisis.
th_kosong <- setdiff(TAHUN, unique(year(hs$acq_date)))
if (length(th_kosong)) {
  stop("Data FIRMS belum ada untuk tahun: ", paste(th_kosong, collapse = ", "),
       ".\nUnduh manual CSV Indonesia (VIIRS S-NPP) tahun tsb dari ",
       "https://firms.modaps.eosdis.nasa.gov/download/ (lihat petunjuk di bagian 2 script),",
       "\nsimpan di folder ", file.path(DIR, "firms"), ", lalu jalankan ulang script ini.")
}

# Buang deteksi ganda (mis. file otomatis & manual tumpang tindih)
kunci <- intersect(c("latitude", "longitude", "acq_date", "acq_time"), names(hs))
hs <- distinct(hs, across(all_of(kunci)), .keep_all = TRUE)

# Buang deteksi kepercayaan rendah (VIIRS: l/n/h; MODIS: 0-100)
if ("confidence" %in% names(hs)) {
  conf_num <- suppressWarnings(as.numeric(hs$confidence))
  hs <- hs[ifelse(is.na(conf_num),
                  !tolower(hs$confidence) %in% c("l", "low"),
                  conf_num >= 30), ]
}
# Hanya kebakaran vegetasi (type 0); buang sumber panas statis/offshore
if ("type" %in% names(hs)) hs <- filter(hs, is.na(type) | type == "0")

hs_kab <- st_as_sf(hs, coords = c("longitude", "latitude"), crs = 4326) |>
  st_join(st_transform(kab["id_kab"], 4326), join = st_within) |>
  st_drop_geometry() |>
  filter(!is.na(id_kab)) |>
  mutate(bulan = floor_date(acq_date, "month")) |>
  count(id_kab, bulan, name = "hotspot")

message("Total titik panas Kalimantan: ", sum(hs_kab$hotspot))
rm(hs); gc()

# ---- 3. Curah hujan bulanan CHIRPS v3.0 (1 file NetCDF dari server IRI) --
# Unduh SEKALI lewat browser (Chrome) link di bawah ini, hasilnya file "data.nc"
# (~35 MB). Biarkan di folder Downloads, atau pindahkan ke folder paper_karhutla.
# https://iridl.ldeo.columbia.edu/SOURCES/.UCSB/.CHIRPS/.v3p0/.monthly/.global/.precipitation/X/108.5/119.5/RANGEEDGES/Y/-5/7.5/RANGEEDGES/T/(Jan%202012)/(Dec%202025)/RANGEEDGES/data.nc
kandidat_nc <- c("data.nc", file.path(DIR, "data.nc"),
                 file.path(Sys.getenv("USERPROFILE"), "Downloads", "data.nc"))
f_nc <- kandidat_nc[file.exists(kandidat_nc)][1]
if (is.na(f_nc)) {
  stop("File CHIRPS 'data.nc' tidak ditemukan.\n",
       "Unduh lewat Chrome dari link di bagian 3 script ini, lalu taruh di folder ",
       normalizePath(getwd()), " (atau biarkan di Downloads), kemudian jalankan ulang.")
}
message("Membaca CHIRPS dari: ", normalizePath(f_nc))

# Sumbu waktu dibaca langsung dari NetCDF (IRI memakai "months since 1960-01-01")
nc <- ncdf4::nc_open(f_nc)
nm_t <- intersect(c("T", "time"), names(nc$dim))[1]
if (is.na(nm_t)) { ncdf4::nc_close(nc); stop("Sumbu waktu tidak ditemukan di ", f_nc) }
tv <- nc$dim[[nm_t]]$vals
tu <- nc$dim[[nm_t]]$units
ncdf4::nc_close(nc)

asal <- as.Date(sub(".*since\\s+([0-9]{4}-[0-9]{1,2}-[0-9]{1,2}).*", "\\1", tu))
bulan_nc <- if (grepl("^months", tu)) {
  asal %m+% months(as.integer(floor(tv)))
} else if (grepl("^days", tu)) {
  floor_date(asal + floor(tv), "month")
} else stop("Satuan waktu NetCDF tidak dikenali: ", tu)

r_nc <- terra::rast(f_nc)
if (terra::nlyr(r_nc) != length(bulan_nc)) {
  stop("Jumlah layer (", terra::nlyr(r_nc), ") tidak sama dengan jumlah waktu (",
       length(bulan_nc), ") di file NetCDF.")
}
if (terra::crs(r_nc) == "") terra::crs(r_nc) <- "EPSG:4326"
r_nc <- terra::classify(r_nc, cbind(-Inf, 0, NA))        # nilai < 0 = missing
names(r_nc) <- format(bulan_nc)
r_nc <- r_nc[[which(year(bulan_nc) %in% TAHUN)]]
bulan_nc <- bulan_nc[year(bulan_nc) %in% TAHUN]

bulan_harap <- seq(as.Date(sprintf("%d-01-01", min(TAHUN))),
                   as.Date(sprintf("%d-12-01", max(TAHUN))), by = "month")
kurang <- setdiff(format(bulan_harap), format(bulan_nc))
if (length(kurang)) {
  message("  Perhatian: CHIRPS belum tersedia untuk ", length(kurang), " bulan: ",
          paste(substr(kurang, 1, 7), collapse = ", "),
          " -> hujan_mm bulan tsb akan NA.")
}
message("  CHIRPS: ", length(bulan_nc), " bulan (", format(min(bulan_nc), "%Y-%m"),
        " s.d. ", format(max(bulan_nc), "%Y-%m"), ")")

ch <- terra::extract(r_nc, terra::vect(kab), fun = mean, na.rm = TRUE,
                     touches = TRUE, ID = TRUE)
names(ch)[1] <- "id_kab"

hujan <- ch |>
  pivot_longer(-id_kab, names_to = "bulan", values_to = "hujan_mm") |>
  mutate(bulan = as.Date(bulan),
         hujan_mm = ifelse(is.nan(hujan_mm), NA_real_, hujan_mm))

# ---- 4. Indeks ENSO (Nino 3.4) dan IOD (DMI) dari NOAA PSL -------------
baca_psl <- function(url, nama) {
  x <- tryCatch(readLines(url, warn = FALSE), error = function(e) NULL)
  if (is.null(x)) { message("  ! Gagal unduh ", nama); return(NULL) }
  tok <- strsplit(trimws(x), "\\s+")
  ok  <- vapply(tok, function(t) length(t) == 13 && grepl("^\\d{4}$", t[1]), logical(1))
  m   <- do.call(rbind, lapply(tok[ok], as.numeric))
  df  <- data.frame(tahun = rep(m[, 1], each = 12),
                    bl    = rep(1:12, nrow(m)),
                    nilai = as.vector(t(m[, -1])))
  df$nilai[df$nilai < -90] <- NA                       # kode missing
  df |>
    transmute(bulan = as.Date(sprintf("%d-%02d-01", tahun, bl)),
              !!nama := nilai)
}
enso <- baca_psl("https://psl.noaa.gov/data/correlation/nina34.anom.data", "nino34")
iod  <- baca_psl("https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/dmi.had.long.data", "dmi")

# ---- 5. (Opsional) proporsi lahan gambut -------------------------------
# Taruh peta gambut (shapefile/gpkg) di data_kalimantan/gambut.gpkg
# (mis. peta gambut BBSDLP Kementan atau layer "Indonesia peat lands" Global Forest Watch)
f_gambut <- list.files(DIR, pattern = "^gambut\\.(gpkg|shp)$", full.names = TRUE)
if (length(f_gambut)) {
  g <- st_read(f_gambut[1], quiet = TRUE) |>
    st_transform(6933) |> st_make_valid() |> st_union()
  kab_ea <- st_transform(kab, 6933)
  luas_g <- st_intersection(kab_ea, g) |>
    mutate(luas_gambut_km2 = as.numeric(st_area(geometry)) / 1e6) |>
    st_drop_geometry() |>
    select(id_kab, luas_gambut_km2)
  kab <- kab |>
    left_join(luas_g, by = "id_kab") |>
    mutate(prop_gambut = coalesce(luas_gambut_km2 / luas_km2, 0)) |>
    select(-luas_gambut_km2)
  message("Proporsi gambut ditambahkan.")
} else {
  message("Peta gambut belum ada -> prop_gambut dilewati (bisa ditambah nanti).")
}

# ---- 6. Susun panel kabupaten x bulan ----------------------------------
bulan_seq <- bulan_harap

panel <- expand_grid(id_kab = kab$id_kab, bulan = bulan_seq) |>
  left_join(hs_kab, by = c("id_kab", "bulan")) |>
  mutate(hotspot = replace_na(hotspot, 0L)) |>
  left_join(hujan, by = c("id_kab", "bulan")) |>
  left_join(st_drop_geometry(kab), by = "id_kab")
if (!is.null(enso)) panel <- left_join(panel, enso, by = "bulan")
if (!is.null(iod))  panel <- left_join(panel, iod,  by = "bulan")

panel <- panel |>
  mutate(bl = month(bulan), tahun = year(bulan)) |>
  group_by(id_kab, bl) |>
  mutate(hujan_anom = (hujan_mm - mean(hujan_mm, na.rm = TRUE)) /
                       sd(hujan_mm, na.rm = TRUE)) |>
  ungroup() |>
  arrange(id_kab, bulan) |>
  mutate(t = as.integer(factor(bulan)))

# ---- 7. Graph ketetanggaan untuk BYM2 (R-INLA) -------------------------
nb <- spdep::poly2nb(kab, queen = TRUE)
n_pulau <- sum(spdep::card(nb) == 0)
if (n_pulau > 0) message("Perhatian: ", n_pulau,
                         " kabupaten tanpa tetangga (pulau). Perlu ditangani di model.")
spdep::nb2INLA(file.path(DIR, "kalimantan.graph"), nb)

# ---- 8. Simpan ---------------------------------------------------------
saveRDS(panel, file.path(DIR, "panel_kalimantan.rds"))
readr::write_csv(panel, file.path(DIR, "panel_kalimantan.csv"))
st_write(kab, file.path(DIR, "kab_kalimantan.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# ---- 9. Ringkasan deskriptif cepat -------------------------------------
cat("\n================ RINGKASAN ================\n")
cat("Kabupaten          :", n_distinct(panel$id_kab), "\n")
cat("Periode            :", format(min(panel$bulan)), "s.d.", format(max(panel$bulan)), "\n")
cat("Jumlah observasi   :", nrow(panel), "\n")
cat("Proporsi nol       :", round(mean(panel$hotspot == 0), 3), "\n")
cat("Rata-rata / varians:", round(mean(panel$hotspot), 2), "/",
    round(var(panel$hotspot), 2), "(overdispersi jika varians >> rata-rata)\n")
cat("Missing hujan      :", sum(is.na(panel$hujan_mm)), "\n")

p1 <- panel |>
  group_by(bulan) |>
  summarise(hotspot = sum(hotspot), hujan = mean(hujan_mm, na.rm = TRUE)) |>
  ggplot(aes(bulan, hotspot)) +
  geom_col(fill = "firebrick") +
  labs(x = NULL, y = "Jumlah titik panas", title = "Titik panas bulanan Kalimantan (VIIRS)") +
  theme_minimal()

p2 <- kab |>
  left_join(panel |> group_by(id_kab) |> summarise(total = sum(hotspot)), by = "id_kab") |>
  ggplot() +
  geom_sf(aes(fill = log1p(total)), colour = "grey40", linewidth = 0.1) +
  scale_fill_viridis_c(option = "inferno", name = "log(1+total)") +
  labs(title = "Total titik panas per kabupaten") +
  theme_minimal()

ggsave(file.path(DIR, "gambar_deret_waktu.png"), p1, width = 9, height = 4, dpi = 200)
ggsave(file.path(DIR, "gambar_peta.png"),        p2, width = 7, height = 7, dpi = 200)
print(p1); print(p2)

message("Selesai. Hasil tersimpan di folder: ", normalizePath(DIR))
