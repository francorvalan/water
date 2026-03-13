library(rstac)
library(sf)
library(dplyr)
library(purrr)
library(curl)

# -------------------------
# 1) Crear items_sf con item_id robusto
# -------------------------
prep_items_sf <- function(items, strict = TRUE) {
  items_sf <- items_as_sf(items)

  # Caso 1: items_as_sf ya trae id
  if ("id" %in% names(items_sf)) {
    items_sf$item_id <- as.character(items_sf$id)
    items_sf$feature_index <- seq_len(nrow(items_sf))
    return(items_sf)
  }

  # Caso 2: ya existe item_id
  if ("item_id" %in% names(items_sf)) {
    if (!"feature_index" %in% names(items_sf)) {
      items_sf$feature_index <- seq_len(nrow(items_sf))
    }
    return(items_sf)
  }

  # Caso 3: inyectar desde items$features (alineación por orden)
  n_feat <- length(items$features)
  n_row  <- nrow(items_sf)

  if (n_row == n_feat) {
    items_sf$item_id <- vapply(items$features, function(f) f$id, character(1))
    items_sf$feature_index <- seq_len(n_feat)
    return(items_sf)
  }

  msg <- paste0(
    "No puedo crear item_id: nrow(items_sf)=", n_row,
    " != length(items$features)=", n_feat,
    ". (items_as_sf no está 1:1 con features)."
  )
  if (strict) stop(msg) else warning(msg)

  items_sf$item_id <- NA_character_
  items_sf$feature_index <- NA_integer_
  items_sf
}

# -------------------------
# 2) Extraer nubosidad y datetime (robusto a nombres)
# -------------------------
get_cloud <- function(items_sf) {
  cand <- c(
    "eo:cloud_cover", "eo.cloud_cover", "cloud_cover", "landsat:cloud_cover",
    "properties.eo:cloud_cover", "properties.cloud_cover", "properties.eo.cloud_cover"
  )

  hit <- cand[cand %in% names(items_sf)]
  if (length(hit) >= 1) return(as.numeric(items_sf[[hit[1]]]))

  # fallback: cualquier columna que contenga cloud_cover
  i <- grep("cloud[_\\.:]?cover", names(items_sf), ignore.case = TRUE)
  if (length(i) == 0) return(rep(NA_real_, nrow(items_sf)))

  as.numeric(items_sf[[i[1]]])
}

get_datetime <- function(items_sf) {
  cand <- c("datetime", "properties.datetime", "properties.datetime_utc")
  hit <- cand[cand %in% names(items_sf)]
  if (length(hit) >= 1) return(as.POSIXct(items_sf[[hit[1]]], tz = "UTC"))

  # fallback: columnas con "datetime"
  i <- grep("datetime", names(items_sf), ignore.case = TRUE)
  if (length(i) == 0) stop("No encuentro campo datetime en items_sf.")
  as.POSIXct(items_sf[[i[1]]], tz = "UTC")
}
# 
# -------------------------
# 3) Rankear escenas: score = nube + overlap + recencia
# -------------------------
rank_scenes <- function(items_sf, roi_wgs84,
                        w_cloud = 0.6, w_overlap = 0.3, w_recent = 0.1) {

  if (!"item_id" %in% names(items_sf)) stop("Falta 'item_id'. Usá prep_items_sf().")

  items_sf  <- st_transform(items_sf, 4326)
  roi_wgs84 <- st_transform(roi_wgs84, 4326)

  # Para evitar problemas con list-columns: solo ID + geometry
  items_min <- items_sf %>% dplyr::select(item_id, geometry)

  suppressWarnings({
    inter <- st_intersection(st_make_valid(items_min), st_make_valid(roi_wgs84))
  })
  if (nrow(inter) == 0) stop("No hay intersección ROI vs escenas (tras bbox).")

  roi_area <- as.numeric(st_area(st_union(roi_wgs84)))
  inter$.__area <- as.numeric(st_area(inter))

  overlap_by_id <- inter %>%
    st_drop_geometry() %>%
    group_by(item_id) %>%
    summarise(overlap_area = sum(.__area, na.rm = TRUE), .groups = "drop")

  cloud <- get_cloud(items_sf)
  dt    <- get_datetime(items_sf)

  df <- items_sf %>%
    mutate(
      cloud = cloud,
      dt    = dt
    ) %>%
    st_drop_geometry() %>%
    left_join(overlap_by_id, by = "item_id") %>%
    mutate(
      overlap_area = ifelse(is.na(overlap_area), 0, overlap_area),
      overlap_frac = pmin(1, pmax(0, overlap_area / roi_area)),
      cloud_frac   = pmin(1, pmax(0, cloud / 100)),
      recent_frac  = {
        dtn <- as.numeric(dt)
        rng <- range(dtn, na.rm = TRUE)
        if (diff(rng) == 0) rep(1, length(dtn)) else (dtn - rng[1]) / diff(rng)
      },
      score = (1 - cloud_frac) * w_cloud + overlap_frac * w_overlap + recent_frac * w_recent
    ) %>%
    arrange(desc(score), cloud, desc(dt))

  df
}

# -------------------------
# 4) Descargar assets del item elegido
# -------------------------
# download_scene_assets <- function(items, scene_id, out_dir,
#                                   assets = c("SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","QA_PIXEL"),
#                                   overwrite = FALSE) {
# 
#   feat <- purrr::keep(items$features, ~ .x$id == scene_id)
#   if (length(feat) != 1) stop("scene_id no único / no encontrado.")
#   feat <- feat[[1]]
# 
#   dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# 
#   available <- names(feat$assets)
#   assets_to_get <- intersect(assets, available)
# 
#   if (length(assets_to_get) == 0) {
#     stop("Ninguno de los assets solicitados está disponible en este item.")
#   }
# 
#   paths <- purrr::map_chr(assets_to_get, function(a) {
#     href <- feat$assets[[a]]$href
#     dest <- file.path(out_dir, paste0(scene_id, "_", a, ".tif"))
#     if (!file.exists(dest) || overwrite) {
#       curl_download(href, destfile = dest, quiet = FALSE)
#     }
#     dest
#   })
# 
#   names(paths) <- assets_to_get
#   paths
# }

# download_landsat_sr_best <- function(roi,
#                                      start_date,
#                                      end_date,
#                                      cloud_max = 20,
#                                      platforms = c("landsat-8","landsat-9"),
#                                      out_dir = "landsat_download") {
# 
#   roi_wgs84 <- st_transform(roi, 4326)
#   bb <- st_bbox(roi_wgs84)
# 
#   stac_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")
# 
#   items <- stac_obj %>%
#     stac_search(
#       collections = "landsat-c2-l2",
#       bbox = c(bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"]),
#       datetime = paste0(start_date, "/", end_date),
#       limit = 500
#     ) %>%
#     post_request() %>%
#     items_fetch()
# 
#   items <- items_sign_planetary_computer(items)
# 
#   items_sf <- prep_items_sf(items, strict = TRUE)
# 
#   # Filtro plataforma (si existe)
#   if ("platform" %in% names(items_sf)) {
#     items_sf <- dplyr::filter(items_sf, platform %in% platforms)
#   }
# 
#   # Filtro nubosidad (si existe)
#   cloud <- get_cloud(items_sf)
#   if (!all(is.na(cloud))) {
#     items_sf <- items_sf %>%
#       mutate(.cloud = cloud) %>%
#       filter(is.na(.cloud) | .cloud <= cloud_max)
#   }
# 
#   if (nrow(items_sf) == 0) stop("0 escenas luego de filtros (fecha/ROI/nube/plataforma).")
# 
#   ranked  <- rank_scenes(items_sf, roi_wgs84)
#   best_id <- ranked$item_id[1]
# 
#   message(
#     "Mejor escena: ", best_id,
#     " | cloud=", ifelse(is.na(ranked$cloud[1]), "NA", round(ranked$cloud[1], 1)),
#     " | overlap=", round(100 * ranked$overlap_frac[1], 1), "%"
#   )
# 
#   scene_dir <- file.path(out_dir, best_id)
#   paths <- download_scene_assets(items, best_id, scene_dir)
# 
#   list(best_id = best_id, ranked = ranked, downloaded = paths)
# }
#------------------------------------------------------------------------------
library(purrr)
library(dplyr)
library(curl)

# --- alias de nombres posibles (PC suele usar common names: red, nir08, swir16, etc.) ---
landsat_asset_aliases <- function() {
  list(
    # Reflectancia (L8/9)
    SR_B2    = c("SR_B2", "sr_b2", "B2", "blue"),
    SR_B3    = c("SR_B3", "sr_b3", "B3", "green"),
    SR_B4    = c("SR_B4", "sr_b4", "B4", "red"),
    SR_B5    = c("SR_B5", "sr_b5", "B5", "nir08", "nir"),
    SR_B6    = c("SR_B6", "sr_b6", "B6", "swir16"),
    SR_B7    = c("SR_B7", "sr_b7", "B7", "swir22"),
    
    # QA
    QA_PIXEL  = c("QA_PIXEL", "qa_pixel", "pixel_qa"),
    QA_RADSAT = c("QA_RADSAT", "qa_radsat"),
    
    # Temperatura (si está expuesto como asset)
    ST_B10   = c("ST_B10", "st_b10", "lwir11", "thermal", "st")
  )
}

# normaliza strings para comparar
.norm_key <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
}

# resuelve 1 "requested" a un key real disponible en el item (case-insensitive, tolera separadores)
resolve_asset_key <- function(available_keys, requested_key, aliases = landsat_asset_aliases()) {
  av0 <- available_keys
  av  <- .norm_key(available_keys)
  
  # 1) match directo
  req_norm <- .norm_key(requested_key)
  hit <- which(av == req_norm)
  if (length(hit) >= 1) return(av0[hit[1]])
  
  # 2) match por alias
  cand <- aliases[[requested_key]]
  if (!is.null(cand)) {
    cand_norm <- .norm_key(cand)
    hit <- which(av %in% cand_norm)
    if (length(hit) >= 1) return(av0[hit[1]])
  }
  
  # 3) match por contiene (fallback)
  hit <- which(grepl(req_norm, av, fixed = TRUE))
  if (length(hit) >= 1) return(av0[hit[1]])
  
  NA_character_
}

download_scene_assets <- function(items, scene_id, out_dir,
                                  assets = c("SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","QA_PIXEL"),
                                  overwrite = FALSE,
                                  verbose = TRUE) {
  
  feat <- purrr::keep(items$features, ~ .x$id == scene_id)
  if (length(feat) != 1) stop("scene_id no único / no encontrado.")
  feat <- feat[[1]]
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  available <- names(feat$assets)
  
  # Resolver assets solicitados a keys reales
  resolved <- vapply(assets, \(a) resolve_asset_key(available, a), character(1))
  ok <- !is.na(resolved)
  
  if (verbose) {
    cat("Assets disponibles (primeros 40):\n",
        paste(head(available, 40), collapse = ", "), "\n\n")
    cat("Solicitados -> resueltos:\n")
    print(data.frame(requested = assets, resolved = resolved))
  }
  
  if (!any(ok)) {
    stop("Ninguno de los assets solicitados está disponible en este item (ni por alias).")
  }
  
  # Descargar
  paths <- purrr::map_chr(which(ok), function(i) {
    req <- assets[i]
    key <- resolved[i]
    href <- feat$assets[[key]]$href
    
    dest <- file.path(out_dir, paste0(scene_id, "_", req, ".tif"))
    if (!file.exists(dest) || overwrite) {
      curl_download(href, destfile = dest, quiet = FALSE)
    }
    dest
  })
  
  names(paths) <- assets[ok]
  paths
}

# download_landsat_sr_best <- function(roi,
#                                      start_date,
#                                      end_date,
#                                      cloud_max = 20,
#                                      platforms = c("landsat-8","landsat-9"),
#                                      out_dir = "landsat_download",
#                                      assets = c("blue","green","red","nir08","swir16","swir22","qa_pixel")) {
#   
#   # ... todo tu código igual hasta best_id ...
#   
#   scene_dir <- file.path(out_dir, best_id)
#   paths <- download_scene_assets(items, best_id, scene_dir, assets = assets)
#   
#   list(best_id = best_id, ranked = ranked, downloaded = paths)
# }

#-------------------------------------------------------------------------------

download_landsat_sr_best <- function(roi,
                                     start_date,
                                     end_date,
                                     cloud_max = 20,
                                     platforms = c("landsat-8","landsat-9"),
                                     out_dir = "landsat_download",
                                     assets = c("SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","QA_PIXEL"),
                                     limit = 500) {
  
  roi_wgs84 <- sf::st_transform(roi, 4326)
  bb <- sf::st_bbox(roi_wgs84)
  
  stac_obj <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1")
  
  items <- stac_obj %>%
    rstac::stac_search(
      collections = "landsat-c2-l2",
      bbox = c(bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"]),
      datetime = paste0(start_date, "/", end_date),
      limit = limit
    ) %>%
    rstac::post_request() %>%
    rstac::items_fetch()
  
  # Firmar URLs (Planetary Computer)
  items <- rstac::items_sign_planetary_computer(items)
  
  # SF + item_id robusto
  items_sf <- prep_items_sf(items, strict = TRUE)
  
  # --- filtros ---
  # Plataforma (si existe)
  if ("platform" %in% names(items_sf)) {
    items_sf <- dplyr::filter(items_sf, .data$platform %in% platforms)
  }
  
  # Nubosidad (si existe)
  cloud <- get_cloud(items_sf)
  if (!all(is.na(cloud))) {
    items_sf <- items_sf %>%
      dplyr::mutate(.cloud = cloud) %>%
      dplyr::filter(is.na(.cloud) | .cloud <= cloud_max)
  }
  
  if (nrow(items_sf) == 0) stop("0 escenas luego de filtros (fecha/ROI/nube/plataforma).")
  
  # --- ranking + best_id ---
  ranked <- rank_scenes(items_sf, roi_wgs84)
  if (nrow(ranked) == 0) stop("No pude rankear escenas (rank_scenes devolvió 0 filas).")
  
  best_id <- ranked$item_id[1]
  
  message(
    "Mejor escena: ", best_id,
    " | cloud=", ifelse(is.na(ranked$cloud[1]), "NA", round(ranked$cloud[1], 1)),
    " | overlap=", round(100 * ranked$overlap_frac[1], 1), "%"
  )
  
  # --- descarga (con alias resolver) ---
  scene_dir <- file.path(out_dir, best_id)
  paths <- download_scene_assets(items, best_id, scene_dir, assets = assets, verbose = TRUE)
  
  list(best_id = best_id, ranked = ranked, downloaded = paths)
}

library(rstac)
library(sf)
library(dplyr)
library(purrr)
library(curl)
library(sf)

WeatherStation_site<- sf::st_point(c( -68.86469,-33.00513)) %>% st_sfc(crs=4326)

roi <- WeatherStation_site %>% st_buffer(500) %>% st_bbox() %>% st_as_sfc()#
#roi <-st_as_sfc(st_bbox(c(xmin=-69.2, ymin=-33.2, xmax=-68.9, ymax=-33.0), crs=4326)) |> st_as_sf()
#res <- download_landsat_sr_best(roi, "2016-03-01", "2016-03-20", cloud_max = 100, out_dir = "./")

res <- download_landsat_sr_best(
  roi,
  start_date = "2016-03-01",
  end_date   = "2016-04-15",
  cloud_max  = 15,
  out_dir    = "Landsat_SR"
)


