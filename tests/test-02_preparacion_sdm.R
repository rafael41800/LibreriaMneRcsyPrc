# Archivo: tests/testthat/test-02_preparacion_sdm.R

library(testthat)
library(sf)
library(terra)
library(LibreriaMneRcsyPrc)

# --- 1. Generacion de Datos Sinteticos ---
set.seed(42)

# Dataframe de puntos con presencias
df_puntos <- data.frame(
  longitud = c(-99.1, -99.2, -99.3, -99.4, -99.5),
  latitud = c(19.1, 19.2, 19.3, 19.4, 19.5),
  especievalida = c("Sp_A", "Sp_A", "Sp_B", "Sp_B", "Sp_B")
)
puntos_sf <- st_as_sf(df_puntos, coords = c("longitud", "latitud"), crs = 4326)

# Raster ambiental continuo sintetico
r_env <- rast(nrows=20, ncols=20, ext=ext(-100, -98, 18, 20), crs="EPSG:4326")
values(r_env) <- runif(ncell(r_env))
names(r_env) <- "Bio1"

# --- 2. Pruebas Unitarias ---

test_that("PresencesSfNoNaListDf extrae variables vectorizadamente y genera lista de data.frames", {
  lista_limpia <- PresencesSfNoNaListDf(
    DataSf = puntos_sf,
    RasterReferencia = r_env,
    EspecieValida = "especievalida",
    verbose = FALSE
  )

  # Validacion estructural de la lista y particion[cite: 3]
  expect_type(lista_limpia, "list")
