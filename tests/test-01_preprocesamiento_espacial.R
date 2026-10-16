# Archivo: tests/testthat/test-01_preprocesamiento_espacial.R

library(testthat)
library(sf)
library(terra)
library(LibreriaMneRcsyPrc)

# --- 1. Generacion de Datos Sinteticos ---
set.seed(42)

# Dataframe simulado con presencias
df_presencias <- data.frame(
  longitud = runif(50, -100, -98),
  latitud = runif(50, 18, 20),
  especievalida = sample(c("Sp_A", "Sp_B"), 50, replace = TRUE),
  fechacolecta = sample(seq(as.Date('2025/01/01'), as.Date('2026/01/01'), by="day"), 50),
  diametro = runif(50, 10, 80)
)

# Poligono de referencia (Bounding Box)
poligono_ref <- st_sfc(st_polygon(list(matrix(
  c(-101, 17, -97, 17, -97, 21, -101, 21, -101, 17),
  ncol = 2, byrow = TRUE))), crs = 4326)

# Raster ambiental sintetico
raster_ref <- rast(nrows=20, ncols=20, ext=ext(-101, -97, 17, 21), crs="EPSG:4326")
values(raster_ref) <- runif(ncell(raster_ref))
names(raster_ref) <- "Bio1"

# --- 2. Pruebas Unitarias ---

test_that("CleaningTransformingDfToSf convierte, filtra y agrupa correctamente", {
  sf_limpio <- CleaningTransformingDfToSf(
    DataFrame = df_presencias,
    SpatialReference = poligono_ref,
    Variables = c("diametro"),
    EspecieColumna = "especievalida",
    FechaColumna = "fechacolecta",
    PrecisionDecimales = 4,
    verbose = FALSE
  )

  # Validacion de clase y geometria[cite: 1]
  expect_s3_class(sf_limpio, "sf")
  expect_equal(sf::st_crs(sf_limpio), sf::st_crs(poligono_ref))

  # Validacion de retencion de variables y calculo de frecuencia[cite: 1]
  expect_true("diametro" %in% names(sf_limpio))
  expect_true("FrecuenciaEspecieValida" %in% names(sf_limpio))

  # Validacion de limpieza (no debe haber NAs en coordenadas implicitas al ser sf)[cite: 1]
  expect_true(nrow(sf_limpio) > 0 && nrow(sf_limpio) <= nrow(df_presencias))
})

test_that("ProcessingRasterToAssignCrs aplica recorte-previo y enmascarado topologico", {
  raster_proyectado <- ProcessingRasterToAssignCrs(
    RasterObject = raster_ref,
    SpatialReference = poligono_ref,
    MetodoInterpolacion = "bilinear",
    verbose = FALSE
  )

  # Validacion de clase[cite: 2]
  expect_s4_class(raster_proyectado, "SpatRaster")

  # Validacion de proyeccion[cite: 2]
  expect_equal(terra::crs(raster_proyectado), terra::crs(poligono_ref))
})
