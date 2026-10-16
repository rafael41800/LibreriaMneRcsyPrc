# Archivo: tests/testthat/test-03_modelado_rf.R

library(testthat)
library(randomForest)
library(terra)
library(LibreriaMneRcsyPrc)

# --- 1. Generacion de Datos Sinteticos ---
set.seed(42)

# Dataframe de modelado (emulando la salida plana y balanceada de la Funcion 4)
df_modelo <- data.frame(
  Presence = sample(c(0, 1), 60, replace = TRUE),
  Bio1 = runif(60, 10, 30),
  Bio2 = runif(60, 50, 100),
  X = runif(60, -100, -98),
  Y = runif(60, 18, 20)
)
lista_datos <- list("Sp_A" = df_modelo)

# Raster sintetico para probar la prediccion espacial
r_env <- rast(nrows=10, ncols=10, ext=ext(0, 10, 0, 10))
values(r_env) <- runif(100)
names(r_env) <- c("Bio1")
r_env$Bio2 <- r_env$Bio1 * runif(100, 1, 2)

# --- 2. Pruebas Unitarias ---

test_that("RfMdeMultiEspecies ajusta el modelo, particiona datos y evalua metricas", {
  res <- RfMdeMultiEspecies(
    ResultadosList = lista_datos,
    Variables = c("Bio1", "Bio2"),
    VariablePresencia = "Presence",
    RegistrosMinimos = 10,
    StackRasterReproyectado = r_env,
    RutaAlmacenamiento = NULL, # Evita escribir en disco durante los tests automatizados
    TuneMTry = FALSE, # Desactivado para acelerar el test
    GuardarDatos = "both",
    verbose = FALSE
  )

  # Validacion de la estructura de salida vectorizada[cite: 5]
  expect_type(res, "list")
  expect_true("Sp_A" %in% names(res))

  res_sp <- res[["Sp_A"]]

  # Validacion del ajuste algoritmico del Random Forest[cite: 5]
  expect_s3_class(res_
