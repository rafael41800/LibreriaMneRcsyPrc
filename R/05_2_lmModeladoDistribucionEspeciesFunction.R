#' Modelado de Distribucion de Especies Multiples con Modelos Lineales (LM)
#'
#' @description
#' Ejecuta modelos lineales (regresion lineal de probabilidad) para multiples especies
#' a partir de datos de presencia/ausencia. Incluye particion de datos, evaluacion de
#' metricas (R2, RMSE, AUC), diagnostico de residuales y prediccion espacial.
#'
#' @param ResultadosList Lista de objetos `sf` o `data.frame` (salida de la Funcion 4).
#' @param Variables Vector de caracteres con los nombres de las variables predictoras.
#' @param VariablePresencia Nombre de la columna de respuesta. Por defecto `"Presence"`.
#' @param RegistrosMinimos Minimo de datos necesarios para intentar el modelo. Por defecto `10`.
#' @param StackRasterReproyectado Objeto `SpatRaster` para la prediccion espacial.
#' @param RutaAlmacenamiento Carpeta para guardar `.rds`, `.tif` y graficos. Si es `NULL`, no guarda.
#' @param SetSeed Semilla para reproducibilidad.
#' @param GuardarDatos Opcion para guardar los datos usados en el modelo:
#'        `"none"` (default), `"train"`, `"valid"` o `"both"`.
#' @param verbose Logico. Muestra el progreso del modelado.
#'
#' @return
#' Una lista nombrada por especie que contiene: modelo, coeficientes, metricas
#' de ajuste y validacion (R2, RMSE, AUC), raster de prediccion y los conjuntos de datos.
#'
#' @importFrom caret createDataPartition
#' @importFrom pROC roc auc
#' @importFrom terra predict writeRaster
#' @importFrom dplyr select %>% all_of
#' @importFrom stats as.formula na.omit lm predict
#' @importFrom grDevices png dev.off
#' @export
LmMdeMultiEspecies <- function(ResultadosList,
                               Variables,
                               VariablePresencia = "Presence",
                               RegistrosMinimos = 10,
                               StackRasterReproyectado = NULL,
                               RutaAlmacenamiento = NULL,
                               SetSeed = NULL,
                               GuardarDatos = c("none", "train", "valid", "both"),
                               verbose = TRUE) {

  GuardarDatos <- match.arg(GuardarDatos)

  if(!is.null(SetSeed)) set.seed(SetSeed)
  modelos_especies <- list()

  VariablePresencia <- as.character(VariablePresencia)
  Variables <- as.character(Variables)

  if(verbose) message("\n--- [INICIO] Modelado Lineal (LM) Multiespecie (75% train / 25% valid) ---")

  NombresEspecies <- names(ResultadosList)

  for (idx in seq_along(NombresEspecies)) {
    m <- NombresEspecies[idx]

    datos_sp <- ResultadosList[[m]]
    if(inherits(datos_sp, "sf")) {
      datos_sp <- as.data.frame(sf::st_drop_geometry(datos_sp))
    } else if(!inherits(datos_sp, "data.frame")) {
      datos_sp <- as.data.frame(datos_sp)
    }

    DatosParaModelado <- datos_sp %>%
      dplyr::select(dplyr::all_of(VariablePresencia), dplyr::all_of(Variables)) %>%
      stats::na.omit()

    n_registros <- nrow(DatosParaModelado)
    if(n_registros < RegistrosMinimos) {
      if(verbose) message(sprintf("  ! %s: Registros insuficientes (%d). Saltando...", m, n_registros))
      next
    }

    # Para LM la variable respuesta debe mantenerse numerica (0 y 1)
    DatosParaModelado[[VariablePresencia]] <- as.numeric(as.character(DatosParaModelado[[VariablePresencia]]))

    if(length(unique(DatosParaModelado[[VariablePresencia]])) < 2) {
      if(verbose) message(sprintf("  ! %s: Solo se encontro una clase (Presencia o Ausencia). Saltando...", m))
      next
    }

    if(verbose) message(sprintf("-> Procesando [%d/%d]: %s", idx, length(NombresEspecies), m))

    # Particion estratificada basada en la variable respuesta original
    idx_train <- caret::createDataPartition(as.factor(DatosParaModelado[[VariablePresencia]]), p = 0.75, list = FALSE)
    DatosTrain <- DatosParaModelado[idx_train, ]
    DatosValid <- DatosParaModelado[-idx_train, ]

    if(length(unique(DatosValid[[VariablePresencia]])) < 2) {
      if(verbose) message(sprintf("  ! %s: Conjunto de validacion con una sola clase. Saltando...", m))
      next
    }

    # Ajuste del Modelo Lineal
    formula_lm <- stats::as.formula(paste(VariablePresencia, "~ ."))
    lm_final <- stats::lm(formula_lm, data = DatosTrain)
    summary_lm <- summary(lm_final)

    # Predicciones y Metricas en Validacion
    pred_valid <- stats::predict(lm_final, DatosValid)

    # Metricas clasicas de regresion
    rmse_val <- sqrt(mean((DatosValid[[VariablePresencia]] - pred_valid)^2))
    r2_val <- stats::cor(DatosValid[[VariablePresencia]], pred_valid)^2

    # Metricas de discriminacion (Truncando predicciones al rango [0,1] para pROC)
    pred_valid_trunc <- pmax(0, pmin(1, pred_valid))
    roc_obj <- pROC::roc(DatosValid[[VariablePresencia]], pred_valid_trunc, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))

    # Prediccion Espacial Raster
    rast_pred <- NULL
    if(!is.null(StackRasterReproyectado)) {
      # terra::predict mapea de forma nativa objetos lm
      rast_pred <- terra::predict(StackRasterReproyectado, lm_final, na.rm = TRUE)
      # Truncar el mapa resultante para evitar probabilidades negativas o mayores a 1
      rast_pred <- terra::clamp(rast_pred, lower = 0, upper = 1)
      names(rast_pred) <- paste0("pred_lm_", m)
    }

    # Exportacion de archivos (Corregido: iniciando archivos con el nombre de la especie)
    if(!is.null(RutaAlmacenamiento)) {
      dir.create(RutaAlmacenamiento, showWarnings = FALSE, recursive = TRUE)

      saveRDS(lm_final, file.path(RutaAlmacenamiento, paste0(m, "_LM_Model.rds")))

      if(GuardarDatos %in% c("train", "both")) {
        saveRDS(DatosTrain, file.path(RutaAlmacenamiento, paste0(m, "_Train_Data.rds")))
      }
      if(GuardarDatos %in% c("valid", "both")) {
        saveRDS(DatosValid, file.path(RutaAlmacenamiento, paste0(m, "_Valid_Data.rds")))
      }

      # Grafico de diagnostico basico (Residuales vs Ajustados)
      grDevices::png(file.path(RutaAlmacenamiento, paste0(m, "_Diagnostic_Plot.png")), width = 800, height = 600)
      plot(lm_final, which = 1, main = paste("Residuales vs Ajustados -", m))
      grDevices::dev.off()

      if(!is.null(rast_pred)) {
        terra::writeRaster(rast_pred, file.path(RutaAlmacenamiento, paste0(m, "_Map.tif")), overwrite = TRUE)
      }
    }

    # Estructura del Output
    modelos_especies[[m]] <- list(
      modelo = lm_final,
      coefficients = stats::coefficients(lm_final),
      r2_train = summary_lm$r.squared,
      r2_valid = r2_val,
      rmse_valid = rmse_val,
      auc = auc_val,
      raster_prediction = rast_pred,
      train_data = if(GuardarDatos %in% c("train", "both")) DatosTrain else NULL,
      valid_data = if(GuardarDatos %in% c("valid", "both")) DatosValid else NULL
    )
  }

  if(verbose) message("\n--- [FIN] Modelado LM completado para todas las especies ---")
  return(modelos_especies)
}
