#' Modelado de Distribucion de Especies Multiples con Random Forest
#'
#' @description
#' Ejecuta modelos Random Forest para multiples especies a partir de datos de
#' presencia/ausencia. Incluye validacion cruzada, optimizacion de hiperparametros
#' (mtry), evaluacion de metricas y prediccion espacial.
#'
#' @param ResultadosList Lista de objetos `sf` o `data.frame` (salida de la Funcion 4).
#' @param Variables Vector de caracteres con los nombres de las variables predictoras.
#' @param VariablePresencia Nombre de la columna de respuesta. Por defecto `"Presence"`.
#' @param RegistrosMinimos Minimo de datos necesarios para intentar el modelo. Por defecto `10`.
#' @param StackRasterReproyectado Objeto `SpatRaster` para la prediccion espacial.
#' @param RutaAlmacenamiento Carpeta para guardar `.rds`, `.tif` y graficos. Si es `NULL`, no guarda.
#' @param NTree Numero de arboles. Por defecto `100`.
#' @param SetSeed Semilla para reproducibilidad.
#' @param CrossValidationParticiones Numero de folds para k-fold CV. Por defecto `5`.
#' @param TuneMTry Logico. ¿Optimizar `mtry` automaticamente? Por defecto `TRUE`.
#' @param GuardarDatos Opcion para guardar los datos usados en el modelo:
#'        `"none"` (default), `"train"`, `"valid"` o `"both"`.
#' @param verbose Logico. Muestra el progreso del modelado.
#'
#' @return
#' Una lista nombrada por especie que contiene: modelo, importancia, metricas
#' (AUC, Kappa, F1), raster de prediccion (si se proporciono el stack) y,
#' opcionalmente, los conjuntos de entrenamiento y validacion.
#'
#' @details
#' La funcion realiza una particion de datos: 75% para entrenamiento y 25% para validacion.
#' Si `TuneMTry` es `TRUE`, utiliza `tuneRF` para optimizar el parametro `mtry`.
#'
#' @section Advertencias:
#' \itemize{
#'   \item Requiere ambas clases (0 y 1) para funcionar.
#'   \item La calidad del mapa depende de la resolucion del raster.
#' }
#'
#' @seealso
#' \code{\link[randomForest]{randomForest}}, \code{\link[caret]{confusionMatrix}}
#'
#' @keywords random forest sdm modelado
#'
#' @importFrom randomForest randomForest tuneRF varImpPlot
#' @importFrom caret createDataPartition createFolds confusionMatrix
#' @importFrom pROC roc auc
#' @importFrom terra predict writeRaster
#' @importFrom dplyr select %>% all_of
#' @importFrom stats as.formula na.omit predict
#' @importFrom grDevices png dev.off
#' @export
RfMdeMultiEspecies <- function(ResultadosList,
                               Variables,
                               VariablePresencia = "Presence",
                               RegistrosMinimos = 10,
                               StackRasterReproyectado = NULL,
                               RutaAlmacenamiento = NULL,
                               NTree = 100,
                               SetSeed = NULL,
                               CrossValidationParticiones = 5,
                               TuneMTry = TRUE,
                               GuardarDatos = c("none", "train", "valid", "both"),
                               verbose = TRUE) {

  GuardarDatos <- match.arg(GuardarDatos)

  if(!is.null(SetSeed)) set.seed(SetSeed)
  modelos_especies <- list()

  VariablePresencia <- as.character(VariablePresencia)
  Variables <- as.character(Variables)

  if(verbose) message("\n--- [INICIO] Modelado Random Forest Multiespecie (75% train / 25% valid) ---")

  NombresEspecies <- names(ResultadosList)

  for (idx in seq_along(NombresEspecies)) {
    m <- NombresEspecies[idx]

    datos_sp <- ResultadosList[[m]]
    if(inherits(datos_sp, "sf")) {
      datos_sp <- as.data.frame(datos_sp)
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

    DatosParaModelado[[VariablePresencia]] <- as.factor(DatosParaModelado[[VariablePresencia]])

    if(length(unique(DatosParaModelado[[VariablePresencia]])) < 2) {
      if(verbose) message(sprintf("  ! %s: Solo se encontro una clase (Presencia o Ausencia). Saltando...", m))
      next
    }

    if(verbose) message(sprintf("-> Procesando [%d/%d]: %s", idx, length(NombresEspecies), m))

    idx_train <- caret::createDataPartition(DatosParaModelado[[VariablePresencia]], p = 0.75, list = FALSE)
    DatosTrain <- DatosParaModelado[idx_train, ]
    DatosValid <- DatosParaModelado[-idx_train, ]

    if(length(unique(DatosValid[[VariablePresencia]])) < 2) {
      if(verbose) message(sprintf("  ! %s: Conjunto de validacion con una sola clase. Saltando...", m))
      next
    }

    mtry_opt <- floor(sqrt(length(Variables)))
    if(TuneMTry) {
      try({
        tuner <- randomForest::tuneRF(x = DatosTrain[, Variables],
                                      y = DatosTrain[[VariablePresencia]],
                                      ntreeTry = 50, stepFactor = 1.5, trace = FALSE, plot = FALSE)
        mtry_opt <- tuner[which.min(tuner[,2]), 1]
      }, silent = TRUE)
    }

    cv_resumen <- NULL
    if(CrossValidationParticiones > 1 && nrow(DatosTrain) > (CrossValidationParticiones * 2)) {
      folds <- caret::createFolds(DatosTrain[[VariablePresencia]], k = CrossValidationParticiones)
      cv_metrics <- lapply(folds, function(idx_f){
        m_cv <- randomForest::randomForest(x = DatosTrain[-idx_f, Variables],
                                           y = DatosTrain[-idx_f, VariablePresencia],
                                           ntree = NTree, mtry = mtry_opt)
        pred_cv <- stats::predict(m_cv, DatosTrain[idx_f, ])
        cm <- caret::confusionMatrix(pred_cv, DatosTrain[[VariablePresencia]][idx_f])
        return(c(acc = cm$overall["Accuracy"], kappa = cm$overall["Kappa"]))
      })
      cv_resumen <- as.data.frame(do.call(rbind, cv_metrics))
    }

    formula_rf <- stats::as.formula(paste(VariablePresencia, "~ ."))
    rf_final <- randomForest::randomForest(formula_rf, data = DatosTrain,
                                           importance = TRUE, ntree = NTree, mtry = mtry_opt)

    pred_valid <- stats::predict(rf_final, DatosValid)
    cm_valid <- caret::confusionMatrix(pred_valid, DatosValid[[VariablePresencia]])

    if(length(unique(DatosValid[[VariablePresencia]])) > 1) {
      prob_valid <- stats::predict(rf_final, DatosValid, type = "prob")[, "1"]
      roc_obj <- pROC::roc(DatosValid[[VariablePresencia]], prob_valid, quiet = TRUE)
      auc_val <- as.numeric(pROC::auc(roc_obj))
      f1_score <- 2 * (cm_valid$byClass["Sensitivity"] * cm_valid$byClass["Pos Pred Value"]) /
        (cm_valid$byClass["Sensitivity"] + cm_valid$byClass["Pos Pred Value"])
    } else {
      auc_val <- NA
      f1_score <- NA
    }

    rast_pred <- NULL
    if(!is.null(StackRasterReproyectado)) {
      rast_pred <- terra::predict(StackRasterReproyectado, rf_final, type = "prob", na.rm = TRUE)[[2]]
      names(rast_pred) <- paste0("prob_", m)
    }

    if(!is.null(RutaAlmacenamiento)) {
      dir.create(RutaAlmacenamiento, showWarnings = FALSE, recursive = TRUE)

      saveRDS(rf_final, file.path(RutaAlmacenamiento, paste0("RF_Model_", m, ".rds")))

      if(GuardarDatos %in% c("train", "both")) {
        saveRDS(DatosTrain, file.path(RutaAlmacenamiento, paste0("Train_Data_", m, ".rds")))
      }
      if(GuardarDatos %in% c("valid", "both")) {
        saveRDS(DatosValid, file.path(RutaAlmacenamiento, paste0("Valid_Data_", m, ".rds")))
      }

      grDevices::png(file.path(RutaAlmacenamiento, paste0("Imp_Plot_", m, ".png")), width = 800, height = 600)
      randomForest::varImpPlot(rf_final, main = paste("Importancia de Variables -", m))
      grDevices::dev.off()

      if(!is.null(rast_pred)) {
        terra::writeRaster(rast_pred, file.path(RutaAlmacenamiento, paste0("Map_", m, ".tif")), overwrite = TRUE)
      }
    }

    modelos_especies[[m]] <- list(
      modelo = rf_final,
      importance = rf_final$importance,
      cv_summary = cv_resumen,
      accuracy_valid = cm_valid$overall["Accuracy"],
      kappa_valid = cm_valid$overall["Kappa"],
      auc = auc_val,
      f1_score = f1_score,
      raster_prediction = rast_pred,
      train_data = if(GuardarDatos %in% c("train", "both")) DatosTrain else NULL,
      valid_data = if(GuardarDatos %in% c("valid", "both")) DatosValid else NULL
    )
  }

  if(verbose) message("\n--- [FIN] Modelado completado para todas las especies ---")
  return(modelos_especies)
}
