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
#' @param verbose Logico. Muestra el progreso del modelado.
#'
#' @return
#' Una lista nombrada por especie que contiene: modelo, importancia, metricas
#' (AUC, Kappa, F1) y el raster de prediccion (si se proporciono el stack).
#'
#' @details
#' La funcion realiza una particion de datos: 10% para validacion externa,
#' y del resto, un 77% para entrenamiento y 23% para prueba. Si `TuneMTry` es `TRUE`,
#' utiliza `tuneRF` para optimizar el parametro `mtry`.
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
RfMdeMultiEspeciesRobusta <- function(ResultadosList,
                               Variables,
                               VariablePresencia = "Presence",
                               RegistrosMinimos = 10,
                               StackRasterReproyectado = NULL,
                               RutaAlmacenamiento = NULL,
                               NTree = 100,
                               SetSeed = NULL,
                               CrossValidationParticiones = 5,
                               TuneMTry = TRUE,
                               verbose = TRUE) {

  if(!is.null(SetSeed)) set.seed(SetSeed)
  modelos_especies <- list()

  VariablePresencia <- as.character(VariablePresencia)
  Variables <- as.character(Variables)

  if(verbose) message("\n--- [INICIO] Modelado Random Forest Multiespecie ---")

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

    idx_val <- caret::createDataPartition(DatosParaModelado[[VariablePresencia]], p = 0.1, list = FALSE)
    DatosVal <- DatosParaModelado[idx_val, ]
    Resto <- DatosParaModelado[-idx_val, ]

    idx_train <- caret::createDataPartition(Resto[[VariablePresencia]], p = 0.77, list = FALSE)
    DatosTrain <- Resto[idx_train, ]
    DatosTest <- Resto[-idx_train, ]

    if(length(unique(DatosTest[[VariablePresencia]])) < 2) {
      if(verbose) message(sprintf("  ! %s: Conjunto de prueba con una sola clase. Saltando...", m))
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

    pred_test <- stats::predict(rf_final, DatosTest)
    cm_test <- caret::confusionMatrix(pred_test, DatosTest[[VariablePresencia]])

    if(length(unique(DatosTest[[VariablePresencia]])) > 1) {
      prob_test <- stats::predict(rf_final, DatosTest, type = "prob")[, "1"]
      roc_obj <- pROC::roc(DatosTest[[VariablePresencia]], prob_test, quiet = TRUE)
      auc_val <- as.numeric(pROC::auc(roc_obj))
      f1_score <- 2 * (cm_test$byClass["Sensitivity"] * cm_test$byClass["Pos Pred Value"]) /
        (cm_test$byClass["Sensitivity"] + cm_test$byClass["Pos Pred Value"])
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
      saveRDS(DatosVal, file.path(RutaAlmacenamiento, paste0("Val_Data_", m, ".rds")))

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
      accuracy_test = cm_test$overall["Accuracy"],
      kappa = cm_test$overall["Kappa"],
      auc = auc_val,
      f1_score = f1_score,
      raster_prediction = rast_pred
    )
  }

  if(verbose) message("\n--- [FIN] Modelado completado para todas las especies ---")
  return(modelos_especies)
}
