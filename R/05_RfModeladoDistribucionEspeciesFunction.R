#' Modelado optimizado de Distribucion de Especies Multiples con Random Forest
#'
#' @description
#' Ejecuta modelos Random Forest para multiples especies a partir de datos tabulares de
#' presencia/ausencia. La funcion esta optimizada mediante iteracion funcional (`lapply`)
#' e indexacion matricial en R base, eliminando dependencias de Tidyverse para acelerar
#' la validacion cruzada, la optimizacion de hiperparametros y la prediccion espacial.
#'
#' @param ResultadosList Lista de objetos `data.frame` (salida tabular de la Funcion 4).
#' @param Variables Vector de caracteres con los nombres de las variables predictoras ambientales.
#' @param VariablePresencia Caracter. Nombre de la columna de respuesta. Por defecto `"Presence"`.
#' @param RegistrosMinimos Numerico. Minimo de datos validos necesarios para intentar ajustar el modelo. Por defecto `10`.
#' @param StackRasterReproyectado Objeto `SpatRaster` opcional utilizado para la proyeccion espacial continua del modelo.
#' @param RutaAlmacenamiento Caracter. Directorio para guardar `.rds`, `.tif` y graficos de importancia. Si es `NULL`, la funcion opera solo en memoria.
#' @param NTree Numerico. Numero de arboles a construir en el bosque. Por defecto `100`.
#' @param SetSeed Numerico. Semilla para garantizar la reproducibilidad de las particiones y el remuestreo.
#' @param CrossValidationParticiones Numerico. Numero de particiones (folds) para k-fold CV sobre el set de entrenamiento. Por defecto `5`.
#' @param TuneMTry Logico. Si es `TRUE`, utiliza \code{tuneRF} para encontrar el hiperparametro espacial optimo automaticamente.
#' @param GuardarDatos Caracter. Define que subconjuntos de datos exportar: `"none"` (default), `"train"`, `"valid"` o `"both"`.
#' @param verbose Logico. Si es `TRUE`, imprime el progreso en consola.
#'
#' @return
#' Una lista estructurada y nombrada por especie. Cada elemento contiene el modelo ajustado,
#' la importancia de variables, las metricas de validacion interna (CV) y externa (AUC, Kappa, F1),
#' el mapa de idoneidad (si se configuro) y los datos particionados.
#'
#' @details
#' El flujo de trabajo interno sigue estos pasos independientes por especie:
#' \enumerate{
#'   \item \strong{Limpieza Matricial}: Extrae las columnas estrictamente necesarias usando R base y elimina los registros con `NA`.
#'   \item \strong{Particion de Datos}: Divide el set aleatoriamente en 75\% para entrenamiento y 25\% para validacion externa.
#'   \item \strong{Optimizacion}: Si se solicita, busca el valor optimo de variables candidatas por division (`mtry`).
#'   \item \strong{Validacion Interna (K-Fold CV)}: Estima la precision y el estadistico Kappa en datos no vistos del set de entrenamiento.
#'   \item \strong{Ajuste y Metricas}: Entrena el modelo Random Forest final y calcula el Area Bajo la Curva (AUC) y F1-Score con el set del 25\%.
#'   \item \strong{Espacializacion}: Proyecta la probabilidad de presencia sobre el raster de referencia.
#'   \item \strong{Exportacion Organizada}: Guarda los resultados en disco estructurando los nombres de archivo con la especie como prefijo (ej. `Especie_RF_Model.rds`).
#' }
#'
#' @note
#' Esta funcion esta diseñada para integrarse en un flujo asincrono. Si cuenta con multiples nucleos,
#' puede sustituir `lapply` internamente por `future.apply::future_lapply` para modelar especies en paralelo.
#'
#' @section Advertencias:
#' \itemize{
#'   \item \strong{Balance de Clases}: El modelo omitira la especie si los conjuntos particionados terminan con solo presencias o solo ausencias.
#'   \item \strong{Rendimiento}: La prediccion espacial sobre rasters de alta resolucion es el paso mas costoso en RAM; ajuste sus recursos en consecuencia.
#' }
#'
#' @section Visualizacion de Resultados:
#' Para revisar los mapas de idoneidad espacial generados (`.tif`), los cuales por defecto
#' almacenan valores continuos de probabilidad (0 a 1), se recomienda visualizarlos utilizando
#' una paleta de colores en formato "semaforo". En esta representacion, el color verde
#' indica las areas mas probables (idoneas) y el rojo las menos probables.
#'
#' @examples
#' \dontrun{
#' library(randomForest)
#' library(terra)
#'
#' # 1. Crear datos simulados
#' df_sp <- data.frame(
#'    Presence = sample(c(0, 1), 50, replace = TRUE),
#'    Bio1 = runif(50, 10, 30),
#'    Bio2 = runif(50, 50, 100)
#' )
#' lista_datos <- list("Sp_A" = df_sp)
#'
#' # 2. Raster sintetico
#' r_env <- rast(nrows=10, ncols=10, ext=ext(0, 10, 0, 10))
#' values(r_env) <- runif(100)
#' names(r_env) <- c("Bio1")
#' add(r_env) <- r_env * 2
#' names(r_env)[2] <- "Bio2"
#'
#' # 3. Correr modelo
#' res <- RfMdeMultiEspecies(
#'    ResultadosList = lista_datos,
#'    Variables = c("Bio1", "Bio2"),
#'    VariablePresencia = "Presence",
#'    StackRasterReproyectado = r_env,
#'    RutaAlmacenamiento = tempdir()
#' )
#'
#' # 4. Revisar metricas y visualizar el mapa con paleta semaforo
#' cat("Exactitud:", round(res[["Sp_A"]]$accuracy_valid, 3), "\n")
#' terra::plot(res[["Sp_A"]]$raster_prediction,
#'             main = "Probabilidad de Distribucion - Sp_A",
#'             col = colorRampPalette(c("red", "yellow", "green"))(100))
#' }
#'
#' @seealso
#' \code{\link[randomForest]{randomForest}}, \code{\link[caret]{createDataPartition}}, \code{\link[pROC]{roc}}
#'
#' @keywords random forest sdm validacion prediccion
#'
#' @importFrom randomForest randomForest tuneRF varImpPlot
#' @importFrom caret createDataPartition createFolds confusionMatrix
#' @importFrom pROC roc auc
#' @importFrom terra predict writeRaster
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

  VariablePresencia <- as.character(VariablePresencia)
  Variables <- setdiff(as.character(Variables), VariablePresencia)

  if(verbose) message("\n--- [INICIO] Modelado Random Forest Multiespecie (75% train / 25% valid) ---")

  NombresEspecies <- names(ResultadosList)
  TotalEspecies <- length(NombresEspecies)

  # Optimizacion 1: Sustitucion del bucle for por lapply para facilitar procesamiento y futura paralelización
  modelos_especies <- lapply(seq_along(NombresEspecies), function(idx) {
    m <- NombresEspecies[idx]

    datos_sp <- ResultadosList[[m]]
    if(!inherits(datos_sp, "data.frame")) {
      datos_sp <- as.data.frame(datos_sp)
    }

    # Optimizacion 2: Indexación matricial en R base en lugar de dplyr (Elimina NSE)
    ColumnasRequeridas <- c(VariablePresencia, Variables)
    DatosParaModelado <- stats::na.omit(datos_sp[, ColumnasRequeridas, drop = FALSE])

    n_registros <- nrow(DatosParaModelado)
    if(n_registros < RegistrosMinimos) {
      if(verbose) message(sprintf("  ! %s: Registros insuficientes (%d). Saltando...", m, n_registros))
      return(NULL)
    }

    DatosParaModelado[[VariablePresencia]] <- as.factor(DatosParaModelado[[VariablePresencia]])

    if(length(unique(DatosParaModelado[[VariablePresencia]])) < 2) {
      if(verbose) message(sprintf("  ! %s: Solo se encontro una clase (Presencia o Ausencia). Saltando...", m))
      return(NULL)
    }

    if(verbose) message(sprintf("-> Procesando [%d/%d]: %s", idx, TotalEspecies, m))

    idx_train <- caret::createDataPartition(DatosParaModelado[[VariablePresencia]], p = 0.75, list = FALSE)
    DatosTrain <- DatosParaModelado[idx_train, , drop = FALSE]
    DatosValid <- DatosParaModelado[-idx_train, , drop = FALSE]

    if(length(unique(DatosValid[[VariablePresencia]])) < 2) {
      if(verbose) message(sprintf("  ! %s: Conjunto de validacion con una sola clase. Saltando...", m))
      return(NULL)
    }

    mtry_opt <- floor(sqrt(length(Variables)))
    if(TuneMTry) {
      try({
        tuner <- randomForest::tuneRF(x = DatosTrain[, Variables, drop = FALSE],
                                      y = DatosTrain[[VariablePresencia]],
                                      ntreeTry = 50, stepFactor = 1.5, trace = FALSE, plot = FALSE)
        mtry_opt <- tuner[which.min(tuner[,2]), 1]
      }, silent = TRUE)
    }

    cv_resumen <- NULL
    if(CrossValidationParticiones > 1 && nrow(DatosTrain) > (CrossValidationParticiones * 2)) {
      folds <- caret::createFolds(DatosTrain[[VariablePresencia]], k = CrossValidationParticiones)
      cv_metrics <- lapply(folds, function(idx_f){
        m_cv <- randomForest::randomForest(x = DatosTrain[-idx_f, Variables, drop = FALSE],
                                           y = DatosTrain[-idx_f, VariablePresencia],
                                           ntree = NTree, mtry = mtry_opt)
        pred_cv <- stats::predict(m_cv, DatosTrain[idx_f, , drop = FALSE])
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

    # Optimizacion 3: Concatenación ajustada para que empiece con el nombre de la especie
    if(!is.null(RutaAlmacenamiento)) {
      if(!dir.exists(RutaAlmacenamiento)) dir.create(RutaAlmacenamiento, showWarnings = FALSE, recursive = TRUE)

      saveRDS(rf_final, file.path(RutaAlmacenamiento, paste0(m, "_RF_Model.rds")))

      if(GuardarDatos %in% c("train", "both")) {
        saveRDS(DatosTrain, file.path(RutaAlmacenamiento, paste0(m, "_Train_Data.rds")))
      }
      if(GuardarDatos %in% c("valid", "both")) {
        saveRDS(DatosValid, file.path(RutaAlmacenamiento, paste0(m, "_Valid_Data.rds")))
      }

      grDevices::png(file.path(RutaAlmacenamiento, paste0(m, "_Imp_Plot.png")), width = 800, height = 600)
      randomForest::varImpPlot(rf_final, main = paste("Importancia de Variables -", m))
      grDevices::dev.off()

      if(!is.null(rast_pred)) {
        terra::writeRaster(rast_pred, file.path(RutaAlmacenamiento, paste0(m, "_Map.tif")), overwrite = TRUE)
      }
    }

    return(list(
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
    ))
  })

  names(modelos_especies) <- NombresEspecies
  modelos_especies <- Filter(Negate(is.null), modelos_especies)

  if(verbose) message("\n--- [FIN] Modelado completado para todas las especies ---")
  return(modelos_especies)
}
