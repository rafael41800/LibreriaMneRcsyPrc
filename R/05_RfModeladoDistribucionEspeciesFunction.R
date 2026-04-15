#' @title Modelado de Distribucion de Especies Multiples con Random Forest
#' / Multiple Species Distribution Modeling with Random Forest
#'
#' @description
#' Ejecuta modelos Random Forest para multiples especies a partir de datos de
#' presencia/ausencia, incluyendo validacion cruzada, evaluacion de metricas y
#' prediccion espacial. Genera modelos individuales por especie cuando se cumplen
#' los requisitos minimos de registros.
#' /
#' Executes Random Forest models for multiple species from presence/absence data,
#' including cross-validation, metric evaluation, and spatial prediction.
#' Generates individual models per species when minimum record requirements are met.
#'
#' @param ResultadosList Lista nombrada donde cada elemento corresponde a una
#' especie y contiene un data.frame con columnas de presencia/ausencia y
#' variables predictoras. Debe incluir la columna especificada en VariablePresencia
#' (factor) y las variables especificadas en Variables.
#' / Named list where each element corresponds to a species and contains a
#' data.frame with presence/absence columns and predictor variables. Must
#' include the column specified in VariablePresencia (factor) and variables
#' specified in Variables.
#'
#' @param Variables Vector de caracteres con los nombres de las variables
#' predictoras a utilizar en el modelado. Estas columnas deben estar presentes
#' en los data.frames de ResultadosList. La funcion automaticamente incluye
#' la columna especificada en VariablePresencia.
#' / Character vector with names of predictor variables to use in modeling.
#' These columns must be present in the data.frames of ResultadosList.
#' The function automatically includes the column specified in VariablePresencia.
#'
#' @param VariablePresencia Nombre de la columna que contiene los datos de
#' presencia/ausencia (0/1). Por defecto es "Presence". Debe ser un factor con
#' niveles 0 (ausencia) y 1 (presencia).
#' / Name of the column containing presence/absence data (0/1). Default is
#' "Presence". Must be a factor with levels 0 (absence) and 1 (presence).
#'
#' @param RegistrosMinimos Numero minimo de registros (presencia + ausencia)
#' requeridos para modelar una especie. Por defecto es 10.
#' / Minimum number of records (presence + absence) required to model a species.
#' Default is 10.
#'
#' @param StackRasterReproyectado Objeto raster stack (de terra o raster)
#' con las variables predictoras proyectadas espacialmente. Si se proporciona,
#' se generan mapas de prediccion de probabilidad.
#' / Raster stack object (from terra or raster) with spatially projected
#' predictor variables. If provided, probability prediction maps are generated.
#'
#' @param RutaAlmacenamiento Ruta de directorio donde guardar los resultados
#' (modelos, datos de validacion, graficos, raster). Si es NULL, no se guardan
#' archivos.
#' / Directory path to save results (models, validation data, plots, raster).
#' If NULL, no files are saved.
#'
#' @param NTree Numero de arboles en cada modelo Random Forest. Por defecto 100.
#' / Number of trees in each Random Forest model. Default 100.
#'
#' @param SetSeed Semilla para reproducibilidad. Si es NULL, no se fija semilla.
#' / Seed for reproducibility. If NULL, no seed is set.
#'
#' @param CrossValidationParticiones Numero de folds para validacion cruzada. Si es 1 o menos,
#' no se realiza validacion cruzada. Por defecto 5.
#' / Number of folds for cross-validation. If 1 or less, no cross-validation is
#' performed. Default 5.
#'
#' @param TuneMTry Logico. Si TRUE, se ajusta automaticamente el parametro
#' mtry (numero de variables muestreadas en cada split). Si FALSE, se usa mtry=2.
#' / Logical. If TRUE, automatically tunes the mtry parameter (number of variables
#' sampled at each split). If FALSE, uses mtry=2.
#'
#' @return
#' Una lista nombrada (por especie) con los siguientes componentes para cada
#' modelo exitoso:
#' \itemize{
#' \item modelo: Objeto Random Forest entrenado.
#' \item importance: Importancia de variables.
#' \item ResumenModeloTemporalCrossValidation: Resumen de metricas de validacion cruzada (si aplica).
#' \item ExactitudDatosEntrenamiento, ExactitudDatosPrueba: Precision en entrenamiento y prueba.
#' \item sensitivity, specificity: Sensibilidad y especificidad.
#' \item precision, recall, f1_score: Metricas de clasificacion.
#' \item kappa: Coeficiente Kappa.
#' \item auc: area bajo la curva ROC (si aplica).
#' \item MatrizConfusionDatosEntrenamiento, MatrizConfusionDatosPrueba, MatrizConfusionDatosValidacion: Matrices de confusion.
#' \item RocCurve: Objeto curva ROC (si aplica).
#' \item mtry_used, ntree_used: Parametros usados.
#' \item VariablePresencia_used: Nombre de la variable de presencia utilizada.
#' \item raster_prediction: Raster de prediccion (si se proporciono stack).
#' }
#' /
#' A named list (by species) with the following components for each successful model:
#' \itemize{
#' \item modelo: Trained Random Forest object.
#' \item importance: Variable importance.
#' \item ResumenModeloTemporalCrossValidation: Cross-validation metrics summary (if applicable).
#' \item ExactitudDatosEntrenamiento, ExactitudDatosPrueba: Accuracy on training and test sets.
#' \item sensitivity, specificity: Sensitivity and specificity.
#' \item precision, recall, f1_score: Classification metrics.
#' \item kappa: Kappa coefficient.
#' \item auc: Area under ROC curve (if applicable).
#' \item MatrizConfusionDatosEntrenamiento, MatrizConfusionDatosPrueba, MatrizConfusionDatosValidacion: Confusion matrices.
#' \item RocCurve: ROC curve object (if applicable).
#' \item mtry_used, ntree_used: Parameters used.
#' \item VariablePresencia_used: Name of the presence variable used.
#' \item raster_prediction: Prediction raster (if stack provided).
#' }
#'
#' @details
#' La funcion sigue este flujo por cada especie:
#' 1. Filtrado: Solo especies con al menos RegistrosMinimos.
#' 2. Preprocesamiento: Limpieza de NA, factorizacion de VariablePresencia.
#' 3. Particion: 10% validacion final, del resto 70% entrenamiento y 30% prueba.
#' 4. Ajuste mtry: Si TuneMTry=TRUE, optimiza con tuneRF.
#' 5. Validacion cruzada: Si CrossValidationParticiones>1, calcula metricas por fold.
#' 6. Modelado: Entrena Random Forest con parametros optimos.
#' 7. Evaluacion: Calcula multiples metricas en entrenamiento, prueba y validacion.
#' 8. Prediccion espacial: Si hay stack raster, genera mapa de probabilidad.
#' 9. Almacenamiento: Si hay ruta, guarda modelos, graficos y rasters.
#' /
#' The function follows this workflow per species:
#' 1. Filtering: Only species with at least RegistrosMinimos.
#' 2. Preprocessing: NA cleaning, VariablePresencia factorization.
#' 3. Splitting: 10% final validation, from remainder 70% training and 30% test.
#' 4. mtry tuning: If TuneMTry=TRUE, optimizes with tuneRF.
#' 5. Cross-validation: If CrossValidationParticiones>1, calculates metrics per fold.
#' 6. Modeling: Trains Random Forest with optimal parameters.
#' 7. Evaluation: Calculates multiple metrics on training, test, and validation.
#' 8. Spatial prediction: If raster stack exists, generates probability map.
#' 9. Storage: If path provided, saves models, plots, and rasters.
#'
#' @note
#' * Se requieren al menos 2 clases (presencia y ausencia) para modelar.
#' * La validacion cruzada solo se realiza si hay suficientes datos.
#' * Las graficas de importancia de variables y curvas ROC se guardan automaticamente.
#' * La funcion maneja errores internamente y omite especies problematicas.
#' * Se asume que VariablePresencia es factor con niveles 0 (ausencia) y 1 (presencia).
#' * La columna VariablePresencia se incluye automaticamente en los datos de modelado.
#' /
#' * At least 2 classes (presence and absence) are required for modeling.
#' * Cross-validation is only performed if enough data exists.
#' * Variable importance plots and ROC curves are saved automatically.
#' * The function handles errors internally and skips problematic species.
#' * Assumes VariablePresencia column is factor with levels 0 (absence) and 1 (presence).
#' * The VariablePresencia column is automatically included in the modeling data.
#'
#' @section Advertencias/Warnings:
#' * Con muchas especies o variables, el tiempo de ejecucion puede ser alto.
#' * La memoria requerida depende del tamaño del raster y numero de especies.
#' * La calidad del modelo depende criticamente de la calidad de los datos de entrada.
#' * El umbral de RegistrosMinimos debe ajustarse segun la realidad de los datos.
#' /
#' * With many species or variables, execution time may be high.
#' * Required memory depends on raster size and number of species.
#' * Model quality critically depends on input data quality.
#' * RegistrosMinimos threshold should be adjusted based on data reality.
#'
#' @examples
#' \dontrun{
#' # Cargar datos de ejemplo / Load example data
#' data(presence_absence_data) # Lista con datos por especie / List with data per species
#' data(predictor_stack) # Stack raster con variables / Raster stack with variables
#'
#' # Ejecutar modelos con nombre por defecto / Run models with default name
#' modelos <- RfMdeMultiEspecies(
#' ResultadosList = presence_absence_data,
#' Variables = c("temp_mean", "precip", "elevation", "soil_type"),
#' RegistrosMinimos = 15,
#' StackRasterReproyectado = predictor_stack,
#' RutaAlmacenamiento = "./resultados_sdm",
#' NTree = 500,
#' SetSeed = 123,
#' CrossValidationParticiones = 5,
#' TuneMTry = TRUE
#' )
#'
#' # Ejecutar modelos con nombre personalizado / Run models with custom name
#' modelos2 <- RfMdeMultiEspecies(
#' ResultadosList = presence_absence_data,
#' Variables = c("temp_mean", "precip", "elevation", "soil_type"),
#' VariablePresencia = "Ocurrencia",
#' RegistrosMinimos = 15,
#' StackRasterReproyectado = predictor_stack,
#' RutaAlmacenamiento = "./resultados_sdm2",
#' NTree = 500,
#' SetSeed = 123,
#' CrossValidationParticiones = 5,
#' TuneMTry = TRUE
#' )
#'
#' # Explorar resultados / Explore results
#' names(modelos)
#' modelos[["Quercus_robur"]]$ExactitudDatosPrueba
#' modelos[["Quercus_robur"]]$auc
#' modelos[["Quercus_robur"]]$VariablePresencia_used
#'
#' # Visualizar importancia de variables / Visualize variable importance
#' barplot(sort(modelos[["Pinus_sylvestris"]]$importance[, "MeanDecreaseGini"],
#' decreasing = TRUE))
#' }
#'
#' @seealso
#' * \code{\link[randomForest]{randomForest}} para detalles del algoritmo.
#' * \code{\link[caret]{confusionMatrix}} para metricas de clasificacion.
#' * \code{\link[pROC]{roc}} para curvas ROC.
#' * \code{\link[terra]{predict}} para prediccion espacial.
#' /
#' * \code{\link[randomForest]{randomForest}} for algorithm details.
#' * \code{\link[caret]{confusionMatrix}} for classification metrics.
#' * \code{\link[pROC]{roc}} for ROC curves.
#' * \code{\link[terra]{predict}} for spatial prediction.
#'
#' @keywords
#' modelado de distribucion de especies, SDM, random forest, aprendizaje automatico,
#' validacion cruzada, evaluacion de modelos, prediccion espacial,
#' species distribution modeling, SDM, random forest, machine learning,
#' cross-validation, model evaluation, spatial prediction
#'
#' @family
#' funciones_modelado_especies, funciones_random_forest,
#' funciones_distribucion_especies /
#' species_modeling_functions, random_forest_functions,
#' species_distribution_functions
#'
#' @importFrom randomForest randomForest tuneRF varImpPlot
#' @importFrom caret createDataPartition createFolds confusionMatrix
#' @importFrom pROC roc auc
#' @importFrom terra predict writeRaster
#' @importFrom dplyr select %>% all_of
#' @importFrom stats as.formula complete.cases predict
#' @importFrom grDevices png dev.off
#' @export
RfMdeMultiEspecies <- function(ResultadosList, Variables, VariablePresencia = "Presence", RegistrosMinimos = 10, StackRasterReproyectado = NULL, RutaAlmacenamiento = NULL, NTree = 100, SetSeed = NULL, CrossValidationParticiones = 5, TuneMTry = TRUE) {
  modelos_especies <- list()
  for (m in names(ResultadosList)) {
    if(!is.null(ResultadosList[[m]]) && nrow(as.data.frame(ResultadosList[[m]])) >= RegistrosMinimos) {
      if(!is.null(SetSeed)) set.seed(SetSeed)
      DatosPresenciasAusenciasJuntosDf <- as.data.frame(ResultadosList[[m]])
      DatosPresenciasAusenciasJuntosDfClean <- DatosPresenciasAusenciasJuntosDf %>%
        dplyr::select(all_of(VariablePresencia), all_of(Variables))  # MODIFICADO
      DatosParaModelado <- as.data.frame(DatosPresenciasAusenciasJuntosDfClean)
      DatosParaModelado <- DatosParaModelado[stats::complete.cases(DatosParaModelado), ]
      if(nrow(DatosParaModelado) < RegistrosMinimos) next
      DatosParaModelado[[VariablePresencia]] <- as.factor(DatosParaModelado[[VariablePresencia]])  # MODIFICADO
      if(length(unique(DatosParaModelado[[VariablePresencia]])) < 2) next  # MODIFICADO
      IndicesValidacion <- caret::createDataPartition(DatosParaModelado[[VariablePresencia]], p = 0.1, list = FALSE)  # MODIFICADO
      DatosParaValidacion <- DatosParaModelado[IndicesValidacion, ]
      DatosParaModelacion <- DatosParaModelado[-IndicesValidacion, ]
      IndicesEntrenamiento <- caret::createDataPartition(DatosParaModelacion[[VariablePresencia]], p = 0.7, list = FALSE)  # MODIFICADO
      DatosParaEntrenamiento <- DatosParaModelacion[IndicesEntrenamiento, ]
      DatosParaTesteo <- DatosParaModelacion[-IndicesEntrenamiento, ]
      MTryOptimal <- 2
      if(TuneMTry && nrow(DatosParaEntrenamiento) >= RegistrosMinimos) {
        tryCatch({
          TuneRfResult <- randomForest::tuneRF(
            x = DatosParaEntrenamiento[, !names(DatosParaEntrenamiento) %in% VariablePresencia],  # MODIFICADO
            y = DatosParaEntrenamiento[[VariablePresencia]],  # MODIFICADO
            ntreeTry = min(100, NTree),
            stepFactor = 1.5,
            improve = 0.01,
            trace = FALSE,
            plot = FALSE
          )
          if(!is.null(TuneRfResult) && nrow(TuneRfResult) > 0) {
            MTryOptimal <- TuneRfResult[which.min(TuneRfResult[,2]), 1]
          }
        }, error = function(e) {})
      }
      if(CrossValidationParticiones > 1 && nrow(DatosParaEntrenamiento) > CrossValidationParticiones * 5) {
        IndicesCrossValidation <- caret::createFolds(DatosParaEntrenamiento[[VariablePresencia]], k = CrossValidationParticiones)  # MODIFICADO
        MetricasCrossValidation <- list()
        for(fold in 1:CrossValidationParticiones) {
          IndicesParticionEntrenamiento <- unlist(IndicesCrossValidation[-fold])
          IndicesParticionValidacion <- IndicesCrossValidation[[fold]]
          ParticionesParaEntrenamiento <- DatosParaEntrenamiento[IndicesParticionEntrenamiento, ]
          ParticionesParaValidacion <- DatosParaEntrenamiento[IndicesParticionValidacion, ]
          FormulaRf <- as.formula(paste(VariablePresencia, "~ ."))  # MODIFICADO
          ModeloTemporalCrossValidation <- randomForest::randomForest(
            FormulaRf,  # MODIFICADO
            data = ParticionesParaEntrenamiento,
            importance = FALSE,
            ntree = NTree,
            mtry = MTryOptimal,
            do.trace = FALSE,
            verbose = FALSE
          )
          PrediccionCrossValidation <- stats::predict(ModeloTemporalCrossValidation, ParticionesParaValidacion)
          MatrizConfusionCrossValidation <- caret::confusionMatrix(PrediccionCrossValidation, ParticionesParaValidacion[[VariablePresencia]], positive = "1")  # MODIFICADO
          MetricasCrossValidation[[fold]] <- list(
            accuracy = MatrizConfusionCrossValidation$overall["Accuracy"],
            sensitivity = MatrizConfusionCrossValidation$byClass["Sensitivity"],
            specificity = MatrizConfusionCrossValidation$byClass["Specificity"],
            kappa = MatrizConfusionCrossValidation$overall["Kappa"]
          )
        }
        ResumenModeloTemporalCrossValidation <- data.frame(
          AccuracyMean = mean(sapply(MetricasCrossValidation, function(x) x$accuracy)),
          SensitivityMean = mean(sapply(MetricasCrossValidation, function(x) x$sensitivity)),
          SpecificityMean = mean(sapply(MetricasCrossValidation, function(x) x$specificity)),
          KappaMean = mean(sapply(MetricasCrossValidation, function(x) x$kappa))
        )
      } else {
        ResumenModeloTemporalCrossValidation <- NULL
      }
      FormulaRf <- as.formula(paste(VariablePresencia, "~ ."))  # MODIFICADO
      RfModelo <- randomForest::randomForest(
        FormulaRf,  # MODIFICADO
        data = DatosParaEntrenamiento,
        importance = TRUE,
        ntree = NTree,
        mtry = MTryOptimal,
        do.trace = FALSE,
        verbose = FALSE
      )
      PrediccionDatosEntrenamiento <- stats::predict(RfModelo, DatosParaEntrenamiento)
      PrediccionDatosPrueba <- stats::predict(RfModelo, DatosParaTesteo)
      MatrizConfusionDatosEntrenamiento <- caret::confusionMatrix(PrediccionDatosEntrenamiento,
                                                                DatosParaEntrenamiento[[VariablePresencia]])  # MODIFICADO
      MatrizConfusionDatosPrueba <- caret::confusionMatrix(PrediccionDatosPrueba,
                                                         DatosParaTesteo[[VariablePresencia]])  # MODIFICADO
      ExactitudDatosEntrenamiento <- MatrizConfusionDatosEntrenamiento$overall["Accuracy"]
      ExactitudDatosPrueba <- MatrizConfusionDatosPrueba$overall["Accuracy"]
      SensibilidadDatosPrueba <- MatrizConfusionDatosPrueba$byClass["Sensitivity"]
      EspecificidadDatosPrueba <- MatrizConfusionDatosPrueba$byClass["Specificity"]
      KappaTest <- MatrizConfusionDatosPrueba$overall["Kappa"]
      PrediccionProbaDatosPrueba <- stats::predict(RfModelo, DatosParaTesteo, type = "prob")
      AucValor <- NA
      RocCurve <- NULL
      if(length(unique(DatosParaTesteo[[VariablePresencia]])) > 1) {  # MODIFICADO
        RocCurve <- pROC::roc(response = DatosParaTesteo[[VariablePresencia]],  # MODIFICADO
                               predictor = PrediccionProbaDatosPrueba[, "1"])
        AucValor <- pROC::auc(RocCurve)
      }
      PrediccionDatosValidacion <- stats::predict(RfModelo, DatosParaValidacion)
      MatrizConfusionDatosValidacion <- caret::confusionMatrix(PrediccionDatosValidacion, DatosParaValidacion[[VariablePresencia]])  # MODIFICADO
      PrecisionTest <- MatrizConfusionDatosPrueba$byClass["Pos Pred Value"]
      RecallTest <- MatrizConfusionDatosPrueba$byClass["Sensitivity"]
      F1Test <- if(!is.na(PrecisionTest) && !is.na(RecallTest) && (PrecisionTest + RecallTest) > 0) {
        2 * (PrecisionTest * RecallTest) / (PrecisionTest + RecallTest)
      } else NA
      InformacionModelo <- list(
        modelo = RfModelo,
        importance = as.data.frame(RfModelo$importance),
        ResumenModeloTemporalCrossValidation = ResumenModeloTemporalCrossValidation,
        ExactitudDatosEntrenamiento = ExactitudDatosEntrenamiento,
        ExactitudDatosPrueba = ExactitudDatosPrueba,
        sensitivity = SensibilidadDatosPrueba,
        specificity = EspecificidadDatosPrueba,
        precision = PrecisionTest,
        recall = RecallTest,
        f1_score = F1Test,
        kappa = KappaTest,
        auc = AucValor,
        MatrizConfusionDatosEntrenamiento = MatrizConfusionDatosEntrenamiento,
        MatrizConfusionDatosPrueba = MatrizConfusionDatosPrueba,
        MatrizConfusionDatosValidacion = MatrizConfusionDatosValidacion,
        RocCurve = RocCurve,
        mtry_used = MTryOptimal,
        ntree_used = NTree,
        VariablePresencia_used = VariablePresencia  # MODIFICADO
      )
      if(!is.null(RutaAlmacenamiento)) {
        dir.create(RutaAlmacenamiento, showWarnings = FALSE, recursive = TRUE)
        saveRDS(RfModelo, file.path(RutaAlmacenamiento, paste0("modelo_randomforest_", m, ".rds")))
        saveRDS(DatosParaValidacion, file.path(RutaAlmacenamiento, paste0("datos_validacion_final_", m, ".rds")))
        png(file.path(RutaAlmacenamiento, paste0("var_importance_", m, ".png")),
            width = 800, height = 600)
        randomForest::varImpPlot(RfModelo)
        dev.off()
        if(!is.null(RocCurve)) {
          png(file.path(RutaAlmacenamiento, paste0("RocCurve_", m, ".png")),
              width = 600, height = 600)
          plot(RocCurve)
          dev.off()
        }
      }
      if(!is.null(StackRasterReproyectado)) {
        RasterDePrediccion <- terra::predict(StackRasterReproyectado, RfModelo,
                                      na.rm = TRUE, type = "prob")
        if(!is.null(RutaAlmacenamiento)) {
          terra::writeRaster(RasterDePrediccion[[2]],
                             file.path(RutaAlmacenamiento, paste0("prediccion_prob_", m, ".tif")),
                             overwrite = TRUE)
        }
        InformacionModelo$raster_prediction <- RasterDePrediccion[[2]]
      }
      modelos_especies[[m]] <- InformacionModelo
    }
  }
  return(modelos_especies)
}









