#' @title Modelado de Distribución de Especies Múltiples con Random Forest
#' / Multiple Species Distribution Modeling with Random Forest
#'
#' @description
#' Ejecuta modelos Random Forest para múltiples especies a partir de datos de
#' presencia/ausencia, incluyendo validación cruzada, evaluación de métricas y
#' predicción espacial. Genera modelos individuales por especie cuando se cumplen
#' los requisitos mínimos de registros.
#' /
#' Executes Random Forest models for multiple species from presence/absence data,
#' including cross-validation, metric evaluation, and spatial prediction.
#' Generates individual models per species when minimum record requirements are met.
#'
#' @param Resultados_List Lista nombrada donde cada elemento corresponde a una
#' especie y contiene un data.frame con columnas de presencia/ausencia y
#' variables predictoras. Debe incluir la columna Presence (factor) y las
#' variables especificadas en Variables.
#' / Named list where each element corresponds to a species and contains a
#' data.frame with presence/absence columns and predictor variables. Must
#' include the Presence column (factor) and variables specified in Variables.
#'
#' @param Variables Vector de caracteres con los nombres de las variables
#' predictoras a utilizar en el modelado. Estas columnas deben estar presentes
#' en los data.frames de Resultados_List.
#' / Character vector with names of predictor variables to use in modeling.
#' These columns must be present in the data.frames of Resultados_List.
#'
#' @param Registros_Minimos Número mínimo de registros (presencia + ausencia)
#' requeridos para modelar una especie. Por defecto es 10.
#' / Minimum number of records (presence + absence) required to model a species.
#' Default is 10.
#'
#' @param Stack_Raster_Reproyectado Objeto raster stack (de terra o raster)
#' con las variables predictoras proyectadas espacialmente. Si se proporciona,
#' se generan mapas de predicción de probabilidad.
#' / Raster stack object (from terra or raster) with spatially projected
#' predictor variables. If provided, probability prediction maps are generated.
#'
#' @param Ruta_Almacenamiento Ruta de directorio donde guardar los resultados
#' (modelos, datos de validación, gráficos, raster). Si es NULL, no se guardan
#' archivos.
#' / Directory path to save results (models, validation data, plots, raster).
#' If NULL, no files are saved.
#'
#' @param N_Tree Número de árboles en cada modelo Random Forest. Por defecto 100.
#' / Number of trees in each Random Forest model. Default 100.
#'
#' @param Set_Seed Semilla para reproducibilidad. Si es NULL, no se fija semilla.
#' / Seed for reproducibility. If NULL, no seed is set.
#'
#' @param Cross_V_Particiones Número de folds para validación cruzada. Si es 1 o menos,
#' no se realiza validación cruzada. Por defecto 5.
#' / Number of folds for cross-validation. If 1 or less, no cross-validation is
#' performed. Default 5.
#'
#' @param Tune_M_Try Lógico. Si TRUE, se ajusta automáticamente el parámetro
#' mtry (número de variables muestreadas en cada split). Si FALSE, se usa mtry=2.
#' / Logical. If TRUE, automatically tunes the mtry parameter (number of variables
#' sampled at each split). If FALSE, uses mtry=2.
#'
#' @return
#' Una lista nombrada (por especie) con los siguientes componentes para cada
#' modelo exitoso:
#' \itemize{
#' \item modelo: Objeto Random Forest entrenado.
#' \item importance: Importancia de variables.
#' \item Resumen_Modelo_Temporal_CrossValidation: Resumen de métricas de validación cruzada (si aplica).
#' \item accuracy_train, accuracy_test: Precisión en entrenamiento y prueba.
#' \item sensitivity, specificity: Sensibilidad y especificidad.
#' \item precision, recall, f1_score: Métricas de clasificación.
#' \item kappa: Coeficiente Kappa.
#' \item auc: Área bajo la curva ROC (si aplica).
#' \item conf_matrix_train, conf_matrix_test, conf_matrix_val: Matrices de confusión.
#' \item roc_curve: Objeto curva ROC (si aplica).
#' \item mtry_used, ntree_used: Parámetros usados.
#' \item raster_prediction: Raster de predicción (si se proporcionó stack).
#' }
#' /
#' A named list (by species) with the following components for each successful model:
#' \itemize{
#' \item modelo: Trained Random Forest object.
#' \item importance: Variable importance.
#' \item Resumen_Modelo_Temporal_CrossValidation: Cross-validation metrics summary (if applicable).
#' \item accuracy_train, accuracy_test: Accuracy on training and test sets.
#' \item sensitivity, specificity: Sensitivity and specificity.
#' \item precision, recall, f1_score: Classification metrics.
#' \item kappa: Kappa coefficient.
#' \item auc: Area under ROC curve (if applicable).
#' \item conf_matrix_train, conf_matrix_test, conf_matrix_val: Confusion matrices.
#' \item roc_curve: ROC curve object (if applicable).
#' \item mtry_used, ntree_used: Parameters used.
#' \item raster_prediction: Prediction raster (if stack provided).
#' }
#'
#' @details
#' La función sigue este flujo por cada especie:
#' 1. Filtrado: Solo especies con al menos Registros_Minimos.
#' 2. Preprocesamiento: Limpieza de NA, factorización de Presence.
#' 3. Partición: 10% validación final, del resto 70% entrenamiento y 30% prueba.
#' 4. Ajuste mtry: Si Tune_M_Try=TRUE, optimiza con tuneRF.
#' 5. Validación cruzada: Si Cross_V_Particiones>1, calcula métricas por fold.
#' 6. Modelado: Entrena Random Forest con parámetros óptimos.
#' 7. Evaluación: Calcula múltiples métricas en entrenamiento, prueba y validación.
#' 8. Predicción espacial: Si hay stack raster, genera mapa de probabilidad.
#' 9. Almacenamiento: Si hay ruta, guarda modelos, gráficos y rasters.
#' /
#' The function follows this workflow per species:
#' 1. Filtering: Only species with at least Registros_Minimos.
#' 2. Preprocessing: NA cleaning, Presence factorization.
#' 3. Splitting: 10% final validation, from remainder 70% training and 30% test.
#' 4. mtry tuning: If Tune_M_Try=TRUE, optimizes with tuneRF.
#' 5. Cross-validation: If Cross_V_Particiones>1, calculates metrics per fold.
#' 6. Modeling: Trains Random Forest with optimal parameters.
#' 7. Evaluation: Calculates multiple metrics on training, test, and validation.
#' 8. Spatial prediction: If raster stack exists, generates probability map.
#' 9. Storage: If path provided, saves models, plots, and rasters.
#'
#' @note
#' * Se requieren al menos 2 clases (presencia y ausencia) para modelar.
#' * La validación cruzada solo se realiza si hay suficientes datos.
#' * Las gráficas de importancia de variables y curvas ROC se guardan automáticamente.
#' * La función maneja errores internamente y omite especies problemáticas.
#' * Se asume que la columna Presence es factor con niveles 0 (ausencia) y 1 (presencia).
#' /
#' * At least 2 classes (presence and absence) are required for modeling.
#' * Cross-validation is only performed if enough data exists.
#' * Variable importance plots and ROC curves are saved automatically.
#' * The function handles errors internally and skips problematic species.
#' * Assumes Presence column is factor with levels 0 (absence) and 1 (presence).
#'
#' @section Advertencias/Warnings:
#' * Con muchas especies o variables, el tiempo de ejecución puede ser alto.
#' * La memoria requerida depende del tamaño del raster y número de especies.
#' * La calidad del modelo depende críticamente de la calidad de los datos de entrada.
#' * El umbral de Registros_Minimos debe ajustarse según la realidad de los datos.
#' /
#' * With many species or variables, execution time may be high.
#' * Required memory depends on raster size and number of species.
#' * Model quality critically depends on input data quality.
#' * Registros_Minimos threshold should be adjusted based on data reality.
#'
#' @examples
#' \dontrun{
#' # Cargar datos de ejemplo / Load example data
#' data(presence_absence_data) # Lista con datos por especie / List with data per species
#' data(predictor_stack) # Stack raster con variables / Raster stack with variables
#'
#' # Ejecutar modelos / Run models
#' modelos <- Rf_SDM_multiple_species(
#' Resultados_List = presence_absence_data,
#' Variables = c("temp_mean", "precip", "elevation", "soil_type"),
#' Registros_Minimos = 15,
#' Stack_Raster_Reproyectado = predictor_stack,
#' Ruta_Almacenamiento = "./resultados_sdm",
#' N_Tree = 500,
#' Set_Seed = 123,
#' Cross_V_Particiones = 5,
#' Tune_M_Try = TRUE
#' )
#'
#' # Explorar resultados / Explore results
#' names(modelos)
#' modelos[["Quercus_robur"]]$accuracy_test
#' modelos[["Quercus_robur"]]$auc
#'
#' # Visualizar importancia de variables / Visualize variable importance
#' barplot(sort(modelos[["Pinus_sylvestris"]]$importance[, "MeanDecreaseGini"],
#' decreasing = TRUE))
#' }
#'
#' @seealso
#' * \code{\link[randomForest]{randomForest}} para detalles del algoritmo.
#' * \code{\link[caret]{confusionMatrix}} para métricas de clasificación.
#' * \code{\link[pROC]{roc}} para curvas ROC.
#' * \code{\link[terra]{predict}} para predicción espacial.
#' /
#' * \code{\link[randomForest]{randomForest}} for algorithm details.
#' * \code{\link[caret]{confusionMatrix}} for classification metrics.
#' * \code{\link[pROC]{roc}} for ROC curves.
#' * \code{\link[terra]{predict}} for spatial prediction.
#'
#' @keywords
#' modelado de distribución de especies, SDM, random forest, aprendizaje automático,
#' validación cruzada, evaluación de modelos, predicción espacial,
#' species distribution modeling, SDM, random forest, machine learning,
#' cross-validation, model evaluation, spatial prediction
#'
#' @family
#' funciones_modelado_especies, funciones_random_forest,
#' funciones_distribucion_especies /
#' species_modeling_functions, random_forest_functions,
#' species_distribution_functions
#'
#' @importFrom randomForest randomForest tuneRF
#' @importFrom caret createDataPartition createFolds confusionMatrix
#' @importFrom pROC roc auc
#' @importFrom terra predict writeRaster
#' @importFrom dplyr select %>% all_of
#' @importFrom stats complete.cases predict
#' @importFrom grDevices png dev.off
#' @importFrom utils globalVariables
#' @export
Rf_SDM_multiple_species <- function(Resultados_List, Variables, Registros_Minimos = 10,
                                    Stack_Raster_Reproyectado = NULL, Ruta_Almacenamiento = NULL,
                                    N_Tree = 100, Set_Seed = NULL,
                                    Cross_V_Particiones = 5, Tune_M_Try = TRUE) {
  modelos_especies <- list()
  for (m in names(Resultados_List)) {
    if(!is.null(Resultados_List[[m]]) && nrow(as.data.frame(Resultados_List[[m]])) >= Registros_Minimos) {
      if(!is.null(Set_Seed)) set.seed(Set_Seed)
      Datos_Presencias_Ausencias_Juntos_Df <- as.data.frame(Resultados_List[[m]])
      Datos_Presencias_Ausencias_Juntos_Df_Clean <- Datos_Presencias_Ausencias_Juntos_Df %>%
        dplyr::select(all_of(Variables))
      Datos_Para_Modelado <- as.data.frame(Datos_Presencias_Ausencias_Juntos_Df_Clean)
      Datos_Para_Modelado <- Datos_Para_Modelado[stats::complete.cases(Datos_Para_Modelado), ]
      if(nrow(Datos_Para_Modelado) < Registros_Minimos) next
      Datos_Para_Modelado$Presence <- as.factor(Datos_Para_Modelado$Presence)
      if(length(unique(Datos_Para_Modelado$Presence)) < 2) next
      Indices_Validacion <- caret::createDataPartition(Datos_Para_Modelado$Presence, p = 0.1, list = FALSE)
      Datos_Para_Validacion <- Datos_Para_Modelado[Indices_Validacion, ]
      Datos_Para_Modelacion <- Datos_Para_Modelado[-Indices_Validacion, ]
      Indices_Entrenamiento <- caret::createDataPartition(Datos_Para_Modelacion$Presence, p = 0.7, list = FALSE)
      Datos_Para_Entrenamiento <- Datos_Para_Modelacion[Indices_Entrenamiento, ]
      Datos_Para_Testeo <- Datos_Para_Modelacion[-Indices_Entrenamiento, ]
      mtry_optimal <- 2
      if(Tune_M_Try && nrow(Datos_Para_Entrenamiento) >= Registros_Minimos) {
        tryCatch({
          tuneRF_result <- randomForest::tuneRF(
            x = Datos_Para_Entrenamiento[, !names(Datos_Para_Entrenamiento) %in% "Presence"],
            y = Datos_Para_Entrenamiento$Presence,
            ntreeTry = min(100, N_Tree),
            stepFactor = 1.5,
            improve = 0.01,
            trace = FALSE,
            plot = FALSE
          )
          if(!is.null(tuneRF_result) && nrow(tuneRF_result) > 0) {
            mtry_optimal <- tuneRF_result[which.min(tuneRF_result[,2]), 1]
          }
        }, error = function(e) {})
      }
      if(Cross_V_Particiones > 1 && nrow(Datos_Para_Entrenamiento) > Cross_V_Particiones * 5) {
        Indices_Cross_Validation <- caret::createFolds(Datos_Para_Entrenamiento$Presence, k = Cross_V_Particiones)
        Metricas_Cross_Validation <- list()
        for(fold in 1:Cross_V_Particiones) {
          Indices_Particion_Entrenamiento <- unlist(Indices_Cross_Validation[-fold])
          Indices_Particion_Validacion <- Indices_Cross_Validation[[fold]]
          Particiones_Para_Entrenamiento <- Datos_Para_Entrenamiento[Indices_Particion_Entrenamiento, ]
          Particiones_Para_Validacion <- Datos_Para_Entrenamiento[Indices_Particion_Validacion, ]
          Modelo_Temporal_CrossValidation <- randomForest::randomForest(
            Presence ~ .,
            data = Particiones_Para_Entrenamiento,
            importance = FALSE,
            ntree = N_Tree,
            mtry = mtry_optimal,
            do.trace = FALSE,
            verbose = FALSE
          )
          Prediccion_CrossValidation <- stats::predict(Modelo_Temporal_CrossValidation, Particiones_Para_Validacion)
          Matriz_Confusion_CrossValidation <- caret::confusionMatrix(Prediccion_CrossValidation, Particiones_Para_Validacion$Presence, positive = "1")
          Metricas_Cross_Validation[[fold]] <- list(
            accuracy = Matriz_Confusion_CrossValidation$overall["Accuracy"],
            sensitivity = Matriz_Confusion_CrossValidation$byClass["Sensitivity"],
            specificity = Matriz_Confusion_CrossValidation$byClass["Specificity"],
            kappa = Matriz_Confusion_CrossValidation$overall["Kappa"]
          )
        }
        Resumen_Modelo_Temporal_CrossValidation <- data.frame(
          accuracy_mean = mean(sapply(Metricas_Cross_Validation, function(x) x$accuracy)),
          sensitivity_mean = mean(sapply(Metricas_Cross_Validation, function(x) x$sensitivity)),
          specificity_mean = mean(sapply(Metricas_Cross_Validation, function(x) x$specificity)),
          kappa_mean = mean(sapply(Metricas_Cross_Validation, function(x) x$kappa))
        )
      } else {
        Resumen_Modelo_Temporal_CrossValidation <- NULL
      }
      RF_model <- randomForest::randomForest(
        Presence ~ .,
        data = Datos_Para_Entrenamiento,
        importance = TRUE,
        ntree = N_Tree,
        mtry = mtry_optimal,
        do.trace = FALSE,
        verbose = FALSE
      )
      Predict.Train.Data <- stats::predict(RF_model, Datos_Para_Entrenamiento)
      Predict.Test.Data <- stats::predict(RF_model, Datos_Para_Testeo)
      conf_matrix_train <- caret::confusionMatrix(Predict.Train.Data,
                                                  Datos_Para_Entrenamiento$Presence)
      conf_matrix_test <- caret::confusionMatrix(Predict.Test.Data,
                                                 Datos_Para_Testeo$Presence)
      accuracy_train <- conf_matrix_train$overall["Accuracy"]
      accuracy_test <- conf_matrix_test$overall["Accuracy"]
      sensitivity_test <- conf_matrix_test$byClass["Sensitivity"]
      specificity_test <- conf_matrix_test$byClass["Specificity"]
      kappa_test <- conf_matrix_test$overall["Kappa"]
      Predict.Test.Prob <- stats::predict(RF_model, Datos_Para_Testeo, type = "prob")
      auc_value <- NA
      roc_curve <- NULL
      if(length(unique(Datos_Para_Testeo$Presence)) > 1) {
        roc_curve <- pROC::roc(response = Datos_Para_Testeo$Presence,
                               predictor = Predict.Test.Prob[, "1"])
        auc_value <- pROC::auc(roc_curve)
      }
      Predict.Validation.Data <- stats::predict(RF_model, Datos_Para_Validacion)
      conf_matrix_val <- caret::confusionMatrix(Predict.Validation.Data, Datos_Para_Validacion$Presence)
      precision_test <- conf_matrix_test$byClass["Pos Pred Value"]
      recall_test <- conf_matrix_test$byClass["Sensitivity"]
      f1_test <- if(!is.na(precision_test) && !is.na(recall_test) && (precision_test + recall_test) > 0) {
        2 * (precision_test * recall_test) / (precision_test + recall_test)
      } else NA
      modelo_info <- list(
        modelo = RF_model,
        importance = as.data.frame(RF_model$importance),
        Resumen_Modelo_Temporal_CrossValidation = Resumen_Modelo_Temporal_CrossValidation,
        accuracy_train = accuracy_train,
        accuracy_test = accuracy_test,
        sensitivity = sensitivity_test,
        specificity = specificity_test,
        precision = precision_test,
        recall = recall_test,
        f1_score = f1_test,
        kappa = kappa_test,
        auc = auc_value,
        conf_matrix_train = conf_matrix_train,
        conf_matrix_test = conf_matrix_test,
        conf_matrix_val = conf_matrix_val,
        roc_curve = roc_curve,
        mtry_used = mtry_optimal,
        ntree_used = N_Tree
      )
      if(!is.null(Ruta_Almacenamiento)) {
        dir.create(Ruta_Almacenamiento, showWarnings = FALSE, recursive = TRUE)
        saveRDS(RF_model, file.path(Ruta_Almacenamiento, paste0("modelo_randomforest_", m, ".rds")))
        saveRDS(Datos_Para_Validacion, file.path(Ruta_Almacenamiento, paste0("datos_validacion_final_", m, ".rds")))
        png(file.path(Ruta_Almacenamiento, paste0("var_importance_", m, ".png")),
            width = 800, height = 600)
        randomForest::varImpPlot(RF_model)
        dev.off()
        if(!is.null(roc_curve)) {
          png(file.path(Ruta_Almacenamiento, paste0("roc_curve_", m, ".png")),
              width = 600, height = 600)
          plot(roc_curve)
          dev.off()
        }
      }
      if(!is.null(Stack_Raster_Reproyectado)) {
        raster_pred <- terra::predict(Stack_Raster_Reproyectado, RF_model,
                                      na.rm = TRUE, type = "prob")
        if(!is.null(Ruta_Almacenamiento)) {
          terra::writeRaster(raster_pred[[2]],
                             file.path(Ruta_Almacenamiento, paste0("prediccion_prob_", m, ".tif")),
                             overwrite = TRUE)
        }
        modelo_info$raster_prediction <- raster_pred[[2]]
      }
      modelos_especies[[m]] <- modelo_info
    }
  }
  return(modelos_especies)
}
