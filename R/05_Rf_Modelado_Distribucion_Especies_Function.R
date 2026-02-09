#' @title Modelado de Distribucion de Especies con Random Forest Multi-Especies / Species Distribution Modeling with Random Forest for Multiple Species
#'
#' @description
#' Entrena modelos Random Forest para modelado de distribucion de especies (SDM) para
#' multiples especies en paralelo. Realiza particion de datos, evaluacion del modelo
#' y almacenamiento de resultados. La funcion incluye validacion cruzada y calculo
#' de metricas de desempeno.
#' /
#' Trains Random Forest models for Species Distribution Modeling (SDM) for
#' multiple species in parallel. Performs data splitting, model evaluation
#' and result storage. The function includes cross-validation and calculation
#' of performance metrics.
#'
#' @param Resultados_List Lista nombrada donde cada elemento corresponde a una especie.
#'   Cada elemento debe ser un data.frame con datos de presencia/ausencia, incluyendo
#'   una columna "Presence" (factor o numerica) y las variables predictoras.
#'   / Named list where each element corresponds to a species.
#'   Each element must be a data.frame with presence/absence data, including
#'   a "Presence" column (factor or numeric) and predictor variables.
#'
#' @param variables Vector de caracteres con los nombres de las variables predictoras
#'   a incluir en el modelo. Deben corresponder a columnas en los data.frames.
#'   / Character vector with names of predictor variables
#'   to include in the model. Must correspond to columns in the data.frames.
#'
#' @param Registros_Minimos Numero minimo de registros requeridos para modelar una especie.
#'   Especies con menos registros son omitidas. Por defecto es 10.
#'   / Minimum number of records required to model a species.
#'   Species with fewer records are omitted. Default is 10.
#'
#' @return
#' Retorna una lista nombrada (por especie) donde cada elemento contiene:
#' * modelo: Objeto Random Forest entrenado
#' * importance: Data.frame con importancia de variables
#' * accuracy_train: Exactitud en datos de entrenamiento
#' * accuracy_test: Exactitud en datos de prueba
#' * sensitivity: Sensibilidad en datos de prueba
#' * specificity: Especificidad en datos de prueba
#' * kappa: Estadistico Kappa en datos de prueba
#' * auc: Area bajo la curva ROC (si aplicable)
#' * conf_matrix_train: Matriz de confusion entrenamiento
#' * conf_matrix_test: Matriz de confusion prueba
#' * conf_matrix_val: Matriz de confusion validacion
#' * roc_curve: Objeto curva ROC (si aplicable)
#' Adicionalmente, guarda en disco los modelos (.rds) y datos de validacion.
#' /
#' Returns a named list (by species) where each element contains:
#' * modelo: Trained Random Forest object
#' * importance: Data.frame with variable importance
#' * accuracy_train: Accuracy on training data
#' * accuracy_test: Accuracy on test data
#' * sensitivity: Sensitivity on test data
#' * specificity: Specificity on test data
#' * kappa: Kappa statistic on test data
#' * auc: Area under ROC curve (if applicable)
#' * conf_matrix_train: Confusion matrix training
#' * conf_matrix_test: Confusion matrix test
#' * conf_matrix_val: Confusion matrix validation
#' * roc_curve: ROC curve object (if applicable)
#' Additionally, saves models (.rds) and validation data to disk.
#'
#' @details
#' El flujo de procesamiento para cada especie es:
#' 1. Filtrado: Solo especies con >= Registros_Minimos registros
#' 2. Limpieza: Seleccion variables y eliminacion de NA
#' 3. Validacion: 10% datos para validacion final
#' 4. Particion: 70% entrenamiento, 30% prueba del 90% restante
#' 5. Modelado: Random Forest con 50 arboles y mtry=2
#' 6. Evaluacion: Metricas en entrenamiento, prueba y validacion
#' 7. Almacenamiento: Modelos y datos guardados en "ResultadosModelacionRf_Todo/"
#' La funcion omite especies con:
#' * Menos de 5 registros despues de limpieza
#' * Una sola clase en la variable Presence
#' /
#' The processing flow for each species is:
#' 1. Filtering: Only species with >= Registros_Minimos records
#' 2. Cleaning: Variable selection and NA removal
#' 3. Validation: 10% data for final validation
#' 4. Splitting: 70% training, 30% testing from remaining 90%
#' 5. Modeling: Random Forest with 50 trees and mtry=2
#' 6. Evaluation: Metrics on training, testing and validation
#' 7. Storage: Models and data saved in "ResultadosModelacionRf_Todo/"
#' The function skips species with:
#' * Fewer than 5 records after cleaning
#' * Only one class in the Presence variable
#'
#' @note
#' * La columna "Presence" debe estar presente en todos los data.frames
#' * Las variables predictoras deben ser numericas o factoriales
#' * La carpeta "ResultadosModelacionRf_Todo" se crea en el directorio de trabajo
#' * Para calcular AUC, se necesitan ambas clases en datos de prueba
#' * Los modelos usan parametros fijos (ntree=50, mtry=2) para reproducibilidad
#' /
#' * The "Presence" column must be present in all data.frames
#' * Predictor variables must be numeric or factor
#' * Folder "ResultadosModelacionRf_Todo" is created in working directory
#' * To calculate AUC, both classes are needed in test data
#' * Models use fixed parameters (ntree=50, mtry=2) for reproducibility
#'
#' @section Advertencias/Warnings:
#' * El proceso puede ser lento con muchas especies o muchos registros
#' * El muestreo estratificado requiere suficientes datos de cada clase
#' * Metricas pueden ser sesgadas con datos desbalanceados
#' * La validacion cruzada simple puede no ser optima para SDM
#' /
#' * Process can be slow with many species or many records
#' * Stratified sampling requires sufficient data from each class
#' * Metrics may be biased with unbalanced data
#' * Simple cross-validation may not be optimal for SDM
#'
#' @examples
#' \dontrun{
#' # Crear datos de ejemplo / Create example data
#' set.seed(123)
#' datos_especie1 <- data.frame(
#'   Presence = as.factor(sample(c(0,1), 50, replace = TRUE)),
#'   temp = rnorm(50, 20, 5),
#'   prec = rnorm(50, 1000, 200),
#'   elev = runif(50, 0, 2000)
#' )
#'
#' datos_especie2 <- data.frame(
#'   Presence = as.factor(sample(c(0,1), 30, replace = TRUE)),
#'   temp = rnorm(30, 18, 4),
#'   prec = rnorm(30, 1200, 150),
#'   elev = runif(30, 500, 2500)
#' )
#'
#' lista_especies <- list(
#'   Especie_A = datos_especie1,
#'   Especie_B = datos_especie2
#' )
#'
#' variables_modelo <- c("temp", "prec", "elev")
#'
#' # Ejecutar modelado / Run modeling
#' resultados <- Rf_SDM_multiple_species(
#'   Resultados_List = lista_especies,
#'   variables = variables_modelo,
#'   Registros_Minimos = 20
#' )
#'
#' # Examinar resultados / Examine results
#' names(resultados)
#' resultados$Especie_A$accuracy_test
#' resultados$Especie_A$importance
#'
#' # Ver modelos guardados / View saved models
#' list.files("ResultadosModelacionRf_Todo/", pattern = ".rds")
#' }
#'
#' @seealso
#' * \code{\link[randomForest]{randomForest}} para detalles del algoritmo
#' * \code{\link[caret]{createDataPartition}} para particion estratificada
#' * \code{\link[pROC]{roc}} para calculo de curvas ROC
#' /
#' * \code{\link[randomForest]{randomForest}} for algorithm details
#' * \code{\link[caret]{createDataPartition}} for stratified splitting
#' * \code{\link[pROC]{roc}} for ROC curve calculation
#'
#' @keywords
#'   modelado distribucion especies, random forest, aprendizaje automatico,
#'   presencia-ausencia, evaluacion modelo, SDM, machine learning,
#'   species distribution modeling, random forest, machine learning,
#'   presence-absence, model evaluation, SDM
#'
#' @family
#'   funciones_modelado_especies, funciones_random_forest,
#'   species_modeling_functions, random_forest_functions
#'
#' @importFrom randomForest randomForest
#' @importFrom caret createDataPartition confusionMatrix
#' @importFrom pROC roc auc
#' @importFrom dplyr select %>%
#' @importFrom stats complete.cases predict
#' @export
Rf_SDM_multiple_species <- function(Resultados_List, variables, Registros_Minimos = 10){
  modelos_especies <- list()
  for (m in names(Resultados_List)) {
    if(!is.null(Resultados_List[[m]]) && nrow(as.data.frame(Resultados_List[[m]])) >= Registros_Minimos) {
      Datos_Presencias_Ausencias_Juntos_Df <- as.data.frame(Resultados_List[[m]])
      Datos_Presencias_Ausencias_Juntos_Df_Clean <- Datos_Presencias_Ausencias_Juntos_Df %>%
        dplyr::select(all_of(variables))
      Datos_Para_Modelado <- as.data.frame(Datos_Presencias_Ausencias_Juntos_Df_Clean)
      Datos_Para_Modelado <- Datos_Para_Modelado[stats::complete.cases(Datos_Para_Modelado), ]
      if(nrow(Datos_Para_Modelado) <= 4) next
      Datos_Para_Modelado$Presence <- as.factor(Datos_Para_Modelado$Presence)
      if(length(unique(Datos_Para_Modelado$Presence)) < 2) next
      Indices_Validacion <- caret::createDataPartition(Datos_Para_Modelado$Presence, p = 0.1, list = FALSE)
      Datos_Para_Validacion <- Datos_Para_Modelado[Indices_Validacion, ]
      Datos_Para_Modelacion <- Datos_Para_Modelado[-Indices_Validacion, ]
      Indices_Entrenamiento <- caret::createDataPartition(Datos_Para_Modelacion$Presence, p = 0.7, list = FALSE)
      Datos_Para_Entrenamiento <- Datos_Para_Modelacion[Indices_Entrenamiento, ]
      Datos_Para_Testeo <- Datos_Para_Modelacion[-Indices_Entrenamiento, ]
      RF_model <- randomForest::randomForest(Presence ~ .,
                               data = Datos_Para_Entrenamiento,
                               importance = TRUE,
                               ntree = 50,
                               mtry = 2,
                               do.trace = FALSE,
                               verbose = FALSE)
      Predict.Train.Data <- stats::predict(RF_model, Datos_Para_Entrenamiento)
      Predict.Test.Data <- stats::predict(RF_model, Datos_Para_Testeo)
      conf_matrix_train <- caret::confusionMatrix(Predict.Train.Data, Datos_Para_Entrenamiento$Presence)
      conf_matrix_test <- caret::confusionMatrix(Predict.Test.Data, Datos_Para_Testeo$Presence)
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
      Predict.Validation.Data <- predict(RF_model, Datos_Para_Validacion)
      conf_matrix_val <- caret::confusionMatrix(Predict.Validation.Data, Datos_Para_Validacion$Presence)
      modelo_info <- list(
        modelo = RF_model,
        importance = as.data.frame(RF_model$importance),
        accuracy_train = accuracy_train,
        accuracy_test = accuracy_test,
        sensitivity = sensitivity_test,
        specificity = specificity_test,
        kappa = kappa_test,
        auc = auc_value,
        conf_matrix_train = conf_matrix_train,
        conf_matrix_test = conf_matrix_test,
        conf_matrix_val = conf_matrix_val,
        roc_curve = roc_curve
      )
      modelos_especies[[m]] <- modelo_info
      ruta_almacenamiento <- file.path(getwd(), "ResultadosModelacionRf_Todo")
      # modelo y datos de validacion en disco
      saveRDS(RF_model, paste0(ruta_almacenamiento, "modelo_randomforest_", m, ".rds"))
      saveRDS(Datos_Para_Validacion, paste0(ruta_almacenamiento, "datos_validacion_final_", m, ".rds"))
    }
  }
  return(modelos_especies)
}
