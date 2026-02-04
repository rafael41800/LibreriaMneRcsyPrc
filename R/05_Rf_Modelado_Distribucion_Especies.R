#' Modelado de distribucion de especies con Random Forest
#'
#' Esta funcion realiza modelado de distribucion de especies (SDM) utilizando
#' el algoritmo Random Forest para multiples especies simultaneamente.
#'
#' @note Esta funcion utiliza los paquetes: caret para particion de datos y matrices
#'   de confusion, randomForest para el algoritmo de clasificacion, pROC para curvas
#'   ROC, y dplyr para manipulacion de datos.
#'
#' @param Data_Df Dataframe con los datos generales (no utilizado directamente).
#' @param variables Vector de caracteres con los nombres de las variables predictoras
#'   que seran utilizadas en el modelo.
#' @param Resultados_List Lista con los datos de presencias y ausencias para cada especie.
#'   Cada elemento de la lista debe ser un dataframe con columnas para las variables
#'   predictoras y la variable respuesta 'Presence' (1 para presencia, 0 para ausencia).
#' @param Raster_Tif_reproyectado Objeto raster con las variables ambientales ya
#'   reproyectadas al sistema de coordenadas de referencia.
#' @param Registros_Minimos Numero minimo de registros que debe tener una especie
#'   para ser modelada (por defecto 5).
#'
#' @returns Retorna una lista donde cada elemento corresponde a una especie y contiene:
#'   \item{modelo}{Objeto randomForest entrenado}
#'   \item{importance}{Importancia de las variables}
#'   \item{accuracy_train}{Accuracy en datos de entrenamiento}
#'   \item{accuracy_test}{Accuracy en datos de prueba}
#'   \item{sensitivity}{Sensibilidad del modelo}
#'   \item{specificity}{Especificidad del modelo}
#'   \item{kappa}{Coeficiente Kappa}
#'   \item{auc}{Area bajo la curva ROC (AUC)}
#'   \item{conf_matrix_train}{Matriz de confusion de entrenamiento}
#'   \item{conf_matrix_test}{Matriz de confusion de prueba}
#'   \item{conf_matrix_val}{Matriz de confusion de validacion}
#'   \item{roc_curve}{Objeto de curva ROC}
#'
#' @importFrom caret createDataPartition confusionMatrix
#' @importFrom randomForest randomForest
#' @importFrom pROC roc auc
#' @importFrom dplyr select %>%
#'
#' @examples
#' \dontrun{
#' # Ejemplo de uso
#' resultados_modelos <- Rf_SDM_multiples_especies(
#'   Data_Df = datos_generales,
#'   variables = c("temperatura", "precipitacion", "altitud"),
#'   Resultados = lista_presencias_ausencias,
#'   Raster_Tif_reproyectado = raster_variables,
#'   Registros_Minimos = 10
#' )
#'
#' # Acceder a resultados de la primera especie
#' resultados_modelos[[1]]$accuracy_test
#' resultados_modelos[[1]]$auc
#' }
#' @export
Rf_SDM_multiples_especies <- function(Data_Df, variables, Resultados_List, Raster_Tif_reproyectado, Registros_Minimos = 5){
  modelos_especies <- list()
  for (m in names(Resultados_List)) {
    if(!is.null(Resultados_List[[m]]) && nrow(as.data.frame(Resultados_List[[m]])) >= Registros_Minimos) {
      Datos_Presencias_Ausencias_Juntos_Df <- as.data.frame(Resultados_List[[m]])
      Datos_Presencias_Ausencias_Juntos_Df_Clean <- Datos_Presencias_Ausencias_Juntos_Df %>%
        dplyr::select(all_of(variables))

      Datos_Para_Modelado <- as.data.frame(Datos_Presencias_Ausencias_Juntos_Df_Clean)
      Datos_Para_Modelado <- Datos_Para_Modelado[complete.cases(Datos_Para_Modelado), ]
      if(nrow(Datos_Para_Modelado) <= 4) next
      Datos_Para_Modelado$Presence <- as.factor(Datos_Para_Modelado$Presence)

      if(length(unique(Datos_Para_Modelado$Presence)) < 2) next

      Indices_Validacion <- caret::createDataPartition(Datos_Para_Modelado$Presence, p = 0.1, list = FALSE)
      Datos_Para_Validacion <- Datos_Para_Modelado[Indices_Validacion, ]
      Datos_Para_Modelacion <- Datos_Para_Modelado[-Indices_Validacion, ]

      Indices_Entrenamiento <- caret::createDataPartition(Datos_Para_Modelacion$Presence, p = 0.7, list = FALSE)
      Datos_Para_Entrenamiento <- Datos_Para_Modelacion[Indices_Entrenamiento, ]
      Datos_Para_Testeo <- Datos_Para_Modelacion[-Indices_Entrenamiento, ]

      RF_model <- randomForest(Presence ~ .,
                               data = Datos_Para_Entrenamiento,
                               importance = TRUE,
                               ntree = 50,
                               mtry = 2,
                               do.trace = FALSE,
                               verbose = FALSE)

      Predict.Train.Data <- predict(RF_model, Datos_Para_Entrenamiento)
      Predict.Test.Data <- predict(RF_model, Datos_Para_Testeo)

      conf_matrix_train <- confusionMatrix(Predict.Train.Data, Datos_Para_Entrenamiento$Presence)
      conf_matrix_test <- confusionMatrix(Predict.Test.Data, Datos_Para_Testeo$Presence)

      accuracy_train <- conf_matrix_train$overall["Accuracy"]
      accuracy_test <- conf_matrix_test$overall["Accuracy"]
      sensitivity_test <- conf_matrix_test$byClass["Sensitivity"]
      specificity_test <- conf_matrix_test$byClass["Specificity"]
      kappa_test <- conf_matrix_test$overall["Kappa"]

      Predict.Test.Prob <- predict(RF_model, Datos_Para_Testeo, type = "prob")

      auc_value <- NA
      roc_curve <- NULL
      if(length(unique(Datos_Para_Testeo$Presence)) > 1) {
        roc_curve <- roc(response = Datos_Para_Testeo$Presence,
                         predictor = Predict.Test.Prob[, "1"])
        auc_value <- auc(roc_curve)
      }

      Predict.Validation.Data <- predict(RF_model, Datos_Para_Validacion)
      conf_matrix_val <- confusionMatrix(Predict.Validation.Data, Datos_Para_Validacion$Presence)

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

      ruta_almacenamiento <- "C:/Users/Admin/Documents/RafaelChavez/Carpeta_inicial_tesis/codigo_R_tesis/Resultados_Rf"
      saveRDS(RF_model, paste0(ruta_almacenamiento, "modelo_randomforest_", m, ".rds"))
      saveRDS(Datos_Para_Validacion, paste0(ruta_almacenamiento, "datos_validacion_final_", m, ".rds"))
    }
  }
  return(modelos_especies)
}
