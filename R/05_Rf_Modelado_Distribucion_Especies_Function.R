#' @title Modelado de distribucion de especies con Random Forest / Species Distribution Modeling with Random Forest
#'
#' @description
#' Realiza modelado de distribucion de especies (SDM, por sus siglas en ingles)
#' utilizando el algoritmo Random Forest para multiples especies simultaneamente.
#' La funcion implementa un pipeline completo de modelado que incluye particion
#' de datos en entrenamiento, prueba y validacion, entrenamiento de modelos y
#' evaluacion de metricas para cada especie.
#' /
#' Performs Species Distribution Modeling (SDM) using the Random Forest algorithm
#' for multiple species simultaneously. The function implements a complete modeling
#' pipeline including data partitioning into training, testing, and validation sets,
#' model training, and metric evaluation for each species.
#'
#' @param Data_Df Dataframe con los datos generales del estudio. Aunque no se
#'   utiliza directamente en los modelos, proporciona contexto y puede ser usado
#'   para validaciones o referencias cruzadas. Debe contener al menos las
#'   variables ambientales utilizadas en el modelado.
#'   /
#'   Dataframe with general study data. Although not directly used in models,
#'   it provides context and can be used for validations or cross-references.
#'   Must contain at least the environmental variables used in modeling.
#'
#' @param variables Vector de caracteres con los nombres de las variables
#'   predictoras que seran utilizadas en los modelos. Estas variables deben:
#'   * Estar presentes en todos los dataframes de `Resultados_List`
#'   * Corresponder a capas en `Raster_Tif_reproyectado`
#'   * No contener valores NA (deben haber sido limpiados previamente)
#'   * Ser numericas para Random Forest
#'   /
#'   Character vector with names of predictor variables to be used in models.
#'   These variables must:
#'   * Be present in all dataframes of `Resultados_List`
#'   * Correspond to layers in `Raster_Tif_reproyectado`
#'   * Not contain NA values (should have been cleaned previously)
#'   * Be numeric for Random Forest
#'
#' @param Resultados_List Lista con los datos de presencias y pseudoausencias
#'   para cada especie. Cada elemento debe ser un dataframe con:
#'   * Una fila por registro (presencia o pseudoausencia)
#'   * Columnas para todas las variables especificadas en `variables`
#'   * Columna `Presence` con valores `1` (presencia) y `0` (pseudoausencia)
#'   * Geometria espacial si es objeto sf (sera convertida a dataframe)
#'   /
#'   List with presence and pseudoabsence data for each species. Each element
#'   must be a dataframe with:
#'   * One row per record (presence or pseudoabsence)
#'   * Columns for all variables specified in `variables`
#'   * `Presence` column with values `1` (presence) and `0` (pseudoabsence)
#'   * Spatial geometry if it's an sf object (will be converted to dataframe)
#'
#' @param Raster_Tif_reproyectado Objeto `SpatRaster` de terra con las variables
#'   ambientales ya reproyectadas al sistema de coordenadas de referencia. Se
#'   utiliza para:
#'   * Validar que las variables en los datos correspondan a capas del raster
#'   * Posiblemente para predicciones espaciales futuras (aunque no en esta funcion)
#'   * Referencia de resolucion y extension espacial
#'   /
#'   `SpatRaster` object from terra with environmental variables already
#'   reprojected to the reference coordinate system. It is used for:
#'   * Validating that variables in data correspond to raster layers
#'   * Possibly for future spatial predictions (though not in this function)
#'   * Reference of spatial resolution and extent
#'
#' @param Registros_Minimos Numero minimo de registros que debe tener una especie
#'   para ser modelada. Incluye tanto presencias como pseudoausencias.
#'   Por defecto es 10. Especies con menos registros seran omitidas.
#'   Este umbral asegura modelos estadisticamente robustos.
#'   /
#'   Minimum number of records a species must have to be modeled. Includes
#'   both presences and pseudoabsences. Default is 10. Species with fewer
#'   records will be omitted. This threshold ensures statistically robust models.
#'
#' @return
#' Retorna una lista nombrada (por especie) donde cada elemento es otra lista
#' con los siguientes componentes:
#'
#' \itemize{
#'   \item{\strong{modelo}: Objeto `randomForest` entrenado para la especie}
#'   \item{\strong{importance}: Dataframe con la importancia de cada variable,
#'         incluyendo MeanDecreaseAccuracy y MeanDecreaseGini}
#'   \item{\strong{accuracy_train}: Accuracy (exactitud) en datos de entrenamiento
#'         (proporcion de predicciones correctas)}
#'   \item{\strong{accuracy_test}: Accuracy en datos de prueba}
#'   \item{\strong{sensitivity}: Sensibilidad (true positive rate) - capacidad
#'         para detectar presencias correctamente}
#'   \item{\strong{specificity}: Especificidad (true negative rate) - capacidad
#'         para detectar ausencias correctamente}
#'   \item{\strong{kappa}: Coeficiente Kappa de Cohen, ajusta el accuracy por
#'         acuerdo por azar (valores >0.6 indican buen acuerdo)}
#'   \item{\strong{auc}: Area bajo la curva ROC (AUC), medida de capacidad
#'         discriminativa (1 = perfecto, 0.5 = aleatorio). Puede ser NA si
#'         solo hay una clase en los datos de prueba.}
#'   \item{\strong{conf_matrix_train}: Matriz de confusion para datos de
#'         entrenamiento (objeto de caret)}
#'   \item{\strong{conf_matrix_test}: Matriz de confusion para datos de prueba}
#'   \item{\strong{conf_matrix_val}: Matriz de confusion para datos de
#'         validacion (10% de los datos originales)}
#'   \item{\strong{roc_curve}: Objeto ROC (de pROC) para analisis detallado
#'         de la curva, o NULL si no se pudo calcular}
#' }
#'
#' Ademas, la funcion guarda en disco dos archivos por especie:
#' 1. El modelo Random Forest entrenado
#' 2. Los datos de validacion utilizados
#' Los archivos se guardan en la ruta: `C:/Users/Admin/Documents/RafaelChavez/Carpeta_inicial_tesis/codigo_R_tesis/Resultados_Rf/`
#' /
#' Returns a named list (by species) where each element is another list
#' with the following components:
#'
#' \itemize{
#'   \item{\strong{modelo}: Trained `randomForest` object for the species}
#'   \item{\strong{importance}: Dataframe with importance of each variable,
#'         including MeanDecreaseAccuracy and MeanDecreaseGini}
#'   \item{\strong{accuracy_train}: Accuracy on training data
#'         (proportion of correct predictions)}
#'   \item{\strong{accuracy_test}: Accuracy on test data}
#'   \item{\strong{sensitivity}: Sensitivity (true positive rate) - ability
#'         to correctly detect presences}
#'   \item{\strong{specificity}: Specificity (true negative rate) - ability
#'         to correctly detect absences}
#'   \item{\strong{kappa}: Cohen's Kappa coefficient, adjusts accuracy for
#'         chance agreement (values >0.6 indicate good agreement)}
#'   \item{\strong{auc}: Area under the ROC curve (AUC), measure of
#'         discriminative ability (1 = perfect, 0.5 = random). May be NA if
#'         only one class is present in test data.}
#'   \item{\strong{conf_matrix_train}: Confusion matrix for training data
#'         (caret object)}
#'   \item{\strong{conf_matrix_test}: Confusion matrix for test data}
#'   \item{\strong{conf_matrix_val}: Confusion matrix for validation data
#'         (10% of original data)}
#'   \item{\strong{roc_curve}: ROC object (from pROC) for detailed curve
#'         analysis, or NULL if not calculable}
#' }
#'
#' Additionally, the function saves two files per species to disk:
#' 1. The trained Random Forest model
#' 2. The validation data used
#' Files are saved to: `C:/Users/Admin/Documents/RafaelChavez/Carpeta_inicial_tesis/codigo_R_tesis/Resultados_Rf/`
#'
#' @details
#' La funcion implementa el siguiente pipeline de modelado para cada especie:
#'
#' 1. **Preprocesamiento y validacion**:
#'    - Filtra especies con suficientes registros (`Registros_Minimos` = 10)
#'    - Convierte objetos `sf` a dataframes si es necesario
#'    - Selecciona solo las variables especificadas
#'    - Elimina filas con valores NA
#'
#' 2. **Particion triple de datos**:
#'    - **Validacion**: 10% de los datos (usando `createDataPartition`)
#'    - **Modelacion**: 90% restante dividido en:
#'      - Entrenamiento: 70% de los datos originales (77.8% de los de modelacion)
#'      - Prueba: 30% de los datos originales (22.2% de los de modelacion)
#'
#' 3. **Entrenamiento del modelo**:
#'    - Configura Random Forest con parametros especificos:
#'      - `ntree = 50` (arboles)
#'      - `mtry = 2` (variables por split)
#'      - `importance = TRUE` (calcula importancia de variables)
#'
#' 4. **Evaluacion del modelo**:
#'    - Predicciones en datos de entrenamiento, prueba y validacion
#'    - Calculo de matrices de confusion para los tres conjuntos
#'    - Calculo de metricas: accuracy, sensibilidad, especificidad, Kappa
#'    - Generacion de curva ROC y calculo de AUC (si hay ambas clases en prueba)
#'
#' 5. **Guardado de resultados**:
#'    - Modelo entrenado guardado como RDS
#'    - Datos de validacion guardados como RDS
#'    - Metricas compiladas en lista de retorno
#'
#' /
#' The function implements the following modeling pipeline for each species:
#'
#' 1. **Preprocessing and validation**:
#'    - Filters species with enough records (`Registros_Minimos` = 10)
#'    - Converts `sf` objects to dataframes if necessary
#'    - Selects only specified variables
#'    - Removes rows with NA values
#'
#' 2. **Triple data partitioning**:
#'    - **Validation**: 10% of data (using `createDataPartition`)
#'    - **Modeling**: Remaining 90% divided into:
#'      - Training: 70% of original data (77.8% of modeling data)
#'      - Testing: 30% of original data (22.2% of modeling data)
#'
#' 3. **Model training**:
#'    - Configures Random Forest with specific parameters:
#'      - `ntree = 50` (trees)
#'      - `mtry = 2` (variables per split)
#'      - `importance = TRUE` (calculates variable importance)
#'
#' 4. **Model evaluation**:
#'    - Predictions on training, testing, and validation data
#'    - Calculation of confusion matrices for all three sets
#'    - Calculation of metrics: accuracy, sensitivity, specificity, Kappa
#'    - Generation of ROC curve and AUC calculation (if both classes in test data)
#'
#' 5. **Results saving**:
#'    - Trained model saved as RDS
#'    - Validation data saved as RDS
#'    - Metrics compiled in return list
#'
#' @note
#' * **Parametros fijos**: Los parametros de Random Forest estan fijos en:
#'   - `ntree = 50` (menos que el default de 500 para velocidad)
#'   - `mtry = 2` (solo 2 variables por split)
#'   Considere ajustar estos valores segun sus necesidades.
#' * **Guardado automatico**: Los modelos y datos se guardan automaticamente
#'   en la ruta especificada en el codigo. Esta ruta es fija y no configurable
#'   desde los parametros.
#' * **Validacion triple**: Se usan tres conjuntos: entrenamiento (63%),
#'   prueba (27%) y validacion (10%).
#' * **Filtrado estricto**: Especies con menos de 10 registros o con solo
#'   una clase despues de limpieza son omitidas sin advertencia.
#' * **Reproducibilidad**: Use `set.seed()` antes de llamar la funcion para
#'   resultados reproducibles en la particion de datos.
#'
#' /
#' * **Fixed parameters**: Random Forest parameters are fixed at:
#'   - `ntree = 50` (less than default 500 for speed)
#'   - `mtry = 2` (only 2 variables per split)
#'   Consider adjusting these values according to your needs.
#' * **Automatic saving**: Models and data are automatically saved to the
#'   path specified in the code. This path is fixed and not configurable
#'   from parameters.
#' * **Triple validation**: Three sets are used: training (63%),
#'   testing (27%), and validation (10%).
#' * **Strict filtering**: Species with fewer than 10 records or with only
#'   one class after cleaning are omitted without warning.
#' * **Reproducibility**: Use `set.seed()` before calling the function for
#'   reproducible results in data partitioning.
#'
#' @section Advertencias/Warnings:
#' * **Ruta fija de guardado**: La ruta de guardado esta hardcodeada y puede
#'   no existir en otros sistemas. Verifique o cree la carpeta antes de usar.
#' * **Parametros no optimizados**: Los parametros de Random Forest no se
#'   optimizan automaticamente. Para mejor rendimiento, considere tuning.
#' * **Sin seleccion de variables**: Todas las variables especificadas se
#'   usan sin evaluacion previa de correlacion o importancia.
#' * **AUC puede faltar**: Si los datos de prueba tienen solo una clase,
#'   no se calcula AUC y el valor sera NA.
#' * **Sin metadatos en retorno**: La lista retornada no incluye metadatos
#'   como numero de registros o variables usadas.
#'
#' /
#' * **Fixed save path**: The save path is hardcoded and may not exist on
#'   other systems. Verify or create the folder before using.
#' * **Non-optimized parameters**: Random Forest parameters are not
#'   automatically optimized. For better performance, consider tuning.
#' * **No variable selection**: All specified variables are used without
#'   prior evaluation of correlation or importance.
#' * **AUC may be missing**: If test data has only one class,
#'   AUC is not calculated and value will be NA.
#' * **No metadata in return**: The returned list does not include metadata
#'   such as number of records or variables used.
#'
#' @examples
#' \dontrun{
#' # Ejemplo 1: Uso basico con datos simulados
#' # Example 1: Basic usage with simulated data
#'
#' # Cargar paquetes requeridos
#' # Load required packages
#' library(randomForest)
#' library(caret)
#' library(pROC)
#' library(dplyr)
#'
#' # Crear datos de ejemplo para multiples especies
#' # Create example data for multiple species
#' set.seed(123)
#' especies <- c("Especie_A", "Especie_B", "Especie_C")
#' lista_datos <- list()
#'
#' for (especie in especies) {
#'   n <- 100
#'   datos <- data.frame(
#'     temperatura = rnorm(n, 20, 5),
#'     precipitacion = rnorm(n, 1000, 200),
#'     altitud = rnorm(n, 500, 100),
#'     humedad = rnorm(n, 70, 10),
#'     Presence = sample(c(0, 1), n, replace = TRUE, prob = c(0.7, 0.3)),
#'     especie = especie
#'   )
#'   lista_datos[[especie]] <- datos
#' }
#'
#' # Crear raster de referencia (solo para validacion)
#' # Create reference raster (for validation only)
#' library(terra)
#' raster_ref <- rast(nrows = 10, ncols = 10,
#'                    xmin = 0, xmax = 1, ymin = 0, ymax = 1)
#' values(raster_ref) <- matrix(1:100, 10, 10)
#' names(raster_ref) <- c("temperatura")
#'
#' # Ejecutar modelado para multiples especies
#' # Run modeling for multiple species
#' resultados_modelos <- Rf_SDM_multiple_species(
#'   Data_Df = data.frame(),  # Datos generales (puede estar vacio para ejemplo)
#'   variables = c("temperatura", "precipitacion", "altitud", "humedad"),
#'   Resultados_List = lista_datos,
#'   Raster_Tif_reproyectado = raster_ref,
#'   Registros_Minimos = 10
#' )
#'
#' # Acceder a resultados de una especie especifica
#' # Access results for a specific species
#' resultados_especie_a <- resultados_modelos[["Especie_A"]]
#'
#' # Ver metricas clave
#' # View key metrics
#' cat("Accuracy test:", resultados_especie_a$accuracy_test, "\n")
#' cat("AUC:", resultados_especie_a$auc, "\n")
#' cat("Kappa:", resultados_especie_a$kappa, "\n")
#'
#' # Ver importancia de variables
#' # View variable importance
#' print(resultados_especie_a$importance)
#'
#' # Los modelos tambien estan guardados en disco:
#' # Models are also saved to disk:
#' # C:/Users/Admin/Documents/RafaelChavez/Carpeta_inicial_tesis/codigo_R_tesis/Resultados_Rf/
#' }
#'
#' @seealso
#'   \itemize{
#'     \item \code{\link[randomForest]{randomForest}} para detalles del algoritmo de Random Forest /
#'           \code{\link[randomForest]{randomForest}} for details on the Random Forest algorithm
#'     \item \code{\link[caret]{confusionMatrix}} para interpretacion de matrices de confusion /
#'           \code{\link[caret]{confusionMatrix}} for interpretation of confusion matrices
#'     \item \code{\link[pROC]{roc}} para analisis detallado de curvas ROC /
#'           \code{\link[pROC]{roc}} for detailed ROC curve analysis
#'   }
#'
#' @keywords
#'   SDM, Random Forest, modelado de distribucion de especies, aprendizaje automatico,
#'   presencia-ausencia, evaluacion de modelos,
#'   SDM, Random Forest, species distribution modeling, machine learning,
#'   presence-absence, model evaluation
#'
#' @family
#'   funciones_SDM, funciones_modelado, funciones_aprendizaje_automatico /
#'   SDM_functions, modeling_functions, machine_learning_functions
#'
#' @importFrom randomForest randomForest
#' @importFrom caret createDataPartition confusionMatrix
#' @importFrom pROC roc auc
#' @importFrom dplyr select %>% all_of
#' @importFrom stats complete.cases predict
#' @export
Rf_SDM_multiple_species <- function(Data_Df, variables, Resultados_List, Raster_Tif_reproyectado, Registros_Minimos = 10){
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
