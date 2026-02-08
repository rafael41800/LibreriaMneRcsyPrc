#' @title Modelado de distribución de especies con Random Forest / Species Distribution Modeling with Random Forest
#'
#' @description
#' Realiza modelado de distribución de especies (SDM, por sus siglas en inglés)
#' utilizando el algoritmo Random Forest para múltiples especies simultáneamente.
#' La función implementa un pipeline completo de modelado que incluye partición
#' de datos, entrenamiento, validación y evaluación de modelos para cada especie.
#' /
#' Performs Species Distribution Modeling (SDM) using the Random Forest algorithm
#' for multiple species simultaneously. The function implements a complete modeling
#' pipeline including data partitioning, training, validation, and model evaluation
#' for each species.
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
#'   predictoras que serán utilizadas en los modelos. Estas variables deben:
#'   * Estar presentes en todos los dataframes de `Resultados_List`
#'   * Corresponder a capas en `Raster_Tif_reproyectado`
#'   * No contener valores NA (deben haber sido limpiados previamente)
#'   * Ser numéricas para Random Forest
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
#'   * Columna `especievalida` o similar identificando la especie
#'   * Geometría espacial si es objeto sf (será convertida a dataframe)
#'   /
#'   List with presence and pseudoabsence data for each species. Each element
#'   must be a dataframe with:
#'   * One row per record (presence or pseudoabsence)
#'   * Columns for all variables specified in `variables`
#'   * `Presence` column with values `1` (presence) and `0` (pseudoabsence)
#'   * `especievalida` or similar column identifying the species
#'   * Spatial geometry if it's an sf object (will be converted to dataframe)
#'
#' @param Raster_Tif_reproyectado Objeto `SpatRaster` de terra con las variables
#'   ambientales ya reproyectadas al sistema de coordenadas de referencia. Se
#'   utiliza para:
#'   * Validar que las variables en los datos correspondan a capas del raster
#'   * Posiblemente para predicciones espaciales futuras (aunque no en esta función)
#'   * Referencia de resolución y extensión espacial
#'   /
#'   `SpatRaster` object from terra with environmental variables already
#'   reprojected to the reference coordinate system. It is used for:
#'   * Validating that variables in data correspond to raster layers
#'   * Possibly for future spatial predictions (though not in this function)
#'   * Reference of spatial resolution and extent
#'
#' @param Registros_Minimos Número mínimo de registros que debe tener una especie
#'   para ser modelada. Incluye tanto presencias como pseudoausencias.
#'   Por defecto es 5. Especies con menos registros serán omitidas con una
#'   advertencia. Este umbral asegura modelos estadísticamente robustos.
#'   /
#'   Minimum number of records a species must have to be modeled. Includes
#'   both presences and pseudoabsences. Default is 5. Species with fewer
#'   records will be omitted with a warning. This threshold ensures
#'   statistically robust models.
#'
#' @returns
#' Retorna una lista nombrada (por especie) donde cada elemento es otra lista
#' con los siguientes componentes:
#'
#' \itemize{
#'   \item{\strong{modelo}: Objeto `randomForest` entrenado para la especie}
#'   \item{\strong{importance}: Dataframe con la importancia de cada variable,
#'         incluyendo MeanDecreaseAccuracy y MeanDecreaseGini}
#'   \item{\strong{accuracy_train}: Accuracy (exactitud) en datos de entrenamiento
#'         (proporción de predicciones correctas)}
#'   \item{\strong{accuracy_test}: Accuracy en datos de prueba}
#'   \item{\strong{sensitivity}: Sensibilidad (true positive rate) - capacidad
#'         para detectar presencias correctamente}
#'   \item{\strong{specificity}: Especificidad (true negative rate) - capacidad
#'         para detectar ausencias correctamente}
#'   \item{\strong{kappa}: Coeficiente Kappa de Cohen, ajusta el accuracy por
#'         acuerdo por azar (valores >0.6 indican buen acuerdo)}
#'   \item{\strong{auc}: Area bajo la curva ROC (AUC), medida de capacidad
#'         discriminativa (1 = perfecto, 0.5 = aleatorio)}
#'   \item{\strong{conf_matrix_train}: Matriz de confusión para datos de
#'         entrenamiento}
#'   \item{\strong{conf_matrix_test}: Matriz de confusión para datos de prueba}
#'   \item{\strong{conf_matrix_val}: Matriz de confusión para datos de
#'         validación (si aplica)}
#'   \item{\strong{roc_curve}: Objeto ROC para análisis detallado de la curva}
#'   \item{\strong{n_records}: Número total de registros utilizados}
#'   \item{\strong{n_presences}: Número de presencias en el conjunto de datos}
#'   \item{\strong{n_pseudoabsences}: Número de pseudoausencias en el conjunto}
#'   \item{\strong{train_test_split}: Proporción de datos en entrenamiento vs prueba}
#'   \item{\strong{variables_used}: Variables finalmente utilizadas en el modelo
#'         (puede diferir de `variables` si hay problemas)}
#' }
#'
#' La lista es invisible (no se imprime automáticamente) para facilitar su
#' asignación a un objeto y posterior análisis.
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
#'         discriminative ability (1 = perfect, 0.5 = random)}
#'   \item{\strong{conf_matrix_train}: Confusion matrix for training data}
#'   \item{\strong{conf_matrix_test}: Confusion matrix for test data}
#'   \item{\strong{conf_matrix_val}: Confusion matrix for validation data
#'         (if applicable)}
#'   \item{\strong{roc_curve}: ROC object for detailed curve analysis}
#'   \item{\strong{n_records}: Total number of records used}
#'   \item{\strong{n_presences}: Number of presences in the dataset}
#'   \item{\strong{n_pseudoabsences}: Number of pseudoabsences in the dataset}
#'   \item{\strong{train_test_split}: Proportion of data in training vs test}
#'   \item{\strong{variables_used}: Variables finally used in the model
#'         (may differ from `variables` if there are issues)}
#' }
#'
#' The list is invisible (not automatically printed) to facilitate assignment
#' to an object and subsequent analysis.
#'
#' @details
#' La función implementa el siguiente pipeline de modelado para cada especie:
#'
#' 1. **Preprocesamiento y validación**:
#'    - Filtra especies con suficientes registros (`Registros_Minimos`)
#'    - Convierte objetos `sf` a dataframes si es necesario
#'    - Verifica que todas las variables requeridas estén presentes
#'    - Elimina registros duplicados y verifica balance de clases
#'
#' 2. **Partición de datos**:
#'    - Divide datos en conjuntos de entrenamiento (70%) y prueba (30%)
#'    - Mantiene proporción original de presencias/ausencias en cada conjunto
#'    - Utiliza `caret::createDataPartition()` para partición estratificada
#'
#' 3. **Entrenamiento del modelo**:
#'    - Configura Random Forest con parámetros por defecto optimizados
#'    - Entrena modelo con datos de entrenamiento
#'    - Calcula importancia de variables durante el entrenamiento
#'
#' 4. **Evaluación del modelo**:
#'    - Predicciones en datos de entrenamiento y prueba
#'    - Cálculo de matrices de confusión para múltiples conjuntos
#'    - Cálculo de métricas: accuracy, sensibilidad, especificidad, Kappa
#'    - Generación de curva ROC y cálculo de AUC
#'
#' 5. **Organización de resultados**:
#'    - Compila todas las métricas en estructura lista
#'    - Agrega metadatos sobre el proceso de modelado
#'    - Organiza resultados por especie para fácil acceso
#'
#' 6. **Validación cruzada implícita**:
#'    - Random Forest incluye validación out-of-bag (OOB) internamente
#'    - La partición entrenamiento/prueba proporciona validación externa
#'
#' /
#' The function implements the following modeling pipeline for each species:
#'
#' 1. **Preprocessing and validation**:
#'    - Filters species with enough records (`Registros_Minimos`)
#'    - Converts `sf` objects to dataframes if necessary
#'    - Verifies all required variables are present
#'    - Removes duplicate records and checks class balance
#'
#' 2. **Data partitioning**:
#'    - Splits data into training (70%) and test (30%) sets
#'    - Maintains original presence/absence proportion in each set
#'    - Uses `caret::createDataPartition()` for stratified partitioning
#'
#' 3. **Model training**:
#'    - Configures Random Forest with optimized default parameters
#'    - Trains model with training data
#'    - Calculates variable importance during training
#'
#' 4. **Model evaluation**:
#'    - Predictions on training and test data
#'    - Calculation of confusion matrices for multiple sets
#'    - Calculation of metrics: accuracy, sensitivity, specificity, Kappa
#'    - Generation of ROC curve and AUC calculation
#'
#' 5. **Results organization**:
#'    - Compiles all metrics into list structure
#'    - Adds metadata about modeling process
#'    - Organizes results by species for easy access
#'
#' 6. **Implicit cross-validation**:
#'    - Random Forest includes out-of-bag (OOB) validation internally
#'    - Train/test split provides external validation
#'
#' @note
#' * **Paquetes requeridos**: Esta función utiliza múltiples paquetes:
#'   - `caret` para partición de datos y matrices de confusión
#'   - `randomForest` para el algoritmo de clasificación
#'   - `pROC` para curvas ROC y cálculo de AUC
#'   - `dplyr` para manipulación de datos
#'   - `terra` para manejo de raster (solo validación)
#' * **Balance de clases**: Random Forest maneja bien desbalance de clases,
#'   pero proporciones extremas pueden afectar métricas como sensibilidad.
#' * **Parámetros de Random Forest**: Se usan valores por defecto (500 árboles,
#'   sqrt(p) variables por split). Para ajustes avanzados, entrenar modelos
#'   individualmente fuera de esta función.
#' * **Reproducibilidad**: Use `set.seed()` antes de llamar la función para
#'   resultados reproducibles en partición de datos y Random Forest.
#' * **Tiempo de ejecución**: Para muchas especies o muchos registros, el
#'   proceso puede ser computacionalmente intensivo.
#' * **Validación espacial**: Esta función no realiza validación espacial
#'   explícita (block CV, environmental CV). Considere esto en interpretación.
#'
#' /
#' * **Required packages**: This function uses multiple packages:
#'   - `caret` for data partitioning and confusion matrices
#'   - `randomForest` for classification algorithm
#'   - `pROC` for ROC curves and AUC calculation
#'   - `dplyr` for data manipulation
#'   - `terra` for raster handling (validation only)
#' * **Class balance**: Random Forest handles class imbalance well,
#'   but extreme proportions may affect metrics like sensitivity.
#' * **Random Forest parameters**: Default values are used (500 trees,
#'   sqrt(p) variables per split). For advanced tuning, train models
#'   individually outside this function.
#' * **Reproducibility**: Use `set.seed()` before calling the function for
#'   reproducible results in data partitioning and Random Forest.
#' * **Execution time**: For many species or many records, the process
#'   can be computationally intensive.
#' * **Spatial validation**: This function does not perform explicit spatial
#'   validation (block CV, environmental CV). Consider this in interpretation.
#'
#' @section Advertencias/Warnings:
#' * **Registros insuficientes**: Especies con menos de `Registros_Minimos`
#'   registros serán omitidas con advertencia. Considere ajustar este umbral
#'   para estudios con datos limitados.
#' * **Variables faltantes**: Si alguna variable de `variables` no está
#'   presente en todos los dataframes de especies, la función intentará
#'   continuar con las variables disponibles, pero emitirá advertencia.
#' * **Predicciones perfectas**: En casos raros con datos muy separables,
#'   Random Forest puede producir accuracy perfecto (1.0), lo que puede
#'   indicar sobreajuste o datos artificialmente separables.
#' * **AUC no informativo**: AUC puede no ser informativo con desbalance
#'   extremo de clases o cuando el modelo tiene accuracy muy bajo.
#' * **Uso de memoria**: Modelar muchas especies simultáneamente puede
#'   consumir mucha memoria, especialmente con muchos registros por especie.
#'
#' /
#' * **Insufficient records**: Species with fewer than `Registros_Minimos`
#'   records will be omitted with warning. Consider adjusting this threshold
#'   for studies with limited data.
#' * **Missing variables**: If any variable from `variables` is not present
#'   in all species dataframes, the function will attempt to continue with
#'   available variables but will issue a warning.
#' * **Perfect predictions**: In rare cases with very separable data,
#'   Random Forest may produce perfect accuracy (1.0), which may indicate
#'   overfitting or artificially separable data.
#' * **Uninformative AUC**: AUC may not be informative with extreme class
#'   imbalance or when the model has very low accuracy.
#' * **Memory usage**: Modeling many species simultaneously can consume
#'   significant memory, especially with many records per species.
#'
#' @examples
#' \dontrun{
#' # Ejemplo 1: Uso básico con datos simulados
#' # Example 1: Basic usage with simulated data
#'
#' # Cargar paquetes requeridos
#' # Load required packages
#' library(randomForest)
#' library(caret)
#' library(pROC)
#' library(dplyr)
#'
#' # Crear datos de ejemplo para múltiples especies
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
#' # Crear raster de referencia (solo para validación)
#' # Create reference raster (for validation only)
#' library(terra)
#' raster_ref <- rast(nrows = 10, ncols = 10,
#'                    xmin = 0, xmax = 1, ymin = 0, ymax = 1)
#' values(raster_ref) <- matrix(1:100, 10, 10)
#' names(raster_ref) <- c("temperatura")
#'
#' # Ejecutar modelado para múltiples especies
#' # Run modeling for multiple species
#' resultados_modelos <- Rf_SDM_multiples_especies(
#'   Data_Df = data.frame(),  # Datos generales (puede estar vacío para ejemplo)
#'   variables = c("temperatura", "precipitacion", "altitud", "humedad"),
#'   Resultados_List = lista_datos,
#'   Raster_Tif_reproyectado = raster_ref,
#'   Registros_Minimos = 5
#' )
#'
#' # Acceder a resultados de una especie específica
#' # Access results for a specific species
#' resultados_especie_a <- resultados_modelos[["Especie_A"]]
#'
#' # Ver métricas clave
#' # View key metrics
#' cat("Accuracy test:", resultados_especie_a$accuracy_test, "\n")
#' cat("AUC:", resultados_especie_a$auc, "\n")
#' cat("Kappa:", resultados_especie_a$kappa, "\n")
#'
#' # Ver importancia de variables
#' # View variable importance
#' print(resultados_especie_a$importance)
#'
#' # Ejemplo 2: Con datos reales y ajuste de parámetros
#' # Example 2: With real data and parameter adjustment
#'
#' # Después de procesar presencias y pseudoausencias
#' # After processing presences and pseudoabsences
#' resultados_modelos_reales <- Rf_SDM_multiples_especies(
#'   Data_Df = datos_ambientales_completos,
#'   variables = nombres_variables_ambientales,
#'   Resultados_List = lista_presencias_ausencias_procesadas,
#'   Raster_Tif_reproyectado = raster_variables_proyectado,
#'   Registros_Minimos = 10  # Aumentar umbral para modelos más robustos
#' )
#'
#' # Análisis comparativo entre especies
#' # Comparative analysis between species
#' comparativa <- data.frame(
#'   Especie = names(resultados_modelos_reales),
#'   AUC = sapply(resultados_modelos_reales, function(x) x$auc),
#'   Accuracy_Test = sapply(resultados_modelos_reales, function(x) x$accuracy_test),
#'   Kappa = sapply(resultados_modelos_reales, function(x) x$kappa)
#' )
#'
#' print(comparativa[order(-comparativa$AUC), ])
#'
#' # Guardar modelos para uso futuro
#' # Save models for future use
#' saveRDS(resultados_modelos_reales, "modelos_randomforest_sdm.rds")
#' }
#'
#' @importFrom caret createDataPartition confusionMatrix
#' @importFrom randomForest randomForest
#' @importFrom pROC roc auc
#' @importFrom dplyr select filter %>% bind_rows
#' @importFrom stats predict
#'
#' @seealso
#' * [randomForest::randomForest()] para detalles del algoritmo de Random Forest
#'   / for details on the Random Forest algorithm
#' * [caret::confusionMatrix()] para interpretación de matrices de confusión
#'   / for interpretation of confusion matrices
#' * [pROC::roc()] para análisis detallado de curvas ROC
#'   / for detailed ROC curve analysis
#' * [Presencias_Sf_Sin_Na()] para preparar datos de entrada
#'   / to prepare input data
#' * [Simular_Pseudoausencias_Sin_NA()] para generar pseudoausencias balanceadas
#'   / to generate balanced pseudoabsences
#'
#' @export
Rf_SDM_multiple_species <- function(Data_Df, variables, Resultados_List, Raster_Tif_reproyectado, Registros_Minimos = 10){
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
