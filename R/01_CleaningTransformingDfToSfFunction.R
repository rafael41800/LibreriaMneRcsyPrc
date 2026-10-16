#' Generación optimizada de Pseudo-ausencias para Múltiples Especies
#'
#' @description
#' Genera puntos de pseudo-ausencia (background points) para una lista de especies,
#' evitando áreas cercanas a las presencias conocidas mediante un radio de exclusión.
#' La función implementa vectorización en R base (`lapply`) y filtrado de matrices nativo
#' para maximizar la eficiencia computacional sin dependencias de evaluación no estándar (NSE).
#'
#' @param Resultados Lista de `data.frame`. Contiene los datos de presencia por especie (típicamente la salida de la extracción de presencias).
#' @param RasterCrsReferencia Objeto `SpatRaster`. Define la extensión espacial, la resolución de las celdas y contiene las variables ambientales.
#' @param ObjetoCrsReferencia Sistema de Referencia de Coordenadas (CRS). Destino final para las coordenadas transformadas.
#' @param RadioGradosExclusion Numerico. Radio de exclusión a aplicar si el CRS de las presencias es geográfico (lat/lon). Por defecto es `0.001`.
#' @param RadioMetrosExclusion Numerico. Radio de exclusión a aplicar si el CRS de las presencias es proyectado (ej. UTM). Por defecto es `100`.
#' @param NombreVariablesAmbientales Vector de caracteres. Nombres exactos de las capas predictoras en `RasterCrsReferencia`.
#' @param verbose Logico. Si es `TRUE` (por defecto), imprime el progreso por especie en consola.
#' @param seed Numerico o `NULL`. Semilla aleatoria para garantizar la reproducibilidad en el sorteo de celdas.
#'
#' @return Una lista de objetos `data.frame` tabulares. Cada elemento corresponde a una especie
#'   e incluye los datos de presencia (`Presence = 1`), pseudo-ausencia (`Presence = 0`),
#'   los valores ambientales sin NAs, y las coordenadas transformadas (`X`, `Y`).
#'
#' @details
#' El flujo de trabajo interno sigue estos pasos vectorizados:
#' \enumerate{
#'   \item \strong{Espacialización}: Convierte las coordenadas del `data.frame` a un objeto `sf` y lo proyecta al CRS del raster ambiental.
#'   \item \strong{Zona de Exclusión}: Construye un buffer topológico (`st_buffer` + `st_union`) alrededor de las presencias usando el radio estipulado.
#'   \item \strong{Muestreo Estratificado}: Identifica las celdas del raster que contienen datos válidos (sin NAs) y descarta aquellas donde ya existe una presencia. Sortea aleatoriamente las celdas candidatas.
#'   \item \strong{Filtro Espacial}: Intersecta los candidatos con el polígono de exclusión para retener solo los puntos verdaderamente alejados (`st_intersects`).
#'   \item \strong{Extracción Ambiental}: Obtiene los valores climáticos o topográficos para los puntos retenidos mediante `terra::extract`.
#'   \item \strong{Ensamble Tabular}: Une presencias y pseudo-ausencias usando indexación matricial rápida (R base), elimina la geometría compleja (`st_drop_geometry`) y extrae las coordenadas finales en formato numérico plano para evitar errores en modelos de *Machine Learning*.
#' }
#'
#' @note
#' Para obtener una guía detallada sobre la calibración del radio de exclusión,
#' consulte la viñeta de la librería: \code{vignette("spatial_cleaning_guide", package = "LibreriaMneRcsyPrc")}.
#'
#' @section Advertencias:
#' \itemize{
#'   \item \strong{Extensión del Raster}: Si el área de estudio es muy pequeña o el radio de exclusión es excesivamente grande, la función podría no encontrar celdas disponibles y la especie será omitida con un `warning`.
#'   \item \strong{Compatibilidad con Algoritmos}: La destrucción intencional de la geometría espacial al final del proceso garantiza que paquetes como `randomForest` o `gbm` no colapsen al ingerir matrices espaciales.
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' # 1. Crear presencias simuladas
#' df_sp <- data.frame(
#'    longitud = c(-99.1, -99.2, -99.3),
#'    latitud = c(19.4, 19.5, 19.6),
#'    especievalida = "Sp_A"
#' )
#' lista_presencias <- list("Sp_A" = df_sp)
#'
#' # 2. Crear raster ambiental
#' r_env <- rast(nrows=100, ncols=100, ext=ext(-100, -98, 18, 20), crs="EPSG:4326")
#' values(r_env) <- runif(ncell(r_env))
#' names(r_env) <- "Bio1"
#'
#' # 3. Ejecutar simulación
#' pseudo_abs <- SimulatePseudoabsencesNoNAListDf(
#'    Resultados = lista_presencias,
#'    RasterCrsReferencia = r_env,
#'    ObjetoCrsReferencia = 4326,
#'    RadioGradosExclusion = 0.05,
#'    NombreVariablesAmbientales = c("Bio1"),
#'    verbose = TRUE,
#'    seed = 42
#' )
#' }
#'
#' @seealso
#' \code{\link[sf]{st_buffer}}, \code{\link[terra]{cellFromXY}}, \code{\link[terra]{extract}}
#'
#' @keywords sdm pseudo-ausencias background-points espacial
#'
#' @importFrom sf st_as_sf st_transform st_crs st_buffer st_union st_is_longlat st_coordinates st_intersects st_drop_geometry
#' @importFrom terra cellFromXY values xyFromCell extract vect
#' @importFrom stats complete.cases
#' @export
CleaningTransformingDfToSf <- function(DataFrame,
                                       SpatialReference,
                                       Latitude = 0,
                                       Variables,
                                       EspecieColumna = "especievalida",
                                       FechaColumna = "fechacolecta",
                                       SCR = 4326,
                                       PrecisionDecimales = 4,
                                       verbose = TRUE) {

  nombres_req <- c("longitud", "latitud", EspecieColumna, FechaColumna)
  faltantes <- setdiff(nombres_req, names(DataFrame))
  if(length(faltantes) > 0) {
    stop(paste("Error algoritmico: El DataFrame carece de las columnas:", paste(faltantes, collapse = ", ")))
  }
  if(verbose) message("\n--- Inicia procesamiento de datos geoespaciales ---")
  ColumnasNecesarias <- unique(c(Variables, EspecieColumna, FechaColumna, "longitud", "latitud"))
  DFNorth <- DataFrame %>%
    dplyr::filter(!is.na(.data$longitud) & !is.na(.data$latitud)) %>%
    dplyr::filter(.data$latitud >= Latitude) %>%
    dplyr::select(dplyr::any_of(ColumnasNecesarias))
  if(verbose) message("1/7 Filtrado inicial: ", nrow(DFNorth), " registros retenidos.")
  DFNorth <- sf::st_as_sf(DFNorth, coords = c("longitud", "latitud"), crs = SCR)
  transformation <- sf::st_transform(DFNorth, sf::st_crs(SpatialReference))
  if(verbose) message("2/7 & 3/7 Transformacion a CRS objetivo: OK")
  mask <- sf::st_filter(transformation, SpatialReference, .predicate = sf::st_intersects)
  if(verbose) message("4/7 Filtro topologico espacial: ", nrow(mask), " puntos retenidos.")
  coords <- sf::st_coordinates(mask)
  DFNorthClean <- mask %>%
    dplyr::mutate(
      lon_round = round(coords[,1], digits = PrecisionDecimales),
      lat_round = round(coords[,2], digits = PrecisionDecimales)
    ) %>%
    dplyr::distinct(.data$lon_round, .data$lat_round, .data[[EspecieColumna]], .data[[FechaColumna]], .keep_all = TRUE) %>%
    dplyr::select(-dplyr::any_of(c("lon_round", "lat_round")))
  if(verbose) {
    eliminados <- nrow(mask) - nrow(DFNorthClean)
    message("5/7 Eliminacion de duplicados (resolucion ", PrecisionDecimales, " decimales): ",
            nrow(DFNorthClean), " unicos (", eliminados, " eliminados).")
  }
  DFNorthClean <- DFNorthClean %>%
    dplyr::add_count(.data[[EspecieColumna]], name = "FrecuenciaEspecieValida")
  if(verbose) message("6/7 Calculo de frecuencias: OK")
  if(verbose) message("7/7 Procesamiento completado con exito.")
  return(DFNorthClean)
}
