#' Preparación y limpieza de datos para análisis geoespacial
#'
#' @description
#' Transforma un data.frame con coordenadas en un objeto `sf` (Simple Features)
#' listo para análisis espacial, aplicando filtros, transformaciones de CRS,
#' recorte y eliminación de duplicados.
#'
#' Transforms a data.frame with coordinates into an `sf` (Simple Features) object
#' ready for spatial analysis, applying filters, CRS transformations,
#' clipping, and duplicate removal.
#'
#' @param DataFrame `data.frame` con los datos a procesar. Debe contener las columnas
#'   `longitud`, `latitud` y la columna especificada en `EspecieColumna`.
#'   `data.frame` with data to process. Must contain columns `longitud`, `latitud`,
#'   and the column specified in `EspecieColumna`.
#' @param SpatialReference Objeto espacial (`sf` o `sfc`) que define la referencia
#'   espacial de destino y área de recorte.
#'   Spatial object (`sf` or `sfc`) defining the target spatial reference
#'   and clipping area.
#' @param Latitude Valor numérico que define el límite mínimo de latitud para
#'   filtrar observaciones. Por defecto es `0` (ecuador).
#'   Numeric value defining the minimum latitude threshold for filtering observations.
#'   Default is `0` (equator).
#' @param Variables Vector de caracteres con los nombres de las columnas a
#'   conservar en el resultado.
#'   Character vector with column names to retain in the output.
#' @param EspecieColumna Nombre de la columna que contiene la información de especie.
#'   Por defecto es `"especie"`.
#'   Name of the column containing species information. Default is `"especie"`.
#' @param verbose Logical. Si `TRUE` (por defecto), muestra mensajes de progreso
#'   en la consola. Si `FALSE`, ejecuta silenciosamente.
#'   Logical. If `TRUE` (default), displays progress messages in console.
#'   If `FALSE`, runs silently.
#'
#' @return
#' Un objeto `sf` con las siguientes características:
#' \itemize{
#'   \item Geometrías en el sistema de referencia de `SpatialReference`
#'   \item Variables seleccionadas mediante el parámetro `Variables`
#'   \item Columna adicional `FrecuenciaEspecieValida` con el conteo de observaciones por especie
#'   \item Sin duplicados espaciales (mismas coordenadas y misma especie)
#'   \item Recortado al área definida por `SpatialReference`
#'   \item Observaciones al norte de `Latitude`
#' }
#'
#' An `sf` object with the following features:
#' \itemize{
#'   \item Geometries in the coordinate reference system of `SpatialReference`
#'   \item Variables selected through the `Variables` parameter
#'   \item Additional column `FrecuenciaEspecieValida` with observation count per species
#'   \item No spatial duplicates (same coordinates and same species)
#'   \item Clipped to the area defined by `SpatialReference`
#'   \item Observations north of `Latitude`
#' }
#'
#' @details
#' \strong{Español:} La función ejecuta el siguiente pipeline de procesamiento:
#' \enumerate{
#'   \item \strong{Filtrado}: Conserva observaciones con `latitud >= Latitude`
#'   \item \strong{Selección}: Mantiene solo las variables especificadas en `Variables`
#'   \item \strong{Conversión a `sf`}: Crea geometrías puntuales desde coordenadas
#'         geográficas (CRS: WGS84, EPSG:4326)
#'   \item \strong{Transformación CRS}: Proyecta al sistema de coordenadas de `SpatialReference`
#'   \item \strong{Recorte}: Intersecta con el área de `SpatialReference`
#'   \item \strong{Deduplicación}: Elimina puntos con mismas coordenadas y especie
#'   \item \strong{Conteo}: Agrega columna `FrecuenciaEspecieValida` con frecuencia por especie
#' }
#'
#' \strong{English:} The function executes the following processing pipeline:
#' \enumerate{
#'   \item \strong{Filtering}: Keeps observations with `latitud >= Latitude`
#'   \item \strong{Selection}: Retains only variables specified in `Variables`
#'   \item \strong{Conversion to `sf`}: Creates point geometries from geographic
#'         coordinates (CRS: WGS84, EPSG:4326)
#'   \item \strong{CRS Transformation}: Projects to the coordinate system of `SpatialReference`
#'   \item \strong{Clipping}: Intersects with the area of `SpatialReference`
#'   \item \strong{Deduplication}: Removes points with same coordinates and species
#'   \item \strong{Counting}: Adds column `FrecuenciaEspecieValida` with frequency per species
#' }
#'
#' @note
#' \strong{Español:}
#' \itemize{
#'   \item La columna especificada en `EspecieColumna` debe estar presente aunque no se incluya en `Variables`
#'   \item La deduplicación se basa en coordenadas transformadas (no originales)
#'   \item Se crea un directorio temporal `~/temp` para manejo de datos raster si es necesario
#'   \item La función asume que las coordenadas están en grados decimales (WGS84)
#'   \item Los mensajes de progreso pueden suprimirse estableciendo `verbose = FALSE`
#' }
#'
#' \strong{English:}
#' \itemize{
#'   \item The column specified in `EspecieColumna` must be present even if not included in `Variables`
#'   \item Deduplication is based on transformed coordinates (not originals)
#'   \item A temporary directory `~/temp` is created for raster data handling if needed
#'   \item The function assumes coordinates are in decimal degrees (WGS84)
#'   \item Progress messages can be suppressed by setting `verbose = FALSE`
#' }
#'
#' @section Warning/Advertencia:
#' \itemize{
#'   \item \strong{Español}: Si `Latitude` es mayor que el límite norte de `SpatialReference`, el resultado podría ser vacío
#'   \item \strong{English}: If `Latitude` is greater than the northern limit of `SpatialReference`, the result could be empty
#'   \item \strong{Español}: La proyección podría distorsionar distancias si se usa entre hemisferios
#'   \item \strong{English}: Projection may distort distances if used across hemispheres
#'   \item \strong{Español}: El proceso de intersección podría ser lento con muchos polígonos complejos
#'   \item \strong{English}: Intersection process could be slow with many complex polygons
#' }
#'
#' @examples
#' \dontrun{
#' # Crear datos de ejemplo
#' set.seed(123)
#' datos_ejemplo <- data.frame(
#'   longitud = runif(100, -100, -80),
#'   latitud = runif(100, 15, 30),
#'   especie = sample(c("Quercus", "Pinus", "Abies"), 100, replace = TRUE),
#'   altura = runif(100, 5, 25),
#'   diametro = runif(100, 10, 50)
#' )
#'
#' # Crear polígono de referencia (Zona UTM 14N)
#' library(sf)
#' bbox <- sf::st_bbox(c(xmin = -95, xmax = -85, ymin = 20, ymax = 25), crs = 4326)
#' referencia <- sf::st_as_sfc(bbox) |> sf::st_transform(32614)
#'
#' # Procesar datos
#' datos_espaciales <- CleaningTransformingDfToSf(
#'   DataFrame = datos_ejemplo,
#'   SpatialReference = referencia,
#'   Latitude = 18,
#'   Variables = c("especie", "altura", "diametro")
#' )
#' }
#'
#' @seealso
#' \code{\link[sf]{st_as_sf}}, \code{\link[sf]{st_transform}}, \code{\link[sf]{st_intersection}}
#'
#' @keywords data_cleaning spatial_transformation gis
#'
#' @importFrom sf st_as_sf st_transform st_intersection st_coordinates st_crs
#' @importFrom dplyr filter select mutate distinct add_count all_of
#' @importFrom terra terraOptions
#' @importFrom rlang .data
#' @export
CleaningTransformingDfToSf <- function(DataFrame, SpatialReference, Latitude=0, Variables, EspecieColumna="especie", verbose = TRUE) {
  if(!dir.exists("~/temp")) dir.create("~/temp")
  terra::terraOptions(tempdir=("~/temp"), memfrac=0.8)
  if(verbose) cat("\n--- Procesando datos geoespaciales ---\n")
  ColumnasNecesarias <- unique(c(Variables, EspecieColumna, "longitud", "latitud"))
  DFNorth <- DataFrame %>%
    dplyr::filter(.data$latitud >= Latitude) %>%
    dplyr::select(dplyr::all_of(ColumnasNecesarias))
  if(verbose) cat("1/7 Filtrado: ", nrow(DFNorth), "/", nrow(DataFrame), " registros\n")
  DFNorth <- DFNorth %>%
    sf::st_as_sf(coords = c("longitud", "latitud"), crs = 4326)
  if(verbose) cat("2/7 Conversion a sf: OK\n")
  transformation <- sf::st_transform(DFNorth, sf::st_crs(SpatialReference))
  if(verbose) cat("3/7 Transformacion CRS: OK\n")
  mask <- sf::st_intersection(transformation, SpatialReference)
  if(verbose) cat("4/7 Recorte espacial: ", nrow(mask), " puntos retenidos\n")
  coords <- sf::st_coordinates(mask)
  DFNorthClean <- mask %>%
    dplyr::mutate(XY = paste0(coords[,1], coords[,2], .data[[EspecieColumna]])) %>%
    dplyr::distinct(.data$XY, .keep_all = TRUE) %>%
    dplyr::select(-.data$XY)
  if(verbose) cat("5/7 Eliminacion duplicados: ", nrow(DFNorthClean),
                  " unicos (", nrow(mask) - nrow(DFNorthClean), " eliminados)\n")
  DFNorthClean <- DFNorthClean %>%
    dplyr::add_count(.data[[EspecieColumna]], name = "FrecuenciaEspecieValida")
  if(verbose) cat("6/7 Frecuencia por especie: OK\n")
  if(verbose) cat("7/7 Procesamiento completado!\n")
  return(DFNorthClean)
}
