#' @title Preparacion y limpieza de datos para analisis geoespacial / Data Preparation and Cleaning for Geospatial Analysis
#'
#' @description
#' Transforma un data.frame con coordenadas en un objeto `sf` (Simple Features)
#' listo para analisis espacial, aplicando filtros, transformaciones de CRS,
#' recorte y eliminacion de duplicados.
#' /
#' Transforms a data.frame with coordinates into an `sf` (Simple Features) object
#' ready for spatial analysis, applying filters, CRS transformations,
#' clipping, and duplicate removal.
#'
#' @param DataFrame `data.frame` con los datos a procesar. Debe contener las columnas
#'   `longitud`, `latitud` y `especie`.
#'   / `data.frame` with data to process. Must contain columns
#'   `longitud`, `latitud`, and `especie`.
#'
#' @param SpatialReference Objeto espacial (`sf` o `sfc`) que define la referencia
#'   espacial de destino y area de recorte.
#'   / Spatial object (`sf` or `sfc`) defining the target spatial reference
#'   and clipping area.
#'
#' @param latitude Valor numerico que define el limite minimo de latitud para
#'   filtrar observaciones. Por defecto es 0 (ecuador).
#'   / Numeric value defining the minimum latitude threshold for
#'   filtering observations. Default is 0 (equator).
#'
#' @param variables Vector de caracteres con los nombres de las columnas a
#'   conservar en el resultado.
#'   / Character vector with column names to
#'   retain in the output.
#'
#' @return
#' Un objeto `sf` con las siguientes caracteristicas:
#' * Geometrias en el sistema de referencia de `SpatialReference`
#' * Variables seleccionadas mediante el parametro `variables`
#' * Columna adicional `espcvld` con el conteo de observaciones por especie
#' * Sin duplicados espaciales (mismas coordenadas y misma especie)
#' * Recortado al area definida por `SpatialReference`
#' * Observaciones al norte de `latitude`
#' * /
#'   * An `sf` object with the following features:
#'   * Geometries in the coordinate reference system of `SpatialReference`
#' * Variables selected through the `variables` parameter
#' * Additional column `espcvld` with observation count per species
#' * No spatial duplicates (same coordinates and same species)
#' * Clipped to the area defined by `SpatialReference`
#' * Observations north of `latitude`
#' *
#' @details
#' La funcion ejecuta el siguiente pipeline de procesamiento:
#' 1. **Filtrado**: Conserva observaciones con `latitud >= latitude`
#' 2. **Seleccion**: Mantiene solo las variables especificadas en `variables`
#' 3. **Conversion a `sf`**: Crea geometrias puntuales desde coordenadas
#'    geograficas (CRS: WGS84, EPSG:4326)
#' 4. **Transformacion CRS**: Proyecta al sistema de coordenadas de
#'    `SpatialReference`
#' 5. **Recorte**: Intersecta con el area de `SpatialReference`
#' 6. **Deduplicacion**: Elimina puntos con mismas coordenadas y especie
#' 7. **Conteo**: Agrega columna `espcvld` con frecuencia por especie
#' /
#' The function executes the following processing pipeline:
#' 1. **Filtering**: Keeps observations with `latitud >= latitude`
#' 2. **Selection**: Retains only variables specified in `variables`
#' 3. **Conversion to `sf`**: Creates point geometries from geographic
#'    coordinates (CRS: WGS84, EPSG:4326)
#' 4. **CRS Transformation**: Projects to the coordinate system of
#'    `SpatialReference`
#' 5. **Clipping**: Intersects with the area of `SpatialReference`
#' 6. **Deduplication**: Removes points with same coordinates and species
#' 7. **Counting**: Adds column `espcvld` with frequency per species
#'
#' @note
#' * La columna `especie` debe estar presente aunque no se incluya en `variables`
#' * La deduplicacion se basa en coordenadas transformadas (no originales)
#' * Se crea un directorio temporal `~/temp` para manejo de datos raster si es necesario
#' * La funcion asume que las coordenadas estan en grados decimales (WGS84)
#' /
#' * The `especie` column must be present even if not included in `variables`
#' * Deduplication is based on transformed coordinates (not originals)
#' * A temporary directory `~/temp` is created for raster data handling if needed
#' * The function assumes coordinates are in decimal degrees (WGS84)
#'
#' @section Advertencias/Warnings:
#' * Si `latitude` es mayor que el limite norte de `SpatialReference`, el resultado
#'   podria ser vacio
#' * La proyeccion podria distorsionar distancias si se usa entre hemisferios
#' * El proceso de interseccion podria ser lento con muchos poligonos complejos
#' /
#' * If `latitude` is greater than the northern limit of `SpatialReference`, the result
#'   could be empty
#' * Projection may distort distances if used across hemispheres
#' * Intersection process could be slow with many complex polygons
#'
#' @examples
#' \dontrun{
#' # Crear datos de ejemplo / Create example data
#' set.seed(123)
#' datos_ejemplo <- data.frame(
#'   longitud = runif(100, -100, -80),
#'   latitud = runif(100, 15, 30),
#'   especie = sample(c("Quercus", "Pinus", "Abies"), 100, replace = TRUE),
#'   altura = runif(100, 5, 25),
#'   diametro = runif(100, 10, 50)
#' )
#'
#' # Crear poligono de referencia (Zona UTM 14N) / Create reference polygon (UTM Zone 14N)
#' library(sf)
#' bbox <- st_bbox(c(xmin = -95, xmax = -85, ymin = 20, ymax = 25), crs = 4326)
#' referencia <- st_as_sfc(bbox) %>% st_transform(32614)
#'
#' # Procesar datos / Process data
#' datos_espaciales <- CleaningTransformingDfToSf(
#'   DataFrame = datos_ejemplo,
#'   SpatialReference = referencia,
#'   latitude = 18,
#'   variables = c("especie", "altura", "diametro")
#' )
#'
#' # Inspeccionar resultados / Inspect results
#' class(datos_espaciales)
#' nrow(datos_espaciales)
#' table(datos_espaciales$especie, datos_espaciales$espcvld)
#'
#' # Visualizacion basica / Basic visualization
#' plot(st_geometry(datos_espaciales), pch = 16,
#'      col = factor(datos_espaciales$especie))
#' }
#'
#' @seealso
#' * \code{\link[sf]{st_as_sf}} para conversion a objetos espaciales / for conversion to spatial objects
#' * \code{\link[sf]{st_transform}} para transformaciones de CRS / for CRS transformations
#' * \code{\link[sf]{st_intersection}} para operaciones de recorte / for clipping operations
#'
#' @keywords
#'   limpieza de datos, transformacion espacial, GIS, datos geoespaciales,
#'   data cleaning, spatial transformation, GIS, geospatial data
#'
#' @family
#'   funciones_espaciales, funciones_limpieza_datos /
#'   spatial_functions, data_cleaning_functions
#'
#'
#' @importFrom sf st_as_sf st_transform st_intersection st_coordinates st_crs
#' @importFrom dplyr filter select mutate distinct %>% all_of
#' @importFrom terra terraOptions
#' @export
CleaningTransformingDfToSf <- function(DataFrame, SpatialReference, latitude=0, variables) {
  if(!dir.exists("~/temp")) dir.create("~/temp")
  terra::terraOptions(tempdir=("~/temp"), memfrac=0.8)
  DF.north <- DataFrame %>%
    dplyr::filter(latitud >= latitude) %>%
    dplyr::select(all_of(variables)) %>%
    sf::st_as_sf(coords = c("longitud", "latitud"), crs = 4326)
  transformation <- sf::st_transform(DF.north, sf::st_crs(SpatialReference))
  mask <- sf::st_intersection(transformation, SpatialReference)
  DF.north.clean <- mask %>%
    dplyr::mutate(XY = paste0(st_coordinates(.)[,1], st_coordinates(.)[,2], especie)) %>%
    dplyr::distinct(XY, .keep_all = TRUE) %>%
    dplyr::select(-XY)
  DF.north.clean <- DF.north.clean %>%
    dplyr::mutate(espcvld = ave(seq_len(nrow(.)), especie, FUN = length))
  return(DF.north.clean)
}
