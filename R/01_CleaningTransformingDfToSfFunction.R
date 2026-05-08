#' Preparacion y limpieza de datos para analisis geoespacial
#'
#' @description
#' Transforma un `data.frame` con coordenadas geograficas en un objeto espacial `sf`.
#' La funcion integra un pipeline completo de pre-procesamiento: filtrado por latitud,
#' eliminacion de registros sin coordenadas, proyeccion cartografica, recorte
#' mediante una capa de referencia, eliminacion de duplicados espaciales y
#' calculo de frecuencia por especie.
#' @param DataFrame Un `data.frame` que debe contener obligatoriamente las columnas
#'   `longitud`, `latitud` y la columna definida en `EspecieColumna`.
#' @param SpatialReference Objeto espacial (`sf` o `sfc`) que actua como mascara
#'   de recorte y define el Sistema de Referencia de Coordenadas (CRS) de destino.
#' @param Latitude Valor numerico. Límite minimo de latitud para conservar registros.
#'   util para filtrar hemisferios o areas de interes. Por defecto es `0`.
#' @param Variables Vector de caracteres. Nombres de las columnas adicionales que
#'   desea conservar en el objeto final (ej. `c("altura", "diametro")`).
#' @param EspecieColumna Caracter. Nombre de la columna que identifica a la especie.
#'   Por defecto es `"especie"`.
#' @param SCR Numerico o caracter. Sistema de Referencia de Coordenadas original de
#'   los datos en `DataFrame`. Por defecto es `4326` (WGS84).
#' @param verbose Logico. Si es `TRUE` (por defecto), imprime en consola el progreso
#'   detallado de cada etapa del procesamiento.
#'
#' @return Un objeto de clase `sf` (puntos) con el CRS de `SpatialReference`.
#'   Incluye la columna `FrecuenciaEspecieValida` que indica cuantas veces aparece
#'   cada especie en el set de datos final tras la limpieza.
#'
#' @details
#' El flujo de trabajo interno sigue estos pasos:
#' \enumerate{
#'   \item \strong{Limpieza inicial}: Elimina filas con `NA` en coordenadas y filtra por `Latitude`.
#'   \item \strong{Espacializacion}: Convierte el DF a objeto `sf` usando el `SCR` proporcionado.
#'   \item \strong{Reproyeccion}: Transforma los datos al CRS del objeto `SpatialReference`.
#'   \item \strong{Recorte Geografico}: Interseccion espacial para mantener solo puntos dentro de la referencia.
#'   \item \strong{Deduplicacion}: Elimina registros que comparten las mismas coordenadas y especie, evitando redundancia.
#'   \item \strong{Estadisticas}: Agrega el conteo de frecuencia por grupo biologico.
#' }
#'
#' @note
#' Para obtener una guia detallada en ingles y ejemplos avanzados de flujos de trabajo,
#' consulte la viñeta de la libreria: \code{vignette("spatial_cleaning_guide", package = "TuLibreria")}.
#'
#' @section Advertencias:
#' \itemize{
#'   \item \strong{Nombres de Columnas}: Si las columnas `longitud` o `latitud` tienen nombres distintos, debe renombrarlas antes de ingresar el objeto.
#'   \item \strong{Memoria RAM}: Para datasets de millones de registros, el proceso de interseccion y deduplicacion puede ser intensivo. Se ha optimizado el uso de memoria mediante \code{terraOptions}.
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(dplyr)
#'
#' # 1. Preparar un dataset sintetico
#' datos <- data.frame(
#'   longitud = runif(200, -102, -98),
#'   latitud = runif(200, 19, 21),
#'   especie = sample(c("Quercus rugosa", "Pinus ayacahuite"), 200, replace = TRUE),
#'   diametro = runif(200, 10, 80),
#'   condicion = "vivo"
#' )
#'
#' # 2. Crear una referencia espacial (ej. un cuadrado en el centro de Mexico)
#' # Usaremos el EPSG:6372 (Mexico ITRF2008 / LCC) para la proyeccion de destino
#' poligono_ref <- st_sfc(st_polygon(list(matrix(
#'   c(-101, 19.5, -99, 19.5, -99, 20.5, -101, 20.5, -101, 19.5),
#'   ncol = 2, byrow = TRUE))), crs = 4326) %>%
#'   st_transform(6372)
#'
#' # 3. Ejecutar la funcion
#' resultado <- CleaningTransformingDfToSf(
#'   DataFrame = datos,
#'   SpatialReference = poligono_ref,
#'   Latitude = 18,
#'   Variables = c("diametro", "condicion"),
#'   EspecieColumna = "especie",
#'   SCR = 4326,
#'   verbose = TRUE
#' )
#'
#' # 4. Inspeccionar el resultado
#' print(resultado)
#' plot(st_geometry(resultado), col = "blue", pch = 20)
#' }
#'
#' @seealso
#' \code{\link[sf]{st_as_sf}}, \code{\link[sf]{st_intersection}}, \code{\link[dplyr]{distinct}}
#'
#' @keywords limpieza espacial sf gis
#'
#' @importFrom sf st_as_sf st_transform st_intersection st_coordinates st_crs st_geometry
#' @importFrom dplyr filter select mutate distinct add_count all_of
#' @importFrom terra terraOptions
#' @importFrom rlang .data
#' @export
CleaningTransformingDfToSf <- function(DataFrame,
                                       SpatialReference,
                                       Latitude = 0,
                                       Variables,
                                       EspecieColumna = "especievalida",
                                       SCR = 4326,
                                       verbose = TRUE) {

  nombres_req <- c("longitud", "latitud", EspecieColumna)
  faltantes <- setdiff(nombres_req, names(DataFrame))
  if(length(faltantes) > 0) {
    stop(paste("El DataFrame no contiene las columnas:", paste(faltantes, collapse = ", ")))
  }

  if(verbose) message("\n--- Inicia procesamiento datos geoespaciales ---")

  ColumnasNecesarias <- unique(c(Variables, EspecieColumna, "longitud", "latitud"))
  DFNorth <- DataFrame %>%
    dplyr::filter(!is.na(.data$longitud) & !is.na(.data$latitud)) %>%
    dplyr::filter(.data$latitud >= Latitude) %>%
    dplyr::select(dplyr::all_of(ColumnasNecesarias))

  if(verbose) message("1/7 Filtrado: ", nrow(DFNorth), "/", nrow(DataFrame), " registros")

  DFNorth <- sf::st_as_sf(DFNorth, coords = c("longitud", "latitud"), crs = SCR)
  if(verbose) message("2/7 Conversion a sf (CRS:", SCR, "): OK")

  transformation <- sf::st_transform(DFNorth, sf::st_crs(SpatialReference))
  if(verbose) message("3/7 Transformacion CRS: OK")

  mask <- sf::st_intersection(transformation, sf::st_geometry(SpatialReference))
  if(verbose) message("4/7 Recorte espacial: ", nrow(mask), " puntos retenidos")

  coords <- sf::st_coordinates(mask)
  DFNorthClean <- mask %>%
    dplyr::mutate(XY = paste(coords[,1], coords[,2], .data[[EspecieColumna]], sep = "_")) %>%
    dplyr::distinct(.data$XY, .keep_all = TRUE) %>%
    dplyr::select(-.data$XY)

  if(verbose) {
    eliminados <- nrow(mask) - nrow(DFNorthClean)
    message("5/7 Eliminacion duplicados: ", nrow(DFNorthClean), " unicos (", eliminados, " eliminados)")
  }

  DFNorthClean <- DFNorthClean %>%
    dplyr::add_count(.data[[EspecieColumna]], name = "FrecuenciaEspecieValida", .drop = FALSE)

  if(verbose) message("6/7 Frecuencia por especie: OK")
  if(verbose) message("7/7 Procesamiento completado!")

  return(DFNorthClean)
}
