#' Preparacion y limpieza de datos para analisis geoespacial
#'
#' @description
#' Transforma un `data.frame` con coordenadas geograficas en un objeto espacial `sf`.
#' La funcion integra un pipeline completo de pre-procesamiento vectorizado: filtrado por latitud,
#' eliminacion de registros sin coordenadas, proyeccion cartografica, filtro topologico
#' mediante una capa de referencia, eliminacion de duplicados espaciales (con control de precision) y
#' calculo de frecuencia por especie.
#'
#' @param DataFrame Un `data.frame` que debe contener obligatoriamente las columnas
#'   `longitud`, `latitud`, la definida en `EspecieColumna` y la definida en `FechaColumna`.
#' @param SpatialReference Objeto espacial (`sf` o `sfc`) que actua como mascara
#'   de recorte y define el Sistema de Referencia de Coordenadas (CRS) de destino.
#' @param Latitude Valor numerico. Limite minimo de latitud para conservar registros.
#'   util para filtrar hemisferios o areas de interes. Por defecto es `0`.
#' @param Variables Vector de caracteres. Nombres de las columnas adicionales que
#'   desea conservar. Puede ser un subconjunto propio o impropio; se conservaran solo aquellas que existan en el `DataFrame`.
#' @param EspecieColumna Caracter. Nombre de la columna que identifica a la especie.
#'   Por defecto es `"especievalida"`.
#' @param FechaColumna Caracter. Nombre de la columna que identifica la fecha de colecta.
#'   Por defecto es `"fechacolecta"`.
#' @param SCR Numerico o caracter. Sistema de Referencia de Coordenadas original de
#'   los datos en `DataFrame`. Por defecto es `4326` (WGS84).
#' @param PrecisionDecimales Numerico. Cantidad de decimales a los que se redondearan
#'   las coordenadas antes de buscar duplicados. Modificar este valor aumenta o reduce el tamaño
#'   del area cubierta para la deduplicacion: reducir el numero de decimales (ej. `3`) aumenta
#'   el area de tolerancia en la que los puntos se consideran identicos (~111 metros), agrupando
#'   mas registros. Aumentar el numero de decimales (ej. `5`) reduce el area cubierta, exigiendo
#'   mayor proximidad espacial para considerarlos duplicados (~1.1 metros). Un valor de `4`
#'   mitiga errores de coma flotante y equivale a una resolucion base de ~11 metros. Por defecto es `4`.
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
#'   \item \strong{Filtro Topologico}: Utiliza `st_filter` para retener de forma eficiente solo los puntos que intersectan con la referencia espacial.
#'   \item \strong{Deduplicacion}: Redondea las coordenadas segun `PrecisionDecimales` para evitar falsos unicos por artefactos de coma flotante y elimina registros identicos en espacio, especie y fecha.
#'   \item \strong{Estadisticas}: Agrega el conteo de frecuencia por grupo biologico.
#' }
#'
#' @note
#' Para obtener una guia detallada en ingles y ejemplos avanzados de flujos de trabajo,
#' consulte la viñeta de la libreria: \code{vignette("spatial_cleaning_guide", package = "LibreriaMneRcsyPrc")}.
#'
#' @section Advertencias:
#' \itemize{
#'   \item \strong{Nombres de Columnas}: Si las columnas `longitud` o `latitud` tienen nombres distintos, debe renombrarlas antes de ingresar el objeto.
#'   \item \strong{Rendimiento Optimo}: Al utilizar `st_filter` en lugar de intersecciones geometricas duras y evaluar duplicados vectorizadamente, esta funcion esta disenada para ser segura en memoria (RAM) incluso con datasets masivos.
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(dplyr)
#'
#' # 1. Preparar un dataset sintetico
#' datos <- data.frame(
#'    longitud = runif(200, -102, -98),
#'    latitud = runif(200, 19, 21),
#'    especievalida = sample(c("Quercus rugosa", "Pinus ayacahuite"), 200, replace = TRUE),
#'    fechacolecta = sample(seq(as.Date('2020/01/01'), as.Date('2023/01/01'), by="day"), 200),
#'    diametro = runif(200, 10, 80),
#'    condicion = "vivo"
#' )
#'
#' # 2. Crear una referencia espacial
#' poligono_ref <- st_sfc(st_polygon(list(matrix(
#'    c(-101, 19.5, -99, 19.5, -99, 20.5, -101, 20.5, -101, 19.5),
#'    ncol = 2, byrow = TRUE))), crs = 4326) %>%
#'    st_transform(6372)
#'
#' # 3. Ejecutar la funcion
#' resultado <- CleaningTransformingDfToSf(
#'    DataFrame = datos,
#'    SpatialReference = poligono_ref,
#'    Variables = c("diametro", "condicion", "variable_inexistente"),
#'    EspecieColumna = "especievalida",
#'    FechaColumna = "fechacolecta",
#'    PrecisionDecimales = 4 # Modificar a 3 para agrupar puntos en un area mayor (~111m)
#' )
#' }
#'
#' @seealso
#' \code{\link[sf]{st_as_sf}}, \code{\link[sf]{st_filter}}, \code{\link[dplyr]{distinct}}
#'
#' @keywords limpieza espacial sf gis
#'
#' @importFrom sf st_as_sf st_transform st_filter st_coordinates st_crs st_intersects
#' @importFrom dplyr filter select mutate distinct add_count any_of %>%
#' @importFrom rlang .data
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
