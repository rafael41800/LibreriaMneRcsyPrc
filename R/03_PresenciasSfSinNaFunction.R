#' Extraccion vectorizada de variables ambientales para puntos de presencia
#'
#' @description
#' Extrae valores de variables ambientales desde un objeto `SpatRaster` para
#' cada punto de presencia en un objeto `sf`. Implementa un particionamiento
#' vectorizado nativo en C para segmentar la matriz por especie y eliminar registros
#' incompletos sin penalizaciones de memoria RAM.
#'
#' @param DataSf Objeto `sf` con geometrias de puntos (presencias).
#' @param RasterReferencia Objeto `SpatRaster` que contiene las variables ambientales continuas o categoricas.
#' @param EspecieValida Caracter. Nombre de la columna en `DataSf` que identifica a la especie.
#'   Por defecto es `"especievalida"`.
#' @param verbose Logico. Si es `TRUE` (por defecto), imprime en consola el progreso
#'   detallado del particionamiento espacial y limpieza.
#'
#' @return Una lista invisible de objetos `data.frame`. Cada elemento de la lista
#'   corresponde a una especie y contiene solo registros espaciales sin valores `NA`.
#'
#' @details
#' El flujo de trabajo interno optimizado sigue estos pasos:
#' \enumerate{
#'   \item \strong{Extraccion Espacial}: Ejecuta \code{terra::extract} de forma global sobre todos los puntos, maximizando el rendimiento.
#'   \item \strong{Estructuracion tabular}: Ensambla una matriz unica desvinculando la geometria espacial para acelerar el procesamiento logico.
#'   \item \strong{Particionamiento ($O(N)$)}: Utiliza \code{split} nativo para segmentar el dataset en submatrices independientes por especie en un solo escaneo.
#'   \item \strong{Filtro Algebraico}: Aplica logica booleana vectorizada (\code{vapply}) para descartar variables (columnas) nulas para especies especificas.
#'   \item \strong{Limpieza de Tuplas}: Intersecta el espacio de datos con \code{complete.cases} para retener estrictamente registros utiles para el modelado SDM.
#' }
#'
#' @note
#' Para ejemplos avanzados de ensamble de modelos de distribucion, consulte la viñeta
#' de la libreria: \code{vignette("spatial_cleaning_guide", package = "LibreriaMneRcsyPrc")}.
#'
#' @section Advertencias:
#' \itemize{
#'   \item \strong{Alineacion Topologica}: `DataSf` y `RasterReferencia` deben compartir el mismo CRS previo a la ejecucion.
#'   \item \strong{NAs Geograficos}: Puntos que caigan en el oceano o fuera del *Bounding Box* del raster seran eliminados automaticamente.
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' # 1. Crear datos de presencias simulados
#' puntos_sf <- st_as_sf(data.frame(
#'   lon = runif(100, -100, -90),
#'   lat = runif(100, 15, 25),
#'   especievalida = sample(c("SpA", "SpB"), 100, replace = TRUE)
#' ), coords = c("lon", "lat"), crs = 4326)
#'
#' # 2. Crear raster ambiental simulado
#' r_env <- rast(nrows=50, ncols=50, crs="EPSG:4326", ext=ext(-105, -85, 10, 30))
#' values(r_env) <- runif(ncell(r_env))
#' names(r_env) <- c("Bio1")
#'
#' # 3. Ejecutar extraccion optimizada
#' lista_presencias <- PresencesSfNoNaListDf(
#'   DataSf = puntos_sf,
#'   RasterReferencia = r_env,
#'   EspecieValida = "especievalida",
#'   verbose = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link[terra]{extract}}, \code{\link[stats]{complete.cases}}, \code{\link[base]{split}}
#'
#' @keywords sdm ambiental extraccion vectorizacion
#'
#' @importFrom terra extract
#' @importFrom sf st_coordinates st_drop_geometry
#' @importFrom stats complete.cases
#' @export
PresencesSfNoNaListDf <- function(DataSf,
                                  RasterReferencia,
                                  EspecieValida = "especievalida",
                                  verbose = TRUE) {
  if(verbose) message("\n--- Iniciando extraccion de variables ambientales ---")
  if(!inherits(DataSf, "sf")) stop("DataSf debe ser un objeto sf.")
  if(!inherits(RasterReferencia, "SpatRaster")) stop("RasterReferencia debe ser un SpatRaster.")

  if(verbose) message("-> Extrayendo valores del raster para todos los puntos...")
  ValoresRaster <- terra::extract(RasterReferencia, DataSf)

  Coords <- sf::st_coordinates(DataSf)
  Atributos <- sf::st_drop_geometry(DataSf)

  FullDf <- data.frame(
    longitud = Coords[, 1],
    latitud = Coords[, 2],
    Atributos,
    ValoresRaster[, names(ValoresRaster) != "ID", drop = FALSE],
    Presencia = 1
  )

  if(verbose) message("-> Procesando y limpiando datos por especie...")
  ListaSplit <- split(FullDf, FullDf[[EspecieValida]])

  Resultados <- lapply(names(ListaSplit), function(sp) {
    DfSp <- ListaSplit[[sp]]

    valido <- vapply(DfSp, function(x) !all(is.na(x)), logical(1))
    DfCleanCols <- DfSp[, valido, drop = FALSE]

    if(verbose) {
      eliminadas <- names(DfSp)[!valido]
      if(length(eliminadas) > 0) {
        message("  ! ", sp, " - Columnas todo NA eliminadas: ", paste(eliminadas, collapse = ", "))
      }
    }

    DfFinal <- DfCleanCols[stats::complete.cases(DfCleanCols), , drop = FALSE]

    if(verbose) {
      message("  > ", sp, " - Registros finales: ", nrow(DfFinal), " (", nrow(DfSp) - nrow(DfFinal), " eliminados)")
    }

    return(DfFinal)
  })

  names(Resultados) <- names(ListaSplit)

  if(verbose) message("--- Extraccion completada exitosamente ---\n")
  return(invisible(Resultados))
}
PresencesSfNoNaListDf <- function(DataSf, RasterReferencia, EspecieValida = "especievalida", verbose = TRUE) {
  if(verbose) message("\n--- Iniciando extraccion de variables ambientales ---")
  if(!inherits(DataSf, "sf")) stop("DataSf debe ser un objeto sf.")
  if(!inherits(RasterReferencia, "SpatRaster")) stop("RasterReferencia debe ser un SpatRaster.")

  if(verbose) message("-> Extrayendo valores del raster para todos los puntos...")
  ValoresRaster <- terra::extract(RasterReferencia, DataSf)

  Coords <- sf::st_coordinates(DataSf)
  Atributos <- sf::st_drop_geometry(DataSf)

  FullDf <- data.frame(
    longitud = Coords[, 1],
    latitud = Coords[, 2],
    Atributos,
    ValoresRaster[, names(ValoresRaster) != "ID", drop = FALSE],
    Presencia = 1
  )

  if(verbose) message("-> Procesando y limpiando datos por especie...")
  ListaSplit <- split(FullDf, FullDf[[EspecieValida]])

  Resultados <- lapply(names(ListaSplit), function(sp) {
    DfSp <- ListaSplit[[sp]]

    valido <- vapply(DfSp, function(x) !all(is.na(x)), logical(1))
    DfCleanCols <- DfSp[, valido, drop = FALSE]

    if(verbose) {
      eliminadas <- names(DfSp)[!valido]
      if(length(eliminadas) > 0) {
        message("  ! ", sp, " - Columnas todo NA eliminadas: ", paste(eliminadas, collapse = ", "))
      }
    }

    DfFinal <- DfCleanCols[stats::complete.cases(DfCleanCols), , drop = FALSE]

    if(verbose) {
      message("  > ", sp, " - Registros finales: ", nrow(DfFinal), " (", nrow(DfSp) - nrow(DfFinal), " eliminados)")
    }

    return(DfFinal)
  })

  names(Resultados) <- names(ListaSplit)

  if(verbose) message("--- Extraccion completada exitosamente ---\n")
  return(invisible(Resultados))
}
