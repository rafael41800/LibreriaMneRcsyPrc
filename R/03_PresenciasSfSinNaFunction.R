#' @title Presencias desde un sf sin entradas NA / Presence Data from sf without NA Entries
#'
#' @description
#' Extrae y procesa datos de presencia de especies desde un objeto `sf`,
#' combinándolos con valores de variables ambientales de un raster y
#' eliminando todos los registros incompletos (con valores NA).
#' /
#' Extracts and processes species presence data from an `sf` object,
#' combining them with environmental variable values from a raster and
#' removing all incomplete records (with NA values).
#'
#' @param DataSf Objeto `sf` con datos de presencia de especies. Debe contener
#'   geometrías puntuales y la columna especificada en `EspecieValida`.
#'   / `sf` object with species presence data. Must contain point geometries
#'   and the column specified in `EspecieValida`.
#'
#' @param RasterReferenciaReproyectado Objeto `SpatRaster` proyectado con
#'   las variables ambientales.
#'   / Projected `SpatRaster` object with environmental variables.
#'
#' @param EspecieValida Nombre de la columna que contiene los identificadores
#'   de especie. Por defecto `"especievalida"`.
#'   / Name of the column containing species identifiers. Default is `"especievalida"`.
#'
#' @param verbose Lógico. Si `TRUE` muestra mensajes de progreso. Por defecto `TRUE`.
#'   / Logical. If `TRUE` shows progress messages. Default is `TRUE`.
#'
#' @return
#' Una lista (retornada invisiblemente) con un dataframe por cada especie.
#' Cada dataframe contiene coordenadas, atributos originales, valores del raster
#' y una columna `Presencia = 1`. Los dataframes han sido limpiados eliminando
#' columnas completamente NA y filas con cualquier NA.
#' /
#' A list (returned invisibly) with one dataframe per species.
#' Each dataframe contains coordinates, original attributes, raster values,
#' and a `Presencia = 1` column. Dataframes are cleaned by removing
#' completely NA columns and rows with any NA.
#'
#' @details
#' La función procesa cada especie: filtra registros, extrae coordenadas y valores
#' raster, elimina columnas completamente NA, luego filas con NA, y añade
#' `Presencia = 1`. El orden de limpieza (columnas luego filas) maximiza la
#' retención de datos válidos.
#' /
#' The function processes each species: filters records, extracts coordinates
#' and raster values, removes completely NA columns, then rows with NA, and adds
#' `Presencia = 1`. The cleaning order (columns then rows) maximizes retention
#' of valid data.
#'
#' @note
#' * `DataSf` y `RasterReferenciaReproyectado` deben estar en el mismo CRS.
#' * Las geometrías deben ser puntos (`POINT`).
#' * La columna `EspecieValida` debe existir en `DataSf`.
#' /
#' * `DataSf` and `RasterReferenciaReproyectado` must be in the same CRS.
#' * Geometries must be points (`POINT`).
#' * The `EspecieValida` column must exist in `DataSf`.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' # Ejemplo básico / Basic example
#' datos <- st_as_sf(data.frame(
#'   especie = c("A", "B"),
#'   x = c(-100, -101),
#'   y = c(20, 21)
#' ), coords = c("x", "y"), crs = 4326)
#'
#' raster_ejemplo <- rast(matrix(1:100, 10, 10))
#' crs(raster_ejemplo) <- "EPSG:4326"
#'
#' resultados <- PresencesSfNoNaListDf(datos, raster_ejemplo, "especie")
#' }
#'
#' @seealso
#' \code{\link[terra]{extract}}, \code{\link[sf]{st_coordinates}},
#' \code{\link[stats]{complete.cases}}
#'
#' @importFrom terra extract
#' @importFrom sf st_coordinates st_drop_geometry
#' @importFrom dplyr select %>% where
#' @export
PresencesSfNoNaListDf <- function(DataSf, RasterReferenciaReproyectado, EspecieValida = "especievalida", verbose = TRUE) {
  if(verbose) cat("Procesando presencias...\n")
  ListaEspecies <- unique(DataSf[[EspecieValida]])
  resultados <- list()
  for (l in ListaEspecies) {
    if(verbose) cat("  Procesando especie:", as.character(l), "\n")
    SubsetSfMamOaxEspecie <- DataSf[DataSf[[EspecieValida]] == as.character(l), ]
    VariablesStackRaster <- terra::extract(RasterReferenciaReproyectado, SubsetSfMamOaxEspecie)
    coords <- sf::st_coordinates(SubsetSfMamOaxEspecie)
    atributos <- sf::st_drop_geometry(SubsetSfMamOaxEspecie)
    DfPresencias <- data.frame(
      longitud = coords[, "X"],
      latitud = coords[, "Y"],
      atributos,
      VariablesStackRaster[, !names(VariablesStackRaster) == "ID", drop = FALSE],
      Presencia = 1
    )
    PresenciasParcialNa <- DfPresencias %>%
      dplyr::select(!where(~ all(is.na(.))))
    VariablesPuroNA <- names(DfPresencias)[!names(DfPresencias) %in% names(PresenciasParcialNa)]
    DfPresenciasSinNa <- PresenciasParcialNa[complete.cases(PresenciasParcialNa), ]
    resultados[[as.character(l)]] <- DfPresenciasSinNa
    if(verbose && length(VariablesPuroNA) > 0) {
      cat("    Eliminadas columnas:", paste(VariablesPuroNA, collapse = ", "), "\n")
    }
  }
  if(verbose) cat("  Procesamiento completado\n")
  return(invisible(resultados))
}
