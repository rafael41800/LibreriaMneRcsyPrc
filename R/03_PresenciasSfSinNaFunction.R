#' Extraccion de variables ambientales para puntos de presencia
#'
#' @description
#' Extrae valores de variables ambientales desde un objeto `SpatRaster` para
#' cada punto de presencia en un objeto `sf`. Organiza los datos por especie
#' y elimina registros incompletos.
#'
#' @param DataSf Objeto `sf` con geometrias de puntos.
#' @param RasterReferencia Objeto `SpatRaster` con las variables ambientales.
#' @param EspecieValida Nombre de la columna de especie. Por defecto `"especievalida"`.
#' @param verbose Logico. Si es `TRUE`, muestra mensajes de progreso.
#'
#' @return Una lista invisible de `data.frame` segmentados por especie.
#'
#' @details
#' \strong{Limpieza de datos}: La funcion elimina variables (columnas) que son
#' totalmente `NA` para una especie y luego aplica \code{complete.cases} para
#' asegurar que no existan registros vacios antes del modelado.
#'
#' @section Advertencias:
#' \itemize{
#'   \item \strong{Alineacion}: Los objetos deben compartir el mismo CRS; de lo contrario, la extraccion fallara o sera erronea.
#'   \item \strong{NAs Geograficos}: Puntos fuera de la extension del raster seran eliminados.
#' }
#'
#' @seealso
#' \code{\link[terra]{extract}}, \code{\link[stats]{complete.cases}}
#'
#' @keywords sdm ambiental extraccion
#'
#' @importFrom terra extract
#' @importFrom sf st_coordinates st_drop_geometry
#' @importFrom dplyr select mutate where
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
    ValoresRaster[, !names(ValoresRaster) == "ID", drop = FALSE],
    Presencia = 1
  )

  ListaEspecies <- unique(as.character(FullDf[[EspecieValida]]))
  Resultados <- list()

  for (sp in ListaEspecies) {
    if(verbose) message("  - Procesando: ", sp)

    DfSp <- FullDf[FullDf[[EspecieValida]] == sp, ]

    DfCleanCols <- DfSp %>% dplyr::select(where(~ !all(is.na(.))))

    cols_eliminadas <- setdiff(names(DfSp), names(DfCleanCols))
    if(verbose && length(cols_eliminadas) > 0) {
      message("    ! Columnas eliminadas por ser todo NA: ", paste(cols_eliminadas, collapse = ", "))
    }

    DfFinal <- DfCleanCols[stats::complete.cases(DfCleanCols), ]

    if(verbose) {
      perdidos <- nrow(DfSp) - nrow(DfFinal)
      message("    > Registros finales: ", nrow(DfFinal), " (", perdidos, " eliminados por NAs)")
    }

    Resultados[[sp]] <- DfFinal
  }

  if(verbose) message("--- Extraccion completada exitosamente ---\n")

  return(invisible(Resultados))
}
