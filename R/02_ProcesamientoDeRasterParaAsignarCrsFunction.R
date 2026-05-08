#' Asignar CRS y procesar objetos Raster
#'
#' @description
#' Proyecta un objeto `SpatRaster` a un sistema de referencia de coordenadas (CRS)
#' objetivo. La funcion permite estandarizar capas ambientales, asegurando que
#' coincidan con el CRS y la extension geografica de un area de estudio definida.
#'
#' @param RasterObject Objeto de clase `SpatRaster` que se desea procesar.
#' @param SpatialReference Definicion del destino. Puede ser un objeto espacial
#'   (`sf`, `sfc`, `SpatVector`) o una definicion directa de CRS (EPSG, WKT).
#' @param verbose Logico. Si es `TRUE` (por defecto), imprime progreso detallado.
#'
#' @return Un objeto `SpatRaster` proyectado y, opcionalmente, recortado.
#'
#' @details
#' El flujo de trabajo interno sigue estos pasos:
#' \enumerate{
#'   \item \strong{Configuracion}: Ajusta opciones de `terra` para optimizar el uso de RAM.
#'   \item \strong{Proyeccion}: Reproyecta el raster al CRS de destino usando interpolacion bilineal.
#'   \item \strong{Recorte}: Si se provee geometria, realiza el recorte (`crop`) y enmascarado (`mask`).
#' }
#'
#' @note
#' Para detalles sobre el manejo de memoria en rasters globales, consulte la viñeta:
#' \code{vignette("spatial_cleaning_guide", package = "TuLibreria")}.
#'
#' @section Advertencias:
#' \itemize{
#'   \item \strong{Memoria RAM}: Proyectar rasters de muy alta resolucion puede ser intensivo.
#'   \item \strong{Geometria}: El recorte solo se activa si `SpatialReference` posee una geometria valida.
#' }
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r <- rast(nrows=10, ncols=10, crs="EPSG:4326")
#' values(r) <- 1:ncell(r)
#' r_final <- ProcessingRasterToAssignCrs(r, "EPSG:32614")
#' }
#'
#' @seealso
#' \code{\link[terra]{project}}, \code{\link[terra]{crop}}, \code{\link[terra]{mask}}
#'
#' @keywords spatial raster projection
#'
#' @importFrom terra terraOptions project crop vect crs
#' @export
ProcessingRasterToAssignCrs <- function(RasterObject,
                                        SpatialReference,
                                        verbose = TRUE) {

  if(!inherits(RasterObject, "SpatRaster")) {
    stop("El objeto proporcionado en 'RasterObject' debe ser de clase SpatRaster (paquete terra).")
  }

  if(verbose) message("\n--- Iniciando procesamiento de Raster ---")

  is_spatial_obj <- inherits(SpatialReference, c("sf", "sfc", "SpatVector"))

  if(is_spatial_obj) {
    target_crs <- terra::crs(SpatialReference)
  } else if(is.character(SpatialReference) || is.numeric(SpatialReference)) {
    target_crs <- SpatialReference
  } else {
    stop("SpatialReference debe ser un objeto espacial (sf/SpatVector) o una definicion de CRS valida.")
  }

  if(verbose) message("-> Proyectando raster al CRS destino...")
  raster_projected <- terra::project(RasterObject, target_crs)

  if(is_spatial_obj) {
    if(verbose) message("-> Recortando y enmascarando raster segun la referencia...")

    if(!inherits(SpatialReference, "SpatVector")) {
      vect_ref <- terra::vect(SpatialReference)
    } else {
      vect_ref <- SpatialReference
    }

    raster_projected <- terra::crop(raster_projected, vect_ref, mask = TRUE)
  }

  if(verbose) message("--- Procesamiento raster completado ---\n")

  return(raster_projected)
}
