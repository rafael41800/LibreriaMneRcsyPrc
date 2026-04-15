#' @title Asignar CRS y procesar raster / Assign CRS and Process Raster
#'
#' @description
#' Proyecta un raster a un sistema de referencia de coordenadas (CRS) objetivo
#' y opcionalmente lo recorta a un área poligonal.
#' /
#' Projects a raster to a target Coordinate Reference System (CRS)
#' and optionally clips it to a polygon area.
#'
#' @param RasterObject Objeto SpatRaster a ser transformado.
#'   / SpatRaster object to be transformed.
#'
#' @param SpatialReference Puede ser: un objeto espacial (sf/sfc/SpatVector)
#'   para extraer el CRS y opcionalmente recortar, o una definición de CRS
#'   (carácter/WKT/código EPSG).
#'   / Either: a spatial object (sf/sfc/SpatVector) to extract CRS and optionally clip,
#'   or a CRS definition (character/WKT/EPSG code).
#'
#' @param verbose Lógico. Si `TRUE` (por defecto), muestra mensajes de progreso.
#'   / Logical. If `TRUE` (default), displays progress messages.
#'
#' @return
#' Un objeto SpatRaster proyectado (y opcionalmente recortado).
#' / A projected (and optionally clipped) SpatRaster object.
#'
#' @details
#' La función realiza las siguientes operaciones:
#' \enumerate{
#'   \item Configura directorio temporal para manejo de memoria
#'   \item Valida que el objeto de entrada sea un SpatRaster
#'   \item Detecta el tipo de referencia espacial
#'   \item Proyecta el raster al CRS objetivo
#'   \item Si se proporciona un objeto espacial, recorta el raster al área del polígono
#' }
#' /
#' The function performs the following operations:
#' \enumerate{
#'   \item Sets temporary directory for memory management
#'   \item Validates that input is a SpatRaster object
#'   \item Detects the type of spatial reference
#'   \item Projects the raster to the target CRS
#'   \item If a spatial object is provided, clips the raster to the polygon area
#' }
#'
#' @note
#' \itemize{
#'   \item El directorio temporal `~/temp` se crea automáticamente si no existe
#'   \item Para objetos sf/sfc, se convierten automáticamente a SpatVector
#'   \item El recorte utiliza `mask = TRUE` para enmascarar celdas fuera del polígono
#' }
#' /
#' \itemize{
#'   \item The temporary directory `~/temp` is automatically created if it does not exist
#'   \item For sf/sfc objects, they are automatically converted to SpatVector
#'   \item Clipping uses `mask = TRUE` to mask cells outside the polygon
#' }
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r <- rast(nrows=10, ncols=10, xmin=-10, xmax=10, ymin=-10, ymax=10)
#' values(r) <- 1:ncell(r)
#'
#' # Proyectar a un CRS específico / Project to a specific CRS
#' r_proj <- ProcessingRasterToAssignCrs(r, "EPSG:32614")
#'
#' # Proyectar y recortar con un polígono / Project and clip with a polygon
#' library(sf)
#' poly <- st_as_sfc("POLYGON((-5 -5, 5 -5, 5 5, -5 5, -5 -5))", crs = 4326)
#' r_clip <- ProcessingRasterToAssignCrs(r, poly)
#' }
#'
#' @seealso
#' \code{\link[terra]{project}}, \code{\link[terra]{crop}}, \code{\link[terra]{crs}}
#'
#' @importFrom terra terraOptions project crop crs vect
#' @export
ProcessingRasterToAssignCrs <- function(RasterObject, SpatialReference, verbose = TRUE) {
  if(!dir.exists("~/temp")) dir.create("~/temp")
  terra::terraOptions(tempdir = "~/temp", memfrac = 0.8)
  if(!inherits(RasterObject, "SpatRaster")) {
    stop("Esta funcion solo acepta objetos SpatRaster")
  }
  if(verbose) cat("Procesando raster...\n")
  if(inherits(SpatialReference, c("sf", "sfc", "SpatVector"))) {
    target_crs <- terra::crs(SpatialReference)
  } else if(is.character(SpatialReference) || is.numeric(SpatialReference)) {
    target_crs <- SpatialReference
  } else {
    stop("SpatialReference debe ser un objeto espacial o CRS valido")
  }
  if(verbose) cat("  Proyectando raster al CRS de referencia...\n")
  raster_projected <- terra::project(RasterObject, target_crs)
  if(inherits(SpatialReference, c("sf", "sfc", "SpatVector"))) {
    if(verbose) cat("  Recortando raster al area del poligono...\n")
    vect_reference <- terra::vect(SpatialReference)
    raster_masked <- terra::crop(raster_projected, vect_reference, mask = TRUE)
    if(verbose) cat("  Raster recortado exitosamente\n")
    return(raster_masked)
  }
  if(verbose) cat("  Procesamiento completado\n")
  return(raster_projected)
}
