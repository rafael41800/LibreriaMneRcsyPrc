#' Asignar CRS y procesar raster
#' Assign CRS and process raster
#'
#' @param RasterObject Objeto raster al que asignaremos coordenadas de referencia del sistema / Raster object to which we will assign coordinate reference system.
#' @param SpatialReference Referencia espacial de la que tomaremos el CRS. Puede ser: / Spatial reference from which we will take the CRS. Can be:
#' \itemize{
#'   \item{Objeto espacial (sf, sfc, SpatVector) - Se usara su CRS y se recortara el raster al area del poligono / Spatial object - Its CRS will be used and the raster will be cropped to the polygon area}
#'   \item{CRS como caracter o numerico (ej: "EPSG:4326" o 4326) - Solo se proyectara el raster / CRS as character or numeric (e.g., "EPSG:4326" or 4326) - Only raster projection will be applied}
#' }
#'
#' @returns Un objeto `SpatRaster` procesado que puede ser: / A processed `SpatRaster` object that can be:
#' \itemize{
#'   \item{Proyectado al CRS especificado / Projected to the specified CRS}
#'   \item{Proyectado y recortado al area del poligono (si `SpatialReference` es un objeto espacial) / Projected and cropped to the polygon area (if `SpatialReference` is a spatial object)}
#' }
#'
#' @details
#' Esta funcion realiza las siguientes operaciones:
#' 1. Configura el directorio temporal para manejo eficiente de memoria
#' 2. Valida que el objeto de entrada sea un `SpatRaster`
#' 3. Detecta el tipo de referencia espacial proporcionada
#' 4. Proyecta el raster al CRS de referencia
#' 5. Si la referencia es un objeto espacial (poligono), recorta el raster al area del poligono usando `mask = TRUE`
#'
#' @note
#' La funcion utiliza `terra::crop()` con `mask = TRUE` cuando se proporciona un poligono,
#' lo que asegura que solo las celdas dentro del poligono sean retenidas, estableciendo
#' las celdas externas a NA. Esto es mas eficiente que hacer crop y mask por separado.
#'
#' @importFrom terra terraOptions project crop crs vect
#'
#' @examples
#' \dontrun{
#' # Ejemplo 1: Proyectar raster usando un CRS numerico
#' # Example 1: Project raster using a numeric CRS
#' library(terra)
#'
#' # Crear un raster de ejemplo
#' # Create example raster
#' r <- rast(nrows=100, ncols=100, xmin=0, xmax=10, ymin=0, ymax=10)
#' values(r) <- runif(ncell(r))
#'
#' # Proyectar a WGS84 (EPSG:4326)
#' # Project to WGS84 (EPSG:4326)
#' raster_proyectado <- ProcessingRasterToAssigneCrs(r, 4326)
#' print(raster_proyectado)
#'
#' # Ejemplo 2: Proyectar y recortar usando un poligono
#' # Example 2: Project and crop using a polygon
#' library(sf)
#'
#' # Crear un poligono de ejemplo
#' # Create example polygon
#' poligono <- st_as_sfc("POLYGON((2 2, 8 2, 8 8, 2 8, 2 2))", crs = 4326)
#'
#' # Proyectar y recortar el raster
#' # Project and crop the raster
#' raster_recortado <- ProcessingRasterToAssigneCrs(r, poligono)
#' print(raster_recortado)
#'
#' # Visualizar resultados
#' # Visualize results
#' plot(raster_recortado)
#' plot(st_geometry(poligono), add = TRUE, border = "red", lwd = 2)
#'
#' # Ejemplo 3: Usar un objeto SpatVector como referencia
#' # Example 3: Use a SpatVector object as reference
#' vect_ref <- as(poligono, "SpatVector")
#' raster_vect <- ProcessingRasterToAssigneCrs(r, vect_ref)
#' print(raster_vect)
#'
#' # Ejemplo 4: Usar CRS como texto
#' # Example 4: Use CRS as text
#' raster_text <- ProcessingRasterToAssigneCrs(r, "EPSG:4326")
#' print(raster_text)
#' }
#'
#' @export
ProcessingRasterToAssigneCrs <- function(RasterObject, SpatialReference){
  if(!dir.exists("~/temp")) dir.create("~/temp")
  terra::terraOptions(tempdir = "~/temp", memfrac = 0.8)
  if(!inherits(RasterObject, "SpatRaster")) {
    stop("Esta funcion solo acepta objetos SpatRaster")
  }
  cat("Procesando raster... \n")
  if(inherits(SpatialReference, c("sf", "sfc", "SpatVector"))) {
    target_crs <- terra::crs(SpatialReference)
  } else if(is.character(SpatialReference) || is.numeric(SpatialReference)) {
    target_crs <- SpatialReference
  } else {
    stop("SpatialReference debe ser un objeto espacial o CRS valido")
  }
  cat("  Proyectando raster al CRS de referencia... \n")
  raster_projected <- terra::project(RasterObject, target_crs)

  if(inherits(SpatialReference, c("sf", "sfc", "SpatVector"))) {
    cat("  Recortando raster al area del poligono... \n")
    vect_reference <- terra::vect(SpatialReference)
    raster_masked <- terra::crop(raster_projected, vect_reference, mask = TRUE)
    cat("Raster recortado al area del poligono \n")
    return(raster_masked)
  }
  return(raster_projected)
}
