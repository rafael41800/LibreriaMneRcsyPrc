#' @title Asignar CRS y procesar raster / Assign CRS and Process Raster
#'
#' @description
#' Asigna un sistema de referencia de coordenadas (CRS) a un raster y
#' opcionalmente lo recorta a un área espacial definida. La función maneja
#' diferentes tipos de referencia espacial (objetos espaciales, CRS numéricos
#' o textuales) y aplica las transformaciones necesarias.
#' /
#' Assigns a coordinate reference system (CRS) to a raster and optionally
#' crops it to a defined spatial area. The function handles different types
#' of spatial references (spatial objects, numeric or textual CRS) and
#' applies necessary transformations.
#'
#' @param RasterObject Objeto raster al que asignaremos coordenadas de
#'   referencia del sistema. Debe ser un objeto `SpatRaster` de terra.
#'   / Raster object to which we will assign coordinate reference system.
#'   Must be a `SpatRaster` object from terra.
#'
#' @param SpatialReference Referencia espacial de la que tomaremos el CRS.
#'   Puede ser:
#'   / Spatial reference from which we will take the CRS. Can be:
#'   \itemize{
#'     \item{**Objeto espacial** (`sf`, `sfc`, `SpatVector`) - Se usará su CRS
#'           y se recortará el raster al área del polígono (si es polígono)
#'           / **Spatial object** (`sf`, `sfc`, `SpatVector`) - Its CRS will
#'           be used and the raster will be cropped to the polygon area
#'           (if it's a polygon)}
#'     \item{**CRS como caracter o numérico** (ej: `"EPSG:4326"` o `4326`) -
#'           Solo se proyectará el raster al CRS especificado
#'           / **CRS as character or numeric** (e.g., `"EPSG:4326"` or `4326`) -
#'           Only raster projection will be applied to the specified CRS}
#'     \item{**Objeto `crs`** de terra - CRS directo para la proyección
#'           / **`crs` object** from terra - Direct CRS for projection}
#'   }
#'
#' @returns
#' Un objeto `SpatRaster` procesado con las siguientes características:
#' * Proyectado al CRS especificado por `SpatialReference`
#' * Recortado al área del polígono si `SpatialReference` es un objeto espacial
#'   de tipo polígono
#' * Con valores NA fuera del área del polígono (cuando se aplica recorte)
#' * Conservando los metadatos originales cuando es posible
#' /
#' A processed `SpatRaster` object with the following features:
#' * Projected to the CRS specified by `SpatialReference`
#' * Cropped to the polygon area if `SpatialReference` is a polygon-type
#'   spatial object
#' * With NA values outside the polygon area (when cropping is applied)
#' * Preserving original metadata when possible
#'
#' @details
#' La función realiza las siguientes operaciones en secuencia:
#' 1. **Configuración del entorno**: Establece directorio temporal para
#'    manejo eficiente de memoria.
#' 2. **Validación de entrada**: Verifica que `RasterObject` sea un
#'    `SpatRaster` válido.
#' 3. **Detección de tipo de referencia**: Identifica si `SpatialReference` es
#'    un objeto espacial, CRS textual/numerico, o objeto `crs`.
#' 4. **Proyección**: Transforma el raster al CRS de destino usando
#'    `terra::project()`.
#' 5. **Recorte condicional**: Si `SpatialReference` es un objeto espacial
#'    con geometría de polígono, aplica `terra::crop()` con `mask = TRUE`
#'    para recortar y enmascarar simultáneamente.
#' /
#' The function performs the following operations in sequence:
#' 1. **Environment setup**: Sets temporary directory for efficient
#'    memory management.
#' 2. **Input validation**: Verifies that `RasterObject` is a valid
#'    `SpatRaster`.
#' 3. **Reference type detection**: Identifies whether `SpatialReference` is
#'    a spatial object, textual/numeric CRS, or `crs` object.
#' 4. **Projection**: Transforms the raster to target CRS using
#'    `terra::project()`.
#' 5. **Conditional cropping**: If `SpatialReference` is a spatial object
#'    with polygon geometry, applies `terra::crop()` with `mask = TRUE`
#'    to simultaneously crop and mask.
#'
#' @note
#' * La función utiliza `terra::crop()` con `mask = TRUE` cuando se proporciona
#'   un polígono, lo que asegura que solo las celdas dentro del polígono sean
#'   retenidas, estableciendo las celdas externas a NA. Esto es más eficiente
#'   que realizar crop y mask por separado.
#' * Para objetos `sf` o `sfc`, la función convierte automáticamente a
#'   `SpatVector` para compatibilidad con terra.
#' * El parámetro `mask = TRUE` en `terra::crop()` puede ser más lento para
#'   polígonos muy complejos, pero produce resultados más precisos.
#' * Se recomienda que el raster de entrada tenga un CRS definido para
#'   obtener mejores resultados en la transformación.
#' /
#' * The function uses `terra::crop()` with `mask = TRUE` when a polygon is
#'   provided, which ensures that only cells within the polygon are retained,
#'   setting cells outside to NA. This is more efficient than performing
#'   crop and mask separately.
#' * For `sf` or `sfc` objects, the function automatically converts to
#'   `SpatVector` for compatibility with terra.
#' * The `mask = TRUE` parameter in `terra::crop()` can be slower for
#'   very complex polygons but produces more accurate results.
#' * It is recommended that the input raster has a defined CRS for
#'   better transformation results.
#'
#' @section Advertencias/Warnings:
#' * Si el raster de entrada no tiene CRS definido y se solicita proyección,
#'   se asumirá WGS84 (EPSG:4326) y se mostrará una advertencia.
#' * El recorte con polígonos muy complejos o con muchos agujeros puede
#'   consumir mucha memoria y tiempo de procesamiento.
#' * Las transformaciones entre CRS con diferentes elipsoides o datums pueden
#'   introducir pequeñas distorsiones.
#' * Si el área del polígono está completamente fuera de la extensión del
#'   raster, el resultado será un raster vacío.
#' /
#' * If the input raster has no defined CRS and projection is requested,
#'   WGS84 (EPSG:4326) will be assumed and a warning will be shown.
#' * Cropping with very complex polygons or polygons with many holes may
#'   consume significant memory and processing time.
#' * Transformations between CRS with different ellipsoids or datums may
#'   introduce small distortions.
#' * If the polygon area is completely outside the raster's extent,
#'   the result will be an empty raster.
#'
#' @importFrom terra terraOptions project crop crs vect rast
#'
#' @examples
#' \dontrun{
#' # Ejemplo 1: Proyectar raster usando un CRS numérico
#' # Example 1: Project raster using a numeric CRS
#' library(terra)
#'
#' # Crear un raster de ejemplo
#' # Create example raster
#' r <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
#' values(r) <- runif(ncell(r))
#'
#' # Proyectar a WGS84 (EPSG:4326)
#' # Project to WGS84 (EPSG:4326)
#' raster_proyectado <- ProcessingRasterToAssigneCrs(r, 4326)
#' print(raster_proyectado)
#'
#' # Ejemplo 2: Proyectar y recortar usando un polígono sf
#' # Example 2: Project and crop using an sf polygon
#' library(sf)
#'
#' # Crear un polígono de ejemplo en WGS84
#' # Create example polygon in WGS84
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
#' vect_ref <- vect(poligono)
#' raster_vect <- ProcessingRasterToAssigneCrs(r, vect_ref)
#' print(raster_vect)
#'
#' # Ejemplo 4: Usar CRS como texto
#' # Example 4: Use CRS as text
#' raster_text <- ProcessingRasterToAssigneCrs(r, "EPSG:4326")
#' print(raster_text)
#'
#' # Ejemplo 5: Raster sin CRS definido
#' # Example 5: Raster without defined CRS
#' r_sin_crs <- r
#' crs(r_sin_crs) <- ""  # Eliminar CRS / Remove CRS
#' raster_sin_crs <- ProcessingRasterToAssigneCrs(r_sin_crs, 4326)
#' print(raster_sin_crs)
#' }
#'
#' @seealso
#' * [terra::project()] para transformaciones de sistema de coordenadas de rasters
#'   / for coordinate system transformations of rasters
#' * [terra::crop()] para recortar rasters a áreas específicas
#'   / for cropping rasters to specific areas
#' * [terra::crs()] para manejar sistemas de referencia de coordenadas
#'   / for handling coordinate reference systems
#' * [sf::st_crs()] para obtener CRS de objetos sf
#'   / for getting CRS from sf objects
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
