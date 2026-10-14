#' Asignar CRS y procesar objetos Raster de forma optimizada
#'
#' @description
#' Proyecta un objeto `SpatRaster` a un sistema de referencia de coordenadas (CRS)
#' objetivo. La funcion implementa una arquitectura optimizada de "recorte previo"
#' (crop-first) que reduce drasticamente el uso de memoria RAM al estandarizar
#' capas ambientales globales, asegurando ademas el empalme perfecto de bordes.
#'
#' @param RasterObject Objeto de clase `SpatRaster` que se desea procesar.
#' @param SpatialReference Definicion del destino. Puede ser un objeto espacial
#'   (`sf`, `sfc`, `SpatVector`) o una definicion directa de CRS (EPSG, WKT).
#' @param MetodoInterpolacion Caracter. Metodo matematico para recalcular los valores de los
#'   pixeles al proyectar. Use `"bilinear"` (por defecto) para variables continuas
#'   (temperatura, elevacion) y `"near"` (vecino mas cercano) para variables categoricas
#'   (cobertura de suelo, tipos de vegetacion).
#' @param verbose Logico. Si es `TRUE` (por defecto), imprime en consola el progreso
#'   detallado de cada etapa.
#'
#' @return Un objeto `SpatRaster` proyectado y adaptado a la extension exacta de la referencia.
#'
#' @details
#' El flujo de trabajo interno resuelve problemas comunes de los SIG mediante estos pasos:
#' \enumerate{
#'   \item \strong{Traduccion de Geometria}: Extrae el CRS objetivo y transforma temporalmente la referencia al CRS original del raster.
#'   \item \strong{Recorte Previo (Optimizacion)}: Realiza un recorte burdo (Bounding Box) del raster global usando un margen de seguridad (`snap = "out"`). Esto evita procesar millones de pixeles innecesarios en la proyeccion.
#'   \item \strong{Reproyeccion}: Aplica la transformacion matematica (proyeccion) exclusivamente sobre el area recortada, usando el `MetodoInterpolacion` elegido.
#'   \item \strong{Enmascarado Fino}: Realiza un recorte exacto (`mask`) sobre el raster ya proyectado para garantizar que los bordes empaten perfectamente sin generar NAs.
#' }
#'
#' @note
#' Para obtener mas detalles sobre la integracion de variables ambientales, consulte la viñeta:
#' \code{vignette("spatial_cleaning_guide", package = "LibreriaMneRcsyPrc")}.
#'
#' @section Advertencias:
#' \itemize{
#'   \item \strong{Variables Categoricas}: Si intenta usar `"bilinear"` en rasters de clases (ej. uso de suelo), los datos se corromperan al inventar clases intermedias. Cambie siempre a `"near"`.
#'   \item \strong{Eficiencia de RAM}: Gracias a la estrategia "Crop-First", el consumo de RAM se minimiza automaticamente sin necesidad de forzar configuraciones complejas en `terraOptions`.
#' }
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(sf)
#'
#' # 1. Crear raster global sintetico (variable continua)
#' r_temp <- rast(nrows=180, ncols=360, crs="EPSG:4326")
#' values(r_temp) <- runif(ncell(r_temp), 10, 35)
#'
#' # 2. Definir area de estudio (ej. un estado en coordenadas proyectadas)
#' area_estudio <- st_sfc(st_polygon(list(matrix(
#'    c(-100, 20, -99, 20, -99, 21, -100, 21, -100, 20),
#'    ncol = 2, byrow = TRUE))), crs = 4326) %>%
#'    st_transform(32614)
#'
#' # 3. Procesar raster optimizadamente
#' r_final <- ProcessingRasterToAssignCrs(
#'    RasterObject = r_temp,
#'    SpatialReference = area_estudio,
#'    MetodoInterpolacion = "bilinear",
#'    verbose = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link[terra]{project}}, \code{\link[terra]{crop}}, \code{\link[terra]{mask}}
#'
#' @keywords spatial raster projection
#'
#' @importFrom terra project crop vect crs
#' @export
ProcessingRasterToAssignCrs <- function(RasterObject,
                                        SpatialReference,
                                        MetodoInterpolacion = "bilinear",
                                        verbose = TRUE) {

  if(!inherits(RasterObject, "SpatRaster")) {
    stop("Error: 'RasterObject' debe ser de clase SpatRaster (paquete terra).")
  }
  if(verbose) message("\n--- Iniciando procesamiento de Raster ---")
  is_spatial_obj <- inherits(SpatialReference, c("sf", "sfc", "SpatVector"))
  if(is_spatial_obj) {
    if(!inherits(SpatialReference, "SpatVector")) {
      vect_ref <- terra::vect(SpatialReference)
    } else {
      vect_ref <- SpatialReference
    }
    target_crs <- terra::crs(vect_ref)
  } else if(is.character(SpatialReference) || is.numeric(SpatialReference)) {
    target_crs <- SpatialReference
  } else {
    stop("Error: SpatialReference debe ser un objeto espacial (sf/SpatVector) o una definicion de CRS valida.")
  }
  if(is_spatial_obj) {
    if(verbose) message("1/3 Recortando raster original para optimizar memoria...")
    orig_crs <- terra::crs(RasterObject)
    vect_ref_orig_crs <- terra::project(vect_ref, orig_crs)
    RasterObject <- terra::crop(RasterObject, vect_ref_orig_crs, snap = "out")
  }
  if(verbose) message(paste0("2/3 Proyectando raster al CRS destino (Metodo: '", MetodoInterpolacion, "')..."))
  raster_projected <- terra::project(RasterObject, target_crs, method = MetodoInterpolacion)
  if(is_spatial_obj) {
    if(verbose) message("3/3 Ajustando extension exacta y enmascarando...")
    raster_projected <- terra::crop(raster_projected, vect_ref, mask = TRUE)
  }
  if(verbose) message("--- Procesamiento raster completado con exito ---\n")
  return(raster_projected)
}
