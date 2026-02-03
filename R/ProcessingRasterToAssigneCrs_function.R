ProcessingRasterToAssigneCrs <- function(RasterObject, SpatialReference){
  if(!dir.exists("~/temp")) dir.create("~/temp")
  terra::terraOptions(tempdir = "~/temp", memfrac = 0.8)
  if(!inherits(RasterObject, "SpatRaster")) {
    stop("Esta función solo acepta objetos SpatRaster")
  }
  cat("✓ Procesando raster... \n")
  if(inherits(SpatialReference, c("sf", "sfc", "SpatVector"))) {
    target_crs <- terra::crs(SpatialReference)
  } else if(is.character(SpatialReference) || is.numeric(SpatialReference)) {
    target_crs <- SpatialReference
  } else {
    stop("SpatialReference debe ser un objeto espacial o CRS válido")
  }
  cat("  Proyectando raster al CRS de referencia... \n")
  raster_projected <- terra::project(RasterObject, target_crs)

  if(inherits(SpatialReference, c("sf", "sfc", "SpatVector"))) {
    cat("  Recortando raster al área del polígono... \n")
    vect_reference <- terra::vect(SpatialReference)
    raster_cropped <- terra::crop(raster_projected, vect_reference) #, mask = TRUE)
    raster_masked <- terra::mask(raster_cropped, vect_reference) # se omitiria esta linea con lo anterir
    cat("✓ Raster recortado al área del polígono \n")
    return(raster_masked)
  }
  return(raster_projected)
}
