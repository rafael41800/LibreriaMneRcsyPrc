#' Preparación de datos para hacerlos datos geoespaciales
#'
#' @param DataFrame Conjunto de datos a procesar.
#' @param SpatialReference Referencia espacial a asignar.
#' @param latitude Punto a partir del cual se considerará la latitud, por default es 0.
#' @param variables Variables que son de interes, en forma de vector, con entradas de tipo string.
#'
#' @returns Conjunto de datos listos para usarse como datos geoespaciales
#'
#' @importFrom terra terraOptions
#' @importFrom dplyr filter select mutate distinct
#' @importFrom sf st_as_sf st_transform st_crop st_intersection
#'
#' @examples
#'
#' @export
CleaningTransformingDf <- function(DataFrame, SpatialReference, latitude=0, variables) {
  if(!dir.exists("~/temp")) dir.create("~/temp")
  terra::terraOptions(tempdir=("~/temp"), memfrac=0.8)
  DF.north <- DataFrame %>%
    dplyr::filter(latitud >= latitude) %>%
    dplyr::select(variables) %>%
    sf::st_as_sf(coords = c("longitud", "latitud"), crs = 4326)
  transformation <- sf::st_transform(DF.north, st_crs(SpatialReference))
  crop <- sf::st_crop(transformation, SpatialReference)
  mask <- sf::st_intersection(crop, SpatialReference)
  DF.north.clean <- mask %>%
    dplyr::mutate(XY = paste0(st_coordinates(.)[,1], st_coordinates(.)[,2], especie)) %>%
    dplyr::distinct(XY, .keep_all = TRUE) %>%
    dplyr::select(-XY)
  return(DF.north.clean)
}
