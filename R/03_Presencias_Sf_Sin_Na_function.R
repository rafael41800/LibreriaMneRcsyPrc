#' Presencias desde un sf sin entradas NA's
#' Presence data from sf without NA entries
#'
#' @param Data_Sf Objeto `sf` con datos de presencia de especies / `sf` object with species presence data.
#' @param raster_referencia_reproyectado Objeto `SpatRaster` proyectado con las variables ambientales / Projected `SpatRaster` object with environmental variables.
#' @param Especie_Valida Nombre de la columna que contiene las especies / Name of the column containing species.
#'
#' @returns Una lista con dataframes para cada especie que contienen: / A list with dataframes for each species containing:
#' \itemize{
#'   \item{Coordenadas (longitud, latitud) / Coordinates (longitude, latitude)}
#'   \item{Atributos originales del objeto sf / Original attributes from the sf object}
#'   \item{Valores extraidos del raster / Values extracted from the raster}
#'   \item{Variable "Presencia" con valor 1 / "Presence" variable with value 1}
#' }
#' Todos los dataframes tienen filas sin valores NA / All dataframes have rows without NA values.
#'
#' @details
#' Esta funcion realiza las siguientes operaciones:
#' 1. Extrae la lista unica de especies de la columna especificada
#' 2. Para cada especie, filtra los registros correspondientes
#' 3. Extrae valores de las variables ambientales del raster en cada ubicacion
#' 4. Crea un dataframe combinando coordenadas, atributos y valores del raster
#' 5. Elimina columnas que son completamente NA
#' 6. Elimina filas que tienen cualquier valor NA en las columnas restantes
#' 7. Retorna una lista con dataframes limpios para cada especie
#'
#' @note
#' La funcion requiere que el objeto sf y el raster esten en el mismo sistema de coordenadas.
#' Las columnas completamente NA son removidas primero, luego las filas con NA son eliminadas.
#' Esto asegura que solo se retengan registros completos para el analisis.
#'
#' @importFrom terra extract
#' @importFrom sf st_coordinates st_drop_geometry
#' @importFrom dplyr select
#'
#' @examples
#' \dontrun{
#' # Ejemplo 1: Procesar datos de mamiferos de Oaxaca
#' # Example 1: Process mammal data from Oaxaca
#' library(sf)
#' library(dplyr)
#' library(terra)
#'
#' # Cargar datos de presencia
#' # Load presence data
#' Mamiferos_Oax_Sf <- st_read("datos/mamiferos.csv") %>%
#'   filter(estadomapa == "OAXACA")
#'
#' # Cargar y proyectar raster de variables ambientales
#' # Load and project environmental variables raster
#' raster_c1 <- rast("datos/temperatura.tif")
#' raster_c2 <- rast("datos/precipitacion.tif")
#' raster_tif <- c(raster_c1, raster_c2)
#' raster_proyectado <- terra::project(raster_tif, "EPSG:4326")
#'
#' # Procesar datos sin NA
#' # Process data without NA
#' resultados <- Presencias_Sf_Sin_Na(Mamiferos_Oax_Sf, raster_proyectado, "especie")
#'
#' # Ver resultados para una especie
#' # Check results for one species
#' if(length(resultados) > 0) {
#'   print(names(resultados)[1])
#'   print(head(resultados[[1]]))
#' }
#'
#' # Ejemplo 2: Usar diferente nombre de columna de especie
#' # Example 2: Use different species column name
#' resultados2 <- Presencias_Sf_Sin_Na(Mamiferos_Oax_Sf, raster_proyectado, "nombre_cientifico")
#' }
#'
#' @export
Presencias_Sf_Sin_Na <- function(Data_Sf, raster_referencia_reproyectado, Especie_Valida = "especievalida"){
  ListaEspecies <- unique(Data_Sf[[Especie_Valida]])
  resultados <- list()
  for (l in ListaEspecies) {
    Subset_Sf_Mam_Oax_especie <- Data_Sf[Data_Sf[[Especie_Valida]] == as.character(l), ]
    Variables_Stack_Raster <- terra::extract(raster_referencia_reproyectado, Subset_Sf_Mam_Oax_especie)
    coords <- sf::st_coordinates(Subset_Sf_Mam_Oax_especie)
    atributos <- sf::st_drop_geometry(Subset_Sf_Mam_Oax_especie)
    Df_Presencias <- data.frame(
      longitud = coords[, "X"],
      latitud = coords[, "Y"],
      atributos,
      Variables_Stack_Raster[, !names(Variables_Stack_Raster) == "ID", drop = FALSE],
      Presencia = 1
    )
    Presencias_Parcial_Na <- Df_Presencias %>%
      dplyr::select(!where(~ all(is.na(.))))
    Variables_Puro_NA <- names(Df_Presencias)[!names(Df_Presencias) %in% names(Presencias_Parcial_Na)]
    Df_Presencias_Sin_Na <- Presencias_Parcial_Na[complete.cases(Presencias_Parcial_Na), ]
    resultados[[as.character(l)]] <- Df_Presencias_Sin_Na
    print(paste("Se han eliminado las columnas:", paste(Variables_Puro_NA, collapse = ", ")))
  }
  return(invisible(resultados))
}
