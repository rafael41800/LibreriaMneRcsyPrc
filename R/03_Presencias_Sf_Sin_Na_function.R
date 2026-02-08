#' @title Presencias desde un sf sin entradas NA's / Presence Data from sf without NA Entries
#'
#' @description
#' Extrae y procesa datos de presencia de especies desde un objeto `sf`,
#' combinandolos con valores de variables ambientales de un raster y
#' eliminando todos los registros incompletos (con valores NA). La funcion
#' genera dataframes limpios y listos para analisis de distribucion de especies.
#' /
#' Extracts and processes species presence data from an `sf` object,
#' combining them with environmental variable values from a raster and
#' removing all incomplete records (with NA values). The function produces
#' clean dataframes ready for species distribution analysis.
#'
#' @param Data_Sf Objeto `sf` con datos de presencia de especies. Debe contener
#'   geometrias puntuales y la columna especificada en `Especie_Valida`.
#'   / `sf` object with species presence data. Must contain point geometries
#'   and the column specified in `Especie_Valida`.
#'
#' @param raster_referencia_reproyectado Objeto `SpatRaster` proyectado con
#'   las variables ambientales. Debe estar en el mismo sistema de coordenadas
#'   que `Data_Sf` para una extraccion correcta.
#'   / Projected `SpatRaster` object with environmental variables. Must be in
#'   the same coordinate system as `Data_Sf` for correct extraction.
#'
#' @param Especie_Valida Nombre de la columna que contiene los identificadores
#'   de especie. Puede ser nombre cientifico, codigo de especie u otro
#'   identificador unico.
#'   / Name of the column containing species identifiers. Can be scientific
#'   name, species code, or other unique identifier.
#'
#' @return
#' Una lista con los siguientes elementos:
#' * Un dataframe por cada especie unica encontrada en `Data_Sf`
#' * Cada dataframe contiene:
#'   - Coordenadas (`longitud`, `latitud`) extraidas de las geometrias
#'   - Todos los atributos originales del objeto `sf`
#'   - Valores extraidos de todas las capas del raster (`layer1`, `layer2`, ...)
#'   - Variable `Presencia` con valor constante `1` para todos los registros
#' * Todos los dataframes han sido limpiados eliminando:
#'   - Columnas que son completamente NA (sin valores utiles)
#'   - Filas con cualquier valor NA en las columnas restantes
#' * Los nombres de la lista corresponden a las especies unicas
#' /
#' A list with the following elements:
#' * One dataframe per unique species found in `Data_Sf`
#' * Each dataframe contains:
#'   - Coordinates (`longitud`, `latitud`) extracted from geometries
#'   - All original attributes from the `sf` object
#'   - Values extracted from all raster layers (`layer1`, `layer2`, ...)
#'   - `Presencia` variable with constant value `1` for all records
#' * All dataframes have been cleaned by removing:
#'   - Columns that are completely NA (no useful values)
#'   - Rows with any NA values in the remaining columns
#' * List names correspond to unique species
#'
#' @details
#' La funcion ejecuta el siguiente flujo de procesamiento para cada especie:
#' 1. **Filtrado por especie**: Separa los registros para cada especie unica
#' 2. **Extraccion de coordenadas**: Convierte geometrias sf a coordenadas X,Y
#' 3. **Extraccion de valores raster**: Obtiene valores de todas las bandas del
#'    raster en cada ubicacion puntual
#' 4. **Combinacion de datos**: Crea un dataframe integrando coordenadas,
#'    atributos originales y valores raster
#' 5. **Limpieza de columnas**: Elimina columnas que contienen solo valores NA
#' 6. **Limpieza de filas**: Elimina filas que contienen cualquier valor NA
#'    en las columnas restantes
#' 7. **Adicion de variable respuesta**: Agrega columna `Presencia = 1` para
#'    modelado de distribucion de especies
#'
#' La limpieza en dos etapas (columnas luego filas) optimiza la retencion de
#' datos utiles y garantiza datasets completos para analisis posteriores.
#' /
#' The function executes the following processing flow for each species:
#' 1. **Species filtering**: Separates records for each unique species
#' 2. **Coordinate extraction**: Converts sf geometries to X,Y coordinates
#' 3. **Raster value extraction**: Obtains values from all raster bands at
#'    each point location
#' 4. **Data combination**: Creates a dataframe integrating coordinates,
#'    original attributes, and raster values
#' 5. **Column cleaning**: Removes columns containing only NA values
#' 6. **Row cleaning**: Removes rows containing any NA values in the
#'    remaining columns
#' 7. **Response variable addition**: Adds `Presencia = 1` column for
#'    species distribution modeling
#'
#' The two-stage cleaning (columns then rows) optimizes retention of useful
#' data and ensures complete datasets for subsequent analysis.
#'
#' @note
#' * La funcion requiere que `Data_Sf` y `raster_referencia_reproyectado`
#'   esten en el mismo sistema de coordenadas de referencia (CRS). Si no lo
#'   estan, los valores extraidos del raster pueden ser incorrectos.
#' * Las geometrias deben ser puntos (`POINT`). Otros tipos de geometria
#'   pueden causar errores o resultados inesperados.
#' * El orden de limpieza es importante: primero se eliminan columnas
#'   completamente NA, luego filas con NA. Esto maximiza la retencion de
#'   datos validos.
#' * La columna `Especie_Valida` debe existir en `Data_Sf` y contener
#'   valores no nulos para todas las observaciones.
#' * La variable `Presencia` siempre tiene valor 1, ya que todos los
#'   registros son presencias confirmadas.
#' /
#' * The function requires that `Data_Sf` and `raster_referencia_reproyectado`
#'   are in the same coordinate reference system (CRS). If they are not,
#'   extracted raster values may be incorrect.
#' * Geometries must be points (`POINT`). Other geometry types may cause
#'   errors or unexpected results.
#' * The cleaning order is important: first completely NA columns are removed,
#'   then rows with NA. This maximizes retention of valid data.
#' * The `Especie_Valida` column must exist in `Data_Sf` and contain
#'   non-null values for all observations.
#' * The `Presencia` variable always has value 1, since all records are
#'   confirmed presences.
#'
#' @section Advertencias/Warnings:
#' * Si el raster y el objeto sf tienen diferentes CRS, la funcion no fallara
#'   pero los valores extraidos pueden no corresponder a las ubicaciones reales.
#' * Especies con muy pocos registros pueden resultar en dataframes vacios
#'   despues de la limpieza de NA.
#' * Columnas con nombres duplicados entre el objeto sf y el raster pueden
#'   causar conflictos de nombres en el dataframe resultante.
#' * El proceso de extraccion de valores raster puede ser lento con muchos
#'   puntos (miles o mas).
#' /
#' * If the raster and sf object have different CRS, the function will not fail
#'   but extracted values may not correspond to actual locations.
#' * Species with very few records may result in empty dataframes after
#'   NA cleaning.
#' * Columns with duplicate names between the sf object and raster may cause
#'   naming conflicts in the resulting dataframe.
#' * The raster value extraction process can be slow with many points
#'   (thousands or more).
#'
#' @examples
#' \dontrun{
#' # Ejemplo 1: Procesar datos de mamiferos de Oaxaca
#' # Example 1: Process mammal data from Oaxaca
#' library(sf)
#' library(dplyr)
#' library(terra)
#'
#' # Cargar datos de presencia (formato sf con puntos)
#' # Load presence data (sf format with points)
#' Mamiferos_Oax_Sf <- st_read("datos/mamiferos.csv") %>%
#'   filter(estadomapa == "OAXACA") %>%
#'   st_as_sf(coords = c("longitud", "latitud"), crs = 4326)
#'
#' # Cargar y proyectar raster de variables ambientales
#' # Load and project environmental variables raster
#' raster_c1 <- rast("datos/temperatura.tif")
#' raster_c2 <- rast("datos/precipitacion.tif")
#' raster_tif <- c(raster_c1, raster_c2)
#' raster_proyectado <- terra::project(raster_tif, "EPSG:4326")
#'
#' # Verificar que CRS coinciden
#' # Verify that CRS match
#' print(st_crs(Mamiferos_Oax_Sf))
#' print(crs(raster_proyectado))
#'
#' # Procesar datos sin NA
#' # Process data without NA
#' resultados <- Presencias_Sf_Sin_Na(
#'   Data_Sf = Mamiferos_Oax_Sf,
#'   raster_referencia_reproyectado = raster_proyectado,
#'   Especie_Valida = "especie"
#' )
#'
#' # Ver resultados
#' # Check results
#' print(paste("Numero de especies procesadas:", length(resultados)))
#'
#' if(length(resultados) > 0) {
#'   # Ver primera especie
#'   # Check first species
#'   print(paste("Primera especie:", names(resultados)[1]))
#'   print(paste("Numero de registros:", nrow(resultados[[1]])))
#'   print(head(resultados[[1]]))
#'
#'   # Resumen de todas las especies
#'   # Summary of all species
#'   resumen <- data.frame(
#'     Especie = names(resultados),
#'     Registros = sapply(resultados, nrow)
#'   )
#'   print(resumen)
#' }
#'
#' # Ejemplo 2: Usar diferente nombre de columna de especie
#' # Example 2: Use different species column name
#' resultados2 <- Presencias_Sf_Sin_Na(
#'   Data_Sf = Mamiferos_Oax_Sf,
#'   raster_referencia_reproyectado = raster_proyectado,
#'   Especie_Valida = "nombre_cientifico"
#' )
#'
#' # Ejemplo 3: Datos con diferentes tipos de variables
#' # Example 3: Data with different variable types
#' raster_multi <- rast(system.file("ex/elev.tif", package="terra"))
#' datos_ejemplo <- st_as_sf(
#'   data.frame(
#'     especie = rep(c("A", "B"), each = 10),
#'     longitud = runif(20, 5.5, 6.5),
#'     latitud = runif(20, 49, 50)
#'   ),
#'   coords = c("longitud", "latitud"), crs = 4326
#' )
#'
#' resultados3 <- Presencias_Sf_Sin_Na(
#'   Data_Sf = datos_ejemplo,
#'   raster_referencia_reproyectado = raster_multi,
#'   Especie_Valida = "especie"
#' )
#' }
#'
#' @seealso
#'   \itemize{
#'     \item \code{\link[terra]{extract}} para extraer valores raster en puntos /
#'           \code{\link[terra]{extract}} for extracting raster values at points
#'     \item \code{\link[sf]{st_coordinates}} para extraer coordenadas de objetos sf /
#'           \code{\link[sf]{st_coordinates}} for extracting coordinates from sf objects
#'     \item \code{\link[stats]{complete.cases}} para identificar filas completas sin NA /
#'           \code{\link[stats]{complete.cases}} for identifying complete rows without NA
#'     \item \code{\link[dplyr]{filter}} para filtrar dataframes basados en condiciones /
#'           \code{\link[dplyr]{filter}} for filtering dataframes based on conditions
#'   }
#'
#' @keywords
#'   presencias, limpieza de datos, extraccion de valores raster, SDM,
#'   presencia-ausencia, GIS,
#'   presences, data cleaning, raster value extraction, SDM,
#'   presence-absence, GIS
#'
#' @family
#'   funciones_SDM, funciones_limpieza_datos, funciones_espaciales /
#'   SDM_functions, data_cleaning_functions, spatial_functions
#'
#' @importFrom terra extract
#' @importFrom sf st_coordinates st_drop_geometry
#' @importFrom dplyr select %>% where
#' @importFrom stats complete.cases
#' @export
Presences_Sf_No_Na_List_Df <- function(Data_Sf, raster_referencia_reproyectado, Especie_Valida = "especievalida"){
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
