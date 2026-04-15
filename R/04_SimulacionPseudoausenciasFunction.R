#' @title Simulacion de Pseudoausencias sin NA para Lista de DataFrames /
#'   Simulation of Pseudoabsences without NA for List of DataFrames
#'
#' @description
#'   Esta funcion genera pseudoausencias (ausencias simuladas) para modelado de
#'   distribucion de especies (SDM), excluyendo areas cercanas a presencias conocidas
#'   mediante un radio de exclusion, extrayendo variables ambientales de un raster
#'   de referencia y filtrando registros con valores NA. La funcion procesa una lista
#'   de especies y retorna objetos espaciales combinados (presencias + pseudoausencias).
#'
#'   This function generates pseudoabsences (simulated absences) for species
#'   distribution modeling (SDM), excluding areas near known presences using an
#'   exclusion radius, extracting environmental variables from a reference raster,
#'   and filtering records with NA values. The function processes a list of species
#'   and returns combined spatial objects (presences + pseudoabsences).
#'
#' @param Resultados Lista nombrada donde cada elemento es un dataframe con los datos
#'   de presencia de una especie. Cada dataframe debe contener las columnas:
#'   `especievalida` (nombre de la especie), `longitud` (coordenada X/longitud),
#'   `latitud` (coordenada Y/latitud), y las variables ambientales necesarias
#'   (cuyos nombres deben coincidir con `NombreVariablesAmbientales`).
#'   El numero de pseudoausencias generadas sera igual al numero de presencias
#'   (o al maximo posible si hay menos celdas disponibles).
#'
#'   Named list where each element is a dataframe with presence data for a species.
#'   Each dataframe must contain the columns: `especievalida` (species name),
#'   `longitud` (X/longitude coordinate), `latitud` (Y/latitude coordinate),
#'   and the required environmental variables (whose names must match
#'   `NombreVariablesAmbientales`). The number of generated pseudoabsences
#'   will equal the number of presences (or the maximum possible if fewer cells are available).
#'
#' @param RasterCrsReferencia Objeto SpatRaster de \code{terra} que contiene las variables
#'   ambientales en la proyeccion de referencia para el analisis. Debe tener al menos
#'   una capa. Se usa para extraer valores ambientales en las ubicaciones de
#'   pseudoausencias y para identificar celdas disponibles.
#'
#'   \code{terra} SpatRaster object containing environmental variables in the reference
#'   projection for analysis. Must have at least one layer. Used to extract environmental
#'   values at pseudoabsence locations and to identify available cells.
#'
#' @param ObjetoCrsReferencia CRS objetivo (como objeto \code{sf::crs}) para transformar
#'   las coordenadas finales de los datos. Los resultados se transformaran a este
#'   sistema de referencia antes de ser retornados.
#'
#'   Target CRS (as \code{sf::crs} object) to transform the final data coordinates.
#'   Results will be transformed to this reference system before being returned.
#'
#' @param RadioGradosExclusion Radio de exclusion en grados para coordenadas geograficas
#'   (longitud/latitud). Se usa cuando el CRS de \code{RasterCrsReferencia} es geografico
#'   (longlat). Valor por defecto: \code{0.001} (~111 metros en el ecuador).
#'
#'   Exclusion radius in degrees for geographic coordinates (longitude/latitude).
#'   Used when the CRS of \code{RasterCrsReferencia} is geographic (longlat).
#'   Default: \code{0.001} (~111 meters at the equator).
#'
#' @param RadioMetrosExclusion Radio de exclusion en metros para coordenadas proyectadas.
#'   Se usa cuando el CRS de \code{RasterCrsReferencia} es proyectado (ej. UTM).
#'   Valor por defecto: \code{100} (100 metros).
#'
#'   Exclusion radius in meters for projected coordinates. Used when the CRS of
#'   \code{RasterCrsReferencia} is projected (e.g., UTM). Default: \code{100} (100 meters).
#'
#' @param NombreVariablesAmbientales Vector de caracteres con los nombres de las capas
#'   del raster que se utilizaran como variables ambientales. Debe coincidir exactamente
#'   con los nombres de las capas en \code{RasterCrsReferencia}. Estas variables seran
#'   extraidas en las ubicaciones de pseudoausencias y preservadas en el resultado.
#'   Este parametro es obligatorio.
#'
#'   Character vector with the names of the raster layers to be used as environmental
#'   variables. Must exactly match the layer names in \code{RasterCrsReferencia}.
#'   These variables will be extracted at pseudoabsence locations and preserved in the output.
#'   This parameter is required.
#'
#' @param verbose Logico. Si \code{TRUE} (valor por defecto), muestra mensajes de progreso
#'   en la consola indicando la especie en procesamiento, numero de pseudoausencias
#'   generadas, filtros aplicados y total final de registros.
#'
#'   Logical. If \code{TRUE} (default), displays progress messages in the console
#'   indicating the species being processed, number of generated pseudoabsences,
#'   applied filters, and final total of records.
#'
#' @param seed Entero opcional para fijar la semilla del generador de numeros aleatorios,
#'   garantizando reproducibilidad en la seleccion de pseudoausencias. Si es \code{NULL}
#'   (valor por defecto), no se fija semilla.
#'
#'   Optional integer to set the random number generator seed, ensuring reproducibility
#'   in pseudoabsence selection. If \code{NULL} (default), no seed is set.
#'
#' @return
#'   Una lista invisible con los mismos nombres que la lista de entrada \code{Resultados},
#'   donde cada elemento es un objeto \code{sf} (simple features) con los datos combinados
#'   de presencias reales y pseudoausencias. Solo se incluyen especies que generaron
#'   al menos una pseudoausencia valida. Cada objeto contiene las siguientes columnas:
#'   \itemize{
#'     \item \code{especievalida}: Nombre de la especie (caracter) / Species name (character)
#'     \item \code{...}: Variables ambientales especificadas en \code{NombreVariablesAmbientales}
#'           (numericas) / Specified environmental variables from \code{NombreVariablesAmbientales} (numeric)
#'     \item \code{ID}: Identificador unico por registro (entero) / Unique identifier per record (integer)
#'     \item \code{Presence}: Variable binaria: \code{1} para presencias reales,
#'           \code{0} para pseudoausencias / Binary variable: \code{1} for real presences,
#'           \code{0} for pseudoabsences
#'     \item \code{geometry}: Geometria \code{sf} con las coordenadas en el CRS especificado
#'           en \code{ObjetoCrsReferencia} / \code{sf} geometry with coordinates in the CRS
#'           specified in \code{ObjetoCrsReferencia}
#'   }
#'
#'   An invisible list with the same names as the input list \code{Resultados},
#'   where each element is an \code{sf} (simple features) object with combined real presence
#'   and pseudoabsence data. Only species that generated at least one valid pseudoabsence
#'   are included. Each object contains the following columns:
#'   \itemize{
#'     \item \code{especievalida}: Species name (character)
#'     \item \code{...}: Specified environmental variables from \code{NombreVariablesAmbientales} (numeric)
#'     \item \code{ID}: Unique identifier per record (integer)
#'     \item \code{Presence}: Binary variable: \code{1} for real presences,
#'           \code{0} for pseudoabsences
#'     \item \code{geometry}: \code{sf} geometry with coordinates in the CRS
#'           specified in \code{ObjetoCrsReferencia}
#'   }
#'
#' @details
#'   \strong{Proceso realizado para cada especie / Process performed for each species:}
#'   \enumerate{
#'     \item Conversion de coordenadas de presencias a objetos espaciales \code{sf} con CRS WGS84 (EPSG:4326) /
#'           Conversion of presence coordinates to \code{sf} spatial objects with CRS WGS84 (EPSG:4326)
#'     \item Transformacion al CRS del raster de referencia (\code{RasterCrsReferencia}) /
#'           Transformation to the reference raster CRS (\code{RasterCrsReferencia})
#'     \item Creacion de buffer de exclusion alrededor de presencias:
#'           \itemize{
#'             \item Si el CRS es geografico (longlat): usa \code{RadioGradosExclusion}
#'             \item Si el CRS es proyectado: usa \code{RadioMetrosExclusion}
#'           } /
#'           Creation of exclusion buffer around presences:
#'           \itemize{
#'             \item If CRS is geographic (longlat): uses \code{RadioGradosExclusion}
#'             \item If CRS is projected: uses \code{RadioMetrosExclusion}
#'           }
#'     \item Identificacion de celdas raster disponibles:
#'           \itemize{
#'             \item Celdas sin NA en la primera capa del raster
#'             \item Excluyendo celdas que contienen presencias
#'           } /
#'           Identification of available raster cells:
#'           \itemize{
#'             \item Non-NA cells in the first raster layer
#'             \item Excluding cells containing presences
#'           }
#'     \item Generacion de pseudoausencias: seleccion aleatoria de celdas (sin reemplazo)
#'           en numero igual al de presencias (o al maximo disponible) /
#'           Pseudoabsence generation: random selection of cells (without replacement)
#'           equal in number to presences (or maximum available)
#'     \item Extraccion de variables ambientales en ubicaciones de pseudoausencias /
#'           Extraction of environmental variables at pseudoabsence locations
#'     \item Identificacion de pseudoausencias dentro del buffer de exclusion de presencias /
#'           Identification of pseudoabsences within the presence exclusion buffer
#'     \item Filtrado: eliminacion de pseudoausencias:
#'           \itemize{
#'             \item Dentro del radio de exclusion (marcadas como \code{Presence = NA})
#'             \item Con valores NA en las variables ambientales especificadas
#'           } /
#'           Filtering: removal of pseudoabsences:
#'           \itemize{
#'             \item Within exclusion radius (marked as \code{Presence = NA})
#'             \item With NA values in specified environmental variables
#'           }
#'     \item Combinacion de presencias (con \code{Presence = 1}) y pseudoausencias validas
#'           (con \code{Presence = 0}) en un solo dataframe /
#'           Combination of presences (with \code{Presence = 1}) and valid pseudoabsences
#'           (with \code{Presence = 0}) into a single dataframe
#'     \item Transformacion al CRS objetivo (\code{ObjetoCrsReferencia}) /
#'           Transformation to target CRS (\code{ObjetoCrsReferencia})
#'   }
#'
#' @note
#'   \itemize{
#'     \item La funcion se detiene con error si no hay celdas disponibles para pseudoausencias
#'           (todas las celdas contienen presencias o son NA) para alguna especie. /
#'           The function stops with an error if no cells are available for pseudoabsences
#'           (all cells contain presences or are NA) for any species.
#'     \item El parametro \code{NombreVariablesAmbientales} es obligatorio y debe coincidir
#'           con los nombres de las capas del raster. /
#'           The parameter \code{NombreVariablesAmbientales} is required and must match
#'           the raster layer names.
#'     \item El resultado se retorna con \code{return(invisible(...))}, permitiendo asignacion
#'           silenciosa o impresion directa. /
#'           The result is returned with \code{return(invisible(...))}, allowing silent
#'           assignment or direct printing.
#'     \item Si para una especie no se generan pseudoausencias validas despues del filtrado,
#'           se emite una advertencia y se omite dicha especie del resultado. /
#'           If no valid pseudoabsences are generated for a species after filtering,
#'           a warning is issued and that species is omitted from the result.
#'     \item Las presencias originales deben incluir las variables ambientales en el dataframe;
#'           estas se conservan en el resultado final. /
#'           Original presences must include environmental variables in the dataframe;
#'           these are preserved in the final output.
#'     \item La funcion utiliza la primera capa del raster (\code{RasterCrsReferencia[[1]]})
#'           para identificar celdas con presencia y celdas sin NA. Asegurese de que esta
#'           capa sea representativa de la extension y mascara de todas las capas. /
#'           The function uses the first raster layer (\code{RasterCrsReferencia[[1]]})
#'           to identify cells with presence and non-NA cells. Ensure this layer is
#'           representative of the extent and mask of all layers.
#'   }
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(sf)
#' library(dplyr)
#'
#' # Datos de ejemplo / Example data
#' datos_especie <- data.frame(
#'   especievalida = "Especie_ejemplo",
#'   longitud = c(-70.5, -70.6, -70.7),
#'   latitud = c(-33.4, -33.5, -33.6),
#'   var1 = c(100, 120, 110),
#'   var2 = c(800, 850, 820)
#' )
#'
#' lista_especies <- list(especie1 = datos_especie)
#'
#' # Crear raster de ejemplo / Create example raster
#' r <- rast(nrows = 100, ncols = 100,
#'           xmin = -71, xmax = -70,
#'           ymin = -34, ymax = -33)
#' r[[1]] <- setValues(r, runif(10000, 80, 150))
#' r[[2]] <- setValues(r, runif(10000, 700, 900))
#' names(r) <- c("var1", "var2")
#' crs(r) <- "EPSG:4326"
#'
#' # CRS objetivo (ej. UTM zona 19 sur para Chile) / Target CRS (e.g., UTM zone 19 south for Chile)
#' crs_objetivo <- st_crs("EPSG:32719")
#'
#' # Ejecutar funcion / Run function
#' resultado <- SimulatePseudoabsencesNoNAListDf(
#'   Resultados = lista_especies,
#'   RasterCrsReferencia = r,
#'   ObjetoCrsReferencia = crs_objetivo,
#'   NombreVariablesAmbientales = c("var1", "var2"),
#'   verbose = TRUE,
#'   seed = 123
#' )
#'
#' # Acceder a resultados / Access results
#' especie_procesada <- resultado$especie1
#' head(especie_procesada)
#' table(especie_procesada$Presence)
#' }
#'
#' @seealso
#'   \code{\link[terra]{rast}} para manejo de rasters,
#'   \code{\link[sf]{st_as_sf}} para creacion de objetos espaciales,
#'   \code{\link[sf]{st_buffer}} para creacion de buffers de exclusion,
#'   \code{\link[dplyr]{filter}} para operaciones de filtrado,
#'   \code{\link{sample}} para muestreo aleatorio de celdas
#'
#' @keywords SDM pseudoausencias modelado distribucion especies pseudoabsences
#'
#' @family funciones_SDM
#'
#' @importFrom sf st_as_sf st_transform st_is_longlat st_buffer st_union
#'   st_intersects st_coordinates st_crs
#' @importFrom terra cellFromXY values xyFromCell crs vect extract
#' @importFrom dplyr %>% mutate select filter if_all everything bind_rows row_number all_of
#' @importFrom rlang .data
#' @export
SimulatePseudoabsencesNoNAListDf <- function(Resultados, RasterCrsReferencia, ObjetoCrsReferencia, RadioGradosExclusion = 0.001, RadioMetrosExclusion = 100, NombreVariablesAmbientales, verbose = TRUE, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  if(verbose) cat("Simulando pseudo-ausencias...\n")
  ResultadosFinalesList <- list()
  TotalEspecies <- length(names(Resultados))
  for (idx in seq_along(names(Resultados))) {
    i <- names(Resultados)[idx]
    if(verbose) cat(sprintf("  %s [%d/%d]\n", i, idx, TotalEspecies))
    EspecieEnProcesamientoDf <- Resultados[[i]]
    PresenciasSinNaSf <- EspecieEnProcesamientoDf %>%
      sf::st_as_sf(coords = c("longitud", "latitud"), crs = 4326) %>%
      sf::st_transform(sf::st_crs(RasterCrsReferencia))
    if (sf::st_is_longlat(PresenciasSinNaSf)) {
      PoligonoRadio <- PresenciasSinNaSf %>%
        sf::st_buffer(dist = RadioGradosExclusion) %>%
        sf::st_union() %>%
        sf::st_as_sf()
    } else {
      PoligonoRadio <- PresenciasSinNaSf %>%
        sf::st_buffer(dist = RadioMetrosExclusion) %>%
        sf::st_union() %>%
        sf::st_as_sf()
    }
    CeldaRasterConPresencia <- terra::cellFromXY(RasterCrsReferencia[[1]], sf::st_coordinates(PresenciasSinNaSf))
    CeldaRasterSinNa <- which(!is.na(terra::values(RasterCrsReferencia[[1]])))
    CeldasDisponiblesParaPseudoausencias <- setdiff(CeldaRasterSinNa, CeldaRasterConPresencia)
    if (length(CeldasDisponiblesParaPseudoausencias) == 0) {
      stop("ERROR: No hay celdas disponibles (todas tienen presencias o son NA)")
    }
    NPresencias <- nrow(PresenciasSinNaSf)
    NAusenciasPosibles <- min(NPresencias, length(CeldasDisponiblesParaPseudoausencias))
    if(verbose) cat(sprintf("    Generando %d pseudo-ausencias\n", NAusenciasPosibles))
    PseudoausenciasAleatoriasEnCdpp <- sample(CeldasDisponiblesParaPseudoausencias, NAusenciasPosibles, replace = FALSE)
    CoordenadasRaster <- terra::xyFromCell(RasterCrsReferencia[[1]], PseudoausenciasAleatoriasEnCdpp)
    CoordenadasRasterSf <- sf::st_as_sf(
      data.frame(x = CoordenadasRaster[,1], y = CoordenadasRaster[,2]),
      coords = c("x", "y"),
      crs = terra::crs(RasterCrsReferencia)
    )
    SpavectorVariablesAmbientalesPseudoausencias <- terra::vect(CoordenadasRasterSf)
    VariablesAmbientalesPseudoausencias <- terra::extract( RasterCrsReferencia, SpavectorVariablesAmbientalesPseudoausencias, method = 'simple', ID = FALSE, bind = FALSE)
    ExclusionRadioPresencias <- sf::st_intersects(CoordenadasRasterSf, PoligonoRadio, sparse = FALSE)
    VariablesAmbientalesPseudoausencias$ExclusionRadioPresencias <- as.numeric(ExclusionRadioPresencias[,1])
    VariablesAmbientalesPseudoausenciasDf <- as.data.frame(VariablesAmbientalesPseudoausencias)
    VariablesAmbientalesPseudoausenciasDf$Presence <- ifelse(VariablesAmbientalesPseudoausenciasDf$ExclusionRadioPresencias == 1, NA, 0)
    NombreEspecieValida <- unique(EspecieEnProcesamientoDf$especievalida)
    CoordsMat <- sf::st_coordinates(CoordenadasRasterSf)
    VariablesAmbientalesPseudoausenciasDf <- VariablesAmbientalesPseudoausenciasDf %>%
      dplyr::mutate( especievalida = NombreEspecieValida, x = CoordsMat[, "X"], y = CoordsMat[, "Y"]) %>%
      dplyr::select(.data$especievalida, dplyr::everything())
    VariablesAmbientalesPseudoausenciasDfFueraRadio <- VariablesAmbientalesPseudoausenciasDf %>%
      dplyr::filter(.data$Presence == 0)
    RasterNombreVariables <- NombreVariablesAmbientales
    VariablesAmbientalesPseudoausenciasSinNa <- VariablesAmbientalesPseudoausenciasDfFueraRadio %>%
      dplyr::filter(
        dplyr::if_all(
          dplyr::all_of(RasterNombreVariables), ~ !is.na(.)
        )
      )
    if(verbose) cat(sprintf("    Pseudo-ausencias validas despues de filtrado: %d\n", nrow(VariablesAmbientalesPseudoausenciasSinNa)))
    if (nrow(VariablesAmbientalesPseudoausenciasSinNa) == 0) {
      warning(sprintf("Especie %s: No se generaron pseudo-ausencias validas. Usando solo presencias.", i))
      next
    }
    if (inherits(EspecieEnProcesamientoDf, "SpatialPointsDataFrame")) { EspecieEnProcesamientoDf <- as.data.frame(EspecieEnProcesamientoDf)}
    DatosPseudoausenciasParaFusion <- EspecieEnProcesamientoDf %>%
      dplyr::select(.data$especievalida, dplyr::all_of(RasterNombreVariables), .data$longitud, .data$latitud) %>%
      dplyr::mutate( ID = dplyr::row_number(), Presence = 1, x = .data$longitud, y = .data$latitud) %>%
      dplyr::select(-.data$longitud, -.data$latitud)
    DatosAusenciasParaFusion <- VariablesAmbientalesPseudoausenciasSinNa %>%
      dplyr::select(.data$especievalida, dplyr::all_of(RasterNombreVariables), .data$x, .data$y) %>%
      dplyr::mutate( ID = dplyr::row_number() + nrow(DatosPseudoausenciasParaFusion), Presence = 0)
    DatosPresenciasAusenciasJuntosDf <- dplyr::bind_rows(DatosPseudoausenciasParaFusion, DatosAusenciasParaFusion)
    DatosPresenciasAusenciasJuntosDfToSf <- DatosPresenciasAusenciasJuntosDf %>%
      sf::st_as_sf(coords = c("x", "y"), crs = terra::crs(RasterCrsReferencia)) %>%
      sf::st_transform(sf::st_crs(ObjetoCrsReferencia))
    ResultadosFinalesList[[i]] <- DatosPresenciasAusenciasJuntosDfToSf
    if(verbose) cat(sprintf("    Total final para %s: %d registros (presencias + ausencias)\n", i, nrow(DatosPresenciasAusenciasJuntosDf)))
  }
  if(verbose) cat("Procesamiento completado\n")
  return(invisible(ResultadosFinalesList))
}
