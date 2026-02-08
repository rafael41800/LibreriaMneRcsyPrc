#' @title Simulacion de Pseudoausencias para Modelado de Distribucion de Especies /
#'   Simulation of Pseudoabsences for Species Distribution Modeling
#'
#' @description
#'   Esta funcion genera pseudoausencias (ausencias simuladas) para modelado de
#'   distribucion de especies (SDM), excluyendo areas cercanas a presencias conocidas
#'   y asegurando la calidad de los datos ambientales.
#'
#'   This function generates pseudoabsences (simulated absences) for species
#'   distribution modeling (SDM), excluding areas near known presences and ensuring
#'   environmental data quality.
#'
#' @param Resultados Lista nombrada donde cada elemento es un dataframe con los datos
#'   de presencia de una especie. Cada dataframe debe contener las columnas:
#'   `especievalida`, `longitud`, `latitud`, y las variables ambientales necesarias.
#'
#'   Named list where each element is a dataframe with presence data for a species.
#'   Each dataframe must contain the columns: `especievalida`, `longitud`, `latitud`,
#'   and the necessary environmental variables.
#'
#' @param Raster_CRS_Referencia Objeto SpatRaster de terra que contiene las variables
#'   ambientales en la proyeccion de referencia para el analisis.
#'
#'   terra SpatRaster object containing environmental variables in the reference
#'   projection for analysis.
#'
#' @param Objeto_CRS_referencia CRS objetivo (como objeto sf::crs) para transformar
#'   las coordenadas finales de los datos.
#'
#'   Target CRS (as sf::crs object) to transform the final data coordinates.
#'
#' @return
#'   Una lista invisible con el mismo nombre que la lista de entrada, donde cada
#'   elemento es un objeto sf con los datos combinados de presencias y pseudoausencias.
#'   Cada objeto contiene las columnas:
#'   \itemize{
#'     \item \code{especievalida}: Nombre de la especie / Species name
#'     \item \code{CarbonStock_05cmdepth_Mg_x_ha}: Variable ambiental 1 / Environmental variable 1
#'     \item \code{EvapoTransp}: Variable ambiental 2 / Environmental variable 2
#'     \item \code{ID}: Identificador unico por registro / Unique identifier per record
#'     \item \code{Presence}: 1 para presencias reales, 0 para pseudoausencias /
#'           1 for real presences, 0 for pseudoabsences
#'     \item \code{geometry}: Geometria sf con las coordenadas / sf geometry with coordinates
#'   }
#'
#'   An invisible list with the same names as the input list, where each element
#'   is an sf object with combined presence and pseudoabsence data. Each object
#'   contains the columns described above.
#'
#' @details
#'   \strong{Proceso realizado / Process performed:}
#'   \enumerate{
#'     \item Conversion de coordenadas de presencias a objetos espaciales /
#'           Conversion of presence coordinates to spatial objects
#'     \item Creacion de un buffer de exclusion alrededor de las presencias (0.1° para
#'           coordenadas geograficas, 10000 metros para coordenadas proyectadas) /
#'           Creation of an exclusion buffer around presences (0.1° for geographic
#'           coordinates, 10000 meters for projected coordinates)
#'     \item Identificacion de celdas raster disponibles para pseudoausencias /
#'           Identification of available raster cells for pseudoabsences
#'     \item Generacion de pseudoausencias en numero igual al de presencias (o menor
#'           si no hay suficientes celdas disponibles) /
#'           Generation of pseudoabsences equal in number to presences (or fewer
#'           if not enough cells are available)
#'     \item Extraccion de variables ambientales en las ubicaciones de pseudoausencias /
#'           Extraction of environmental variables at pseudoabsence locations
#'     \item Exclusion de pseudoausencias dentro del buffer de presencias /
#'           Exclusion of pseudoabsences within the presence buffer
#'     \item Filtrado de registros sin valores NA en variables ambientales /
#'           Filtering of records without NA values in environmental variables
#'     \item Combinacion de presencias y pseudoausencias en un solo dataset /
#'           Combination of presences and pseudoabsences into a single dataset
#'     \item Transformacion al sistema de coordenadas objetivo /
#'           Transformation to the target coordinate system
#'   }
#'
#' @note
#'   \strong{Notas importantes / Important notes:}
#'   \itemize{
#'     \item La funcion asume que los rasters contienen las variables especificas
#'           `CarbonStock_05cmdepth_Mg_x_ha` y `EvapoTransp`. Modifica el vector
#'           `Raster_Nombre_Variables` dentro de la funcion si necesitas otras variables. /
#'           The function assumes that rasters contain the specific variables
#'           `CarbonStock_05cmdepth_Mg_x_ha` and `EvapoTransp`. Modify the
#'           `Raster_Nombre_Variables` vector inside the function if you need other variables.
#'     \item El buffer de exclusion varia segun el sistema de coordenadas para mantener
#'           distancias biologicamente significativas. /
#'           The exclusion buffer varies according to the coordinate system to maintain
#'           biologically meaningful distances.
#'     \item Se detiene la ejecucion con error si no hay celdas disponibles para
#'           pseudoausencias. /
#'           Execution stops with an error if no cells are available for pseudoabsences.
#'   }
#'
#' @examples
#' \dontrun{
#' # Ejemplo en espanol / Example in Spanish
#' # Cargar librerias necesarias / Load necessary libraries
#' library(terra)
#' library(sf)
#' library(dplyr)
#'
#' # Crear datos de ejemplo / Create example data
#' datos_especie <- data.frame(
#'   especievalida = "Especie_ejemplo",
#'   longitud = c(-70.5, -70.6, -70.7),
#'   latitud = c(-33.4, -33.5, -33.6),
#'   CarbonStock_05cmdepth_Mg_x_ha = c(100, 120, 110),
#'   EvapoTransp = c(800, 850, 820)
#' )
#'
#' lista_especies <- list(especie1 = datos_especie)
#'
#' # Crear raster de ejemplo / Create example raster
#' raster_ref <- rast(nrows = 100, ncols = 100, xmin = -71, xmax = -70,
#'                    ymin = -34, ymax = -33)
#' values(raster_ref[[1]]) <- runif(10000, 80, 150)
#' raster_ref[[2]] <- setValues(raster_ref[[1]], runif(10000, 700, 900))
#' names(raster_ref) <- c("CarbonStock_05cmdepth_Mg_x_ha", "EvapoTransp")
#' crs(raster_ref) <- "EPSG:4326"
#'
#' # Definir CRS objetivo / Define target CRS
#' crs_objetivo <- st_crs("EPSG:32719")
#'
#' # Ejecutar funcion / Execute function
#' resultado <- Simulate_Pseudoabsences_No_NA_List_Df(
#'   Resultados = lista_especies,
#'   Raster_CRS_Referencia = raster_ref,
#'   Objeto_CRS_referencia = crs_objetivo
#' )
#'
#' # Inspeccionar resultados / Inspect results
#' str(resultado)
#' length(resultado)
#' }
#'
#'
#' @seealso
#'   \itemize{
#'     \item \code{\link[terra]{rast}} para crear objetos raster /
#'           \code{\link[terra]{rast}} to create raster objects
#'     \item \code{\link[sf]{st_as_sf}} para conversion a objetos espaciales /
#'           \code{\link[sf]{st_as_sf}} for conversion to spatial objects
#'     \item \code{\link[dplyr]{filter}} para filtrado de datos /
#'           \code{\link[dplyr]{filter}} for data filtering
#'   }
#'
#' @keywords
#'   SDM, pseudoausencias, modelado de distribucion de especies, ecologia,
#'   SIG, pseudoabsences, species distribution modeling, ecology, GIS
#'
#' @family
#'   funciones_SDM, funciones_ecologia, funciones_espaciales /
#'   SDM_functions, ecology_functions, spatial_functions
#'
#' @importFrom sf st_as_sf st_transform st_is_longlat st_buffer st_union
#'   st_intersects st_coordinates st_crs
#' @importFrom terra cellFromXY values xyFromCell crs vect extract
#' @importFrom dplyr %>% mutate select filter if_all everything bind_rows row_number
#' @export
Simulate_Pseudoabsences_No_NA_List_Df <- function(Resultados, Raster_CRS_Referencia, Objeto_CRS_referencia){
  resultados_finales_List <- list()
  for (i in names(Resultados)) {
    Especie_En_Procesamiento_Df <- Resultados[[i]]
    Presencias_Sin_Na_sf <- Especie_En_Procesamiento_Df %>%
      sf::st_as_sf(coords = c("longitud", "latitud"),
                   crs = 4326) %>%
      sf::st_transform(sf::st_crs(Raster_CRS_Referencia))
    if (sf::st_is_longlat(Presencias_Sin_Na_sf)) {
      Radio_Presencia <- 0.1
      Poligono_Radio <- Presencias_Sin_Na_sf %>%
        sf::st_buffer(dist = Radio_Presencia) %>%
        sf::st_union() %>%
        sf::st_as_sf()
    } else {
      Radio_Presencia <- 10000
      Poligono_Radio <- Presencias_Sin_Na_sf %>%
        sf::st_buffer(dist = Radio_Presencia) %>%
        sf::st_union() %>%
        sf::st_as_sf()
    }
    Celda_Raster_Con_Presencia <- terra::cellFromXY(Raster_CRS_Referencia[[1]], sf::st_coordinates(Presencias_Sin_Na_sf))
    Celda_Raster_Sin_Na <- which(!is.na(terra::values(Raster_CRS_Referencia[[1]])))
    Celdas_Disponibles_Para_Pseudoausencias <- setdiff(Celda_Raster_Sin_Na, Celda_Raster_Con_Presencia)
    if (length(Celdas_Disponibles_Para_Pseudoausencias) == 0) {
      stop("ERROR: No hay celdas disponibles (todas tienen presencias o son NA)")
    }
    n_Presencias <- nrow(Presencias_Sin_Na_sf)
    n_ausencias_posibles <- min(n_Presencias, length(Celdas_Disponibles_Para_Pseudoausencias))
    Pseudoausencias_Aleatorias_En_CDPP <- sample(Celdas_Disponibles_Para_Pseudoausencias, n_ausencias_posibles, replace = FALSE)
    Coordenadas_Raster <- terra::xyFromCell(Raster_CRS_Referencia[[1]], Pseudoausencias_Aleatorias_En_CDPP)
    Coordenadas_Raster_Sf <- sf::st_as_sf(
      data.frame(x = Coordenadas_Raster[,1], y = Coordenadas_Raster[,2]),
      coords = c("x", "y"),
      crs = terra::crs(Raster_CRS_Referencia)
    )
    Spavector_Variables_Ambientales_Pseudoausencias <- terra::vect(Coordenadas_Raster_Sf)
    Variables_Ambientales_Pseudoausencias <- terra::extract(
      Raster_CRS_Referencia,
      Spavector_Variables_Ambientales_Pseudoausencias,
      method = 'simple',
      ID = FALSE,
      bind = FALSE
    )
    Exclusion_Radio_Presencias <- sf::st_intersects(Coordenadas_Raster_Sf, Poligono_Radio, sparse = FALSE)
    Variables_Ambientales_Pseudoausencias$Exclusion_Radio_Presencias <- as.numeric(Exclusion_Radio_Presencias[,1])
    Variables_Ambientales_Pseudoausencias_Df <- as.data.frame(Variables_Ambientales_Pseudoausencias)
    Variables_Ambientales_Pseudoausencias_Df$Presence <- ifelse(Variables_Ambientales_Pseudoausencias_Df$Exclusion_Radio_Presencias == 1, NA, 0)
    Nombre_EspecieValida <- unique(Especie_En_Procesamiento_Df$especievalida)
    Variables_Ambientales_Pseudoausencias_Df <- Variables_Ambientales_Pseudoausencias_Df %>%
      dplyr::mutate(
        especievalida = Nombre_EspecieValida,
        x = sf::st_coordinates(Coordenadas_Raster_Sf)[, "X"],
        y = sf::st_coordinates(Coordenadas_Raster_Sf)[, "Y"]
      ) %>%
      dplyr::select(especievalida, everything())
    Variables_Ambientales_Pseudoausencias_Df_Fuera_Radio <- Variables_Ambientales_Pseudoausencias_Df %>%
      dplyr::filter(Presence == 0)
    Raster_Nombre_Variables <- c("CarbonStock_05cmdepth_Mg_x_ha", "EvapoTransp")
    Variables_Ambientales_Pseudoausencias_Sin_Na <- Variables_Ambientales_Pseudoausencias_Df_Fuera_Radio %>%
      dplyr::filter(
        dplyr::if_all(
          dplyr::all_of(Raster_Nombre_Variables),
          ~ !is.na(.)
        )
      )
    if (inherits(Especie_En_Procesamiento_Df, "SpatialPointsDataFrame")) {
      Especie_En_Procesamiento_Df <- as.data.frame(Especie_En_Procesamiento_Df)
    }
    Datos_Pseudoausencias_Para_Fusion <- Especie_En_Procesamiento_Df %>%
      dplyr::select(especievalida, all_of(Raster_Nombre_Variables), longitud, latitud) %>%
      dplyr::mutate(
        ID = dplyr::row_number(),
        Presence = 1,
        x = longitud,
        y = latitud
      ) %>%
      dplyr::select(-longitud, -latitud)
    Datos_Ausencias_Para_Fusion <- Variables_Ambientales_Pseudoausencias_Sin_Na %>%
      dplyr::select(especievalida, all_of(Raster_Nombre_Variables), x, y) %>%
      dplyr::mutate(
        ID = dplyr::row_number() + nrow(Datos_Pseudoausencias_Para_Fusion),
        Presence = 0
      )
    Datos_Presencias_Ausencias_Juntos_Df <- dplyr::bind_rows(Datos_Pseudoausencias_Para_Fusion, Datos_Ausencias_Para_Fusion)
    Datos_Presencias_Ausencias_Juntos_Df_To_Sf <- Datos_Presencias_Ausencias_Juntos_Df %>%
      sf::st_as_sf(coords = c("x", "y"), crs = terra::crs(Raster_CRS_Referencia)) %>%
      sf::st_transform(sf::st_crs(Objeto_CRS_referencia))
    resultados_finales_List[[i]] <- Datos_Presencias_Ausencias_Juntos_Df_To_Sf
  }
  return(invisible(resultados_finales_List))
}
