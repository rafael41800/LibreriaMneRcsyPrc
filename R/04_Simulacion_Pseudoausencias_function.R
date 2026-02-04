#' Pseudoausencias en base a las presencias, con radio de exclusión
#'
#' Genera pseudoausencias balanceadas (1:1) para modelos de nicho ecológico,
#' excluyendo celdas dentro de un radio alrededor de las presencias reales.
#' Las presencias de entrada deben estar limpias de valores NA.
#'
#' @param Resultados Lista de data frames con presencias por especie.
#'   Cada data frame debe contener las columnas: longitud, latitud, especievalida,
#'   y las variables ambientales extraídas del raster.
#' @param Raster_CRS_Referencia Objeto SpatRaster de terra con las variables ambientales
#'   ya proyectadas al CRS de trabajo. Define la resolución espacial para la selección
#'   de pseudoausencias.
#' @param Objeto_CRS_referencia Objeto espacial (sf, sfc o SpatVector) que define
#'   el sistema de coordenadas de referencia final para la salida.
#'
#' @returns Lista invisible de objetos sf, uno por especie, cada uno con:
#'   - Todas las variables ambientales
#'   - Columna 'Presence' (1 para presencias, 0 para ausencias)
#'   - Columna 'especievalida'
#'   - Geometría en el CRS especificado por Objeto_CRS_referencia
#'   Las pseudoausencias se generan en igual número que las presencias,
#'   excluyendo celdas dentro del radio de exclusión (0.1 grados o 10 km).
#'
#' @examples
#' \dontrun{
#' # Después de obtener presencias con Presencias_Sf_Sin_Na()
#' pseudoausencias <- Simular_Pseudoausencias_Sin_NA(
#'   Resultados = Subsets_Df_Presencias_Sin_Na,
#'   Raster_CRS_Referencia = Raster_Tif_reproyectado,
#'   Objeto_CRS_referencia = Oax_sf
#' )
#' }
#' @export
Simular_Pseudoausencias_Sin_NA_List_Df <- function(Resultados, Raster_CRS_Referencia, Objeto_CRS_referencia){
  resultados_finales <- list()
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
    resultados_finales[[i]] <- Datos_Presencias_Ausencias_Juntos_Df_To_Sf
  }
  return(invisible(resultados_finales))
}
