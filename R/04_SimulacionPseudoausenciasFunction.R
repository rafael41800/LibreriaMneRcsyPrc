#' Generación de Pseudo-ausencias para Múltiples Especies
#'
#' @description
#' Genera puntos de pseudo-ausencia para una lista de especies, evitando áreas
#' cercanas a las presencias conocidas mediante un radio de exclusión.
#'
#' @param Resultados Lista de data.frames con los datos de presencia por especie.
#' @param RasterCrsReferencia Objeto `SpatRaster` que define la extensión y resolución.
#' @param ObjetoCrsReferencia Sistema de referencia de coordenadas (CRS) de destino.
#' @param RadioGradosExclusion Radio de exclusión si el CRS es geográfico. Por defecto `0.001`.
#' @param RadioMetrosExclusion Radio de exclusión si el CRS es proyectado. Por defecto `100`.
#' @param NombreVariablesAmbientales Vector con los nombres de los predictores.
#' @param verbose Logico. Muestra el progreso en consola.
#' @param seed Semilla para reproducibilidad.
#'
#' @return Una lista de objetos `sf` con datos de presencia (1) y pseudo-ausencia (0).
#'
#' @importFrom sf st_as_sf st_transform st_crs st_buffer st_union st_is_longlat st_coordinates st_intersects
#' @importFrom terra cellFromXY values xyFromCell crs extract vect
#' @importFrom dplyr select mutate rename filter if_all all_of bind_rows row_number %>%
#' @importFrom stats setNames
#' @export
SimulatePseudoabsencesNoNAListDf <- function(Resultados,
                                             RasterCrsReferencia,
                                             ObjetoCrsReferencia,
                                             RadioGradosExclusion = 0.001,
                                             RadioMetrosExclusion = 100,
                                             NombreVariablesAmbientales,
                                             verbose = TRUE,
                                             seed = NULL) {

  # Prevenir notas de R CMD check
  especievalida <- longitud <- latitud <- X <- Y <- Presence <- FueraRadio <- NULL

  if(!is.null(seed)) set.seed(seed)
  if(verbose) message("\n--- Iniciando simulacion de pseudo-ausencias ---")

  ResultadosFinalesList <- list()
  NombresEspecies <- names(Resultados)
  TotalEspecies <- length(NombresEspecies)

  for (idx in seq_along(NombresEspecies)) {
    i <- NombresEspecies[idx]
    if(verbose) message(sprintf("-> Procesando [%d/%d]: %s", idx, TotalEspecies, i))

    EspecieEnProcesamientoDf <- Resultados[[i]]

    # 1. Convertir presencias a sf y proyectar al CRS del Raster
    PresenciasSinNaSf <- EspecieEnProcesamientoDf %>%
      sf::st_as_sf(coords = c("longitud", "latitud"), crs = 4326) %>%
      sf::st_transform(sf::st_crs(RasterCrsReferencia))

    # 2. Crear buffer de exclusión
    dist_buffer <- if(sf::st_is_longlat(PresenciasSinNaSf)) RadioGradosExclusion else RadioMetrosExclusion
    PoligonoExclusion <- PresenciasSinNaSf %>%
      sf::st_buffer(dist = dist_buffer) %>%
      sf::st_union()

    # 3. Identificar celdas disponibles (sin NA y sin presencias)
    CeldaConPresencia <- terra::cellFromXY(RasterCrsReferencia[[1]], sf::st_coordinates(PresenciasSinNaSf))
    CeldaSinNa <- which(!is.na(terra::values(RasterCrsReferencia[[1]])))
    CeldasDisponibles <- setdiff(CeldaSinNa, CeldaConPresencia)

    if (length(CeldasDisponibles) == 0) {
      warning(paste("No hay celdas disponibles para la especie:", i))
      next
    }

    # 4. Sortear candidatos a pseudo-ausencias
    NAusenciasPedidas <- nrow(PresenciasSinNaSf)
    NAusenciasPosibles <- min(NAusenciasPedidas, length(CeldasDisponibles))

    IndicesSorteados <- sample(CeldasDisponibles, NAusenciasPosibles, replace = FALSE)
    CoordsPseudo <- terra::xyFromCell(RasterCrsReferencia[[1]], IndicesSorteados)

    # 5. Filtrar por Radio de Exclusión
    PseudoSf <- sf::st_as_sf(
      as.data.frame(CoordsPseudo),
      coords = c("x", "y"),
      crs = sf::st_crs(RasterCrsReferencia)
    )

    Interseccion <- sf::st_intersects(PseudoSf, PoligonoExclusion, sparse = FALSE)
    PseudoSf$FueraRadio <- !Interseccion[,1]

    # 6. Extraer valores ambientales y limpiar
    ValoresEnv <- terra::extract(RasterCrsReferencia, terra::vect(PseudoSf), ID = FALSE)

    PseudoFinalDf <- cbind(sf::st_coordinates(PseudoSf), ValoresEnv, FueraRadio = PseudoSf$FueraRadio) %>%
      as.data.frame() %>%
      dplyr::filter(FueraRadio == TRUE) %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(NombreVariablesAmbientales), ~ !is.na(.))) %>%
      dplyr::mutate(especievalida = i, Presence = 0) %>%
      dplyr::rename(X = X, Y = Y) # Garantizar nombres X e Y

    if (nrow(PseudoFinalDf) == 0) {
      warning(sprintf("Especie %s: 0 pseudo-ausencias validas tras filtrado. Se omite.", i))
      next
    }

    # 7. Unir Presencias y Pseudo-ausencias
    PresenciasParaUnir <- EspecieEnProcesamientoDf %>%
      dplyr::select(especievalida, dplyr::all_of(NombreVariablesAmbientales), longitud, latitud) %>%
      dplyr::mutate(Presence = 1) %>%
      dplyr::rename(X = longitud, Y = latitud)

    PseudoParaUnir <- PseudoFinalDf %>%
      dplyr::select(especievalida, dplyr::all_of(NombreVariablesAmbientales), X, Y, Presence)

    UnionFinal <- dplyr::bind_rows(PresenciasParaUnir, PseudoParaUnir) %>%
      dplyr::mutate(ID = dplyr::row_number()) %>%
      sf::st_as_sf(coords = c("X", "Y"), crs = sf::st_crs(RasterCrsReferencia)) %>%
      sf::st_transform(sf::st_crs(ObjetoCrsReferencia))

    ResultadosFinalesList[[i]] <- UnionFinal

    if(verbose) message(sprintf("    > Finalizado: %d presencias + %d pseudo-ausencias",
                                nrow(PresenciasParaUnir), nrow(PseudoParaUnir)))
  }

  if(verbose) message("\n--- Simulacion completada ---\n")
  return(ResultadosFinalesList)
}
