#' Generación optimizada de Pseudo-ausencias para Múltiples Especies
#'
#' @description
#' Genera puntos de pseudo-ausencia (background points) para una lista de especies,
#' evitando áreas cercanas a las presencias conocidas mediante un radio de exclusión.
#' La función implementa vectorización en R base (`lapply`) y filtrado de matrices nativo
#' para maximizar la eficiencia computacional sin dependencias de evaluación no estándar (NSE).
#'
#' @param Resultados Lista de `data.frame`. Contiene los datos de presencia por especie (típicamente la salida de la extracción de presencias).
#' @param RasterCrsReferencia Objeto `SpatRaster`. Define la extensión espacial, la resolución de las celdas y contiene las variables ambientales.
#' @param ObjetoCrsReferencia Sistema de Referencia de Coordenadas (CRS). Destino final para las coordenadas transformadas.
#' @param RadioGradosExclusion Numerico. Radio de exclusión a aplicar si el CRS de las presencias es geográfico (lat/lon). Por defecto es `0.001`.
#' @param RadioMetrosExclusion Numerico. Radio de exclusión a aplicar si el CRS de las presencias es proyectado (ej. UTM). Por defecto es `100`.
#' @param NombreVariablesAmbientales Vector de caracteres. Nombres exactos de las capas predictoras en `RasterCrsReferencia`.
#' @param verbose Logico. Si es `TRUE` (por defecto), imprime el progreso por especie en consola.
#' @param seed Numerico o `NULL`. Semilla aleatoria para garantizar la reproducibilidad en el sorteo de celdas.
#'
#' @return Una lista de objetos `data.frame` tabulares. Cada elemento corresponde a una especie
#'   e incluye los datos de presencia (`Presence = 1`), pseudo-ausencia (`Presence = 0`),
#'   los valores ambientales sin NAs, y las coordenadas transformadas (`X`, `Y`).
#'
#' @details
#' El flujo de trabajo interno sigue estos pasos vectorizados:
#' \enumerate{
#'   \item \strong{Espacialización}: Convierte las coordenadas del `data.frame` a un objeto `sf` y lo proyecta al CRS del raster ambiental.
#'   \item \strong{Zona de Exclusión}: Construye un buffer topológico (`st_buffer` + `st_union`) alrededor de las presencias usando el radio estipulado.
#'   \item \strong{Muestreo Estratificado}: Identifica las celdas del raster que contienen datos válidos (sin NAs) y descarta aquellas donde ya existe una presencia. Sortea aleatoriamente las celdas candidatas.
#'   \item \strong{Filtro Espacial}: Intersecta los candidatos con el polígono de exclusión para retener solo los puntos verdaderamente alejados (`st_intersects`).
#'   \item \strong{Extracción Ambiental}: Obtiene los valores climáticos o topográficos para los puntos retenidos mediante `terra::extract`.
#'   \item \strong{Ensamble Tabular}: Une presencias y pseudo-ausencias usando indexación matricial rápida (R base), elimina la geometría compleja (`st_drop_geometry`) y extrae las coordenadas finales en formato numérico plano para evitar errores en modelos de *Machine Learning*.
#' }
#'
#' @note
#' Para obtener una guía detallada sobre la calibración del radio de exclusión,
#' consulte la viñeta de la librería: \code{vignette("spatial_cleaning_guide", package = "LibreriaMneRcsyPrc")}.
#'
#' @section Advertencias:
#' \itemize{
#'   \item \strong{Extensión del Raster}: Si el área de estudio es muy pequeña o el radio de exclusión es excesivamente grande, la función podría no encontrar celdas disponibles y la especie será omitida con un `warning`.
#'   \item \strong{Compatibilidad con Algoritmos}: La destrucción intencional de la geometría espacial al final del proceso garantiza que paquetes como `randomForest` o `gbm` no colapsen al ingerir matrices espaciales.
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' # 1. Crear presencias simuladas
#' df_sp <- data.frame(
#'    longitud = c(-99.1, -99.2, -99.3),
#'    latitud = c(19.4, 19.5, 19.6),
#'    especievalida = "Sp_A"
#' )
#' lista_presencias <- list("Sp_A" = df_sp)
#'
#' # 2. Crear raster ambiental
#' r_env <- rast(nrows=100, ncols=100, ext=ext(-100, -98, 18, 20), crs="EPSG:4326")
#' values(r_env) <- runif(ncell(r_env))
#' names(r_env) <- "Bio1"
#'
#' # 3. Ejecutar simulación
#' pseudo_abs <- SimulatePseudoabsencesNoNAListDf(
#'    Resultados = lista_presencias,
#'    RasterCrsReferencia = r_env,
#'    ObjetoCrsReferencia = 4326,
#'    RadioGradosExclusion = 0.05,
#'    NombreVariablesAmbientales = c("Bio1"),
#'    verbose = TRUE,
#'    seed = 42
#' )
#' }
#'
#' @seealso
#' \code{\link[sf]{st_buffer}}, \code{\link[terra]{cellFromXY}}, \code{\link[terra]{extract}}
#'
#' @keywords sdm pseudo-ausencias background-points espacial
#'
#' @importFrom sf st_as_sf st_transform st_crs st_buffer st_union st_is_longlat st_coordinates st_intersects st_drop_geometry
#' @importFrom terra cellFromXY values xyFromCell extract vect
#' @importFrom stats complete.cases
#' @export
SimulatePseudoabsencesNoNAListDf <- function(Resultados,
                                             RasterCrsReferencia,
                                             ObjetoCrsReferencia,
                                             RadioGradosExclusion = 0.001,
                                             RadioMetrosExclusion = 100,
                                             NombreVariablesAmbientales,
                                             verbose = TRUE,
                                             seed = NULL) {

  if(!is.null(seed)) set.seed(seed)
  if(verbose) message("\n--- Iniciando simulacion de pseudo-ausencias ---")

  NombresEspecies <- names(Resultados)
  TotalEspecies <- length(NombresEspecies)

  # Optimizacion 1: Sustitucion de ciclo for por lapply para vectorizacion a nivel de lista
  ResultadosFinalesList <- lapply(seq_along(NombresEspecies), function(idx) {
    i <- NombresEspecies[idx]
    if(verbose) message(sprintf("-> Procesando [%d/%d]: %s", idx, TotalEspecies, i))

    EspecieEnProcesamientoDf <- Resultados[[i]]

    PresenciasSinNaSf <- sf::st_as_sf(EspecieEnProcesamientoDf, coords = c("longitud", "latitud"), crs = 4326)
    PresenciasSinNaSf <- sf::st_transform(PresenciasSinNaSf, sf::st_crs(RasterCrsReferencia))

    dist_buffer <- if(sf::st_is_longlat(PresenciasSinNaSf)) RadioGradosExclusion else RadioMetrosExclusion
    PoligonoExclusion <- sf::st_union(sf::st_buffer(PresenciasSinNaSf, dist = dist_buffer))

    CeldaConPresencia <- terra::cellFromXY(RasterCrsReferencia[[1]], sf::st_coordinates(PresenciasSinNaSf))
    CeldaSinNa <- which(!is.na(terra::values(RasterCrsReferencia[[1]])))
    CeldasDisponibles <- setdiff(CeldaSinNa, CeldaConPresencia)

    if (length(CeldasDisponibles) == 0) {
      warning(paste("No hay celdas disponibles para la especie:", i))
      return(NULL)
    }

    NAusenciasPedidas <- nrow(PresenciasSinNaSf)
    NAusenciasPosibles <- min(NAusenciasPedidas, length(CeldasDisponibles))

    IndicesSorteados <- sample(CeldasDisponibles, NAusenciasPosibles, replace = FALSE)
    CoordsPseudo <- terra::xyFromCell(RasterCrsReferencia[[1]], IndicesSorteados)

    PseudoSf <- sf::st_as_sf(as.data.frame(CoordsPseudo), coords = c("x", "y"), crs = sf::st_crs(RasterCrsReferencia))
    Interseccion <- sf::st_intersects(PseudoSf, PoligonoExclusion, sparse = FALSE)

    ValoresEnv <- terra::extract(RasterCrsReferencia, terra::vect(PseudoSf), ID = FALSE)

    # Optimizacion 2: Reemplazo de dplyr por indexacion matricial nativa
    PseudoFinalDf <- as.data.frame(cbind(sf::st_coordinates(PseudoSf), ValoresEnv))
    PseudoFinalDf$FueraRadio <- !Interseccion[,1]

    # Optimizacion 3: Filtrado booleano vectorizado
    filtro_valido <- PseudoFinalDf$FueraRadio == TRUE & stats::complete.cases(PseudoFinalDf[, NombreVariablesAmbientales, drop = FALSE])
    PseudoFinalDf <- PseudoFinalDf[filtro_valido, , drop = FALSE]

    if (nrow(PseudoFinalDf) == 0) {
      warning(sprintf("Especie %s: 0 pseudo-ausencias validas tras filtrado. Se omite.", i))
      return(NULL)
    }

    PseudoFinalDf$especievalida <- i
    PseudoFinalDf$Presence <- 0
    names(PseudoFinalDf)[names(PseudoFinalDf) %in% c("X", "Y")] <- c("X", "Y")

    PresenciasParaUnir <- EspecieEnProcesamientoDf[, c("especievalida", NombreVariablesAmbientales, "longitud", "latitud")]
    names(PresenciasParaUnir)[names(PresenciasParaUnir) %in% c("longitud", "latitud")] <- c("X", "Y")
    PresenciasParaUnir$Presence <- 1

    PseudoParaUnir <- PseudoFinalDf[, c("especievalida", NombreVariablesAmbientales, "X", "Y", "Presence")]

    UnionFinal <- rbind(PresenciasParaUnir, PseudoParaUnir)
    UnionFinal$ID <- seq_len(nrow(UnionFinal))

    UnionFinalSf <- sf::st_as_sf(UnionFinal, coords = c("X", "Y"), crs = sf::st_crs(RasterCrsReferencia))
    UnionFinalSf <- sf::st_transform(UnionFinalSf, sf::st_crs(ObjetoCrsReferencia))

    CoordsTrans <- sf::st_coordinates(UnionFinalSf)
    UnionFinalDf <- sf::st_drop_geometry(UnionFinalSf)
    UnionFinalDf$X <- CoordsTrans[, 1]
    UnionFinalDf$Y <- CoordsTrans[, 2]

    if(verbose) message(sprintf("    > Finalizado: %d presencias + %d pseudo-ausencias",
                                nrow(PresenciasParaUnir), nrow(PseudoParaUnir)))

    return(UnionFinalDf)
  })

  names(ResultadosFinalesList) <- NombresEspecies
  # Filtrar elementos nulos en caso de que alguna especie fallara
  ResultadosFinalesList <- Filter(Negate(is.null), ResultadosFinalesList)

  if(verbose) message("\n--- Simulacion completada ---\n")
  return(ResultadosFinalesList)
}
