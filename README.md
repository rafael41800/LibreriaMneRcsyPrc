# LibreriaMneRcsyPrc 🌍🌲

**LibreriaMneRcsyPrc** es un paquete de R diseñado para automatizar y optimizar el pipeline completo de **Modelado de Nicho Ecológico (SDM)** y distribución de especies. Construido con un enfoque en la eficiencia computacional, el paquete permite procesar limpiezas espaciales masivas, generar pseudo-ausencias y ajustar modelos de *Random Forest* para múltiples especies de forma iterativa y segura para la memoria RAM.

## 🚀 Características Principales

* **Vectorización en R Base:** Erradicación de cuellos de botella de evaluación no estándar (Tidyverse) en rutinas iterativas, utilizando indexación matricial y `lapply` para máxima velocidad.
* **Procesamiento Espacial Seguro:** Manejo inteligente de geometrías con `sf` y `terra`, incluyendo la destrucción controlada de geometrías (`st_drop_geometry`) antes del modelado estadístico para evitar el colapso de algoritmos de Machine Learning.
* **Modelado Multiespecie:** Estructura asíncrona lista para evaluar decenas o cientos de especies iterativamente con validación cruzada (k-fold), optimización de hiperparámetros (`mtry`) y cálculo automático de métricas (AUC, Kappa, F1-Score).
* **Exportación Organizada:** Almacenamiento automatizado de rasters de predicción (`.tif`), objetos `.rds` y gráficos de importancia de variables, estructurados por el nombre de cada especie.

---

## ⚠️ Requisitos Previos del Sistema

Debido a que este paquete realiza geoprocesamiento topológico avanzado y ensamblaje de matrices espaciales, depende fuertemente de librerías compiladas en C++ (como `sf` y `terra`). 

Para usuarios de **Windows**, es **ESTRICTAMENTE NECESARIO** tener instaladas las herramientas de compilación antes de instalar el paquete:

* **RTools (Versión 4.5 o superior)**
  * Puedes descargarlo directamente desde CRAN: [RTools 4.5 para Windows](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html)
  * *Nota:* Durante la instalación de RTools, asegúrate de marcar la casilla para agregarlo al PATH del sistema si el instalador lo solicita. El paquete detectará automáticamente la falta de esta herramienta al cargarse.

Usuarios de **macOS** deben tener instalado `GDAL` y `GEOS` (fácilmente instalables vía Homebrew: `brew install gdal proj geos`). Usuarios de **Linux** (Ubuntu/Debian) deben ejecutar: `sudo apt-get install libgdal-dev libgeos-dev libproj-dev`.

---

## 📦 Instalación

Una vez configurado RTools (o las dependencias en Mac/Linux), puedes instalar la versión de desarrollo directamente desde GitHub utilizando el paquete `devtools`:

```R
# Instalar devtools si no está disponible
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Instalar LibreriaMneRcsyPrc
devtools::install_github("rafael41800/LibreriaMneRcsyPrc")
