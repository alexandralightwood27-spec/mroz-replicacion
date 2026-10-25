# =============================================================================
# 00_master.R  —  Replicacion de Mroz (1987)
# Econometria I - UNMSM 2026-II
#
# QUE ES ESTE ARCHIVO
# Es el "director de orquesta". No hace calculos por si mismo: manda a correr
# los otros scripts en orden. Si alguna vez quieres rehacer TODO el trabajo
# desde cero, abres este archivo y le das a "Source". Nada mas.
#
# COMO SE CORRE
# Abre el proyecto (doble clic en mroz-replicacion.Rproj), abre este archivo,
# y pulsa el boton "Source" arriba a la derecha del editor.
# =============================================================================


# -----------------------------------------------------------------------------
# PASO 1. Instalar los paquetes (SOLO LA PRIMERA VEZ en cada computadora)
# -----------------------------------------------------------------------------
# Un "paquete" es codigo que otra gente escribio y regalo. R base no sabe
# estimar un Heckman; el paquete sampleSelection si. Instalar = descargarlo a
# tu disco (una sola vez). Cargar con library() = ponerlo sobre la mesa de
# trabajo (cada vez que abres R).
#
# Quita el # de las tres lineas de abajo, corre, y vuelve a poner el #.
# Puede tardar 10-20 minutos la primera vez. Es normal.

# paquetes <- c("here", "wooldridge", "dplyr", "AER", "sampleSelection",
#               "censReg", "truncreg", "sandwich", "lmtest", "modelsummary")
# install.packages(paquetes)


# -----------------------------------------------------------------------------
# PASO 2. Cargar los paquetes que se usan en toda la sesion
# -----------------------------------------------------------------------------
library(here)   # encuentra los archivos sin importar en que carpeta este cada uno

# here() apunta a la carpeta donde vive el .Rproj. Por eso el script funciona
# igual en la laptop de los tres integrantes, aunque sus rutas sean distintas.
cat("Carpeta raiz del proyecto:", here(), "\n\n")


# -----------------------------------------------------------------------------
# PASO 3. Correr los scripts en orden
# -----------------------------------------------------------------------------
# source() = "corre este otro archivo completo, de arriba a abajo".

source(here("R", "01_datos.R"))       # carga la base y verifica que sea la correcta
source(here("R", "02_validacion.R"))  # escalera de validacion (Fase 1)
source(here("R", "03_tabla4.R"))      # Tabla IV de Mroz (MCO y MC2E)
source(here("R", "04_tabla3.R"))      # Tabla III de Mroz (medias)

cat("\n>>> Master terminado sin errores.\n")
