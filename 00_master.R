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
# PASO 1. Instalar los paquetes que falten (automatico)
# -----------------------------------------------------------------------------
# Un "paquete" es codigo que otra gente escribio y regalo. R base no sabe
# estimar un Heckman; el paquete sampleSelection si. Instalar = descargarlo a
# tu disco (una sola vez). Cargar con library() = ponerlo sobre la mesa de
# trabajo (cada vez que abres R).
#
# No hay que tocar nada: el script revisa que paquetes faltan en ESTA
# computadora e instala solo esos. La primera vez puede tardar varios minutos
# (necesita internet); las siguientes no instala nada y pasa de largo.
#   here            -> rutas del proyecto        wooldridge -> la base mroz
#   AER             -> ivreg() para MC2E         sandwich   -> errores robustos
#   lmtest          -> coeftest()                maxLik     -> Tobit (script 08)
#   sampleSelection -> heckit() de la validacion (script 02)

paquetes <- c("here", "wooldridge", "AER", "sandwich", "lmtest", "maxLik",
              "sampleSelection")
faltan <- paquetes[!sapply(paquetes, requireNamespace, quietly = TRUE)]
if (length(faltan) > 0) {
  cat("Instalando paquetes que faltan:", paste(faltan, collapse = ", "), "\n")
  install.packages(faltan, repos = "https://cloud.r-project.org")
}
# Si alguno no se pudo instalar, se avisa aqui y no a mitad del trabajo.
faltan <- paquetes[!sapply(paquetes, requireNamespace, quietly = TRUE)]
if (length(faltan) > 0) {
  stop("No se pudieron instalar: ", paste(faltan, collapse = ", "),
       ". Revisa la conexion a internet y vuelve a correr el master.")
}


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
source(here("R", "05_tabla6.R"))      # Tabla VI de Mroz (endogeneidad del salario)
source(here("R", "06_tablas7y8.R"))   # Tablas VII y VIII (exogeneidad de ingreso e hijos)
source(here("R", "07_seleccion_multietapa.R"))   # Tablas IX, X y XI (multietapa)
source(here("R", "08_tobit.R"))                  # Tobit, Tabla XI y experimento 4.097 -> 67
source(here("R", "09_tabla4_pruebas.R"))         # Tabla IV: pruebas F y de sobreidentificacion
source(here("R", "10_tabla12.R"))                # Tabla XII (exogeneidad con correccion por seleccion)

cat("\n>>> Master terminado sin errores.\n")
