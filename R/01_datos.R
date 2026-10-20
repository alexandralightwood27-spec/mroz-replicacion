# =============================================================================
# 01_datos.R  —  Cargar la base de Mroz y verificar que es la correcta
#
# QUE HACE ESTE SCRIPT
# 1. Trae los datos de las 753 mujeres casadas del PSID 1975.
# 2. Construye las variables derivadas que faltan.
# 3. Comprueba los numeros de control. Si alguno no cuadra, PARA el script.
#
# POR QUE ESTE PASO EXISTE
# Todo el trabajo se apoya en que la base este bien cargada. Un error aqui
# (una variable mal leida, filas de menos) no da un mensaje de error: da
# resultados plausibles pero falsos, y te enteras tres semanas despues.
# Por eso el script se detiene solo si algo no cuadra.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Traer los datos
# -----------------------------------------------------------------------------
# Hay dos fuentes independientes de la MISMA base, y el script usa la que
# encuentre. Que existan dos no es redundancia inutil: si alguna vez los
# numeros de control no cuadran con una fuente, se prueba con la otra y se
# sabe de inmediato si el problema es de los datos o del codigo.
#
#   Opcion A: paquete 'wooldridge' -> objeto 'mroz', ya con lwage y nwifeinc.
#   Opcion B: paquete 'AER'        -> objeto 'PSID1976', los mismos datos con
#             los nombres en ingles largo. Hay que renombrar las columnas.
#             AER se necesita igual para ivreg(), asi que no cuesta nada.

if (requireNamespace("wooldridge", quietly = TRUE)) {

  data(mroz, package = "wooldridge")
  cat("Datos cargados desde el paquete 'wooldridge'.\n")

} else if (requireNamespace("AER", quietly = TRUE)) {

  data("PSID1976", package = "AER")
  d <- PSID1976
  mroz <- data.frame(
    inlf     = as.integer(d$participation == "yes"),  # 1 si trabaja
    hours    = d$hours,        # horas trabajadas en el anio
    kidslt6  = d$youngkids,    # hijos menores de 6
    kidsge6  = d$oldkids,      # hijos de 6 a 18
    age      = d$age,
    educ     = d$education,    # anios de escolaridad de ella
    wage     = d$wage,         # salario por hora (solo si trabaja)
    repwage  = d$repwage,
    hushrs   = d$hhours,       # horas del esposo
    husage   = d$hage,
    huseduc  = d$heducation,   # escolaridad del esposo
    huswage  = d$hwage,
    faminc   = d$fincome,      # ingreso familiar total
    mtr      = d$tax,
    motheduc = d$meducation,   # escolaridad de la madre  <- instrumento
    fatheduc = d$feducation,   # escolaridad del padre    <- instrumento
    unem     = d$unemp,
    city     = as.integer(d$city == "yes"),
    exper    = d$experience    # experiencia laboral previa
  )
  cat("Datos cargados desde el paquete 'AER' (PSID1976).\n")

} else {

  stop("No hay ninguna fuente de datos. Instala los paquetes: ",
       "install.packages(c('wooldridge','AER'))")

}


# -----------------------------------------------------------------------------
# 2. Construir las variables derivadas
# -----------------------------------------------------------------------------
# El operador <- significa "guarda esto con este nombre". mroz$educ significa
# "la columna educ del data frame mroz". Un data frame es una tabla: filas =
# mujeres, columnas = variables.

# Ingreso no laboral de la esposa, en miles de dolares.
# Es el ingreso familiar menos lo que ella misma gano trabajando.
# Economicamente: es el "efecto ingreso" puro, la plata que le llega
# independientemente de si trabaja o no.
if (!"nwifeinc" %in% names(mroz)) {
  mroz$nwifeinc <- (mroz$faminc - mroz$wage * mroz$hours) / 1000
}

# Experiencia al cuadrado: permite que el retorno a la experiencia sea
# decreciente (los primeros anios de trabajo ensenian mas que los ultimos).
if (!"expersq" %in% names(mroz)) {
  mroz$expersq <- mroz$exper^2
}

# Log del salario. CLAVE: NA para las que NO trabajan.
# Ese NA no es un dato perdido por descuido: es exactamente el problema del
# paper. El salario de mercado de una mujer que no trabaja no se observa.
# Toda la maquinaria de Heckman existe para lidiar con estos NA.
if (!"lwage" %in% names(mroz)) {
  mroz$lwage <- ifelse(mroz$inlf == 1, log(mroz$wage), NA)
}

# Submuestra de las que si trabajan (para los peldanios 1 y 2).
# ATENCION: esta submuestra se usa SOLO para MCO y MC2E de la ecuacion de
# salarios. A heckit() hay que pasarle la base COMPLETA de 753 filas.
trab <- subset(mroz, inlf == 1)


# -----------------------------------------------------------------------------
# 3. Puntos de control
# -----------------------------------------------------------------------------
# stopifnot() = "esto tiene que ser verdad; si no, detente y avisa".
# Es la red de seguridad del proyecto.

cat("\n--- PUNTOS DE CONTROL DE LA BASE ---\n")
cat("Filas totales        :", nrow(mroz), "   (debe ser 753)\n")
cat("Participantes        :", sum(mroz$inlf == 1), "   (debe ser 428)\n")
cat("Salarios observados  :", sum(!is.na(mroz$lwage)), "   (debe ser 428)\n")
cat("Horas medias, todas  :", round(mean(mroz$hours), 1), " (debe ser 740.6)\n")
cat("Horas medias, trabaj.:", round(mean(mroz$hours[mroz$inlf == 1]), 1),
    "(debe ser 1302.9)\n")

# La identidad de Mroz: la media general es la media de las que trabajan
# ponderada por la tasa de participacion. Si esto cuadra, las tres piezas
# (n total, n participantes, horas) son mutuamente consistentes.
identidad <- sum(mroz$inlf == 1) * mean(mroz$hours[mroz$inlf == 1]) / nrow(mroz)
cat("(428 x 1302.9)/753   :", round(identidad, 1), " (debe ser 740.6)\n")

stopifnot(
  nrow(mroz) == 753,
  sum(mroz$inlf == 1) == 428,
  sum(!is.na(mroz$lwage)) == 428,
  abs(mean(mroz$hours) - 740.6) < 0.1
)

cat("\nOK: la base es la correcta.\n")
