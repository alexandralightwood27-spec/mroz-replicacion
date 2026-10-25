# =============================================================================
# 04_tabla3.R  —  Tabla III de Mroz (1987): medias de la muestra
# Autora: Judith (rol 1: datos, integracion y reproducibilidad)
# -----------------------------------------------------------------------------
# Que hace:
#   Calcula la media y la desviacion estandar de las variables de la Tabla III
#   (p. 769 del paper) en dos grupos: las 753 mujeres y las 428 que trabajan.
#   Compara con los valores del paper, guarda salidas/tabla3.csv y se detiene
#   si alguna cifra no cuadra.
# =============================================================================

library(wooldridge)
data("mroz", package = "wooldridge")
# -----------------------------------------------------------------------------
# 1. Variable que falta: ingreso no laboral del HOGAR
# -----------------------------------------------------------------------------
# OJO: no es 'nwifeinc'. 'nwifeinc' es todo lo que entra al hogar menos lo que
# gana ella (incluye el sueldo del esposo); es la que va en la ecuacion de horas.
# La Tabla III reporta el ingreso que no viene del trabajo de NINGUNO de los dos:
#   ingreso familiar - ingreso laboral de ella - ingreso laboral del esposo
# En las que no trabajan el salario es NA, asi que su ingreso laboral es 0.
ingreso_ella   <- ifelse(mroz$inlf == 1, mroz$wage * mroz$hours, 0)
ingreso_esposo <- mroz$huswage * mroz$hushrs
mroz$nolab <- (mroz$faminc - ingreso_ella - ingreso_esposo) / 1000  # en miles
# -----------------------------------------------------------------------------
# 2. Variables de la tabla, en el orden del paper
# -----------------------------------------------------------------------------
variables <- c(
  age     = "Edad de la esposa",
  educ    = "Educacion de la esposa",
  kidslt6 = "Hijos menores de 6",
  kidsge6 = "Hijos de 6 a 18",
  husage  = "Edad del esposo",
  huseduc = "Educacion del esposo",
  wage    = "Salario de la esposa (por hora)",
  huswage = "Salario del esposo (por hora)",
  nolab   = "Ingreso no laboral del hogar (miles)",
  hours   = "Horas de trabajo de la esposa",
  hushrs  = "Horas de trabajo del esposo"
)

trab <- subset(mroz, inlf == 1)   # las 428 que trabajan
# -----------------------------------------------------------------------------
# 3. Medias y desviaciones estandar
# -----------------------------------------------------------------------------
# sapply aplica la funcion a cada nombre de variable y devuelve un vector.
media <- function(datos) sapply(names(variables), function(v) mean(datos[[v]]))
desv  <- function(datos) sapply(names(variables), function(v) sd(datos[[v]]))

tabla3 <- data.frame(
  variable   = unname(variables),
  media_todas = round(media(mroz), 2),
  de_todas    = round(desv(mroz), 2),
  media_trab  = round(media(trab), 2),
  de_trab     = round(desv(trab), 2),
  row.names   = NULL
)

# El salario de la esposa solo existe para las que trabajan: el paper pone "—"
tabla3[tabla3$variable == variables["wage"], c("media_todas", "de_todas")] <- NA
# -----------------------------------------------------------------------------
# 4. Valores del paper (Tabla III, p. 769)
# -----------------------------------------------------------------------------
paper <- data.frame(
  media_todas = c(42.5, 12.3, 0.24, 1.35, 45.1, 12.5,   NA, 7.48, 3.76,  740.6, 2267.3),
  de_todas    = c( 8.1,  2.3, 0.52, 1.32,  8.1,  3.0,   NA, 4.23, 5.90,  871.3,  595.6),
  media_trab  = c(42.0, 12.7, 0.14, 1.35, 44.6, 12.6, 4.18, 7.23, 3.39, 1302.9, 2233.5),
  de_trab     = c( 7.7,  2.3, 0.39, 1.32,  8.0,  3.0, 3.31, 3.57, 6.07,  776.3,  582.9)
)

cat("\n=========== TABLA III: R vs. paper ===========\n")
comparacion <- data.frame(
  variable = tabla3$variable,
  todas_R = tabla3$media_todas, todas_paper = paper$media_todas,
  trab_R  = tabla3$media_trab,  trab_paper  = paper$media_trab
)
print(comparacion)

# -----------------------------------------------------------------------------
# 5. Guardar y verificar
# -----------------------------------------------------------------------------
dir.create(here::here("salidas"), showWarnings = FALSE)
write.csv(tabla3, here::here("salidas", "tabla3.csv"), row.names = FALSE)
cat("\nGuardado en salidas/tabla3.csv\n")

# El paper redondea a 1 o 2 decimales, asi que la diferencia maxima
# aceptable es medio decimal (0.05), con un pequeno margen.
dif <- abs(as.matrix(tabla3[, -1]) - as.matrix(paper))
stopifnot(max(dif, na.rm = TRUE) <= 0.06)
cat("OK: las", sum(!is.na(dif)), "cifras de la Tabla III cuadran con el paper.\n")