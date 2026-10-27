# =============================================================================
# 05_tabla6.R  —  Tabla VI de Mroz (1987): endogeneidad del salario
# Autora: Judith (codigo) · Experta del bloque: Melany
# -----------------------------------------------------------------------------
# Que hace:
#   Seis especificaciones de la ecuacion de horas (todas con los instrumentos
#   comunes B, C, I, F3 y errores robustos HC0) y tres contrastes de
#   diferencias entre especificaciones, con la covarianza de White que usa
#   Mroz en el Apendice 2.
#
#   Muestra completa de las que trabajan (428):
#     1. MCO
#     2. MC2E de referencia (sin experiencia)
#     3. MC2E agregando la experiencia (E) como instrumento
#   Submuestra restringida (326 que trabajaban tambien en la entrevista de 1976):
#     4. MCO
#     5. MC2E agregando el salario reportado en 1976 como instrumento
#     6. MC2E de referencia
#
# Preguntas que responde:
#   - 3 vs 2: ¿la experiencia es un instrumento valido?            (Mroz: no)
#   - 5 vs 4: ¿hay error de medicion en el salario promedio?        (Mroz: si)
#   - 6 vs 5: ¿el salario de 1976 es exogeno?                       (Mroz: no)
# =============================================================================

library(wooldridge)
data("mroz", package = "wooldridge")
source(here::here("R", "funciones_comunes.R"))

trab <- preparar_trab(mroz)              # 428 que trabajan
restr <- subset(trab, repwage > 0)       # 326 con salario reportado en 1976
stopifnot(nrow(trab) == 428, nrow(restr) == 326)

# -----------------------------------------------------------------------------
# 1. Las seis especificaciones
# -----------------------------------------------------------------------------
modelos <- list(
  "1. MCO"                       = estimar_horas(trab),
  "2. MC2E base"                 = estimar_horas(trab,  c(B, F3)),
  "3. MC2E + E"                  = estimar_horas(trab,  c(B, E, F3)),
  "4. MCO submuestra"            = estimar_horas(restr),
  "5. MC2E + ln(wage1976), subm" = estimar_horas(restr, c(B, F3, "lrepwage")),
  "6. MC2E base, submuestra"     = estimar_horas(restr, c(B, F3))
)

res <- t(sapply(modelos, function(m) {
  x <- coef_robustos(m)
  c(x$coef, se = x$se)
}))

# -----------------------------------------------------------------------------
# 2. Contrastes de diferencias (filas inferiores de la tabla)
# -----------------------------------------------------------------------------
pares <- list("3 - 2" = c(3, 2), "5 - 4" = c(5, 4), "6 - 5" = c(6, 5))
difs <- t(sapply(pares, function(p) {
  d <- contraste_diferencia(modelos[[p[1]]], modelos[[p[2]]])
  c(d$dif, se = d$se)
}))

# -----------------------------------------------------------------------------
# 3. Comparacion con el paper (Tabla VI, p. 775)
# -----------------------------------------------------------------------------
paper <- rbind(
  c( -17,  -4.2, -342, -115,   81, 3.1, 131, 29),
  c(  46,  -4.4, -337, -112,  220, 3.3, 131, 30),
  c( 672,  -6.4, -283,  -85,  217, 3.6, 147, 36),
  c( -32,  -4.4, -275,  -98,  126, 3.4, 168, 31),
  c( 348,  -5.9, -301,  -88,  133, 3.6, 171, 31),
  c(-244,  -3.6, -260, -103,  260, 3.6, 171, 31),
  c( 627,  -2.0,   54,   27,  209, 1.8,  69, 19),
  c( 380,  -1.5,  -26,   10,  142, 1.1,  41, 11),
  c(-592,   2.4,   41,  -15,  247, 1.7,  64, 18)
)

todo <- rbind(res, difs)
tabla6 <- data.frame(
  especificacion = rownames(todo),
  lnw_R   = round(todo[, "lwage"]),    lnw_paper  = paper[, 1],
  se_R    = round(todo[, "se.lwage"]), se_paper   = paper[, 5],
  nwif_R  = round(todo[, "nwifeinc"], 1), nwif_paper = paper[, 2],
  jov_R   = round(todo[, "kidslt6"]),  jov_paper  = paper[, 3],
  may_R   = round(todo[, "kidsge6"]),  may_paper  = paper[, 4],
  row.names = NULL
)

cat("\n=========== TABLA VI: R vs. paper ===========\n")
print(tabla6)

z <- difs[, "lwage"] / difs[, "se.lwage"]
cat("\nEstadistico z de la diferencia en ln(w):\n")
print(round(z, 2))
cat("|z| > 1.96 -> se rechaza la especificacion indicada al 5%.\n")
cat("3 - 2: se rechaza que la experiencia sea exogena (Mroz: 3.0).\n")
cat("5 - 4: se rechaza que el salario promedio no tenga error de medicion.\n")
cat("6 - 5: se rechaza que el salario de 1976 sea exogeno.\n")

dir.create(here::here("salidas"), showWarnings = FALSE)
write.csv(tabla6, here::here("salidas", "tabla6.csv"), row.names = FALSE)
cat("\nGuardado en salidas/tabla6.csv\n")

# -----------------------------------------------------------------------------
# 4. Anclas: coeficiente y error estandar de ln(w) a 1 unidad del paper
# -----------------------------------------------------------------------------
stopifnot(all(abs(tabla6$lnw_R - tabla6$lnw_paper) <= 1),
          all(abs(tabla6$se_R  - tabla6$se_paper)  <= 1))
cat("OK: los 9 coeficientes de ln(w) y sus errores estandar cuadran con el paper.\n")
