# =============================================================================
# 09_tabla4_pruebas.R  —  Las dos pruebas con que Mroz elige la fila base
#                         de la Tabla IV (seccion 1.1, notas 10 y 11)
# Bloque del integrante 2 (variables instrumentales)
# -----------------------------------------------------------------------------
# POR QUE EXISTE ESTE SCRIPT
#   03_tabla4.R reproduce los numeros de la Tabla IV, pero el paper ademas
#   JUSTIFICA por que la especificacion base usa el polinomio cubico (F3) y por
#   que desconfia de la experiencia (E). Lo hace con dos pruebas que aqui se
#   replican. No cambia ningun resultado anterior: solo agrega la justificacion.
#
#   PRUEBA 1 (relevancia). Prueba F en la ecuacion de salarios (forma
#     reducida): ¿agregar los terminos de grado p+1 en edad y educacion mejora
#     la explicacion de ln(w)?  H0: los terminos nuevos tienen coeficiente 0.
#       F = [(SCR_r - SCR_nr)/q] / [SCR_nr/(n - k)]
#     Mroz se queda con el modelo mas simple que NO es rechazado frente al
#     siguiente: el cubico.
#
#   PRUEBA 2 (validez). Restricciones de sobreidentificacion, nota 10:
#       a) estimar por MC2E;
#       b) construir  y* = horas - a1 * ln(w);
#       c) regresionar y* sobre TODOS los instrumentos;
#       d) LR = n * ln(SCR_r / SCR_nr), chi-cuadrado con (instrumentos
#          excluidos - 1) grados de libertad.
#     Intuicion: si los instrumentos son validos, los que NO estan en la
#     ecuacion de horas no deben explicar las horas una vez descontado el
#     salario. Si lo hacen, alguno esta correlacionado con el error.
# =============================================================================

library(wooldridge)
data("mroz", package = "wooldridge")
source(here::here("R", "funciones_comunes.R"))
trab <- preparar_trab(mroz)
stopifnot(nrow(trab) == 428)

# Mismas filas 2 a 10 de la Tabla IV (instrumentos excluidos de cada una)
filas <- list(
  "2. E"        = c(B, E),          "3. E,F2"   = c(B, E, F2),
  "4. E,F3"     = c(B, E, F3),      "5. E,F3,H3*" = c(B, E, F3, H3),
  "6. E,F4"     = c(B, E, F4),      "7. F2"     = c(B, F2),
  "8. F3 (base)" = c(B, F3),        "9. F3,H3*" = c(B, F3, H3),
  "10. F4"      = c(B, F4)
)

# -----------------------------------------------------------------------------
# 1. Pruebas F en la ecuacion de salarios (nota 11)
# -----------------------------------------------------------------------------
# anova() compara dos regresiones anidadas y devuelve la F y su valor p.
salario <- function(inst) lm(reformulate(c(EXOG, inst), "lwage"), data = trab)
prueba_F <- function(chico, grande) {
  a <- anova(salario(filas[[chico]]), salario(filas[[grande]]))
  data.frame(restringido = chico, no_restringido = grande,
             F = round(a$F[2], 2), gl1 = a$Df[2], gl2 = a$Res.Df[2],
             p = round(a$`Pr(>F)`[2], 3), rechaza_5pct = a$`Pr(>F)`[2] < 0.05)
}
pruebasF <- rbind(
  prueba_F("2. E",    "4. E,F3"),       # nota 11: la fila 2 se rechaza frente a la 4
  prueba_F("3. E,F2", "4. E,F3"),       #          la fila 3 se rechaza frente a la 4
  prueba_F("4. E,F3", "6. E,F4"),       #          la 6 no mejora a la 4
  prueba_F("7. F2",   "8. F3 (base)"),  #          la fila 7 se rechaza frente a la 8
  prueba_F("8. F3 (base)", "10. F4")    #          la 10 no mejora a la 8
)
pruebasF$rechaza_paper <- c(TRUE, TRUE, FALSE, TRUE, FALSE)

cat("\n=========== TABLA IV, PRUEBA 1: F en la ecuacion de salarios ===========\n")
print(pruebasF, row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. Pruebas de sobreidentificacion (nota 10)
# -----------------------------------------------------------------------------
sobreid <- function(inst) {
  m  <- estimar_horas(trab, inst)
  ye <- trab$hours - coef(m)["lwage"] * trab$lwage          # y* = h - a1 ln(w)
  scr_r  <- deviance(lm(reformulate(EXOG, "ye"), data = cbind(trab, ye = ye)))
  scr_nr <- deviance(lm(reformulate(c(EXOG, inst), "ye"), data = cbind(trab, ye = ye)))
  LR <- nrow(trab) * log(scr_r / scr_nr)
  gl <- length(inst) - 1                                    # excluidos - 1 endogena
  c(LR = LR, gl = gl, p = pchisq(LR, gl, lower.tail = FALSE))
}
so <- t(sapply(filas, sobreid))
pruebasS <- data.frame(fila = names(filas), LR = round(so[, "LR"], 1), gl = so[, "gl"],
                       p = round(so[, "p"], 3), rechaza_5pct = so[, "p"] < 0.05,
                       usa_experiencia = grepl("E,|^2\\. E", names(filas)), row.names = NULL)

cat("\n=========== TABLA IV, PRUEBA 2: sobreidentificacion ===========\n")
print(pruebasS, row.names = FALSE)

dir.create(here::here("salidas"), showWarnings = FALSE)
write.csv(pruebasF, here::here("salidas", "tabla4_pruebasF.csv"), row.names = FALSE)
write.csv(pruebasS, here::here("salidas", "tabla4_sobreidentificacion.csv"), row.names = FALSE)
cat("\nGuardado en salidas/tabla4_pruebasF.csv y salidas/tabla4_sobreidentificacion.csv\n")

# -----------------------------------------------------------------------------
# 3. Instrumentos debiles (analisis adicional, NO esta en el paper)
# -----------------------------------------------------------------------------
# La literatura posterior a 1987 (Staiger y Stock, 1997) mide la fuerza de los
# instrumentos con la F de los instrumentos EXCLUIDOS en la primera etapa (no
# con el R2 completo, que incluye edad, educacion, hijos e ingreso). Regla
# practica: F > 10. Si es menor, MC2E puede estar sesgado hacia MCO y sus
# errores estandar no son confiables.
# Para no depender de esa fuerza se usa la prueba de Anderson-Rubin (AR): para
# cada valor a0, se regresa h - a0*ln(w) sobre todos los instrumentos y se
# contrasta (con errores HC0) que los excluidos no la expliquen. El conjunto de
# valores a0 que NO se rechazan al 5 % es un intervalo valido aunque los
# instrumentos sean debiles. Si queda vacio, ningun valor de a0 es compatible
# con los instrumentos: senal de que alguno es invalido.
F_excluidos <- function(inst) {
  r <- lm(reformulate(EXOG, "lwage"), data = trab)
  u <- lm(reformulate(c(EXOG, inst), "lwage"), data = trab)
  clasica <- anova(r, u)
  robusta <- waldtest(r, u, vcov = vcovHC(u, type = "HC0"), test = "F")
  c(F = clasica$F[2], p = clasica$`Pr(>F)`[2], F_robusta = robusta$F[2])
}
p_AR <- function(a0, inst) {
  d <- transform(trab, y0 = hours - a0 * lwage)
  r <- lm(reformulate(EXOG, "y0"), data = d)
  u <- lm(reformulate(c(EXOG, inst), "y0"), data = d)
  waldtest(r, u, vcov = vcovHC(u, type = "HC0"), test = "Chisq")$`Pr(>Chisq)`[2]
}
rejilla <- seq(-2000, 3000, by = 5)
debiles <- do.call(rbind, lapply(c("2. E", "4. E,F3", "7. F2", "8. F3 (base)", "10. F4"), function(f) {
  inst <- filas[[f]]
  pa   <- sapply(rejilla, p_AR, inst = inst)
  ok   <- rejilla[pa >= 0.05]
  data.frame(fila = f, excluidos = length(inst), round(t(F_excluidos(inst)), 3),
             p_AR_a0 = round(p_AR(0, inst), 3),
             # en el borde de la rejilla = el conjunto no esta acotado (-Inf / Inf)
             AR_inf = if (!length(ok)) NA else if (min(ok) == min(rejilla)) -Inf else min(ok),
             AR_sup = if (!length(ok)) NA else if (max(ok) == max(rejilla)) Inf else max(ok),
             row.names = NULL)
}))
cat("\n=========== TABLA IV, PRUEBA 3: fuerza de los instrumentos ===========\n")
cat("(F de los excluidos en la primera etapa; conjunto de Anderson-Rubin al 95 %;\n")
cat(" NA = conjunto vacio; -Inf o Inf = no acotado en la rejilla de -2000 a 3000)\n")
print(debiles, row.names = FALSE)
write.csv(debiles, here::here("salidas", "tabla4_instrumentos_debiles.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Verificacion contra lo que dice el paper (p. 772)
# -----------------------------------------------------------------------------
#  - Las cinco pruebas F dan la misma decision que la nota 11.
#  - "Sin la experiencia (filas 7 a 10) ninguna restriccion se rechaza al 5%".
#  - "Con la experiencia hay cierta evidencia de restricciones invalidas".
stopifnot(all(pruebasF$rechaza_5pct == pruebasF$rechaza_paper),
          !any(pruebasS$rechaza_5pct[!pruebasS$usa_experiencia]),
          any(pruebasS$rechaza_5pct[pruebasS$usa_experiencia]))
cat("OK: las pruebas F eligen el polinomio cubico y la sobreidentificacion\n")
cat("solo se rechaza cuando la experiencia esta entre los instrumentos.\n")
