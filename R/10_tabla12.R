# =============================================================================
# 10_tabla12.R  —  Tabla XII de Mroz (1987): exogeneidad del ingreso no
#                  laboral, los hijos y la experiencia CON correccion por
#                  seleccion (multietapa, media condicional normal)
# Bloque del integrante 4 (seleccion multietapa); usa el contraste del bloque 3
# -----------------------------------------------------------------------------
# PREGUNTA: las Tablas VII y VIII concluyeron, SIN corregir la seleccion, que
# el ingreso no laboral (I) y los hijos (C) se pueden tratar como exogenos.
# ¿Se mantiene la conclusion cuando se corrige la seleccion? (seccion 1.7)
#
# METODO: es el multietapa de 07_seleccion_multietapa.R con dos cambios:
#   - instrumentos comunes B, F3 y H3 (como en las Tablas VII y VIII);
#   - en cada especificacion, las variables que NO se suponen exogenas
#     (I, C o ambas) se instrumentan en la ecuacion de horas y ademas SALEN
#     de la ecuacion de participacion, porque Z1 solo admite exogenas.
# Errores estandar y contrastes: Apendice 2, caso 2 (mismas funciones
# influencia, se_ap2 y dif_ap2 de funciones_comunes.R). Los dos estimadores
# de cada contraste tienen ecuaciones de participacion distintas; la varianza
# de la diferencia incluye la covarianza entre ellas (formula "conjunta").
# Sin esa covarianza los errores se acercan mas a los impresos, pero no es la
# varianza correcta: se reporta solo como sensibilidad (se_sin_cov).
# Contrastes: z de cada coeficiente y Wald CONJUNTO de los cuatro,
# W = d' V^-1 d ~ chi2(4), con precision completa.
#
# REPRODUCCION APROXIMADA: H3 no se puede reconstruir completo (falta la
# educacion de los padres del esposo). Los instrumentos no son exactamente los
# de Mroz, asi que la tabla es una aproximacion con los instrumentos
# disponibles; se comparan decisiones, no cifras exactas.
#
# H3 esta incompleto (falta la educacion de los padres del esposo), asi que,
# como en las Tablas VII y VIII, las magnitudes no salen identicas: lo que se
# verifica es que las CONCLUSIONES de los 6 contrastes coincidan con el paper.
# =============================================================================

library(wooldridge)
data("mroz", package = "wooldridge")
source(here::here("R", "funciones_comunes.R"))

todas <- agregar_terminos(mroz)
n <- nrow(todas); trab <- todas$inlf == 1
VARS <- c("lwage", "nwifeinc", "kidslt6", "kidsge6")
K_HIJOS <- c("kidslt6", "kidsge6"); AE <- c("age", "educ")

# -----------------------------------------------------------------------------
# 1. Una especificacion de la Tabla XII
# -----------------------------------------------------------------------------
#   exo_I, exo_C, exo_E : TRUE si el ingreso, los hijos o la experiencia se
#                         suponen exogenos en esa especificacion.
# multietapa_general() esta en funciones_comunes.R (seccion 6).
multietapa_xii <- function(exo_I, exo_C, exo_E) {
  multietapa_general(todas,
    endog = c("lwage", if (!exo_I) "nwifeinc", if (!exo_C) K_HIJOS),
    exog  = c(if (exo_I) "nwifeinc", if (exo_C) K_HIJOS, AE),
    inst  = c(B, F3, H3, if (exo_E) E))
}

# -----------------------------------------------------------------------------
# 2. Las ocho especificaciones (columna "Additional Instruments" del paper)
# -----------------------------------------------------------------------------
espec <- list(                       #  I      C      E
  XII_1 = c(TRUE,  TRUE,  FALSE),    # C, I
  XII_2 = c(TRUE,  FALSE, FALSE),    # I
  XII_3 = c(FALSE, TRUE,  FALSE),    # C
  XII_4 = c(FALSE, FALSE, FALSE),    # (ninguno)
  XII_5 = c(TRUE,  TRUE,  TRUE),     # E, C, I
  XII_6 = c(TRUE,  FALSE, TRUE),     # E, I
  XII_7 = c(FALSE, TRUE,  TRUE),     # E, C
  XII_8 = c(FALSE, FALSE, TRUE)      # E
)
cat("\nEstimando las 8 especificaciones multietapa de la Tabla XII...\n")
inf <- lapply(espec, function(s) influencia(multietapa_xii(s[1], s[2], s[3]), n, trab))

pares <- list("1 - 2" = c("XII_1", "XII_2"), "1 - 3" = c("XII_1", "XII_3"),
              "1 - 4" = c("XII_1", "XII_4"), "5 - 6" = c("XII_5", "XII_6"),
              "5 - 7" = c("XII_5", "XII_7"), "5 - 8" = c("XII_5", "XII_8"))
dif <- lapply(pares, function(p) dif_ap2(inf[[p[1]]], inf[[p[2]]]))
dif_sin_cov <- lapply(pares, function(p) dif_ap2(inf[[p[1]]], inf[[p[2]]], "sin_cov"))
res <- c(lapply(inf, function(x) list(dif = x$b, se = se_ap2(x))), dif)

# -----------------------------------------------------------------------------
# 3. Comparacion con el paper (Tabla XII, p. 785)
# -----------------------------------------------------------------------------
paper <- rbind(          # ln(w), ingreso, hijos<6, hijos 6-18 | errores estandar
  "XII_1" = c( -67, -0.4, -184, -109,  192,  6.2, 278, 32),
  "XII_2" = c( -70, -2.4, -209,   -2,  209,  7.6, 413, 77),
  "XII_3" = c( -87,  0.8, -158, -111,  195, 10.2, 270, 33),
  "XII_4" = c(-120,  1.6, -135,  -10,  225, 11.5, 428, 97),
  "XII_5" = c(   3,  3.7,   33,  -89,  200,  4.4, 173, 34),
  "XII_6" = c(  48,  3.2, -268,  -81,  203,  4.6, 418, 83),
  "XII_7" = c( -14,  0.0,   23,  -88,  201, 10.2, 179, 36),
  "XII_8" = c(  34, -1.2, -301,  -80,  203, 10.4, 421, 83),
  "1 - 2" = c(   2,  2.0,   26, -107,  136,  6.5, 405, 98),
  "1 - 3" = c(  20, -1.2,  -26,    2,  119, 10.0, 209, 15),
  "1 - 4" = c(  52, -2.0,  -48, -100,  156, 11.3, 426, 93),
  "5 - 6" = c( -45,  0.6,  301,   -9,  170,  3.5, 390, 79),
  "5 - 7" = c(  17,  3.6,    9,   -1,  160,  9.9, 117, 22),
  "5 - 8" = c( -31,  4.9,  334,   -9,  168, 10.1, 393, 78)
)
coefs <- sapply(res, `[[`, "dif"); errs <- sapply(res, `[[`, "se")
tabla12 <- data.frame(
  especificacion = rownames(paper),
  lnw_R  = round(coefs[1, ]),     lnw_paper  = paper[, 1],
  se_R   = round(errs[1, ]),      se_paper   = paper[, 5],
  nwif_R = round(coefs[2, ], 1),  nwif_paper = paper[, 2],
  se_nwif_R = round(errs[2, ], 1), se_nwif_paper = paper[, 6],
  jov_R  = round(coefs[3, ]),     jov_paper  = paper[, 3],
  se_jov_R = round(errs[3, ]),    se_jov_paper = paper[, 7],
  may_R  = round(coefs[4, ]),     may_paper  = paper[, 4],
  se_may_R = round(errs[4, ]),    se_may_paper = paper[, 8],
  se_sin_cov = round(c(rep(NA, 8), sapply(dif_sin_cov, function(r) unname(r$se[1])))),
  row.names = NULL
)
cat("\n=========== TABLA XII: R vs. paper ===========\n")
cat("(Reproduccion aproximada: H3 incompleto, se esperan diferencias de magnitud)\n")
cat("(se_sin_cov: sensibilidad sin la covarianza entre etapas 1, solo diferencias)\n\n")
print(tabla12)

# Contrastes, con precision completa. La hipotesis es CONJUNTA (que ninguno
# de los cuatro coeficientes cambie), asi que la decision se toma con el Wald
# chi2(4): p < 0.05 -> se rechaza. Los z individuales se muestran para ver
# que coeficiente se mueve. (Mirar solo el mayor |z| contra 1,96 no es un
# contraste al 5%.)
zz <- (coefs / errs)[, names(pares)]
contrastes <- data.frame(
  contraste = names(pares),
  z_lnw = round(zz[1, ], 2), z_ingreso = round(zz[2, ], 2),
  z_hijos_men = round(zz[3, ], 2), z_hijos_may = round(zz[4, ], 2),
  wald = round(sapply(dif, `[[`, "wald"), 2),
  p_wald = round(sapply(dif, `[[`, "p_wald"), 3),
  p_wald_sin_cov = round(sapply(dif_sin_cov, `[[`, "p_wald"), 3),
  rechaza_R = sapply(dif, `[[`, "p_wald") < 0.05,
  rechaza_paper = FALSE,
  pregunta = c("hijos exogenos", "ingreso exogeno", "hijos e ingreso exogenos",
               "hijos exogenos (con E)", "ingreso exogeno (con E)",
               "hijos e ingreso exogenos (con E)"),
  row.names = NULL
)
cat("\nContrastes: z por coeficiente y Wald conjunto chi2(4) (p < 0.05 -> se rechaza):\n")
print(contrastes)

dir.create(here::here("salidas"), showWarnings = FALSE)
write.csv(tabla12, here::here("salidas", "tabla12.csv"), row.names = FALSE)
write.csv(contrastes, here::here("salidas", "contrastes12.csv"), row.names = FALSE)
cat("Guardado en salidas/tabla12.csv y salidas/contrastes12.csv\n")

# -----------------------------------------------------------------------------
# 4. Verificacion: las 6 decisiones coinciden con el paper (no se rechaza),
#    tambien sin la covarianza entre etapas 1
# -----------------------------------------------------------------------------
stopifnot(all(contrastes$rechaza_R == contrastes$rechaza_paper),
          all(contrastes$p_wald_sin_cov >= 0.05))
cat("OK: las conclusiones de los 6 contrastes de la Tabla XII coinciden con el paper.\n")
cat("Lectura: tambien corrigiendo la seleccion, no se rechaza la exogeneidad\n")
cat("del ingreso no laboral ni de los hijos (con los instrumentos disponibles).\n")
