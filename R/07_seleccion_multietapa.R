# =============================================================================
# 07_seleccion_multietapa.R  —  Tablas IX y X (especificaciones 3 y 5) y
#                               Tabla XI (filas 3 y 5) de Mroz (1987)
# Autora: Judith (codigo) · Experta del bloque: Yamile
# -----------------------------------------------------------------------------
# PROBLEMA: las Tablas IV a VIII usan solo a las 428 que trabajan. Si quienes
# trabajan no son una muestra al azar, esas estimaciones pueden estar sesgadas
# (sesgo de seleccion). Aqui se corrige con el procedimiento multietapa de
# Mroz (ecuaciones 5 y 6, nota 20):
#
#   Etapa 1. Ecuacion de participacion con las 753: P(trabaja) en funcion de
#            todas las variables exogenas Z1 (edad, educacion, B, C, I, F3 y,
#            en la Tabla X, tambien E).
#   Etapa 2. Para cada mujer que trabaja se calcula la media condicional del
#            error de participacion, E(u | trabaja). Depende de la
#            distribucion supuesta:
#              normal   ->  phi(c) / Phi(c)              (inversa de Mills)
#              logistica-> [1+exp(-c)] ln[1+exp(-c)] + c exp(-c)   (p. 780)
#            donde c = f'Z1 es el indice estimado en la etapa 1.
#   Etapa 3. MC2E de la ecuacion de horas en las 428, agregando esa media
#            condicional como regresor. Su coeficiente estima rho1 + a1*rho2
#            (nota 20), no solo la correlacion de los errores.
#
# ERRORES ESTANDAR: la media condicional es ESTIMADA, y su error se arrastra a
# la etapa 3. La formula robusta simple lo ignora y subestima la varianza.
# Mroz usa formulas analiticas (Lee 1982, Apendice 2); aqui se usa BOOTSTRAP:
# se remuestrean mujeres con reposicion, se repiten las tres etapas y la
# desviacion estandar de las estimaciones es el error estandar. Esto es una
# diferencia de procedimiento declarada en documentacion/discrepancias.csv.
#
# La especificacion 4 (lognormal) no se replica: requiere estimar por maxima
# verosimilitud un modelo binario con distribucion sesgada; el paper muestra
# que da resultados casi iguales a la normal y la logistica.
# =============================================================================

library(wooldridge)
data("mroz", package = "wooldridge")
source(here::here("R", "funciones_comunes.R"))

todas <- agregar_terminos(mroz)       # las 753, con los polinomios
VARS  <- c("lwage", "nwifeinc", "kidslt6", "kidsge6")

# -----------------------------------------------------------------------------
# 1. Media condicional del error de participacion
# -----------------------------------------------------------------------------
media_condicional <- function(c, dist) {
  if (dist == "normal") dnorm(c) / pnorm(c)
  else (1 + exp(-c)) * log(1 + exp(-c)) + c * exp(-c)
}

# -----------------------------------------------------------------------------
# 2. Procedimiento multietapa completo
# -----------------------------------------------------------------------------
#   con_E = FALSE -> Tabla IX (sin experiencia);  TRUE -> Tabla X (con E)
multietapa <- function(datos, con_E, dist) {
  z1 <- c("age", "educ", "nwifeinc", "kidslt6", "kidsge6", B, F3, if (con_E) E)
  f  <- as.formula(paste("inlf ~", paste(z1, collapse = " + ")))
  enlace <- if (dist == "normal") "probit" else "logit"
  etapa1 <- glm(f, family = binomial(link = enlace), data = datos)      # 753
  datos$media_cond <- media_condicional(predict(etapa1, type = "link"), dist)
  estimar_horas(subset(datos, inlf == 1),                                # 428
                inst = c(B, F3, if (con_E) E),
                exog = c(EXOG, "media_cond"))
}
mc2e_base <- function(datos, con_E)
  estimar_horas(subset(datos, inlf == 1), inst = c(B, F3, if (con_E) E))

# -----------------------------------------------------------------------------
# 3. Estimaciones puntuales
# -----------------------------------------------------------------------------
est <- list(
  IX_3 = multietapa(todas, FALSE, "normal"),
  IX_5 = multietapa(todas, FALSE, "logit"),
  IX_6 = mc2e_base(todas, FALSE),
  X_3  = multietapa(todas, TRUE,  "normal"),
  X_5  = multietapa(todas, TRUE,  "logit"),
  X_6  = mc2e_base(todas, TRUE)
)
b <- sapply(est, function(m) coef(m)[VARS])          # 4 x 6

# Diferencias que reporta el paper
d <- cbind(
  "IX 6-3" = b[, "IX_6"] - b[, "IX_3"], "IX 6-5" = b[, "IX_6"] - b[, "IX_5"],
  "X 6-3"  = b[, "X_6"]  - b[, "X_3"],  "X 6-5"  = b[, "X_6"]  - b[, "X_5"],
  "XI 3"   = b[, "X_3"]  - b[, "IX_3"], "XI 5"   = b[, "X_5"]  - b[, "IX_5"]
)

# -----------------------------------------------------------------------------
# 4. Bootstrap: errores estandar de las estimaciones y de las diferencias
# -----------------------------------------------------------------------------
# R = 500 replicas. set.seed() fija los numeros aleatorios para que el
# resultado sea identico cada vez que alguien corre el master.
R_BOOT <- 500
set.seed(1987)
cat("\nBootstrap con", R_BOOT, "replicas (tarda alrededor de un minuto)...\n")
replicas <- replicate(R_BOOT, {
  db <- todas[sample(nrow(todas), replace = TRUE), ]
  bb <- sapply(list(multietapa(db, FALSE, "normal"), multietapa(db, FALSE, "logit"),
                    mc2e_base(db, FALSE),
                    multietapa(db, TRUE, "normal"),  multietapa(db, TRUE, "logit"),
                    mc2e_base(db, TRUE)),
               function(m) coef(m)[VARS])
  c(bb, bb[, 3] - bb[, 1], bb[, 3] - bb[, 2], bb[, 6] - bb[, 4], bb[, 6] - bb[, 5],
    bb[, 4] - bb[, 1], bb[, 5] - bb[, 2])
})
se <- matrix(apply(replicas, 1, sd), nrow = 4)       # 4 x 12
se_est <- se[, 1:6];  se_dif <- se[, 7:12]

# -----------------------------------------------------------------------------
# 5. Comparacion con el paper (Tablas IX, X y XI, pp. 782-783)
# -----------------------------------------------------------------------------
paper <- rbind(          # ln(w), ingreso, hijos<6, hijos 6-18 | errores estandar
  "IX_3"   = c(  64, -1.0, -183, -106,  227,  9.2, 408, 34),
  "IX_5"   = c(  66, -1.1, -184, -105,  226,  8.7, 385, 34),
  "IX_6"   = c(  46, -4.4, -337, -112,  220,  3.3, 131, 30),
  "X_3"    = c( 122,  3.9,   53,  -87,  225,  4.5, 173, 34),
  "X_5"    = c( 157,  3.5,   38,  -86,  217,  4.5, 171, 34),
  "X_6"    = c( 672, -6.4, -283,  -85,  217,  3.6, 147, 36),
  "IX 6-3" = c( -19, -3.4, -153,   -7,   84,  8.8, 375, 14),
  "IX 6-5" = c( -20, -3.4, -153,   -7,   74,  8.2, 352, 14),
  "X 6-3"  = c( 551, 10.2,  337,    2,  218,  3.7, 128, 24),
  "X 6-5"  = c( 516, 10.0,  321,    2,  209,  3.7, 121, 22),
  "XI 3"   = c(  57,  4.8,  237,   18,  163,  8.9, 383, 24),
  "XI 5"   = c(  91,  4.6,  222,   18,  100,  8.4, 360, 24)
)
coefs <- cbind(b, d); errs <- cbind(se_est, se_dif)
tabla <- data.frame(
  especificacion = rownames(paper),
  lnw_R  = round(coefs[1, ]),     lnw_paper  = paper[, 1],
  se_R   = round(errs[1, ]),      se_paper   = paper[, 5],
  nwif_R = round(coefs[2, ], 1),  nwif_paper = paper[, 2],
  jov_R  = round(coefs[3, ]),     jov_paper  = paper[, 3],
  may_R  = round(coefs[4, ]),     may_paper  = paper[, 4],
  row.names = NULL
)
cat("\n=========== TABLAS IX, X y XI (multietapa): R vs. paper ===========\n")
cat("3 = normal, 5 = logistica, 6 = MC2E sin correccion por seleccion\n")
cat("Errores estandar de R por bootstrap (", R_BOOT, " replicas)\n\n", sep = "")
print(tabla)

# Contrastes: |z| > 1.96 -> se rechaza
z <- tabla$lnw_R[7:12] / tabla$se_R[7:12]
contrastes <- data.frame(
  contraste = tabla$especificacion[7:12], z_lnw = round(z, 2),
  rechaza_R = abs(z) > 1.96,
  rechaza_paper = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
  pregunta = c("hay sesgo de seleccion (sin E)", "hay sesgo de seleccion (sin E)",
               "hay sesgo de seleccion (con E)", "hay sesgo de seleccion (con E)",
               "experiencia exogena (normal)",   "experiencia exogena (logistica)"),
  row.names = NULL
)
cat("\nContrastes (|z| > 1.96 -> se rechaza al 5%):\n")
print(contrastes)

dir.create(here::here("salidas"), showWarnings = FALSE)
write.csv(tabla, here::here("salidas", "tablas9_10_11_multietapa.csv"), row.names = FALSE)
cat("Guardado en salidas/tablas9_10_11_multietapa.csv\n")

# -----------------------------------------------------------------------------
# 6. Verificacion
# -----------------------------------------------------------------------------
# Coeficientes de ln(w): a 2 unidades del paper (salvo las filas de la Tabla XI,
# que el paper redondea por separado). Conclusiones: identicas al paper.
stopifnot(all(abs(tabla$lnw_R[1:10] - tabla$lnw_paper[1:10]) <= 2),
          all(contrastes$rechaza_R == contrastes$rechaza_paper))
cat("OK: los coeficientes multietapa y las conclusiones de los 6 contrastes",
    "coinciden con el paper.\n")
cat("Lectura: sin experiencia, corregir la seleccion no cambia el efecto del\n")
cat("salario; con experiencia, si, porque la experiencia sesga el MC2E.\n")
