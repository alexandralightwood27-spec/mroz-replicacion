# =============================================================================
# 07_seleccion_multietapa.R  —  Tablas IX y X (especificaciones 3 a 6) y
#                               Tabla XI (filas 3 a 5) de Mroz (1987)
# Autora: Judith (codigo) · Bloque del integrante 4
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
#            distribucion supuesta (p. 780), con c = f'Z1:
#              normal    ->  phi(c) / Phi(c)              (inversa de Mills)
#              lognormal ->  exp(s^2/2) {N[s - ln(-c)/s] / N[-ln(-c)/s] - 1}
#              logistica ->  [1+exp(-c)] ln[1+exp(-c)] + c exp(-c)
#   Etapa 3. MC2E de la ecuacion de horas en las 428, agregando esa media
#            condicional como regresor. Su coeficiente (kappa) no es una
#            correlacion: con u_j = rho_j u_I + v_j (supuesto 1) y el error
#            estructural e1 = u1 - a1 u2, kappa = rho1 - a1*rho2. La nota 20
#            del paper escribe rho1 + a1*rho2; con sus propias definiciones
#            (ecuacion 3) el signo es menos.
# Las funciones de las tres etapas estan en funciones_comunes.R (seccion 5).
#
# ERRORES ESTANDAR: la media condicional es ESTIMADA, y su error se arrastra a
# la etapa 3. Se usa la formula del Apendice 2 del paper (caso 2):
#     Var = A^-1 B C B' A^-1  +  A^-1 K V K' A^-1
# El primer termino es el error robusto del MC2E; el segundo agrega el error
# de la etapa 1 (K = derivada de las condiciones del MC2E respecto de los
# parametros de la etapa 1, V = su covarianza). Como comprobacion se calcula
# tambien un BOOTSTRAP de 500 replicas para las especificaciones normal y
# logistica: remuestrea mujeres, repite las tres etapas y toma la desviacion
# estandar de las estimaciones.
#
# CONTRASTES: para cada diferencia se reportan (a) el z del coeficiente del
# salario, que es el que usa el paper, y (b) el Wald CONJUNTO de los cuatro
# coeficientes, W = d' V^-1 d ~ chi2(4). Todo se calcula con precision
# completa; solo se redondea al mostrar. En la Tabla XI (dos ecuaciones de
# participacion distintas) la varianza incluye la covarianza entre ellas
# (formula "conjunta"); la version sin esa covarianza, que se acerca mas a
# los errores impresos, va como columna de sensibilidad (se_sin_cov).
# =============================================================================

library(wooldridge)
data("mroz", package = "wooldridge")
source(here::here("R", "funciones_comunes.R"))

todas <- agregar_terminos(mroz)
VARS  <- c("lwage", "nwifeinc", "kidslt6", "kidsge6")
n     <- nrow(todas); trab <- todas$inlf == 1
mc2e_base <- function(datos, con_E)
  estimar_horas(subset(datos, inlf == 1), inst = c(B, F3, if (con_E) E))

# -----------------------------------------------------------------------------
# 1. Estimaciones puntuales y funciones de influencia
# -----------------------------------------------------------------------------
cat("\nEstimando los modelos multietapa (la lognormal, por maxima verosimilitud con derivadas exactas)...\n")
est <- list(
  IX_3 = multietapa(todas, FALSE, "normal"),
  IX_4 = multietapa(todas, FALSE, "lognormal"),
  IX_5 = multietapa(todas, FALSE, "logit"),
  IX_6 = mc2e_base(todas, FALSE),
  X_3  = multietapa(todas, TRUE,  "normal"),
  X_4  = multietapa(todas, TRUE,  "lognormal"),
  X_5  = multietapa(todas, TRUE,  "logit"),
  X_6  = mc2e_base(todas, TRUE)
)
inf <- lapply(est, influencia, n = n, trab = trab)

# Diferencias que reporta el paper
pares <- list("IX 6-3" = c("IX_6", "IX_3"), "IX 6-4" = c("IX_6", "IX_4"),
              "IX 6-5" = c("IX_6", "IX_5"), "X 6-3"  = c("X_6", "X_3"),
              "X 6-4"  = c("X_6", "X_4"),   "X 6-5"  = c("X_6", "X_5"),
              "XI 3"   = c("X_3", "IX_3"),  "XI 4"   = c("X_4", "IX_4"),
              "XI 5"   = c("X_5", "IX_5"))
# Formula conjunta (incluye la covarianza entre etapas 1 distintas); la
# version sin esa covarianza se guarda aparte, solo como sensibilidad.
dif <- lapply(pares, function(p) dif_ap2(inf[[p[1]]], inf[[p[2]]]))
dif_sin_cov <- lapply(pares, function(p) dif_ap2(inf[[p[1]]], inf[[p[2]]], "sin_cov"))
res <- c(lapply(inf, function(x) list(dif = x$b, se = se_ap2(x))), dif)

# -----------------------------------------------------------------------------
# 2. Bootstrap (comprobacion): especificaciones normal, logistica y MC2E
# -----------------------------------------------------------------------------
# R = 500 replicas. set.seed() fija los numeros aleatorios para que el
# resultado sea identico cada vez que alguien corre el master.
# Precision del propio bootstrap: con R replicas, el error de Monte Carlo de
# una desviacion estandar es cerca de 1/sqrt(2(R-1)), unos 3% con R = 500.
# Si una replica falla, replicate() se detiene con error: no hay replicas
# descartadas en silencio. El bootstrap incluye por construccion todas las
# covarianzas (entre etapas y entre especificaciones).
R_BOOT <- 500
set.seed(1987)
cat("Bootstrap con", R_BOOT, "replicas (tarda alrededor de un minuto)...\n")
replicas <- replicate(R_BOOT, {
  db <- todas[sample(n, replace = TRUE), ]
  bb <- sapply(list(multietapa(db, FALSE, "normal", FALSE)$m, multietapa(db, FALSE, "logit", FALSE)$m,
                    mc2e_base(db, FALSE),
                    multietapa(db, TRUE, "normal", FALSE)$m,  multietapa(db, TRUE, "logit", FALSE)$m,
                    mc2e_base(db, TRUE)),
               function(m) coef(m)["lwage"])
  c(bb, bb[3] - bb[1], bb[3] - bb[2], bb[6] - bb[4], bb[6] - bb[5], bb[4] - bb[1], bb[5] - bb[2])
})
se_boot <- setNames(apply(replicas, 1, sd),
                    c("IX_3", "IX_5", "IX_6", "X_3", "X_5", "X_6",
                      "IX 6-3", "IX 6-5", "X 6-3", "X 6-5", "XI 3", "XI 5"))

# -----------------------------------------------------------------------------
# 3. Comparacion con el paper (Tablas IX, X y XI, pp. 782-783)
# -----------------------------------------------------------------------------
paper <- rbind(          # ln(w), ingreso, hijos<6, hijos 6-18 | errores estandar
  "IX_3"   = c(  64, -1.0, -183, -106,  227,  9.2, 408, 34),
  "IX_4"   = c( 110, -0.8, -183, -103,  223, 10.2, 418, 36),
  "IX_5"   = c(  66, -1.1, -184, -105,  226,  8.7, 385, 34),
  "IX_6"   = c(  46, -4.4, -337, -112,  220,  3.3, 131, 30),
  "X_3"    = c( 122,  3.9,   53,  -87,  225,  4.5, 173, 34),
  "X_4"    = c( 161,  3.9,   65,  -80,  226,  4.7, 179, 34),
  "X_5"    = c( 157,  3.5,   38,  -86,  217,  4.5, 171, 34),
  "X_6"    = c( 672, -6.4, -283,  -85,  217,  3.6, 147, 36),
  "IX 6-3" = c( -19, -3.4, -153,   -7,   84,  8.8, 375, 14),
  "IX 6-4" = c( -64, -3.6, -153,   -9,   85,  9.7, 382, 18),
  "IX 6-5" = c( -20, -3.4, -153,   -7,   74,  8.2, 352, 14),
  "X 6-3"  = c( 551, 10.2,  337,    2,  218,  3.7, 128, 24),
  "X 6-4"  = c( 511, 10.4,  349,    5,  209,  3.9, 132, 22),
  "X 6-5"  = c( 516, 10.0,  321,    2,  209,  3.7, 121, 22),
  "XI 3"   = c(  57,  4.8,  237,   18,  163,  8.9, 383, 24),
  "XI 4"   = c(  51,  4.8,  249,   24,  171,  9.9, 390, 26),
  "XI 5"   = c(  91,  4.6,  222,   18,  100,  8.4, 360, 24)
)
coefs <- sapply(res, `[[`, "dif"); errs <- sapply(res, `[[`, "se")
coefs <- coefs[, rownames(paper)]; errs <- errs[, rownames(paper)]
se_sc <- sapply(dif_sin_cov, function(r) unname(r$se[1]))[rownames(paper)]
tabla <- data.frame(
  especificacion = rownames(paper),
  lnw_R  = round(coefs[1, ]),     lnw_paper  = paper[, 1],
  se_R   = round(errs[1, ]),      se_paper   = paper[, 5],
  se_sin_cov = round(ifelse(grepl("^XI", rownames(paper)), se_sc, NA)),
  se_boot = round(se_boot[rownames(paper)]),
  nwif_R = round(coefs[2, ], 1),  nwif_paper = paper[, 2],
  jov_R  = round(coefs[3, ]),     jov_paper  = paper[, 3],
  may_R  = round(coefs[4, ]),     may_paper  = paper[, 4],
  row.names = NULL
)
cat("\n=========== TABLAS IX, X y XI (multietapa): R vs. paper ===========\n")
cat("3 = normal, 4 = lognormal, 5 = logistica, 6 = MC2E sin correccion\n")
cat("se_R: Apendice 2 (formula conjunta); se_sin_cov: sensibilidad sin la covarianza\n")
cat("entre etapas 1 (solo Tabla XI); se_boot: bootstrap de comprobacion (", R_BOOT,
    " replicas; no se calcula para la lognormal)\n\n", sep = "")
print(tabla)

# Contrastes, con precision completa (sin redondear):
#   z_lnw : coeficiente del salario, |z| > 1.96 -> se rechaza al 5% (criterio del paper)
#   Wald  : los cuatro coeficientes a la vez, chi2(4), p < 0.05 -> se rechaza
nd <- names(pares)
z  <- coefs[1, nd] / errs[1, nd]
contrastes <- data.frame(
  contraste = nd, z_lnw = round(z, 2),
  rechaza_R = abs(z) > 1.96,
  rechaza_paper = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
  wald = round(sapply(dif[nd], `[[`, "wald"), 2),
  p_wald = round(sapply(dif[nd], `[[`, "p_wald"), 3),
  p_wald_sin_cov = round(sapply(dif_sin_cov[nd], `[[`, "p_wald"), 3),
  pregunta = rep(c("hay sesgo de seleccion (sin E)", "hay sesgo de seleccion (con E)",
                   "experiencia exogena"), each = 3),
  row.names = NULL
)
cat("\nContrastes: z del salario (|z| > 1.96 -> rechaza) y Wald conjunto chi2(4):\n")
print(contrastes)
dir.create(here::here("salidas"), showWarnings = FALSE)
write.csv(contrastes, here::here("salidas", "contrastes9_10_11.csv"), row.names = FALSE)
write.csv(tabla, here::here("salidas", "tablas9_10_11_multietapa.csv"), row.names = FALSE)
cat("Guardado en salidas/tablas9_10_11_multietapa.csv y salidas/contrastes9_10_11.csv\n")

# -----------------------------------------------------------------------------
# 4. Verificacion
# -----------------------------------------------------------------------------
# Coeficientes de ln(w) de las 8 especificaciones y de las diferencias: a 2
# unidades del paper. Decisiones de los 9 contrastes del salario (criterio
# del paper): identicas al paper.
stopifnot(all(abs(tabla$lnw_R - tabla$lnw_paper) <= 2),
          all(contrastes$rechaza_R == contrastes$rechaza_paper))
cat("OK: los coeficientes multietapa (normal, lognormal y logistica) y las\n")
cat("decisiones de los 9 contrastes del salario coinciden con el paper.\n")
cat("Lectura: sin experiencia, los contrastes no encuentran diferencias\n")
cat("significativas al corregir la seleccion (ni en el salario ni en conjunto).\n")
cat("Con experiencia, el coeficiente del salario si cambia (z cerca de 2,5), pero\n")
cat("el Wald conjunto de los cuatro coeficientes queda en el limite del 5%.\n")
