# =====================================================================
# 03_tabla4.R  —  Tabla IV de Mroz (1987): eleccion de la especificacion
# Autor(es): ______________________   (rol 2: variables instrumentales)
# ---------------------------------------------------------------------
# Que hace:
#   Estima la ecuacion de oferta de horas de las 428 mujeres que trabajan
#     hours = a*ln(w) + b*nwifeinc + c1*kidslt6 + c2*kidsge6 + d1*age + d2*educ + u
#   por MCO (fila 1) y por MC2E con distintos conjuntos de instrumentos
#   (filas 2-10), con errores estandar robustos (HC0, nota 7 del paper).
#
# Instrumentos (Tabla V del paper). Todo conjunto incluye constante, edad
# y educacion de la esposa (que ademas estan en la ecuacion de horas).
#   Comunes (B, C, I): unem, city, motheduc, fatheduc + hijos + nwifeinc
#   E  : exper, expersq
#   F2 : age^2, educ^2, age*educ
#   F3 : F2 + age^3, educ^3, age^2*educ, educ^2*age
#   F4 : F3 + age^4, educ^4, age^3*educ, age^2*educ^2, educ^3*age
#   H3 : terminos del esposo -> NO replicables por completo (falta la
#        educacion de los padres del esposo). Filas 5 y 9 = parciales.
#
# Ancla: fila 1 debe dar ln(w) aprox. -17 y fila 8 aprox. 46.
# =====================================================================

library(wooldridge)
library(AER)        # ivreg
library(sandwich)   # vcovHC
library(lmtest)     # coeftest

data("mroz", package = "wooldridge")

# --- Submuestra de trabajadoras y terminos polinomicos -----------------
trab <- subset(mroz, inlf == 1)
trab <- within(trab, {
  # F2
  age2 <- age^2;  educ2 <- educ^2;  ageeduc <- age * educ
  # F3
  age3 <- age^3;  educ3 <- educ^3
  age2educ <- age^2 * educ;  educ2age <- educ^2 * age
  # F4
  age4 <- age^4;  educ4 <- educ^4
  age3educ <- age^3 * educ;  age2educ2 <- age^2 * educ^2
  educ3age <- educ^3 * age
  # Esposo (version parcial de H2/H3: sin educacion de sus padres)
  husage2 <- husage^2;  huseduc2 <- huseduc^2;  husageeduc <- husage * huseduc
  husage3 <- husage^3;  huseduc3 <- huseduc^3
  husage2educ <- husage^2 * huseduc;  huseduc2age <- huseduc^2 * husage
})
stopifnot(nrow(trab) == 428)

# --- Conjuntos de instrumentos (solo los EXCLUIDOS de la ecuacion) ----
B  <- c("unem", "city", "motheduc", "fatheduc")
E  <- c("exper", "expersq")
F2 <- c("age2", "educ2", "ageeduc")
F3 <- c(F2, "age3", "educ3", "age2educ", "educ2age")
F4 <- c(F3, "age4", "educ4", "age3educ", "age2educ2", "educ3age")
H3 <- c("husage", "huseduc", "husage2", "huseduc2", "husageeduc",
        "husage3", "huseduc3", "husage2educ", "huseduc2age")

exog <- c("nwifeinc", "kidslt6", "kidsge6", "age", "educ")   # regresores exogenos

# --- Funcion que estima una fila ----------------------------------------
estimar_fila <- function(excluidos = NULL) {
  rhs <- paste(c("lwage", exog), collapse = " + ")
  if (is.null(excluidos)) {                               # MCO
    m <- lm(as.formula(paste("hours ~", rhs)), data = trab)
    r2_fr <- 1   # ln(w) es su propio instrumento: el paper reporta R2 = 1.0
  } else {                                                # MC2E
    inst <- paste(c(exog, excluidos), collapse = " + ")
    m <- ivreg(as.formula(paste("hours ~", rhs, "|", inst)), data = trab)
    # R^2 de la forma reducida del salario (columna derecha de la Tabla IV)
    r2_fr <- summary(lm(as.formula(paste("lwage ~", inst)), data = trab))$r.squared
  }
  ct <- coeftest(m, vcov. = vcovHC(m, type = "HC0"))
  v  <- c("lwage", "nwifeinc", "kidslt6", "kidsge6")
  c(coef = ct[v, 1], se = ct[v, 2], R2 = r2_fr)
}

# --- Las 10 filas de la Tabla IV ----------------------------------------
filas <- list(
  "1. MCO"          = NULL,
  "2. E"            = c(B, E),
  "3. E,F2"         = c(B, E, F2),
  "4. E,F3"         = c(B, E, F3),
  "5. E,F3,H3*"     = c(B, E, F3, H3),
  "6. E,F4"         = c(B, E, F4),
  "7. F2"           = c(B, F2),
  "8. F3 (base)"    = c(B, F3),
  "9. F3,H3*"       = c(B, F3, H3),
  "10. F4"          = c(B, F4)
)
res <- t(sapply(filas, estimar_fila))

# --- Valores del paper (Tabla IV, p. 770) -------------------------------
paper <- data.frame(
  lnw  = c(-17, 1282, 831, 672, 482, 638, -182,  46, -30, 129),
  lnw_se = c(81, 461, 312, 217, 171, 197, 355, 220, 174, 201),
  nwif = c(-4.2, -8.3, -7.0, -6.4, -5.8, -6.3, -3.7, -4.4, -4.2, -4.7),
  yng  = c(-342, -235, -271, -283, -300, -287, -356, -337, -338, -330),
  old  = c(-115,  -60,  -78,  -85,  -93,  -87, -122, -112, -113, -108),
  R2   = c(1.0, .17, .18, .21, .23, .22, .15, .18, .20, .19)
)

tabla4 <- data.frame(
  fila        = names(filas),
  lnw_R       = round(res[, "coef.lwage"]),
  lnw_paper   = paper$lnw,
  se_R        = round(res[, "se.lwage"]),
  se_paper    = paper$lnw_se,
  nwif_R      = round(res[, "coef.nwifeinc"], 1),
  nwif_paper  = paper$nwif,
  jov_R       = round(res[, "coef.kidslt6"]),
  jov_paper   = paper$yng,
  may_R       = round(res[, "coef.kidsge6"]),
  may_paper   = paper$old,
  R2_R        = round(res[, "R2"], 2),
  R2_paper    = paper$R2,
  row.names   = NULL
)

cat("\n=========== TABLA IV: R vs. paper ===========\n")
cat("(* filas 5 y 9: H3 parcial, sin educacion de los padres del esposo)\n\n")
print(tabla4)

dir.create(here::here("salidas"), showWarnings = FALSE)
write.csv(tabla4, here::here("salidas", "tabla4.csv"), row.names = FALSE)
cat("\nGuardado en salidas/tabla4.csv\n")

# --- Anclas: si fallan, el script se detiene ---------------------------
cat("\n--- ANCLAS ---\n")
cat("Fila 1 (MCO) ln(w):", tabla4$lnw_R[1], " (paper -17)\n")
cat("Fila 8 (MC2E base) ln(w):", tabla4$lnw_R[8], " (paper 46)\n")
stopifnot(abs(tabla4$lnw_R[1] - (-17)) <= 3,
          abs(tabla4$lnw_R[8] - 46) <= 10)
cat("OK: las dos anclas de la Tabla IV cuadran.\n")
