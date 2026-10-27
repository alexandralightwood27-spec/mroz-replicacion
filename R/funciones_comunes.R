# =============================================================================
# funciones_comunes.R  —  Herramientas que usan las Tablas IV, VI, VII y VIII
# Autora: Judith
# -----------------------------------------------------------------------------
# Este archivo no estima nada por si solo: define objetos y funciones que los
# scripts de las tablas cargan con source(). Asi cada definicion (que
# instrumentos forman F3, como se calcula un error estandar robusto, etc.)
# existe en un solo lugar y todas las tablas usan exactamente la misma.
# =============================================================================

library(AER)        # ivreg()
library(sandwich)   # vcovHC()
library(lmtest)     # coeftest()

# -----------------------------------------------------------------------------
# 1. Submuestra de las que trabajan + terminos polinomicos
# -----------------------------------------------------------------------------
preparar_trab <- function(mroz) {
  trab <- subset(mroz, inlf == 1)
  within(trab, {
    # F2: cuadraticos en edad y educacion
    age2 <- age^2;  educ2 <- educ^2;  ageeduc <- age * educ
    # F3: cubicos
    age3 <- age^3;  educ3 <- educ^3
    age2educ <- age^2 * educ;  educ2age <- educ^2 * age
    # F4: cuarticos
    age4 <- age^4;  educ4 <- educ^4
    age3educ <- age^3 * educ;  age2educ2 <- age^2 * educ^2
    educ3age <- educ^3 * age
    # Esposo (version parcial de H2/H3: falta la educacion de sus padres)
    husage2 <- husage^2;  huseduc2 <- huseduc^2;  husageeduc <- husage * huseduc
    husage3 <- husage^3;  huseduc3 <- huseduc^3
    husage2educ <- husage^2 * huseduc;  huseduc2age <- huseduc^2 * husage
    # Salario reportado en la entrevista de 1976 (en logaritmo; NA si es 0)
    lrepwage <- ifelse(repwage > 0, log(repwage), NA)
  })
}

# -----------------------------------------------------------------------------
# 2. Conjuntos de instrumentos (Tabla V del paper), solo los EXCLUIDOS
# -----------------------------------------------------------------------------
B  <- c("unem", "city", "motheduc", "fatheduc")
E  <- c("exper", "expersq")
F2 <- c("age2", "educ2", "ageeduc")
F3 <- c(F2, "age3", "educ3", "age2educ", "educ2age")
F4 <- c(F3, "age4", "educ4", "age3educ", "age2educ2", "educ3age")
H3 <- c("husage", "huseduc", "husage2", "huseduc2", "husageeduc",
        "husage3", "huseduc3", "husage2educ", "huseduc2age")

# Regresores exogenos de la ecuacion de horas (se instrumentan a si mismos)
EXOG <- c("nwifeinc", "kidslt6", "kidsge6", "age", "educ")

# -----------------------------------------------------------------------------
# 3. Estimar la ecuacion de horas
# -----------------------------------------------------------------------------
#   endog : regresores tratados como endogenos (por defecto, solo lwage)
#   exog  : regresores tratados como exogenos
#   inst  : instrumentos excluidos; NULL = MCO
estimar_horas <- function(datos, inst = NULL,
                          endog = "lwage", exog = EXOG) {
  rhs <- paste(c(endog, exog), collapse = " + ")
  if (is.null(inst)) {
    lm(as.formula(paste("hours ~", rhs)), data = datos)
  } else {
    z <- paste(c(exog, inst), collapse = " + ")
    ivreg(as.formula(paste("hours ~", rhs, "|", z)), data = datos)
  }
}

# Coeficientes y errores estandar robustos HC0 (nota 7 del paper)
coef_robustos <- function(m, vars = c("lwage", "nwifeinc", "kidslt6", "kidsge6")) {
  ct <- coeftest(m, vcov. = vcovHC(m, type = "HC0"))
  list(coef = ct[vars, 1], se = ct[vars, 2])
}

# -----------------------------------------------------------------------------
# 4. Contraste de diferencias entre dos especificaciones (Apendice 2 del paper)
# -----------------------------------------------------------------------------
# IDEA: si dos estimadores son consistentes bajo la hipotesis nula, su
# diferencia debe ser cercana a cero. Para saber si una diferencia es "grande"
# necesitamos su error estandar:
#     Var(b1 - b2) = Var(b1) + Var(b2) - 2 Cov(b1, b2)
# El termino Cov(b1, b2) es crucial: los dos estimadores usan las MISMAS
# mujeres, asi que estan correlacionados. Hausman (1978) evita calcularlo
# suponiendo que uno es eficiente; Mroz no lo supone y usa la formula de
# White (1982), robusta a heterocedasticidad:
#     b_k - beta ~ A_k^{-1} * sum_i xhat_ki * u_ki ,   A_k = Xhat_k' X
#     Cov(b1, b2) = A_1^{-1} [ sum_i u_1i u_2i xhat_1i xhat_2i' ] A_2^{-1}'
# donde xhat son los regresores proyectados en los instrumentos (en MCO,
# xhat = x) y u son los residuos de cada especificacion.
piezas_iv <- function(m) {
  if (inherits(m, "ivreg")) {
    X <- model.matrix(m, component = "regressors")
    Z <- model.matrix(m, component = "instruments")
    # Proyeccion de X en Z. Se usa la descomposicion QR en vez de invertir Z'Z:
    # con polinomios de grado 3 las columnas de Z tienen escalas muy distintas
    # (edad^3 ~ 100.000, city ~ 1) y la inversa directa es numericamente inestable.
    Xh <- qr.fitted(qr(Z), X)
  } else {
    X <- model.matrix(m); Xh <- X
  }
  list(Xh = Xh, Ainv = solve(crossprod(Xh, X)), u = residuals(m))
}

contraste_diferencia <- function(m1, m2, vars = c("lwage", "nwifeinc", "kidslt6", "kidsge6")) {
  p1 <- piezas_iv(m1); p2 <- piezas_iv(m2)
  stopifnot(length(p1$u) == length(p2$u))      # mismas mujeres en ambos
  V <- function(a, b) a$Ainv %*% crossprod(a$Xh * (a$u * b$u), b$Xh) %*% t(b$Ainv)
  V11 <- V(p1, p1); V22 <- V(p2, p2); V12 <- V(p1, p2)
  d  <- coef(m1)[vars] - coef(m2)[vars]
  Vd <- V11[vars, vars] + V22[vars, vars] - V12[vars, vars] - t(V12)[vars, vars]
  list(dif = d, se = sqrt(diag(Vd)),
       chi2 = as.numeric(t(d) %*% solve(Vd) %*% d), gl = length(vars))
}
