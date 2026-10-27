# =============================================================================
# 08_tobit.R  —  Modelos Tobit de Mroz (1987)
#   Tablas IX y X, especificaciones 1 (Tobit) y 2 (Tobit condicional),
#   sus diferencias, Tabla XI filas 1 y 2, y el experimento de la seccion 1.7
#   (experiencia como regresor directo: el salario cae de 4.097 a 67).
# Autora: Judith (codigo) · Experto del bloque: Max
# -----------------------------------------------------------------------------
# EL MODELO (Heckman 1974; restricciones 1 a 3 del paper)
#   Salario de mercado:  ln w  = c'Z1 + u2
#   Horas deseadas:      h*    = a0 + a1 ln w + a2 Y + a3'Z + e
#   Horas observadas:    h     = h*  si h* > 0 (trabaja)
#                        h     = 0   si h* <= 0 (no trabaja; ln w no se observa)
#   (u2, e) normal bivariada con correlacion rho.
#
# LA RESTRICCION CLAVE DEL TOBIT: la misma ecuacion h* decide SI trabaja y
# CUANTO trabaja. Todo lo que aumenta la probabilidad de trabajar aumenta las
# horas en la misma proporcion. Mroz la rechaza (filas "1 y 2").
#
# VEROSIMILITUD
#   Trabaja:     densidad conjunta de (ln w, h) = densidad normal bivariada
#                de los errores (u2, e). El jacobiano de la transformacion es 1.
#   No trabaja:  P(h* <= 0) = Phi( -m / s_u1 ), con
#                m    = a0 + a1 c'Z1 + a2 Y + a3'Z  (forma reducida de h*)
#                s_u1 = desviacion estandar de a1 u2 + e
#   Tobit condicional: solo las 428 que trabajan, condicionando en h* > 0:
#                se resta log Phi( m / s_u1 ) a la contribucion de cada una.
#
# ERRORES ESTANDAR: Mroz usa la matriz de productos externos del gradiente
# (BHHH): Var = [ sum_i g_i g_i' ]^(-1), con g_i el gradiente de la
# contribucion de la mujer i. Con este metodo los errores salen identicos al
# paper. Las diferencias entre especificaciones usan la covarianza cruzada de
# White (1982), como en el Apendice 2.
#
# Tarda unos segundos: son cinco maximizaciones numericas (BFGS + BHHH).
# =============================================================================

library(wooldridge)
library(maxLik)
data("mroz", package = "wooldridge")
source(here::here("R", "funciones_comunes.R"))

todas <- agregar_terminos(mroz)
todas$lw0 <- ifelse(todas$inlf == 1, todas$lwage, 0)   # 0 solo para no tener NA
todas$h   <- todas$hours / 1000   # horas en miles: mejora la estabilidad numerica
VARS <- c("lwage", "nwifeinc", "kidslt6", "kidsge6")

# -----------------------------------------------------------------------------
# 1. Estimacion de un Tobit conjunto (salario + horas) por maxima verosimilitud
# -----------------------------------------------------------------------------
#   con_E          : la experiencia entra en la ecuacion de salario (Tabla X)
#   exper_en_horas : la experiencia entra TAMBIEN en la de horas (experimento 1.7)
#   condicional    : Tobit condicional (solo las 428)
tobit_conjunto <- function(con_E, exper_en_horas = FALSE, condicional = FALSE) {
  d  <- if (condicional) subset(todas, inlf == 1) else todas
  z1 <- c("age", "educ", "nwifeinc", "kidslt6", "kidsge6", B, F3,
          if (con_E || exper_en_horas) E)
  zh <- c("nwifeinc", "kidslt6", "kidsge6", "age", "educ", if (exper_en_horas) E)
  # Las columnas de Z1 se estandarizan (media 0, desviacion 1): los polinomios
  # llegan a 10^5 y sin esto la optimizacion es inestable. No afecta a los
  # coeficientes de la ecuacion de horas, que son los que interesan.
  Z1 <- cbind(1, scale(as.matrix(d[, z1])))
  Zh <- as.matrix(d[, zh])
  w  <- d$inlf == 1; lw <- d$lw0; h <- d$h
  k1 <- ncol(Z1); kh <- ncol(Zh)
  i_c <- 1:k1; i_a0 <- k1 + 1; i_a1 <- k1 + 2; i_a <- (k1 + 3):(k1 + 2 + kh)
  i_s2 <- k1 + kh + 3; i_se <- i_s2 + 1; i_r <- i_se + 1

  loglik <- function(p) {                     # vector: una contribucion por mujer
    cc <- p[i_c]; a0 <- p[i_a0]; a1 <- p[i_a1]; aa <- p[i_a]
    s2 <- exp(p[i_s2]); se <- exp(p[i_se]); r <- tanh(p[i_r])
    mu2 <- drop(Z1 %*% cc); xa <- drop(Zh %*% aa)
    out <- numeric(length(h))
    e2  <- lw[w] - mu2[w]
    eps <- h[w] - a0 - a1 * lw[w] - xa[w]
    q   <- (e2 / s2)^2 - 2 * r * (e2 / s2) * (eps / se) + (eps / se)^2
    out[w] <- -log(2 * pi * s2 * se * sqrt(1 - r^2)) - q / (2 * (1 - r^2))
    su1 <- sqrt(a1^2 * s2^2 + se^2 + 2 * a1 * r * s2 * se)
    m   <- a0 + a1 * mu2 + xa
    if (condicional) out <- out - pnorm(m / su1, log.p = TRUE)
    else out[!w] <- pnorm(-m[!w] / su1, log.p = TRUE)
    out
  }

  # Valores iniciales: MCO del salario y MCO de horas en las que trabajan
  ols_w <- lm(lw[w] ~ Z1[w, ] - 1)
  ols_h <- lm(h[w] ~ lw[w] + Zh[w, ])
  inicio <- c(coef(ols_w), coef(ols_h),
              log(sd(resid(ols_w))), log(sd(resid(ols_h))), 0)

  ajuste <- maxLik(loglik, start = inicio, method = "BFGS",
                   control = list(iterlim = 5000))
  # Paso final con BHHH (usa los gradientes por mujer, como Mroz): converge
  # rapido cerca del maximo y deja el gradiente practicamente en cero.
  ajuste <- maxLik(loglik, start = coef(ajuste), method = "BHHH",
                   control = list(tol = 1e-12, reltol = 1e-14, gradtol = 1e-7, iterlim = 500))
  p <- coef(ajuste)
  G <- numericGradient(loglik, p)             # gradiente por mujer (n x k)
  V <- solve(crossprod(G))                    # BHHH

  sel <- c(i_a1, i_a[1:3])                    # ln w, ingreso, hijos<6, hijos 6-18
  su1 <- sqrt(p[i_a1]^2 * exp(2 * p[i_s2]) + exp(2 * p[i_se]) +
              2 * p[i_a1] * tanh(p[i_r]) * exp(p[i_s2] + p[i_se]))
  list(coef = setNames(p[sel] * 1000, VARS),          # de vuelta a horas
       se   = setNames(sqrt(diag(V))[sel] * 1000, VARS),
       s2_u1 = (su1 * 1000)^2, rho = tanh(p[i_r]),
       converge = max(abs(colSums(G))) < 1e-2,
       G = G, V = V, sel = sel, trabaja = w)
}

cat("\nEstimando cinco modelos Tobit por maxima verosimilitud...\n")
mod <- list(
  IX_1  = tobit_conjunto(con_E = FALSE),
  IX_2  = tobit_conjunto(con_E = FALSE, condicional = TRUE),
  X_1   = tobit_conjunto(con_E = TRUE),
  X_2   = tobit_conjunto(con_E = TRUE,  condicional = TRUE),
  Exper = tobit_conjunto(con_E = TRUE,  exper_en_horas = TRUE)
)
stopifnot(all(sapply(mod, `[[`, "converge")))

# -----------------------------------------------------------------------------
# 2. Diferencias entre especificaciones (covarianza cruzada de White, BHHH)
# -----------------------------------------------------------------------------
# El Tobit condicional solo usa a las 428; para cruzarlo con el Tobit (753) se
# completa su gradiente con ceros en las que no trabajan.
completar <- function(m) {
  if (nrow(m$G) == nrow(todas)) return(m$G)
  G <- matrix(0, nrow(todas), ncol(m$G)); G[todas$inlf == 1, ] <- m$G; G
}
diferencia <- function(a, b) {               # a - b
  Ga <- completar(a); Gb <- completar(b)
  C  <- a$V %*% crossprod(Ga, Gb) %*% b$V
  Vd <- a$V[a$sel, a$sel] + b$V[b$sel, b$sel] - C[a$sel, b$sel] - t(C[a$sel, b$sel])
  d  <- a$coef - b$coef
  list(coef = d, se = setNames(sqrt(diag(Vd)) * 1000, VARS),
       wald = drop(t(d / 1000) %*% solve(Vd) %*% (d / 1000)))
}
dif <- list(
  "IX 2-1" = diferencia(mod$IX_2, mod$IX_1),   # Tobit condicional - Tobit
  "X 1-2"  = diferencia(mod$X_1,  mod$X_2),    # (el paper usa este orden en X)
  "XI 1"   = diferencia(mod$X_1,  mod$IX_1),   # con E - sin E, Tobit
  "XI 2"   = diferencia(mod$X_2,  mod$IX_2)    # con E - sin E, Tobit condicional
)

# -----------------------------------------------------------------------------
# 3. Comparacion con el paper (Tablas IX, X y XI; seccion 1.7)
# -----------------------------------------------------------------------------
paper <- rbind(          # ln(w), ingreso, hijos<6, hijos 6-18 | errores estandar
  "IX_1"   = c( 261, -22.9, -1035,  -97,   357,  5.1, 140,  48),
  "IX_2"   = c( 157,  -7.2,  -561, -166,   395,  6.5, 145,  55),
  "X_1"    = c(4097, -37.6,  -849,   59,  1094, 17.9, 354, 135),
  "X_2"    = c(1903, -13.0,  -403,  -87,   502, 10.9, 207,  78),
  "Exper"  = c(  67,  -9.2,    NA,   NA,   300,  4.9,  NA,  NA),
  "IX 2-1" = c(-106,  15.7,   474,   70,   403,  6.6, 158,  56),
  "X 1-2"  = c(2194,  24.6,   446, -145,   788, 12.8, 237,  95),
  "XI 1"   = c(3837,  14.7,   184,  155,  1153, 16.1, 319, 116),
  "XI 2"   = c(1748,  -5.8,   156,   78,   579,  7.1, 161,  63)
)
todo <- c(mod, dif)
tabla <- data.frame(
  especificacion = rownames(paper),
  lnw_R  = sapply(todo, function(x) round(x$coef["lwage"])),     lnw_paper  = paper[, 1],
  se_R   = sapply(todo, function(x) round(x$se["lwage"])),       se_paper   = paper[, 5],
  nwif_R = sapply(todo, function(x) round(x$coef["nwifeinc"], 1)), nwif_paper = paper[, 2],
  jov_R  = sapply(todo, function(x) round(x$coef["kidslt6"])),   jov_paper  = paper[, 3],
  may_R  = sapply(todo, function(x) round(x$coef["kidsge6"])),   may_paper  = paper[, 4],
  row.names = NULL
)
cat("\n=========== TOBIT (Tablas IX, X, XI y seccion 1.7): R vs. paper ===========\n")
cat("1 = Tobit, 2 = Tobit condicional, Exper = experiencia en la ecuacion de horas\n")
cat("(En varias filas de diferencias el paper invierte el signo de algunas\n")
cat(" columnas; las magnitudes coinciden. Ver documentacion/discrepancias.csv)\n\n")
print(tabla)

cat("\nVarianza de la forma reducida de las horas (sigma^2 u1):\n")
print(data.frame(modelo = c("IX_1", "IX_2", "X_1", "X_2"),
                 R = signif(unname(sapply(mod[1:4], `[[`, "s2_u1")), 2),
                 paper = c(1.6e6, 8.1e5, 1.3e6, 7.1e5), row.names = NULL))

zs <- sapply(dif, function(x) abs(x$coef / x$se))
contrastes <- data.frame(
  contraste = names(dif),
  z_lnw = round(sapply(dif, function(x) x$coef["lwage"] / x$se["lwage"]), 2),
  z_max = round(apply(zs, 2, max), 2),
  wald_p = round(1 - pchisq(sapply(dif, `[[`, "wald"), 4), 3),
  rechaza_paper = c(TRUE, TRUE, TRUE, TRUE),
  pregunta = c("restricciones del Tobit (sin E)", "restricciones del Tobit (con E)",
               "experiencia exogena (Tobit)", "experiencia exogena (Tobit cond.)"),
  row.names = NULL
)
cat("\nContrastes: z de ln(w), mayor |z| de los 4 coeficientes y p del Wald conjunto\n")
print(contrastes)

cat("\nEXPERIMENTO (seccion 1.7): con la experiencia dentro de la ecuacion de horas,\n")
cat("el coeficiente del salario pasa de", round(mod$X_1$coef["lwage"]), "a",
    round(mod$Exper$coef["lwage"]), "y el del ingreso de",
    round(mod$X_1$coef["nwifeinc"], 1), "a", round(mod$Exper$coef["nwifeinc"], 1), "\n")

dir.create(here::here("salidas"), showWarnings = FALSE)
write.csv(tabla, here::here("salidas", "tobit.csv"), row.names = FALSE)
cat("Guardado en salidas/tobit.csv\n")

# -----------------------------------------------------------------------------
# 4. Verificacion: coeficiente y error estandar de ln(w) a 3 unidades del paper
# -----------------------------------------------------------------------------
stopifnot(all(abs(abs(tabla$lnw_R) - abs(tabla$lnw_paper)) <= 3),
          all(abs(tabla$se_R - tabla$se_paper) <= 3))
cat("OK: los 9 coeficientes de ln(w) del Tobit y sus errores estandar",
    "cuadran con el paper.\n")
