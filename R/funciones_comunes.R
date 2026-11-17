# =============================================================================
# funciones_comunes.R  —  Herramientas que usan las Tablas IV a XII
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
preparar_trab <- function(mroz) subset(agregar_terminos(mroz), inlf == 1)

# Agrega los terminos polinomicos a cualquier muestra (las 753 o las 428).
# Los modelos de seleccion los necesitan para TODAS las mujeres, porque la
# ecuacion de participacion se estima con las 753.
agregar_terminos <- function(datos) {
  within(datos, {
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

# -----------------------------------------------------------------------------
# 5. Correccion por seleccion en varias etapas (Tablas IX a XI)
# -----------------------------------------------------------------------------
# Tres etapas (ecuaciones 5 y 6 y nota 20 del paper):
#   1. Ecuacion de participacion con las 753: P(trabaja) = F(f'Z1).
#   2. Media condicional del error de participacion, E(u | trabaja).
#   3. MC2E de la ecuacion de horas con las 428, con esa media como regresor.
# Z1 se estandariza (media 0, desviacion 1). Esto no cambia la media
# condicional ni ninguna estimacion de la etapa 3; solo hace estable el
# calculo numerico (los polinomios llegan a 10^5).


# Media condicional E(u | trabaja), con c = f'Z1 (p. 780 del paper)
media_condicional <- function(c, dist, sigma = NA) {
  switch(dist,
    normal    = dnorm(c) / pnorm(c),                              # inversa de Mills
    logit     = (1 + exp(-c)) * log(1 + exp(-c)) + c * exp(-c),
    lognormal = ifelse(c < 0, exp(sigma^2 / 2) *
                  (pnorm(sigma - log(-c) / sigma) / pnorm(-log(-c) / sigma) - 1), 0))
}

# Etapa 1. Devuelve el indice, los parametros, el gradiente de la verosimilitud
# de cada mujer (n x p) y la inversa de menos el hessiano (Vinf).
#   normal / logit : probit / logit con glm().
#   lognormal      : la mujer trabaja si f'Z1 + W > 0, con W lognormal de
#                    mediana 1 y parametro sigma; P(trabaja) = N(-ln(-c)/sigma)
#                    si c < 0. Se estima a mano por maxima verosimilitud,
#                    con gradiente y hessiano exactos (Newton-Raphson).
primera_etapa <- function(datos, con_E, dist) {
  z1 <- c("age", "educ", "nwifeinc", "kidslt6", "kidsge6", B, F3, if (con_E) E)
  W  <- cbind(1, scale(as.matrix(datos[, z1]))); d <- datos$inlf
  if (dist != "lognormal") {
    g <- glm(d ~ W - 1, family = binomial(link = if (dist == "normal") "probit" else "logit"))
    th <- setNames(coef(g), NULL)
    return(list(W = W, theta = th, scores = unname(sandwich::estfun(g)),
                Vinf = unname(vcov(g)), dist = dist,
                lambda = function(th) media_condicional(drop(W %*% th), dist)))
  }
  # Log-verosimilitud con su gradiente por mujer (n x k) y su hessiano, en
  # forma ANALITICA. Con x = -ln(-c)/s y q = 2d - 1, cada mujer aporta
  # ln N(q x); sus derivadas respecto de c y de ln(s) salen de la regla de la
  # cadena. Las derivadas exactas hacen que el resultado no dependa de la
  # computadora (con derivadas numericas, los errores de esta especificacion
  # cambiaban en el tercer digito de una maquina a otra).
  q <- 2 * d - 1; k <- ncol(W) + 1
  partes <- function(p, derivadas = TRUE) {
    cc <- drop(W %*% p[-k]); s <- exp(p[k])
    ok <- cc < 0                               # si c >= 0, P(trabaja) = 1
    if (any(!ok & d == 0)) return(list(ll = -Inf))
    x <- dxc <- numeric(length(cc))
    x[ok] <- -log(-cc[ok]) / s                 # x = -ln(-c)/s
    dxc[ok] <- -1 / (cc[ok] * s)               # dx/dc
    t <- q * x
    ll <- sum(pnorm(t[ok], log.p = TRUE))
    if (!derivadas) return(list(ll = ll))
    mm <- numeric(length(cc))                  # N'(t)/N(t), calculado en logaritmos
    mm[ok] <- exp(dnorm(t[ok], log = TRUE) - pnorm(t[ok], log.p = TRUE))
    l1 <- q * mm; l2 <- -mm * (t + mm)         # d ln N / dx y d2 ln N / dx2
    G <- cbind(W * (l1 * dxc), -l1 * x)        # dx/d ln(s) = -x
    wff <- (l2 + l1 * s) * dxc^2               # d2x/dc2 = s (dx/dc)^2
    wft <- -(l2 * x + l1) * dxc                # d2x/dc d ln(s) = -dx/dc
    wtt <- l2 * x^2 + l1 * x                   # d2x/d ln(s)2 = x
    H <- rbind(cbind(crossprod(W, W * wff), crossprod(W, wft)),
               c(crossprod(wft, W), sum(wtt)))
    list(ll = ll, G = G, H = H)
  }
  ll <- function(p) partes(p, FALSE)$ll
  # Inicio: el logit reescalado y desplazado para que c < 0 en las que no trabajan
  lg <- glm(d ~ W - 1, family = binomial("logit")); c0 <- drop(W %*% coef(lg))
  p  <- c(coef(lg) / 3, 0); p[1] <- p[1] - (max(c0[d == 0]) + 0.5) / 3
  p  <- optim(p, ll, function(p) colSums(partes(p)$G), method = "BFGS",
              control = list(fnscale = -1, maxit = 5000, reltol = 1e-12))$par
  # Newton-Raphson con derivadas exactas hasta el maximo (gradiente ~ 0)
  for (it in 1:200) {
    P <- partes(p); paso <- solve(-P$H, colSums(P$G)); a <- 1
    while (ll(p + a * paso) < P$ll - 1e-12 && a > 1e-8) a <- a / 2
    p <- p + a * paso
    if (max(abs(a * paso)) < 1e-12) break
  }
  P <- partes(p)
  if (max(abs(colSums(P$G))) > 1e-6) warning("lognormal: el maximo no convergio")
  lam <- function(th) media_condicional(drop(W %*% th[-length(th)]), "lognormal", exp(th[length(th)]))
  list(W = W, theta = unname(p), scores = unname(P$G),
       Vinf = solve(-P$H), dist = dist, lambda = lam)
}

# Etapas 2 y 3. Devuelve el MC2E y K = d^2 F / d b d theta' (Apendice 2, caso 2),
# la derivada de las condiciones de primer orden del MC2E respecto de los
# parametros de la etapa 1, calculada numericamente como hace Mroz.
multietapa <- function(datos, con_E, dist, con_K = TRUE) {
  e1 <- primera_etapa(datos, con_E, dist)
  trab <- datos$inlf == 1
  datos$media_cond <- e1$lambda(e1$theta)
  m <- estimar_horas(subset(datos, inlf == 1), inst = c(B, F3, if (con_E) E),
                     exog = c(EXOG, "media_cond"))
  if (!con_K) return(list(m = m))             # el bootstrap solo necesita el MC2E
  X <- model.matrix(m, component = "regressors"); Z <- model.matrix(m, component = "instruments")
  y <- datos$hours[trab]; b <- coef(m)
  cpo <- function(th) {                       # dF/db evaluado en b estimado
    lam <- e1$lambda(th)[trab]; X[, "media_cond"] <- lam; Z[, "media_cond"] <- lam
    -drop(crossprod(qr.fitted(qr(Z), X), y - X %*% b))
  }
  h <- 1e-6 * pmax(abs(e1$theta), 1)
  K <- sapply(seq_along(e1$theta), function(j) {
    tp <- tm <- e1$theta; tp[j] <- tp[j] + h[j]; tm[j] <- tm[j] - h[j]
    (cpo(tp) - cpo(tm)) / (2 * h[j])
  })
  list(m = m, K = K, e1 = e1, trab = trab)
}

# Funcion de influencia de cada mujer (n = 753 filas), en dos partes:
#   a_i = A^-1 xhat_i u_i              (variacion muestral del MC2E)
#   c_i = -A^-1 K Vinf s_i             (error de haber estimado la etapa 1)
# Para un MC2E sin correccion (especificacion 6), c = 0.
influencia <- function(x, n, trab, vars = c("lwage", "nwifeinc", "kidslt6", "kidsge6")) {
  m <- if (inherits(x, "ivreg")) x else x$m
  p <- piezas_iv(m)
  a <- matrix(0, n, ncol(p$Xh), dimnames = list(NULL, colnames(p$Xh)))
  a[trab, ] <- (p$Xh * p$u) %*% t(p$Ainv)
  cc <- a * 0
  if (!inherits(x, "ivreg")) cc[] <- -(x$e1$scores %*% x$e1$Vinf %*% t(x$K) %*% t(p$Ainv))
  list(a = a[, vars], c = cc[, vars], b = coef(m)[vars])
}

# Covarianza entre dos estimadores (Apendice 2, caso 2):
#   Cov = A1^-1 B1 C B2' A2^-1' + A1^-1 K1 V K2' A2^-1'
# es decir, la parte de MC2E mas la de la etapa 1, sin terminos cruzados.
cov_ap2 <- function(x, y) crossprod(x$a, y$a) + crossprod(x$c, y$c)
se_ap2  <- function(x) sqrt(diag(cov_ap2(x, x)))

# Diferencia d = x - y entre dos estimadores, su covarianza V y el contraste
# de Wald CONJUNTO de los cuatro coeficientes: W = d' V^-1 d ~ chi2(4) bajo H0.
#   formula = "conjunta" (por defecto): V = sum_i (psi_xi - psi_yi)(psi_xi - psi_yi)',
#       con psi = a + c. Incluye la covarianza entre las dos etapas 1 aunque
#       sean ecuaciones de participacion DISTINTAS: se estiman con las mismas
#       753 mujeres y estan correlacionadas. Es la formula que describe el
#       paper (p. 781: "accounting for the cross-specification covariance of
#       the binary choice estimated parameters") y la que usamos para decidir.
#   formula = "sin_cov": suma las varianzas de las dos etapas 1 SIN su
#       covarianza. NO es la varianza correcta; se reporta solo como
#       sensibilidad, porque es la que mas se acerca a los errores impresos en
#       las Tablas XI y XII (XI fila 3: 164 vs 163 del paper; con la conjunta,
#       123). No sabemos por que: el texto del paper dice que si incluye la
#       covarianza. Ninguna decision cambia con una u otra (ver discrepancias).
# Si los dos estimadores comparten la etapa 1 (o uno no tiene, como el MC2E
# sin correccion), las dos formulas coinciden.
# Nota: igual que en el Apendice 2, no hay terminos cruzados entre la parte de
# MC2E (a) y la de la etapa 1 (c): asintoticamente no estan correlacionadas,
# porque E(v | trabaja, Z1) = 0 (supuestos 1 y 2 del paper).
dif_ap2 <- function(x, y, formula = c("conjunta", "sin_cov")) {
  formula <- match.arg(formula)
  V <- crossprod(x$a - y$a) +
       if (formula == "conjunta") crossprod(x$c - y$c) else crossprod(x$c) + crossprod(y$c)
  d <- x$b - y$b
  W <- drop(crossprod(d, solve(V, d)))
  list(dif = d, se = sqrt(diag(V)), V = V, wald = W, gl = length(d),
       p_wald = pchisq(W, df = length(d), lower.tail = FALSE))
}

# -----------------------------------------------------------------------------
# 6. Multietapa general (Tabla XII)
# -----------------------------------------------------------------------------
# Igual que multietapa(), pero se le dice que regresores de la ecuacion de
# horas son endogenos, cuales exogenos y cuales son los instrumentos excluidos.
# La ecuacion de participacion (probit, con las 753) usa todas las exogenas:
# 'exog' + 'inst'. Devuelve lo mismo que multietapa(), para pasarlo a influencia().
multietapa_general <- function(datos, endog, exog, inst) {
  trab <- datos$inlf == 1
  W <- cbind(1, scale(as.matrix(datos[, c(exog, inst)])))
  g <- glm(datos$inlf ~ W - 1, family = binomial(link = "probit"))
  e1 <- list(W = W, theta = unname(coef(g)), scores = unname(sandwich::estfun(g)),
             Vinf = unname(vcov(g)), dist = "normal",
             lambda = function(th) media_condicional(drop(W %*% th), "normal"))
  datos$media_cond <- e1$lambda(e1$theta)
  m <- estimar_horas(subset(datos, inlf == 1), inst = inst, endog = endog,
                     exog = c(exog, "media_cond"))
  X <- model.matrix(m, component = "regressors"); Z <- model.matrix(m, component = "instruments")
  y <- datos$hours[trab]; b <- coef(m)
  cpo <- function(th) {                       # dF/db evaluado en b estimado
    lam <- e1$lambda(th)[trab]; X[, "media_cond"] <- lam; Z[, "media_cond"] <- lam
    -drop(crossprod(qr.fitted(qr(Z), X), y - X %*% b))
  }
  h <- 1e-6 * pmax(abs(e1$theta), 1)
  K <- sapply(seq_along(e1$theta), function(j) {
    tp <- tm <- e1$theta; tp[j] <- tp[j] + h[j]; tm[j] <- tm[j] - h[j]
    (cpo(tp) - cpo(tm)) / (2 * h[j])
  })
  list(m = m, K = K, e1 = e1, trab = trab)
}
