# =============================================================================
# 06_tablas7y8.R  —  Tablas VII y VIII de Mroz (1987)
#                    Exogeneidad del ingreso no laboral y de los hijos
# Autora: Judith (codigo) · Bloque del integrante 3
# -----------------------------------------------------------------------------
# PREGUNTA: el ingreso no laboral (I) y los hijos (C) se tratan como exogenos
# en la ecuacion de horas. ¿Es valido? Podrian no serlo: una mujer con "gusto
# por el trabajo" puede haber decidido tener menos hijos o casarse con alguien
# que gana distinto. Para contrastarlo, se estima la ecuacion tratando esas
# variables como endogenas (instrumentadas) y como exogenas, y se compara.
#
# Todas las especificaciones usan H3 (variables del esposo) como instrumentos.
# En nuestra base H3 esta INCOMPLETO (falta la educacion de los padres del
# esposo), asi que los coeficientes no salen identicos al paper. Lo que se
# verifica es que las CONCLUSIONES de los contrastes coincidan.
# =============================================================================

library(wooldridge)
data("mroz", package = "wooldridge")
source(here::here("R", "funciones_comunes.R"))
trab <- preparar_trab(mroz)

K  <- c("kidslt6", "kidsge6")   # hijos
AE <- c("age", "educ")

# Cada especificacion: instrumentos excluidos, regresores endogenos, exogenos.
# Una variable que NO es instrumento de si misma se trata como endogena.
espec <- list(
  # ---- Tabla VII (comunes: B, C, F3, H3) ----
  VII_1  = list(c(B, F3, H3),    c("lwage", "nwifeinc"),      c(K, AE)),
  VII_2  = list(c(B, F3, H3),    "lwage",                     c("nwifeinc", K, AE)),
  VII_3  = list(c(B, E, F3, H3), c("lwage", "nwifeinc"),      c(K, AE)),
  VII_4  = list(c(B, E, F3, H3), "lwage",                     c("nwifeinc", K, AE)),
  # ---- Tabla VIII (comunes: B, F3, H3) ----
  VIII_2 = list(c(B, F3, H3),    c("lwage", K),               c("nwifeinc", AE)),
  VIII_4 = list(c(B, F3, H3),    c("lwage", "nwifeinc", K),   AE),
  VIII_6 = list(c(B, E, F3, H3), c("lwage", K),               c("nwifeinc", AE))
)
m <- lapply(espec, function(s) estimar_horas(trab, inst = s[[1]], endog = s[[2]], exog = s[[3]]))

# En la Tabla VIII, las filas 1, 3 y 5 son las mismas especificaciones que
# VII_2, VII_1 y VII_4 (mismas variables exogenas y mismos instrumentos).
m$VIII_1 <- m$VII_2; m$VIII_3 <- m$VII_1; m$VIII_5 <- m$VII_4

fila <- function(mod) { x <- coef_robustos(mod); c(x$coef, se = x$se) }
dif  <- function(a, b) { d <- contraste_diferencia(m[[a]], m[[b]]); c(d$dif, se = d$se) }

armar <- function(nombres, pares, paper, rechaza_paper, titulo, archivo) {
  est <- rbind(t(sapply(m[nombres], fila)),
               t(sapply(pares, function(p) dif(p[1], p[2]))))
  rownames(est) <- c(nombres, names(pares))
  tab <- data.frame(
    especificacion = rownames(est),
    lnw_R  = round(est[, "lwage"]),          lnw_paper  = paper[, 1],
    se_R   = round(est[, "se.lwage"]),       se_paper   = paper[, 5],
    nwif_R = round(est[, "nwifeinc"], 1),    nwif_paper = paper[, 2],
    jov_R  = round(est[, "kidslt6"]),        jov_paper  = paper[, 3],
    may_R  = round(est[, "kidsge6"]),        may_paper  = paper[, 4],
    row.names = NULL
  )
  cat("\n===========", titulo, "===========\n")
  cat("(H3 incompleto: se esperan diferencias de magnitud con el paper)\n\n")
  print(tab)

  nd <- length(pares)
  z <- tail(est[, "lwage"] / est[, "se.lwage"], nd)
  decision <- data.frame(contraste = names(pares), z_lnw = round(z, 2),
                         rechaza_R = abs(z) > 1.96, rechaza_paper = rechaza_paper,
                         row.names = NULL)
  cat("\nContrastes (|z| > 1.96 -> se rechaza al 5%):\n")
  print(decision)

  write.csv(tab, here::here("salidas", archivo), row.names = FALSE)
  cat("Guardado en salidas/", archivo, "\n", sep = "")
  stopifnot(all(decision$rechaza_R == decision$rechaza_paper))
  cat("OK: las conclusiones de los", nd, "contrastes coinciden con el paper.\n")
}

dir.create(here::here("salidas"), showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Tabla VII: exogeneidad del ingreso no laboral (p. 777)
# -----------------------------------------------------------------------------
armar(
  nombres = c("VII_1", "VII_2", "VII_3", "VII_4"),
  pares   = list("2 - 1" = c("VII_2", "VII_1"),
                 "4 - 3" = c("VII_4", "VII_3"),
                 "4 - 1" = c("VII_4", "VII_1")),
  paper = rbind(
    c( -29,  -3.1, -341, -117,  174, 8.7, 130, 31),
    c( -29,  -4.2, -334, -116,  174, 3.3, 130, 30),
    c( 405, -18.6, -316,  -84,  174, 8.1, 318, 33),
    c( 482,  -5.8, -300,  -93,  171, 3.4, 138, 33),
    c(   0,  -1.1,    2,    1,    5, 8.1,   6,  8),
    c(  78,  12.7,   16,   10,   57, 7.7,  21,  8),
    c( 512,   2.8,  -42,   23,  148, 8.3,  56, 17)),
  rechaza_paper = c(FALSE, FALSE, TRUE),
  titulo = "TABLA VII: R vs. paper", archivo = "tabla7.csv"
)
cat("Lectura: no se rechaza que el ingreso no laboral sea exogeno (2-1, 4-3);\n")
cat("se rechaza la exogeneidad conjunta de ingreso y experiencia (4-1).\n")

# -----------------------------------------------------------------------------
# Tabla VIII: exogeneidad de los hijos y del ingreso no laboral (p. 777)
# -----------------------------------------------------------------------------
armar(
  nombres = paste0("VIII_", 1:6),
  pares   = list("1 - 2" = c("VIII_1", "VIII_2"),
                 "3 - 4" = c("VIII_3", "VIII_4"),
                 "5 - 6" = c("VIII_5", "VIII_6")),
  paper = rbind(
    c( -30, -4.2, -344, -116,  174, 3.3, 130, 30),
    c( -52, -5.4, -298,   14,  189, 3.6, 380, 85),
    c( -29, -3.1, -341, -117,  174, 8.7, 130, 31),
    c( -53, -3.7, -278,   14,  190, 9.2, 385, 84),
    c( 482, -5.8, -301,  -93,  171, 3.4, 138, 33),
    c( 442, -5.3, -966, -176,  175, 3.7, 343, 87),
    c(  22,  1.7,  -45, -129,   57, 1.2, 344, 79),
    c(  24,  0.7,  -63, -131,   58, 2.3, 351, 78),
    c(  40, -0.6,  663,   83,   82, 1.8, 357, 81)),
  rechaza_paper = c(FALSE, FALSE, FALSE),
  titulo = "TABLA VIII: R vs. paper", archivo = "tabla8.csv"
)
cat("Lectura: en ningun caso se rechaza que los hijos sean exogenos.\n")
