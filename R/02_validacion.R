# =============================================================================
# 02_validacion.R  —  Escalera de validacion (FASE 1)
#
# ESTO TODAVIA NO ES MROZ. Son cuatro regresiones conocidas del libro de
# Wooldridge cuyos coeficientes estan publicados. Si nos salen los numeros
# esperados, sabemos tres cosas:
#   (a) la base esta bien cargada,
#   (b) los paquetes funcionan en las tres computadoras,
#   (c) sabemos escribir la sintaxis de cada estimador.
# Recien entonces tiene sentido pasar a las ecuaciones de HORAS, que son las
# que el paper de Mroz realmente estima.
#
# COMO SE LEE UNA FORMULA EN R
#   y ~ x1 + x2     significa   "y explicada por x1 y x2"
# La tilde ~ se lee "en funcion de". El intercepto se agrega solo.
# =============================================================================

library(AER)              # ivreg() para MC2E, tobit()
library(sampleSelection)  # heckit() y selection()


# =============================================================================
# PELDANIO 1 — MCO de la ecuacion de salarios
# =============================================================================
# Ecuacion de Mincer: log(w) = b0 + b1*educ + b2*exper + b3*exper^2 + u
#
# INTUICION: b1 es el retorno a un anio adicional de escolaridad, en tanto por
# uno. Como la variable dependiente esta en logaritmos, el coeficiente se lee
# directamente como cambio porcentual.
#
# EL PROBLEMA: educ no es exogena. La habilidad no observada (llamemosla A)
# esta en el error u y ademas empuja a estudiar mas. Entonces Cov(educ, u) != 0
# y MCO sobreestima el retorno: le atribuye a la escolaridad parte de lo que en
# realidad es habilidad.
#
# NOTA: se estima sobre 'trab' (las 428 que trabajan) porque para las demas
# no existe lwage. Esto ya es, de hecho, una seleccion muestral -- el peldanio
# 4 la corrige.

mco <- lm(lwage ~ educ + exper + expersq, data = trab)

cat("\n=========== PELDANIO 1: MCO ===========\n")
print(summary(mco))
cat("\nCONTROL -> el coeficiente de educ debe ser aprox. 0.1075\n")
cat("Obtenido:", round(coef(mco)["educ"], 4), "\n")


# =============================================================================
# PELDANIO 2 — MC2E (variables instrumentales)
# =============================================================================
# Instrumentos: educacion de la madre y del padre.
#
# LOS DOS AXIOMAS QUE DEBE CUMPLIR UN INSTRUMENTO:
#   1. RELEVANCIA: Cov(motheduc, educ) != 0. Padres mas educados -> hijas mas
#      educadas. Esto SE CONTRASTA (F de la primera etapa).
#   2. EXCLUSION (o exogeneidad): Cov(motheduc, u) = 0. La educacion de los
#      padres afecta el salario de la hija SOLO a traves de la educacion de
#      ella, no por otro canal. Esto NO se contrasta directamente: se argumenta.
#      Es un supuesto economico, y es discutible (padres educados tambien pasan
#      contactos, habitos, habilidad genetica...).
#
# QUE HACE MC2E: se queda solo con la parte de educ que los instrumentos
# explican -- la parte "limpia", no correlacionada con la habilidad -- y usa
# esa para explicar el salario.
#
# LA SINTAXIS: la barra vertical | separa dos listas.
#   ANTES de la barra: la ecuacion que interesa.
#   DESPUES de la barra: TODOS los instrumentos, incluidas las variables
#   exogenas que ya estaban (exper y expersq se instrumentan a si mismas).
#   Olvidarse de repetirlas es EL error clasico de MC2E en R.

mc2e <- ivreg(lwage ~ educ + exper + expersq |
                motheduc + fatheduc + exper + expersq,
              data = trab)

cat("\n=========== PELDANIO 2: MC2E ===========\n")
print(summary(mc2e, diagnostics = TRUE))
cat("\nCONTROL -> el coeficiente de educ debe caer a aprox. 0.0614\n")
cat("Obtenido:", round(coef(mc2e)["educ"], 4), "\n")
cat("El retorno cae casi a la mitad y el error estandar se duplica.\n")
cat("Ese es el precio de la identificacion por instrumentos.\n")

# En el bloque 'diagnostics = TRUE' aparecen tres contrastes:
#   Weak instruments : F de la primera etapa. Regla practica: F > 10.
#   Wu-Hausman       : contrasta MCO vs MC2E. Si NO se rechaza, educ podria
#                      tratarse como exogena y MCO seria preferible (mas eficiente).
#   Sargan           : solo existe con sobreidentificacion (2 instrumentos para
#                      1 endogena). Contrasta la validez conjunta de las
#                      restricciones de exclusion.


# =============================================================================
# PELDANIO 3 — PROBIT de participacion
# =============================================================================
# Ahora cambiamos de pregunta: ya no "cuanto gana", sino "trabaja o no".
#
# EL MODELO DE VARIABLE LATENTE (es el mismo esqueleto que el Tobit y el
# Heckman, por eso conviene entenderlo bien aqui):
#   Existe una variable no observada  y* = x'b + e , que es el excedente neto
#   de trabajar (salario de mercado menos salario de reserva).
#   Observamos solo el signo:  inlf = 1 si y* > 0 ,  inlf = 0 si y* <= 0.
#   Si e ~ Normal(0,1)  ->  P(inlf=1 | x) = Phi(x'b),  con Phi la normal acumulada.
#
# ADVERTENCIA DE INTERPRETACION: los coeficientes NO son efectos marginales.
# b mide el efecto sobre el indice latente y*, no sobre la probabilidad. El
# efecto sobre la probabilidad depende de donde este parada la persona:
#   dP/dx_j = phi(x'b) * b_j
# Por eso el signo si se lee directo, pero la magnitud no.

probit <- glm(inlf ~ nwifeinc + educ + exper + expersq + age + kidslt6 + kidsge6,
              family = binomial(link = "probit"),
              data = mroz)

cat("\n=========== PELDANIO 3: PROBIT ===========\n")
print(summary(probit))
cat("\nCONTROL -> kidslt6 debe ser aprox. -0.868 y muy significativo\n")
cat("Obtenido:", round(coef(probit)["kidslt6"], 3), "\n")
cat("Un hijo menor de 6 anios es, con diferencia, lo que mas reduce la\n")
cat("participacion. kidsge6 (hijos mayores) no es significativo.\n")

# Efectos marginales promedio, que si son interpretables:
efectos <- mean(dnorm(predict(probit, type = "link"))) * coef(probit)
cat("\nEfecto marginal promedio de kidslt6 sobre la probabilidad:",
    round(efectos["kidslt6"], 3), "\n")


# =============================================================================
# PELDANIO 4 — HECKMAN (correccion de seleccion)
# =============================================================================
# EL PROBLEMA: el salario solo se observa para quien trabaja, y quien trabaja
# no es una muestra aleatoria. Si las mujeres con mayor "gusto por el trabajo"
# (no observado) trabajan mas Y ganan mas, el MCO del peldanio 1 esta sesgado.
#
# LA ESTRUCTURA DE DOS ECUACIONES:
#   Seleccion:  inlf_i = 1[ z_i'g + v_i > 0 ]
#   Resultado:  lwage_i = x_i'b + u_i        (observado solo si inlf_i = 1)
#   Supuesto:   (u, v) ~ Normal bivariada, con correlacion rho.
#
# LA IDEA: al condicionar en trabajar, E[u | inlf=1] ya no es cero, sino
# rho*sigma*lambda(z'g), donde lambda = phi/Phi es la razon inversa de Mills.
# Es decir: el sesgo de seleccion es una VARIABLE OMITIDA. Si la calculamos
# (con el probit del peldanio 3) y la metemos como regresor, el sesgo se va.
#
# EL AXIOMA DE IDENTIFICACION: idealmente z debe tener alguna variable que NO
# este en x (restriccion de exclusion). Aqui: age, kidslt6, kidsge6, nwifeinc
# afectan la decision de participar pero no el salario ofrecido. Sin exclusion
# el modelo se identifica solo por la no-linealidad de lambda, lo que lo hace
# muy fragil.
#
# ERROR QUE HAY QUE EVITAR: a heckit() se le pasa 'mroz' (753 filas, con
# lwage = NA en 325 de ellas), NO 'trab'. Si le pasas la muestra recortada, el
# modelo de seleccion no tiene con que comparar y devuelve numeros plausibles
# pero equivocados.

heck_2p <- heckit(
  selection = inlf ~ nwifeinc + educ + exper + expersq + age + kidslt6 + kidsge6,
  outcome   = lwage ~ educ + exper + expersq,
  data      = mroz,
  method    = "2step"
)

cat("\n=========== PELDANIO 4a: HECKMAN DOS ETAPAS ===========\n")
print(summary(heck_2p))
cat("\nCONTROL -> educ aprox. 0.109 ; el coeficiente de la inversa de Mills\n")
cat("(invMillsRatio) debe salir aprox. 0.032 y NO significativo.\n")

# LECTURA ECONOMICA DEL RESULTADO: la correccion de seleccion casi no mueve
# nada (0.1075 sin corregir -> 0.109 corrigiendo). Esto NO es un fracaso: es
# precisamente el mensaje de Mroz. Los supuestos ESTADISTICOS (normalidad,
# seleccion) importan poco; los supuestos ECONOMICOS (que educ sea exogena,
# que la experiencia sea instrumento valido) lo cambian todo -- compara el
# peldanio 1 con el 2.

# Version de maxima verosimilitud: estima las dos ecuaciones a la vez en lugar
# de en dos pasos. Es mas eficiente SI la normalidad bivariada se cumple, y
# mas fragil si no. Comparar ambas es un ejercicio de robustez del trabajo.
heck_mv <- selection(
  selection = inlf ~ nwifeinc + educ + exper + expersq + age + kidslt6 + kidsge6,
  outcome   = lwage ~ educ + exper + expersq,
  data      = mroz,
  method    = "ml"
)

cat("\n=========== PELDANIO 4b: HECKMAN MAXIMA VEROSIMILITUD ===========\n")
print(summary(heck_mv))
cat("\nAqui rho es la correlacion entre los errores de las dos ecuaciones.\n")
cat("Si rho no es significativo, no hay evidencia de sesgo de seleccion.\n")


# =============================================================================
# RESUMEN DE LA ESCALERA
# =============================================================================
cat("\n\n=========== RESUMEN: retorno a la educacion ===========\n")
# Nota tecnica: en un objeto heckit los nombres de los coeficientes se repiten
# entre las dos ecuaciones (educ esta en ambas), asi que hay que pedir
# explicitamente la parte de resultado. Si tu version del paquete no acepta
# 'part', el numero igual esta impreso en el summary de arriba.
educ_heck <- tryCatch(unname(coef(heck_2p, part = "outcome")["educ"]),
                      error = function(e) NA_real_)

resumen <- data.frame(
  metodo   = c("MCO", "MC2E", "Heckman 2 pasos"),
  educ     = round(c(unname(coef(mco)["educ"]),
                     unname(coef(mc2e)["educ"]),
                     educ_heck), 4),
  esperado = c(0.1075, 0.0614, 0.1091)
)
print(resumen)

cat("\nSi las tres columnas coinciden, la Fase 1 esta cerrada y se puede\n")
cat("pasar a la Fase 2: las ecuaciones de HORAS, que son las de Mroz.\n")
