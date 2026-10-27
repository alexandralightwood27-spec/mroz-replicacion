# Ficha del bloque [1 a 5] — Réplica de Mroz (1987)

**Autor/a:** [Nombre]
**Integrante y bloque:** [Integrante 1 · Pregunta, teoría y datos / Integrante 2 · Variables instrumentales y MC2E / Integrante 3 · Contrastes de especificación / Integrante 4 · Selección multietapa / Integrante 5 · Tobit y el experimento]
**Lectura usada:** [la guía y la parte de la traducción que indica la hoja de ruta para tu bloque]
**Scripts del bloque:** [por ejemplo R/05_tabla6.R y R/06_tablas7y8.R]

> **Cómo usar esta plantilla:** copia todo este texto en un archivo nuevo
> llamado `ficha_bloque1.md`, `ficha_bloque2.md`, ... o `ficha_bloque5.md`
> según tu número, borra las indicaciones entre corchetes y escribe con tus
> propias palabras. Extensión sugerida: dos a tres páginas. No copies la guía:
> explícala como se la explicarías a un compañero que no la leyó.
>
> **Formato básico:** `#` y `##` son títulos, `-` al inicio de línea es una
> viñeta, `**texto**` es negrita y el código va entre tres acentos graves
> (```` ``` ````). Nada más hace falta.

---

## 1. La pregunta

[¿Qué problema resuelve tu bloque? Una o dos oraciones.
Ejemplo: "¿El salario es endógeno en la ecuación de horas, y qué instrumentos sirven para corregirlo?"]

## 2. El método

### En palabras simples

[Explica la idea como si hablaras con alguien de primer ciclo. Usa un ejemplo o una analogía.]

### En términos técnicos

[La ecuación o el procedimiento, paso a paso.]

### Supuestos que necesita

- [Supuesto 1, y qué significa económicamente]
- [Supuesto 2]
- [¿Cuál se puede comprobar con los datos y cuál solo se puede argumentar?]

## 3. El código de mi bloque

[Para cada script de tu bloque, explica qué hace por partes. No hace falta copiar
el script entero: elige los fragmentos importantes, pégalos y explica cada uno.]

### Qué hace el script en general

[Dos o tres oraciones: qué datos usa, qué estima y qué archivo guarda en `salidas/`.]

### Los fragmentos clave

```r
# [pega aquí un fragmento corto del script, de 3 a 10 líneas]
```

[Qué hace cada línea y por qué se hace así. Por ejemplo: qué función de R se usa,
qué significan sus argumentos y a qué paso del método corresponde.]

[Repite con dos o tres fragmentos más.]

### Cómo se comprueba que el resultado es correcto

[Explica el `stopifnot()` de tu script: qué compara contra el paper y con qué tolerancia,
y qué mensaje sale en la consola cuando todo está bien.]

## 4. Resultados

[Copia la tabla de tu bloque desde el archivo de `salidas/` (por ejemplo `salidas/tabla6.csv`),
con nuestros números al lado de los del paper. No escribas los números a mano.]

| Especificación | Nuestro resultado | Paper |
|---|---|---|
| [ ] | [ ] | [ ] |

**Qué significan:** [Interpretación económica. Traduce el coeficiente del salario a horas por dólar
(dividir entre 4,5) o a elasticidad (dividir entre 1.302,9) cuando ayude.]

**Qué contrastes se hicieron y qué concluyen:** [z, si se rechaza o no, y qué supuesto se pone a prueba.]

## 5. Diferencias con el paper

[Revisa `documentacion/discrepancias.csv` y explica las que corresponden a tu bloque:
qué difiere, cuánto y por qué. Si alguna es una errata del paper, dilo.]

- Registro [n]: [explicación]

## 6. Mi conclusión en la exposición

[La frase con la que cierras tu turno y cómo conecta con el mensaje del paper:
¿tu bloque muestra que los supuestos económicos pesan más que los estadísticos? ¿Cómo?
¿Qué deja planteado para el integrante que viene después?]

## 7. Tres preguntas que podría hacer el profesor

1. **[Pregunta]**
   [Respuesta en dos o tres líneas]
2. **[Pregunta sobre tu código]**
   [Respuesta]
3. **[Pregunta]**
   [Respuesta]

---

*Para subirla: en github.com, repositorio mroz-replicacion, carpeta `documentacion`,
botón **Add file → Create new file**, pega el texto, ponle el nombre y presiona
**Commit changes**. Queda registrada con tu nombre.*
