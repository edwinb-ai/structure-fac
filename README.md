# Métodos numéricos para el Factor de Estructura

El factor de estructura de un líquido, sólido o gas tiene propiedades interesantes
numéricamente hablando.

En este repositorio se trata de explorar cuáles son los mejores métodos numéricos
para estudiar esta cantidad fundamental que proviene de la teoría avanzada
de la mecánica estadística y la teoría de líquidos.

Hasta ahora, se han empleado los siguientes métodos:

- Cerradura de Percus-Yevick (esta es la realidad base o _ground truth_)
- Transformada de Fourier (TF) de la función de distribución radial (RDF).

Se proponen los siguientes métodos:

- Integración de Romberg
- Transformada rápida de Fourier (FFT)

## Reproducibilidad y resultados

Para reproducibilidad, se emplea el paquete [`LiquidsStructure.jl`](https://github.com/LANIMFE/LiquidsStructure.jl) del LANIMFE
para los valores reales y cerraduras teóricas.

También se comparten los datos de la RDF para un sistema de esferas duras con fracción de volumen de 0.4, simulado
mediante dinámica browniana con 512 partículas.

Los resultados se presentan en forma de imágenes, pero se está trabajando en un tipo de exposición más adecuada en formato de libretas
o HTML.

## Resultados

- La integral de Riemann numérica (método numérico ingenuo) de la RDF muestra oscilaciones importantes en momento inicial de la integración o aplicación del método. También se puede ver que en el pico más prominente del factor de estructura no se llega al mismo valor empleando este método. Es muy posible que sea numéricamente inestable para valores pequeños o en el caso de números cercanos entre sí. Esto se puede argumentar porque se llega al límite físico real donde el factor de estructura tiende a uno conforme k tiende al infinito.
- No es de sorprenderse que la FFT da el mismo resultado que la integración de Riemann. Esto se debe a que ambas formulaciones son la misma, lo único que cambia es la complejidad algorítmica, mientras que la FFT es O(Nlog(N)), la integración de Riemann es O(N^2).
- Es importante notar que los dos métodos numéricos dan el mismo resultados, pero subestiman el pico principal que se espararía de la solución teórica. Es posible que con mejor estadística (diferentes RDF diferentes provenientes de distintas simulaciones) mejore este resultado considerablemente.