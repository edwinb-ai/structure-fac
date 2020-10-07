using LiquidsStructure
using Plots
using DelimitedFiles
using FFTW

gr(size=(1280, 1024))

"""
Esta función calcula el valor del factor de estructura siguiendo la definición
de la integral, como se ve en https://en.wikipedia.org/wiki/Radial_distribution_function#The_structure_factor

Lo que se hace es tomar cada uno de los valores de onda proveniente de `k` y se
evalúa la integral como una suma de Riemann.
"""
function integral(x, y, η, k)
    dens = 6 * η / π
    dr = x[2]
    n = size(x, 1)

    total_sum = zeros(n)
    sum_loop = 0

    for j = axes(k, 1)
        sum_loop = 0
        @inbounds for i = axes(x, 1)
            sum_loop += x[i] * sin(k[j] * x[i]) * (y[i] - 1) * dr
        end
        total_sum[j] = sum_loop / k[j]
        total_sum[j] *= 2 * π * dens
    end

    return 1 .+ total_sum
end

function sq_fft(r, hr, η, k)
    dens = 6 * η / π
    hk = FFTW.r2r(r .* hr, FFTW.RODFT00)
    normalization = 2 * π * r[1] ./ k
    Sk = 1.0 .+ dens .* (hk .* normalization)

    return Sk
end

# Se crea la función del factor de estructura
η = 0.4
S = StructureFactor(HardSpheres(η), PercusYevick)

# Se calcula usando la definición y la implementación
# ingenua usando la RDF
data = readdlm("gr_BD_040_ih.csv", ',', Float64, '\n')
# Para usar un número para de datos, se recortan los datos originales
new_data = data[2:end-1, :]
display(size(new_data))
# Se generan los datos del número de onda
dk = π / new_data[end, 1]
k = range(1, 510; step=1) .* dk
s_integral = integral(new_data[:, 1], new_data[:, 2], η, k)

# Cálculo del S(q) con la FFT
# La tercera columna del conjunto de datos corresonde a la h(r)
s_fft = sq_fft(new_data[:, 1], new_data[:, 3], η, k)
display(s_fft)

# Se calcula el verdadero factor de estructura
s_true = S.(k)

# Se grafica el factor de estructura
plot(k, s_true, lab="Cerradura de Percus-Yevick")
plot!(k, s_integral, lab="Integral de Riemann de g(r)")
plot!(k, s_fft, lab="FFT de g(r)")
xlims!(0, 50)
xaxis!("k")
yaxis!("S(k)")
savefig("factor_estructura")
gui()
