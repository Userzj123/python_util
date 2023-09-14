using LinearAlgebra
using FFTW

x = LinRange(-pi, pi, 16)[1:end-1]
y = LinRange(-1, 1, 16)

xx = x .+ y'.*0

xx2 = x' .+ y.*0

fft(sin.(xx), 2)