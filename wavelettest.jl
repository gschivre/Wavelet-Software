# all this code is obtian by converting the matlab code provided at http://eeweb.poly.edu/iselesni/WaveletSoftware/index.html
using DSP, Plots, Makie
include("waveletsoft.jl")

af, sf = filtdesign("DWT")

x = randn(64)
lo, hi = afb(x, af)
y = sfb(lo, hi, sf)
err = x .- y
maximum(abs.(err))

w = dwt(x, af, 3)
y = idwt(w, sf)
err = x .- y
maximum(abs.(err))

x = zeros(64)
w = dwt(x, af, 3)
w[3][5] = 1.0
y = idwt(w, sf)
Plots.plot(range(0.0, 1.0; length = 64), y; xlabel = "t", ylabel = "Ψ(t)", legend = false)

x = randn(128, 64)
lo, hi = afb(x, af)
y = sfb(lo, hi, sf)
err = x .- y
maximum(abs.(err[:]))

w = dwt(x, af, 3)
y = idwt(w, sf)
err = x .- y
maximum(abs.(err[:]))

J = 5                    
L = 3 * 2 ^ (J + 1)
N = Int(L / 2 ^ J)
x = zeros(L, 3 * L)
w = dwt(x, af, J)
w[J][1][div(N, 2), div(N, 2) + 0 * N] = 1.0
w[J][2][div(N, 2), div(N, 2) + 1 * N] = 1.0
w[J][3][div(N, 2), div(N, 2) + 2 * N] = 1.0
y = idwt(w, sf)
Plots.heatmap(range(0.0, 1.0; length = 3 * L), range(0.0, 1.0; length = L), y; color = :grays)

x = randn(32, 64, 16)
lo, hi = afb(x, af)
y = sfb(lo, hi, sf)
err = x .- y
maximum(abs.(err[:]))

x = randn(128, 64, 64)
w = dwt(x, af, 3)
y = idwt(w, sf)
err = x .- y
maximum(abs.(err[:]))

J = 4
L = 3 * 2 ^ (J + 1)
N = Int(L / 2 ^ J)
x = zeros(L, L, L)
w = dwt(x, af, J)
w[J][7][div(N, 2), div(N, 2), div(N, 2)] = 1.0
y = idwt(w, sf)
thy1 = zeros(size(y))
thy1[y .> 0.017] .= 1.0
thy2 = zeros(size(y))
thy2[y .< -0.017] .= 1.0
Makie.volume(range(0.0, 1.0; length = L),
			range(0.0, 1.0; length = L),
			range(0.0, 1.0; length = L),
			thy1;
			algorithm = :iso,
			isovalue = 1.0)
Makie.volume!(range(0.0, 1.0; length = L),
			range(0.0, 1.0; length = L),
			range(0.0, 1.0; length = L),
			thy2;
			algorithm = :iso,
			isovalue = 1.0)

Faf, Fsf, af, sf = filtdesign("CWT")

x = randn(512)
w, iw = dtcwt(x, Faf, af, 4)
y = idtcwt(w, iw, Fsf, sf)
err = x .- y
maximum(abs.(err))

x = zeros(256)
w, iw = dtcwt(x, Faf, af, 5)
w[5][4] = 1.0
y1 = idtcwt(w, iw, Fsf, sf)
w, iw = dtcwt(x, Faf, af, 5)
iw[5][4] = 1.0
y2 = idtcwt(w, iw, Fsf, sf)
Plots.plot(range(0.0, 1.0; length = 256), y1; xlabel = "t", ylabel = "Ψ(t)", legend = false)
Plots.plot!(range(0.0, 1.0; length = 256), y2)
Plots.plot!(range(0.0, 1.0; length = 256), sqrt.(y1 .^ 2 .+ y2 .^ 2))

x = randn(256, 256)
w, iw = dtcwt(x, Faf, af, 5)
y = idtcwt(w, iw, Fsf, sf)
err = x .- y
maximum(abs.(err[:]))
w = dtrwt(x, Faf, af, 5)
y = idtrwt(w, Fsf, sf)
err = x .- y
maximum(abs.(err[:]))

J = 4
L = 3 * 2 ^ (J + 1)
N = Int(L / 2 ^ J)
x = zeros(2 * L, 6 * L)
w, iw = dtcwt(x, Faf, af, J)
w[J][5][div(N, 2), div(N, 2) + 0 * N] = 1.0
w[J][3][div(N, 2), div(N, 2) + 1 * N] = 1.0
w[J][4][div(N, 2), div(N, 2) + 2 * N] = 1.0
w[J][1][div(N, 2), div(N, 2) + 3 * N] = 1.0
w[J][6][div(N, 2), div(N, 2) + 4 * N] = 1.0
w[J][2][div(N, 2), div(N, 2) + 5 * N] = 1.0
iw[J][5][div(N, 2) + N, div(N, 2) + 0 * N] = 1.0
iw[J][3][div(N, 2) + N, div(N, 2) + 1 * N] = 1.0
iw[J][4][div(N, 2) + N, div(N, 2) + 2 * N] = 1.0
iw[J][1][div(N, 2) + N, div(N, 2) + 3 * N] = 1.0
iw[J][6][div(N, 2) + N, div(N, 2) + 4 * N] = 1.0
iw[J][2][div(N, 2) + N, div(N, 2) + 5 * N] = 1.0
y = idtcwt(w, iw, Fsf, sf)
y = [y; sqrt.(y[1:L, :] .^ 2 .+ y[L .+ (1:L), :] .^ 2)]
Plots.heatmap(range(0.0, 1.0; length = 6 * L), range(0.0, 1.0; length = 3 * L), y[end:-1:1, :]; color = :grays)

J = 4
L = 3 * 2 ^ (J + 1)
N = Int(L / 2 ^ J)
x = zeros(2 * L, 3 * L)
w = dtrwt(x, Faf, af, J)
w[J][1][div(N, 2), div(N, 2) + 0 * N] = 1.0
w[J][2][div(N, 2), div(N, 2) + 1 * N] = 1.0
w[J][3][div(N, 2), div(N, 2) + 2 * N] = 1.0
w[J][4][div(N, 2) + N, div(N, 2) + 0 * N] = 1.0
w[J][5][div(N, 2) + N, div(N, 2) + 1 * N] = 1.0
w[J][6][div(N, 2) + N, div(N, 2) + 2 * N] = 1.0
y = idtrwt(w, Fsf, sf)
Plots.heatmap(range(0.0, 1.0; length = 3 * L), range(0.0, 1.0; length = 2 * L), y[end:-1:1, :]; color = :grays)

x = randn(64, 64, 64)
J = 3
w, iw = dtcwt(x, Faf, af, J)
y = idtcwt(w, iw, Fsf, sf)
err = x .- y
maximum(abs.(err[:]))
w = dtrwt(x, Faf, af, J)
y = idtrwt(w, Fsf, sf)
err = x .- y
maximum(abs.(err[:]))