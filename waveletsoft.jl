# all this code is obtian by converting the matlab code provided at http://eeweb.poly.edu/iselesni/WaveletSoftware/index.html

# upsample an array by inserting n - 1 '0' in between each element of a
function upsample(a, n)
    out = zeros(eltype(a), n * length(a))
    out[1:n:end] .= a
    return out
end

# downsample an array a by a factor n
function downsample(a, n)
    return a[1:n:end]
end

# upsample by a factor p then aply filter f and downsample by a factor q
function upfirdn(x, f, p = 1, q = 1)
    return downsample(filt(f, upsample(x, p)), q)
end

# plus / minus function for 2D
function pm(a, b)
	return ((a .+ b) ./ sqrt(2), (a .- b) ./ sqrt(2))
end

# plus / minus function for 3D
function pm4(a, b, c, d)
	return ((a .- b .- c .- d) ./ 2, (a .- b .+ c .+ d) ./ 2,
			(a .+ b .- c .+ d) ./ 2, (a .+ b .+ c .- d) ./ 2)	
end

# inverse of the plus / minus function for 3D
function pm4inv(a, b, c, d)
	return ((a .+ b .+ c .+ d) ./ 2, (c .+ d .- a .- b) ./ 2,
			(b .- a .+ d .- c) ./ 2, (b .- a .+ c .- d) ./ 2)
end

# filter design
function filtdesign(type = "DWT")
	if type == "DWT"
		af = [0.0 -0.01122679215254;
    		0.0 0.01122679215254;
    		-0.08838834764832 0.08838834764832;
    		0.08838834764832 0.08838834764832;
    		0.69587998903400 -0.69587998903400;
    		0.69587998903400 0.69587998903400;
    		0.08838834764832 -0.08838834764832;
    		-0.08838834764832 -0.08838834764832;
    		0.01122679215254 0.0;
    		0.01122679215254 0.0]
		sf = af[end:-1:1, :]
		return (af, sf)
	elseif type == "CWT"
		Faf = vcat([[0.0 0.0;
					-0.08838834764832 -0.01122679215254;
					0.08838834764832 0.01122679215254;
					0.69587998903400 0.08838834764832;
					0.69587998903400 0.08838834764832;
					0.08838834764832 -0.69587998903400;
					-0.08838834764832 0.69587998903400;
					0.01122679215254 -0.08838834764832;
					0.01122679215254 -0.08838834764832
					0.0 0.0]],
				[[0.01122679215254 0.0;
				0.01122679215254 0.0;
				-0.08838834764832 -0.08838834764832;
				0.08838834764832 -0.08838834764832;
				0.69587998903400 0.69587998903400;
				0.69587998903400 -0.69587998903400;
				0.08838834764832 0.08838834764832;
				-0.08838834764832 0.08838834764832;
				0.0 0.01122679215254;
				0.0 -0.01122679215254]])
		Fsf = vcat([Faf[1][end:-1:1, :]], [Faf[2][end:-1:1, :]])

		af = vcat([[0.03516384000000 0.0;
					0.0 0.0;
					-0.08832942000000 -0.11430184000000;
					0.23389032000000 0.0;
					0.76027237000000 0.58751830000000;
					0.58751830000000 -0.76027237000000;
					0.0 0.23389032000000;
					-0.11430184000000 0.08832942000000;
					0.0 0.0;
					0.0 -0.03516384000000]],
				[[0.0 -0.03516384000000;
				0.0 0.0;
				-0.11430184000000 0.08832942000000;
				0.0 0.23389032000000;
				0.58751830000000 -0.76027237000000;
				0.76027237000000 0.58751830000000;
				0.23389032000000 0.0;
				-0.08832942000000 -0.11430184000000;
				0.0 0.0;
				0.03516384000000 0.0]])
		sf = vcat([af[1][end:-1:1, :]], [af[2][end:-1:1, :]])
		return (Faf, Fsf, af, sf)
	end
end

# analysis filter bank
begin
	# 1D
    function afb(x::Array{T, 1}, af) where {T <: Number}
        N = length(x)
        L = div(size(af, 1), 2)

    	# lowpass filter
        lo = upfirdn(vcat(circshift(x, -L), zeros(2 * L)), af[:, 1], 1, 2)
        lo[1:L] .= lo[div(N, 2) .+ (1:L)] .+ lo[1:L]

    	# highpass filter
        hi = upfirdn(vcat(circshift(x, -L), zeros(2 * L)), af[:, 2], 1, 2)
        hi[1:L] = hi[div(N, 2) .+ (1:L)] .+ hi[1:L]

        return (lo[1:div(N, 2)], hi[1:div(N, 2)])
	end
	
	# 2D
	function afb(x::Array{T, 2}, af1, af2 = af1) where {T <: Number}
		n, m = size(x)
		hi = Array{eltype(x), ndims(x)}[]

		# filter along columns
		temp = [afb(x[:, i], af1) for i in 1:m]
		L = reduce(hcat, [temp[i][1] for i in 1:m])
		H = reduce(hcat, [temp[i][2] for i in 1:m])

		# filter along rows
		temp = [afb(L[i, :], af2) for i in 1:div(n, 2)]
		lo = reduce(vcat, [permutedims(temp[i][1]) for i in 1:div(n, 2)]) # LL subband
		push!(hi, reduce(vcat, [permutedims(temp[i][2]) for i in 1:div(n, 2)])) # LH subband
		temp = [afb(H[i, :], af2) for i in 1:div(n, 2)]
		push!(hi, reduce(vcat, [permutedims(temp[i][1]) for i in 1:div(n, 2)])) # HL subband
		push!(hi, reduce(vcat, [permutedims(temp[i][2]) for i in 1:div(n, 2)])) # HH subband
		return (lo, hi)
	end

	# 3D
	function afb(x::Array{T, 3}, af1, af2 = af1, af3 = af1) where {T <: Number}
		n, m, z = size(x)
		hi = Array{eltype(x), ndims(x)}[]

		# filter along the first dimension
		temp = [afb(x[:, i, j], af1) for i in 1:m, j in 1:z]
		L = reduce((a, b) -> cat(a, b; dims = 3), [reduce(hcat, [temp[i, j][1] for i in 1:m]) for j in 1:z])
		H = reduce((a, b) -> cat(a, b; dims = 3), [reduce(hcat, [temp[i, j][2] for i in 1:m]) for j in 1:z])

		# filter along the second dimension
		temp = [afb(L[i, :, j], af2) for i in 1:div(n, 2), j in 1:z]
		LL = reduce((a, b) -> cat(a, b; dims = 3), [reduce(vcat, [permutedims(temp[i, j][1]) for i in 1:div(n, 2)]) for j in 1:z])
		LH = reduce((a, b) -> cat(a, b; dims = 3), [reduce(vcat, [permutedims(temp[i, j][2]) for i in 1:div(n, 2)]) for j in 1:z])
		temp = [afb(H[i, :, j], af2) for i in 1:div(n, 2), j in 1:z]
		HL = reduce((a, b) -> cat(a, b; dims = 3), [reduce(vcat, [permutedims(temp[i, j][1]) for i in 1:div(n, 2)]) for j in 1:z])
		HH = reduce((a, b) -> cat(a, b; dims = 3), [reduce(vcat, [permutedims(temp[i, j][2]) for i in 1:div(n, 2)]) for j in 1:z])

		# filter along the third dimension
		temp = [afb(LL[i, j, :], af3) for i in 1:div(n, 2), j in 1:div(m, 2)]
		lo = reduce(hcat, [reduce(vcat, [reshape(temp[i, j][1], 1, 1, div(z, 2)) for i in 1:div(n, 2)]) for j in 1:div(m, 2)]) # LLL subband
		push!(hi, reduce(hcat, [reduce(vcat, [reshape(temp[i, j][2], 1, 1, div(z, 2)) for i in 1:div(n, 2)]) for j in 1:div(m, 2)])) # LLH subband
		temp = [afb(LH[i, j, :], af3) for i in 1:div(n, 2), j in 1:div(m, 2)]
		push!(hi, reduce(hcat, [reduce(vcat, [reshape(temp[i, j][1], 1, 1, div(z, 2)) for i in 1:div(n, 2)]) for j in 1:div(m, 2)])) # LHL subband
		push!(hi, reduce(hcat, [reduce(vcat, [reshape(temp[i, j][2], 1, 1, div(z, 2)) for i in 1:div(n, 2)]) for j in 1:div(m, 2)])) # LHH subband
		temp = [afb(HL[i, j, :], af3) for i in 1:div(n, 2), j in 1:div(m, 2)]
		push!(hi, reduce(hcat, [reduce(vcat, [reshape(temp[i, j][1], 1, 1, div(z, 2)) for i in 1:div(n, 2)]) for j in 1:div(m, 2)])) # HLL subband
		push!(hi, reduce(hcat, [reduce(vcat, [reshape(temp[i, j][2], 1, 1, div(z, 2)) for i in 1:div(n, 2)]) for j in 1:div(m, 2)])) # HLH subband
		temp = [afb(HH[i, j, :], af3) for i in 1:div(n, 2), j in 1:div(m, 2)]
		push!(hi, reduce(hcat, [reduce(vcat, [reshape(temp[i, j][1], 1, 1, div(z, 2)) for i in 1:div(n, 2)]) for j in 1:div(m, 2)])) # HHL subband
		push!(hi, reduce(hcat, [reduce(vcat, [reshape(temp[i, j][2], 1, 1, div(z, 2)) for i in 1:div(n, 2)]) for j in 1:div(m, 2)])) # HHH subband
		return (lo, hi)
	end
end

# synthesis filter bank
begin
	# 1D
	function sfb(lo::Array{T, 1}, hi::Array{T, 1}, sf) where {T <: Number}
    	N = 2 * length(lo)
    	L = size(sf, 1)
    	y = upfirdn(vcat(lo, zeros(div(L - 1, 2))), sf[:, 1], 2, 1) .+ upfirdn(vcat(hi, zeros(div(L - 1, 2))), sf[:, 2], 2, 1)
    	y[1:(L - 2)] .= y[1:(L - 2)] .+ y[N .+ (1:(L - 2))]
    	return circshift(y[1:N], Int(1 - L / 2))
	end

	# 2D
	function sfb(lo::Array{T, 2}, hi::Array{Array{T, 2}, 1}, sf1, sf2 = sf1) where {T <: Number}
		n, m = size(lo)

		# filter along rows
		L = reduce(vcat, [permutedims(sfb(lo[i, :], hi[1][i, :], sf2)) for i in 1:n])
		H = reduce(vcat, [permutedims(sfb(hi[2][i, :], hi[3][i, :], sf2)) for i in 1:n])

		# filter along columns
		return reduce(hcat, [sfb(L[:, i], H[:, i], sf1) for i in 1:(2 * m)])
	end

	# 3D
	function sfb(lo::Array{T, 3}, hi::Array{Array{T, 3}, 1}, sf1, sf2 = sf1, sf3 = sf1) where {T <: Number}
		n, m, z = size(lo)

		# filter along the third dimension
		LL = reduce(hcat, [reduce(vcat, [reshape(sfb(lo[i, j, :], hi[1][i, j, :], sf3), 1, 1, 2 * z) for i in 1:n]) for j in 1:m])
		LH = reduce(hcat, [reduce(vcat, [reshape(sfb(hi[2][i, j, :], hi[3][i, j, :], sf3), 1, 1, 2 * z) for i in 1:n]) for j in 1:m])
		HL = reduce(hcat, [reduce(vcat, [reshape(sfb(hi[4][i, j, :], hi[5][i, j, :], sf3), 1, 1, 2 * z) for i in 1:n]) for j in 1:m])
		HH = reduce(hcat, [reduce(vcat, [reshape(sfb(hi[6][i, j, :], hi[7][i, j, :], sf3), 1, 1, 2 * z) for i in 1:n]) for j in 1:m])

		# filter along the second dimension
		L = reduce((a, b) -> cat(a, b; dims = 3), [reduce(vcat, [permutedims(sfb(LL[i, :, j], LH[i, :, j], sf2)) for i in 1:n]) for j in 1:(2 * z)])
		H = reduce((a, b) -> cat(a, b; dims = 3), [reduce(vcat, [permutedims(sfb(HL[i, :, j], HH[i, :, j], sf2)) for i in 1:n]) for j in 1:(2 * z)])

		# filter along the first dimension
		return reduce((a, b) -> cat(a, b; dims = 3), [reduce(hcat, [sfb(L[:, i, j], H[:, i, j], sf1) for i in 1:(2 * m)]) for j in 1:(2 * z)])
	end
end

# discret wavelet transform
function dwt(x, af, J = floor(Int, log2(minimum(size(x)))))
	x_ = deepcopy(x)
	if ndims(x) > 1 # nD DWT
    	w = Array{Array{Float64, ndims(x)}}[]
    	for k in 1:J
        	x_, hi = afb(x_, af)
        	push!(w, hi)
    	end
		push!(w, [x_])
	else # 1D DWT
		w = Array{Float64, 1}[]
		for k in 1:J
        	x_, hi = afb(x_, af)
        	push!(w, hi)
    	end
		push!(w, x_)
	end
    return w
end

# inverse discret wavelet transform
function idwt(w, sf)
	J = length(w) - 1
	if ndims(w[1][1]) > 1 # nD inverse DWT
		y = w[J + 1][1]
	else # 1D inverse DWT
    	y = w[J + 1]
	end
	for k in J:-1:1
		y = sfb(y, w[k], sf)
	end
    return y
end

# dual-tree complex wavelet transform
function dtcwt(x, Faf, af, J = floor(Int, log2(minimum(size(x)))))
	if ndims(x) > 1 # nD dtℂWT
		w = Array{Array{Float64, ndims(x)}}[]
		iw = Array{Array{Float64, ndims(x)}}[]
		(ndims(x) == 2 ? p = 2 : p = 1)
		for i in Iterators.product(collect([[1; 2] for _ in 1:ndims(x)])...)
			lo, hi = afb(x ./ sqrt(2 ^ ndims(x)), collect([Faf[j] for j in i[end:-1:1]])...)
			if all(i[1:length(i) .!= p] .== 1)
				(i[p] == 1 ? push!(w, hi) : push!(iw, hi))
			else
				(i[p] == 1 ? append!(w[1], hi) : append!(iw[1], hi))
			end
			for k in 2:J
				lo, hi = afb(lo, collect([af[j] for j in i[end:-1:1]])...)
				if all(i[1:length(i) .!= p] .== 1)
					(i[p] == 1 ? push!(w, hi) : push!(iw, hi))
				else
					(i[p] == 1 ? append!(w[k], hi) : append!(iw[k], hi))
				end
			end
			if all(i[1:length(i) .!= p] .== 1)
				(i[p] == 1 ? push!(w, [lo]) : push!(iw, [lo]))
			else
				(i[p] == 1 ? append!(w[J + 1], [lo]) : append!(iw[J + 1], [lo]))
			end
		end
		if ndims(x) == 2
			for k in 1:J
				for i in 1:3
					w[k][i], iw[k][3 + i] = pm(w[k][i], iw[k][3 + i])
					w[k][3 + i], iw[k][i] = pm(w[k][3 + i], iw[k][i])
				end
			end
		else
			for k in 1:J
				for i in 1:7
					w[k][i], w[k][21 + i], iw[k][14 + i], iw[k][7 + i] = pm4(w[k][i], w[k][21 + i], iw[k][14 + i], iw[k][7 + i])
					iw[k][21 + i], iw[k][i], w[k][7 + i], w[k][14 + i] = pm4(iw[k][21 + i], iw[k][i], w[k][7 + i], w[k][14 + i])
				end
			end
		end
	else # 1D dtℂWT
		# tree 1
		w = Array{Float64, 1}[]
		x1, hi = afb(x ./ sqrt(2), Faf[1])
    	push!(w, hi)
		for k in 2:J
        	x1, hi = afb(x1, af[1])
        	push!(w, hi)
    	end
		push!(w, x1)

		# tree 2
		iw = Array{Float64, 1}[]
		x2, hi = afb(x ./ sqrt(2), Faf[2])
    	push!(iw, hi)
		for k in 2:J
        	x2, hi = afb(x2, af[2])
        	push!(iw, hi)
    	end
		push!(iw, x2)
	end
	return (w, iw)
end

# inverse dual-tree complex wavelet transform
function idtcwt(w, iw, Fsf, sf)
	J = length(w) - 1
	w_ = deepcopy(w)
	iw_ = deepcopy(iw)
	if ndims(w[1][1]) > 1 # nD inverse dtℂWT
		if ndims(w[1][1]) == 2
			for k in 1:J
				for i in 1:3
					w_[k][i], iw_[k][3 + i] = pm(w_[k][i], iw_[k][3 + i])
					w_[k][3 + i], iw_[k][i] = pm(w_[k][3 + i], iw_[k][i])
				end
			end
		else
			for k in 1:J
				for i in 1:7
					w_[k][i], w_[k][21 + i], iw_[k][14 + i], iw_[k][7 + i] = pm4inv(w_[k][i], w_[k][21 + i], iw_[k][14 + i], iw_[k][7 + i])
					iw_[k][21 + i], iw_[k][i], w_[k][7 + i], w_[k][14 + i] = pm4inv(iw_[k][21 + i], iw_[k][i], w_[k][7 + i], w_[k][14 + i])
				end
			end
		end
		(ndims(w[1][1]) == 2 ? p = 2 : p = 1)
		o = ndims(w[1][1]) * (ndims(w[1][1]) - 1) + 1
		c = 2 ^ (ndims(w[1][1]) - 1) - 1
		y = zeros(eltype(w[1][1]), 2 .* size(w[1][1]))
		for i in Iterators.product(collect([[1; 2] for _ in 1:ndims(w[1][1])])...)
			if all(i[1:length(i) .!= p] .== 1)
				(i[p] == 1 ? lo = w_[J + 1][1] : lo = iw_[J + 1][1])
			elseif all(i[1:length(i) .!= p] .== 2)
				(i[p] == 1 ? lo = w_[J + 1][c + 1] : lo = iw_[J + 1][c + 1])
			else
				if i[2] == 2
					(i[p] == 1 ? lo = w_[J + 1][2] : lo = iw_[J + 1][2])
				else
					(i[p] == 1 ? lo = w_[J + 1][3] : lo = iw_[J + 1][3])
				end
			end
			for k in J:-1:2
				if all(i[1:length(i) .!= p] .== 1)
					(i[p] == 1 ? lo = sfb(lo, w_[k][1:o], collect([sf[j] for j in i[end:-1:1]])...) : lo = sfb(lo, iw_[k][1:o], collect([sf[j] for j in i[end:-1:1]])...))
				elseif all(i[1:length(i) .!= p] .== 2)
					(i[p] == 1 ? lo = sfb(lo, w_[k][c * o .+ (1:o)], collect([sf[j] for j in i[end:-1:1]])...) : lo = sfb(lo, iw_[k][c * o .+ (1:o)], collect([sf[j] for j in i[end:-1:1]])...))
				else
					if i[2] == 2
						(i[p] == 1 ? lo = sfb(lo, w_[k][(c - 2) * o .+ (1:o)], collect([sf[j] for j in i[end:-1:1]])...) : lo = sfb(lo, iw_[k][(c - 2) * o .+ (1:o)], collect([sf[j] for j in i[end:-1:1]])...))
					else
						(i[p] == 1 ? lo = sfb(lo, w_[k][(c - 1) * o .+ (1:o)], collect([sf[j] for j in i[end:-1:1]])...) : lo = sfb(lo, iw_[k][(c - 1) * o .+ (1:o)], collect([sf[j] for j in i[end:-1:1]])...))
					end
				end
			end
			if all(i[1:length(i) .!= p] .== 1)
				(i[p] == 1 ? y .+= sfb(lo, w_[1][1:o], collect([Fsf[j] for j in i[end:-1:1]])...) : y .+= sfb(lo, iw_[1][1:o], collect([Fsf[j] for j in i[end:-1:1]])...))
			elseif all(i[1:length(i) .!= p] .== 2)
				(i[p] == 1 ? y .+= sfb(lo, w_[1][c * o .+ (1:o)], collect([Fsf[j] for j in i[end:-1:1]])...) : y .+= sfb(lo, iw_[1][c * o .+ (1:o)], collect([Fsf[j] for j in i[end:-1:1]])...))
			else
				if i[2] == 2
					(i[p] == 1 ? y .+= sfb(lo, w_[1][(c - 2) * o .+ (1:o)], collect([Fsf[j] for j in i[end:-1:1]])...) : y .+= sfb(lo, iw_[1][(c - 2) * o .+ (1:o)], collect([Fsf[j] for j in i[end:-1:1]])...))
				else
					(i[p] == 1 ? y .+= sfb(lo, w_[1][(c - 1) * o .+ (1:o)], collect([Fsf[j] for j in i[end:-1:1]])...) : y .+= sfb(lo, iw_[1][(c - 1) * o .+ (1:o)], collect([Fsf[j] for j in i[end:-1:1]])...))
				end
			end
		end
		return y ./ sqrt(2 ^ ndims(w[1][1]))
	else # 1D dtℂWT
		# tree 1
		y1 = w[J + 1]
		for k in J:-1:2
			y1 = sfb(y1, w[k], sf[1])
		end
		y1 = sfb(y1, w[1], Fsf[1])

		# tree 2
		y2 = iw[J + 1]
		for k in J:-1:2
			y2 = sfb(y2, iw[k], sf[2])
		end
		y2 = sfb(y2, iw[1], Fsf[2])

		# normalization
		return (y1 .+ y2) ./ sqrt(2)
	end
end

# dual-tree real wavelet transform
function dtrwt(x::Union{Array{T, 2}, Array{T, 3}}, Faf, af, J = floor(Int, log2(minimum(size(x))))) where {T <: Number}
	w = Array{Array{Float64, ndims(x)}}[]

	if ndims(x) == 2 # 2D dtRWT
		# tree 1
		x1, hi = afb(x ./ sqrt(2), Faf[1])
    	push!(w, hi)
		for k in 2:J
       		x1, hi = afb(x1, af[1])
       		push!(w, hi)
    	end
		push!(w, [x1])

		# tree 2
		x2, hi = afb(x ./ sqrt(2), Faf[2])
    	append!(w[1], hi)
		for k in 2:J
       		x2, hi = afb(x2, af[2])
       		append!(w[k], hi)
    	end
		append!(w[J + 1], [x2])

		# sum and difference
		for k in 1:J
			for i in 1:3
				w[k][i], w[k][3 + i] = pm(w[k][i], w[k][3 + i])
			end
		end
	else # 3D dtRWT
		M = [1 1 1;
			2 2 1;
			2 1 2;
			1 2 2]
		
		for i in 1:4
			xi, hi = afb(x ./ 2, collect([Faf[j] for j in M[i, :]])...)
			(i == 1 ? push!(w, hi) : append!(w[1], hi))
			for k in 2:J
				xi, hi = afb(xi, collect([af[j] for j in M[i, :]])...)
				(i == 1 ? push!(w, hi) : append!(w[k], hi))
			end
			(i == 1 ? push!(w, [xi]) : append!(w[J + 1], [xi]))
		end

		# sum and difference
		for k in 1:J
			for i in 1:7
				w[k][i], w[k][7 + i], w[k][14 + i], w[k][21 + i] = pm4(w[k][i], w[k][7 + i], w[k][14 + i], w[k][21 + i])
			end
		end
	end
	return w
end

# inverse dual-tree real wavelet transform
function idtrwt(w, Fsf, sf)
	J = length(w) - 1
	w_ = deepcopy(w)

	if ndims(w[1][1]) == 2 # 2D idtRWT
		# sum and difference
		for k in 1:J
			for i in 1:3
				w_[k][i], w_[k][3 + i] = pm(w_[k][i], w_[k][3 + i])
			end
		end

		# tree 1
		y1 = w_[J + 1][1]
		for k in J:-1:2
			y1 = sfb(y1, w_[k][1:3], sf[1])
		end
		y1 = sfb(y1, w_[1][1:3], Fsf[1])

		# tree 2
		y2 = w_[J + 1][2]
		for k in J:-1:2
			y2 = sfb(y2, w_[k][4:6], sf[2])
		end
		y2 = sfb(y2, w_[1][4:6], Fsf[2])

		# normalization
		return (y1 .+ y2) ./ sqrt(2)
	else
		# sum and difference
		for k in 1:J
			for i in 1:7
				w[k][i], w[k][7 + i], w[k][14 + i], w[k][21 + i] = pm4inv(w[k][i], w[k][7 + i], w[k][14 + i], w[k][21 + i])
			end
		end

		M = [1 1 1;
			2 2 1;
			2 1 2;
			1 2 2]

		y = zeros(eltype(w[1][1]), 2 .* size(w[1][1]))
		for i in 1:4
			yi = w[J + 1][i]
			for k in J:-1:2
				yi = sfb(yi, w[k][(7 * (i - 1)) .+ (1:7)], collect([sf[j] for j in M[i, :]])...)
			end
			y .+= sfb(yi, w[1][(7 * (i - 1)) .+ (1:7)], collect([Fsf[j] for j in M[i, :]])...)
		end

		# normalization
		return y ./ 2
	end
end