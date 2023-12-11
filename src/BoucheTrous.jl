module BoucheTrous

using SimpleAlgebra, AbstractFFTs,OptimPackNextGen,CartesianBoxes,Zygote

export bouchetrous!

function bouchetrous!(img::AbstractArray{T,N}, badpix::AbstractArray{Bool,N}; pupil_radius=nothing, kwds...) where {T,N}
	
	sz = size(img)

	F = make_function(T,sz[1:2],pupil_radius)

	Threads.@threads for n ∈ axes(img,3)
		data = img[:,:,n]
		S = LinOpSelect(T,badpix[:,:,n])
		H = F * (S' + data)
		cost(x) = H*x
		p = vmlmb( cost, S*data; lower=minimum(data), autodiff=true,kwds...)
		img[:,:,n] .+=  S'*p
	end
end


function make_function(::Type{T},sz::NTuple{2,Int},pupil_radius::AbstractFloat) where T
	pupil = sqrt.(rfftfreq(sz[1]).^2 .+ (fftfreq(sz[2]).^2)') .< pupil_radius
	F = LinOpDFT(T,sz[1:2])
	D = LinOpDiag(T,.!pupil)
	C = CostL2(T,size(pupil),0)
	return C*D*F
end

function make_function(::Type{T},sz::NTuple{2,Int},::Nothing) where T
	G = LinOpGrad(T,sz)
	C = CostL2(T,outputsize(G),0)
	return C*G
end

function get_box(badpix)
	bx = findfirst(n -> !all(view(badpix,n,:,:)), axes(badpix,1))
	ex = findlast(n -> !all(view(badpix,n,:,:)), axes(badpix,1))
	by = findfirst(n -> !all(view(badpix,:,n,:)), axes(badpix,2))
	ey = findlast(n -> !all(view(badpix,:,n,:)), axes(badpix,2))
	return CartesianBox(bx:ex,by:ey)
end


end # module BoucheTrous
