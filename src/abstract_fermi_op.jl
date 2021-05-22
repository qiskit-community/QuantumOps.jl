import Random

abstract type AbstractFermiOp <: AbstractOp end

function Random.rand(rng::Random.AbstractRNG,
                     ::Random.SamplerType{FermiOpT}) where {FermiOpT <: AbstractFermiOp}
    return FermiOpT(rand(rng, 0:4))
end

function Base.show(io::IO, ps::AbstractArray{<:AbstractFermiOp})
    for p in ps
        show(io, p)
    end
end
