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

## We do not need to track phase for FermiOps (for the set of Fermi ops that we support)
accumulate_phase(old_phase_data, op1::T, op2::T) where {T <: AbstractFermiOp} = old_phase_data
compute_phase(::Type{<:AbstractFermiOp}, _) = 1
