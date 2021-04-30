import Random

using ..AbstractOps
using ..AbstractOps: _show_op_plain

abstract type AbstractFermiOp <: AbstractOp end

AbstractOps.rand_ind_range(::Type{<:AbstractFermiOp}) = 0:4

# function Random.rand(rng::Random.AbstractRNG,
#                      ::Random.SamplerType{FermiOpT}) where {FermiOpT <: AbstractFermiOp}
#     return FermiOpT(rand(rng, 0:4))
# end

function Base.show(io::IO, ps::AbstractArray{<:AbstractFermiOp})
    for p in ps
        _show_op_plain(io, p)
    end
end

# AbstractOps.accumulate_phase(old_phase_data, op1::T, op2::T) where {T <: AbstractFermiOp} = old_phase_data
# AbstractOps.compute_phase(::Type{<:AbstractFermiOp}, _) = 1
