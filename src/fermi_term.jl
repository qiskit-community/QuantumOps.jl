using .FermiOps: number_op, raise_op, lower_op, no_op, id_op

#struct FermiTerm{T} <: AbstractTerm
struct FermiTerm{W<:AbstractFermiOp, T<:AbstractVector{W}, V} <: AbstractTerm{W}
    ops::T
    coeff::V
end

op_string(t::FermiTerm) = t.ops

function FermiTerm(::Type{T}, s::AbstractString, coeff=_DEFAULT_COEFF) where T <: AbstractFermiOp
    return FermiTerm(Vector{T}(s), coeff)
end

FermiTerm(s::AbstractString, coeff=_DEFAULT_COEFF) = FermiTerm(FermiOp, s, coeff)

index_to_ops(inds) = index_to_ops(inds...)

function index_to_ops(i, j, k, l)
    if i == j || k == l
        return ()  # unstable, use a barrier or move this out
    end
    if i == k || i == l
        op1 = (op=number_op, ind=i)
    else
        op1 = (op=raise_op, ind=i)
    end
    if j == k || j == l
        op2 = (op=number_op, ind=j)
    else
        op2 = (op=raise_op, ind=j)
    end
    if k != i && k != j
        op3 = (op=lower_op, ind=k)
    else
        op3 = (op=no_op, ind=k)
    end
    if l != i && l != j
        op4 = (op=lower_op, ind=l)
    else
        op4 = (op=no_op, ind=k)
    end
    return (op1, op2, op3, op4)
end

index_phase(inds::NTuple{4}) = iseven(index_perm(inds)) ? 1 : -1

function index_perm(inds::NTuple{4})
    count = 0
    @inbounds i = inds[1]
    @inbounds for j in (2,3,4)
        if i > inds[j]
            count += 1
        end
    end
    @inbounds i = inds[2]
    @inbounds for j in (3,4)
        if i > inds[j]
            count += 1
        end
    end
    @inbounds i = inds[3]
    @inbounds if i > inds[4]
        count += 1
    end
    return count
end

index_to_ops_phase(inds) = (ops=index_to_ops(inds), phase=index_phase(inds))

function ferm_term(inds::NTuple{4, Int}, coeff, n_modes::Integer)
    (ops, phase) = index_to_ops_phase(inds)
    return _ferm_term(ops, phase, coeff, n_modes)
end

## Embed `ops` in a string of width `n_modes`
function _ferm_term(ops::NTuple{4}, phase::Integer, coeff, n_modes::Integer)
    factors = fill(id_op, n_modes)
    for (op, ind) in ops
        if op == no_op
            continue
        end
        factors[ind] = op
    end
    return FermiTerm(factors, phase * coeff)
end

"""
    ferm_terms(h2)

Convert rank-4 tensor `h2` to a `Vector` of `FermiTerm`s.
"""
function ferm_terms(h2)
    n_modes = first(size(h2))
    terms = FermiTerm{FermiOp, Vector{FermiOp}, eltype(h2)}[]
    for ind in CartesianIndices(h2)
        if (ind[1] == ind[2] || ind[3] == ind[4]) || h2[ind] == 0
            continue
        end
        push!(terms, ferm_term(ind.I, h2[ind], n_modes))
    end
    return terms
end
