using .FermiOps: number_op, raise_op, lower_op, no_op, I_op

struct FermiTermA{W<:AbstractFermiOp, T<:AbstractVector{W}, CoeffT} <: AbstractTerm{W, CoeffT}
    ops::T
    coeff::CoeffT
end

const FermiTerms = Union{FermiTermA, OpTerm{<:AbstractFermiOp}}

op_string(t::FermiTermA) = t.ops

# Use OpTerm{FermiOp} instead
# term_type(::Type{<:AbstractFermiOp}) = FermiTermA

strip_typeof(::FermiTermA) = FermiTermA

function Base.promote_rule(::Type{FermiTermA{FermiOp, Vector{FermiOp}, T}},
                           ::Type{FermiTermA{FermiOp, Vector{FermiOp}, V}}) where {T, V}
    return FermiTermA{FermiOp, Vector{FermiOp}, promote_type(T, V)}
end

####
#### Constructors
####

function FermiTermA(::Type{T}, s::AbstractString, coeff=_DEFAULT_COEFF) where T <: AbstractFermiOp
    return FermiTermA(Vector{T}(s), coeff)
end

FermiTermA(s::AbstractString, coeff=_DEFAULT_COEFF) = FermiTermA(FermiOp, s, coeff)

function FermiTermA(inds::NTuple{N, Int}, coeff, n_modes::Integer) where N
    (ops, phase) = index_to_ops_phase(inds)
    (factors, new_coeff) = _fermi_term(ops, phase, coeff, n_modes)
    return FermiTermA(factors, new_coeff)
end

## This is defined just for construction by copying the printed type.
function FermiTermA{W, T, V}(s::AbstractString,
                            coeff=_DEFAULT_COEFF) where {W<:AbstractFermiOp,T<:AbstractVector{W},V}
    return FermiTermA{W,T,V}(Vector{W}(s), coeff)
end

index_to_ops(inds) = index_to_ops(inds...)

function index_to_ops(i, j, k, l)
    ## We filter this elsewhere, but it's dangerous to remove this
    # if i == j || k == l
    #     return ()  # unstable, use a barrier or move this out
    # end
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

function index_to_ops(i, j)
    if i == j
        return ((op=number_op, ind=i), )
    end
    if i < j
        return ((op=raise_op, ind=i), (op=lower_op, ind=j))
    else
        return ((op=lower_op, ind=i), (op=raise_op, ind=j))
    end
end

function count_permutations(inds::NTuple{2})
    if inds[2] < inds[1]
        return 1
    else
        return 0
    end
end

index_phase(inds::NTuple) = iseven(count_permutations(inds)) ? 1 : -1

function count_permutations(inds::NTuple{4})
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

## Embed `ops` in a string of width `n_modes`
function _fermi_term(ops::NTuple, phase::Integer, coeff, n_modes::Integer)
    factors = fill(I_op, n_modes)
    for (op, ind) in ops
        if op == no_op
            continue
        end
        factors[ind] = op
    end
    return (factors, phase * coeff)
end

## Unused ?
# """
#     ferm_terms(two_body)

# Convert rank-4 tensor `two_body` to a `Vector` of `FermiTerm`s.
# """
# function ferm_terms(two_body)
#     n_modes = first(size(two_body))
#     ## FIXME: the following is rigid. Could as well be Vector{Any} ?
#     terms = FermiTerm{FermiOp, Vector{FermiOp}, eltype(two_body)}[]
#     for ind in CartesianIndices(two_body)
#         if (ind[1] == ind[2] || ind[3] == ind[4]) || two_body[ind] == 0
#             continue
#         end
#         push!(terms, FermiTerm(ind.I, two_body[ind], n_modes))
#     end
#     return terms
# end

"""
    rand_fermi_term(::Type{FermiT}=FermiDefault, n::Integer; coeff=_DEFAULT_COEFF) where {FermiT <: AbstractFermiOp}

Return a random `FermiTermA` of `n` tensor factors.
"""
function rand_fermi_term(::Type{FermiT}, n::Integer; coeff=_DEFAULT_COEFF) where {FermiT <: AbstractFermiOp}
    return FermiTermA(rand(FermiT, n), coeff)
end

rand_fermi_term(n::Integer; coeff=_DEFAULT_COEFF) = rand_fermi_term(FermiDefault, n, coeff=coeff)

####
#### Compare / predicates
####

function Base.zero(ft::FermiTerms)
    new_string = similar(op_string(ft))
    fill!(new_string, zero_op)
    return strip_typeof(ft)(new_string, zero(ft.coeff))
end

Base.iszero(ft::FermiTerms) = iszero(ft.coeff) || any(iszero, op_string(ft))

####
#### Algebra / mathematical operations
####

function Base.:*(ft1::FermiTerms, ft2::FermiTerms)
    if length(ft1) != length(ft2)
        throw(DimensionMismatch())
    end
    new_string = similar(op_string(ft1))
    phase_count = 0
    @inbounds for i in eachindex(ft1)
        f1 = op_string(ft1)[i]
        f2 = op_string(ft2)[i]
        _prod = f1 * f2
        if iszero(_prod)  # If any factor is zero, the entire product is zero
            fill!(new_string, zero_op)
            return strip_typeof(ft1)(new_string, zero(ft1.coeff))
        end
        new_string[i] = _prod
        # Tracking phase costs sometimes 5% - 20%, in some cases
#        phase_count += do_phase_count(f1, f2)
    end
    phase_fac = iseven(phase_count) ? 1 : -1
    return strip_typeof(ft1)(new_string, ft1.coeff * ft2.coeff * phase_fac)
end

## I have a more generic method in abstract_term.jl, but it is not
## flexible to do what we need here. I don't see a way to avoid this boiler plate
# function Base.:*(z::Number, term::FermiTerm{W, T, CoeffT}) where {W, T, CoeffT}
#     return FermiTerm{W, T, promote_type(CoeffT, typeof(z))}(op_string(term), term.coeff * z)
# end

function Base.:^(ft::FermiTerms, n::Integer)
    n < 0 && throw(DomainError(n))
    n == 0 && return one(ft)
    n == 1 && return ft
    ## Any +,-,0 sends the entire string to zero
    if any(x -> x === raise_op || x === lower_op || x === zero_op, op_string(ft))
        return zero(ft)
    end
    ## E, N, I, are idempotent
    return strip_typeof(ft)(copy(op_string(ft)), ft.coeff^n)
end

Base.adjoint(ft::FermiTerms) = strip_typeof(ft)(adjoint.(op_string(ft)), conj(ft.coeff))

count_bodies(ft::FermiTerms) = count_bodies(op_string(ft))

function count_bodies(v::Vector{FermiOp})
    n = sum(count_raise_lower, v)
    return n ÷ 2
end
