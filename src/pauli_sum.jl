"""
    struct PauliSum{StringT, CoeffT}

Represents a weighted sum (ie a linear combination) of multi-qubit Pauli strings.

By default `PauliSum`s are constructed and maintained with terms sorted in a canonical order
and with no duplicate Pauli strings.
"""
struct PauliSum{StringT, CoeffT} <: AbstractSum{StringT, CoeffT}
    strings::StringT
    coeffs::CoeffT

    function PauliSum(strings, coeffs, already_sorted=false)
        _abstract_sum_inner_constructor_helper!(strings, coeffs, already_sorted)
        return new{typeof(strings), typeof(coeffs)}(strings, coeffs)
    end
end

function term_type(::Type{T}) where T <: PauliSum
    return PauliTerm
end

function sum_type(::Type{T}) where T <: PauliTerm
    return PauliSum
end

####
#### Constructors
####

function Base.similar(ps::PauliSum{Vector{Vector{PauliT}}, Vector{C}}, n=0) where {PauliT, C}
    m = size(ps, 2)
    strings = [Vector{PauliT}(undef, m) for i in 1:n]
    coeffs = Vector{C}(undef, n)
    return PauliSum(strings, coeffs)
end

"""
    PauliSum(::Type{PauliT}) where PauliT <: AbstractPauli

Return an empty `PauliSum` with Paulis of type `PauliT`
and `Complex{Float64}` coefficients.
"""
PauliSum(::Type{PauliT}) where PauliT <: AbstractPauli = PauliSum(Vector{PauliT}[], Complex{Float64}[])

"""
    PauliSum(strings)

Construct a sum from `strings` with coefficients all equal to one.
"""
PauliSum(strings) = PauliSum(strings, fill(_DEFAULT_COEFF, length(strings)))

"""
        PauliSum(v::AbstractVector{<:PauliTerm}, already_sorted=false)

Construct a sum from an array of `PauliTerm`s.
"""
function PauliSum(v::AbstractVector{<:PauliTerm}, already_sorted=false)
    strings = [x.paulis for x in v]
    coeffs = [x.coeff for x in v]
    return PauliSum(strings, coeffs, already_sorted)
end

PauliSum{T,V}(v, c, already_sorted) where {T, V} = PauliSum(v, c, already_sorted)
PauliSum{T,V}(v, c) where {T, V} = PauliSum(v, c, false)

"""
    PauliSum(v::AbstractMatrix{<:AbstractPauli}, coeffs=fill(_DEFAULT_COEFF, size(v, 1)))

Construct a sum from a matrix of single-qubit Pauli operators. If `size(v) == (m, n)`, then
the the sum has `m` terms with `n` Paulis in each string.
"""
function PauliSum(v::AbstractMatrix{<:AbstractPauli}, coeffs=fill(_DEFAULT_COEFF, size(v, 1)))
    strings = @inbounds [v[i,:] for i in 1:size(v,1)]
    return PauliSum(strings, coeffs)
end

"""
    PauliSum(::Type{PauliT}, matrix::AbstractMatrix{<:Number}; threads=true)

Construct a Pauli decomposition of `matrix`, that is, a `PauliSum` representing `matrix`.
If `thread` is `true`, use a multi-threaded algorithm for increased performance.
"""
function PauliSum(::Type{PauliT}, matrix::AbstractMatrix{<:Number}; threads=true) where PauliT
    if threads
        return pauli_sum_from_matrix_threaded(PauliT, matrix)
    else
        return pauli_sum_from_matrix_one_thread(PauliT, matrix)
    end
end

function pauli_sum_from_matrix_one_thread(::Type{PauliT}, matrix::AbstractMatrix{<:Number}) where PauliT
    nside = LinearAlgebra.checksquare(matrix)
    n_qubits = ILog2.checkispow2(nside)
    denom = 2^n_qubits  # == nside
    s = PauliSum(PauliT)
    for pauli in pauli_basis(PauliT, n_qubits)
        mp = SparseArrays.sparse(pauli)  # Much faster than dense
        coeff = LinearAlgebra.dot(mp, matrix)
        if ! isapprox_zero(coeff)
            push!(s, (pauli.paulis, coeff / denom))  # a bit faster than PauliTerm for small `matrix` (eg 2x2)
        end
    end
    return s
end

function pauli_sum_from_matrix_threaded(::Type{PauliT}, matrix::AbstractMatrix{<:Number}) where PauliT
    nside = LinearAlgebra.checksquare(matrix)
    n_qubits = ILog2.checkispow2(nside)
    denom = 2^n_qubits  # == nside
    ## Create a PauliSum for each thread, for accumulation.
    sums = [PauliSum(PauliT) for i in 1:Threads.nthreads()]
    Threads.@threads for j in 0:(4^n_qubits - 1)
        pauli = PauliTerm(PauliT, j, n_qubits)
        mp = SparseArrays.sparse(pauli)  # Much faster than dense
        coeff = LinearAlgebra.dot(mp, matrix)
        if ! isapprox_zero(coeff)
            push!(sums[Threads.threadid()], (pauli.paulis, coeff / denom))
        end
    end
    for ind in 2:length(sums)  # Collate results from all threads.
        add!(sums[1], sums[ind])
    end
    return sums[1]
end

"""
    rand_pauli_sum(::Type{PauliT}=PauliDefault, n_factors::Integer, n_terms::Integer; coeff_func=nothing) where {PauliT <: AbstractPauli}

Return a `PauliSum` of `n_terms` terms of `n_factors` factors each.

If `coeff_func` is `nothing`, then the coefficients are all equal to one. Otherwise `coeff_func` must be a function that
takes one argument `n_terms`, and returns `n_terms` coefficients.

# Examples
```julia
julia> rand_pauli_sum(4, 3)
(1 + 0im) * IIYX
(1 + 0im) * IIZY
(1 + 0im) * ZXIX

julia> rand_pauli_sum(3, 4; coeff_func=randn)
-0.1929668228083002 * XZY
-0.5661154429051878 * YXX
0.01820960911335266 * YXZ
-0.33717575080125783 * ZII
```
"""
function rand_pauli_sum(::Type{PauliT}, n_factors::Integer, n_terms::Integer; coeff_func=nothing) where {PauliT <: AbstractPauli}
    paulis = [rand(PauliT, n_factors) for i in 1:n_terms]
    if coeff_func == nothing
        coeffs = fill(_DEFAULT_COEFF, n_terms)
    else
        coeffs = coeff_func(n_terms)
    end
    return PauliSum(paulis, coeffs)
end

rand_pauli_sum(n_factors::Integer, n_terms::Integer; coeff_func=nothing) =
    rand_pauli_sum(PauliDefault, n_factors, n_terms; coeff_func=coeff_func)

#####
##### Conversion
#####

"""
    Matrix(ps::PauliSum)

Convert `ps` to a dense `Matrix`.

# Examples

We convert a matrix to a `PauliSum` and then back to a matrix.
```jldoctest
julia> m = [0.1 0.2; 0.3 0.4];

julia> PauliSum(m)
(0.25 + 0.0im) * I
(0.25 + 0.0im) * X
(0.0 - 0.04999999999999999im) * Y
(-0.15000000000000002 + 0.0im) * Z

julia> Matrix(PauliSum(m)) ≈ m
true
```
"""
Base.Matrix(ps::PauliSum) = Matrix(SparseArrays.sparse(ps))

## ThreadsX helps enormously for large sums. 22x faster for 4^8 terms
"""
    sparse(ps::PauliSum)

Convert `ps` to a sparse matrix.
"""
SparseArrays.sparse(ps::PauliSum) = ThreadsX.sum(SparseArrays.sparse(ps[i]) for i in eachindex(ps))

# Using Z4Group0 is 30% faster in many tests, for dense matrices
# Base.Matrix(ps::PauliSum) = ThreadsX.sum(Matrix(Z4Group0, ps[i]) for i in eachindex(ps))

####
#### Algebra / mathematical operations
####

# We use lmul! because that's how LinearAlgebra offers "scaling" of a Matrix (or rmul!)
"""
    lmul!(psum::PauliSum, n)

Left multiplies the coefficient of `psum` by `n` in place.
"""
function LinearAlgebra.lmul!(psum::PauliSum, n)
    @. psum.coeffs = n * psum.coeffs
    return psum
end

####
#### Math
####

"""
    cis(p::PauliTerm)::PauliSum

Compute ``\\exp(i p)``
"""
function Base.cis(pt::PauliTerm)
    return PauliSum([one(pt).paulis, copy(pt.paulis)], [cos(pt.coeff), im * sin(pt.coeff)], true)
end

"""
    exp(p::PauliTerm)::PauliSum

Compute ``\\exp(p)``
"""
function Base.exp(pt::PauliTerm)
    coeff = im * pt.coeff
    return PauliSum([one(pt).paulis, copy(pt.paulis)], [cos(coeff), -im * sin(coeff)], true)
end

"""
    numeric_function(pt::PauliTerm, f)::PauliSum

Compute `f(pt)` by decomposing `f` into odd and even functions.
"""
function numeric_function(pt::PauliTerm, f)
    c = pt.coeff
    fe = (f(c) + f(-c)) / 2
    fo = (f(c) - f(-c)) / 2
    strings = [one(pt).paulis, copy(pt.paulis)]
    coeffs = [fe, fo]
    already_sorted = true # else sorting takes 30x longer
    return PauliSum(strings, coeffs, already_sorted)
end

# Julia 1.5 does not have cispi
for f in (:cos, :sin, :tan, :sqrt, :sind, :sinpi, :cospi, :sinh, :tanh,
          :acos, :asin, :atan, :sec, :csc, :cot, :log, :log2, :log10,
          :log1p)
    @eval Base.$f(pt::PauliTerm) = numeric_function(pt, $f)
end
