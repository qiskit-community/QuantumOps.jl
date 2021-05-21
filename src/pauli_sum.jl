"""
    struct PauliSum{StringT, CoeffT}

Represents a weighted sum (ie a linear combination) of multi-qubit Pauli strings.

By default `PauliSum`s are constructed and maintained with terms sorted in a canonical order
and with no duplicate Pauli strings.
"""
struct PauliSum{StringT, CoeffT}
    strings::StringT
    coeffs::CoeffT

    function PauliSum(strings, coeffs, already_sorted=false)
        if length(strings) != length(coeffs)
            throw(DimensionMismatch("bad dims"))
        end
        if ! isempty(strings)
            n = length(first(strings))
            if ! all(x -> length(x) == n, strings)
                throw(DimensionMismatch("Pauli strings are of differing lengths."))
            end
            if ! already_sorted
                sort_and_sum_duplicates!(strings, coeffs)
            end
        end
        return new{typeof(strings), typeof(coeffs)}(strings, coeffs)
    end
end

####
#### Constructors
####

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

function Base.copy(ps::PauliSum)
    (new_strings, new_coeffs) = (copy.(ps.strings), copy(ps.coeffs))
    already_sorted = true
    return PauliSum(new_strings, new_coeffs, already_sorted)
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

####
#### Canonicalization / sorting
####

### These are helpers for constructors

function sort_and_sum_duplicates!(psum::PauliSum)
    sort_and_sum_duplicates!(psum.strings, psum.coeffs)
    return psum
end

# function sort_and_sum_duplicates!(terms, coeffs)
#     sort_sums!(terms, coeffs)
#     sum_duplicates!(terms, coeffs)
#     remove_zeros!(terms, coeffs)
#     return nothing
# end

"""
    remove_zeros!(psum::PauliSum)
    remove_zeros!(terms, coeffs)

Remove terms from `psum` with coefficient (approximately) equal to zero.
If `terms` and `coeffs` are supplied, then elements are deleted from both `terms`
and `coeffs` at indices corresponding to vanishing elements of `coeff`.
"""
function remove_zeros!(psum::PauliSum)
    remove_zeros!(psum.strings, psum.coeffs)
    return psum
end

# function remove_zeros!(terms, coeffs)
#     # ThreadsX is very slow for small arrays. We need to discriminate
#     # inds = ThreadsX.findall(isapprox_zero, coeffs)
#     # The following is 500ns for two non-zero floats. What is wrong?
#     # Appears to be this: iszero.(array) is taking almost all the time.
#     # The following is what Base does, but writing it out is faster. A bug.
#     inds = findall(isapprox_zero.(coeffs))
#     deleteat!(coeffs, inds)
#     deleteat!(terms, inds)
#     return nothing
# end

## This is 10x faster than the sorting step, even though we don't check if all
## strings are already unique.
sum_duplicates!(psum::PauliSum) = sum_duplicates!(psum.strings, psum.coeffs)

## Modeled on code in unique!
"""
    sum_duplicates!(paulis, coeffs)

Find groups of terms whose members differ only in the coefficient.
Replace each group by one term carrying the sum of the coefficients
in that group. This routine assumes that `paulis` are sorted.
"""
function sum_duplicates!(paulis, coeffs)
    last_pauli::eltype(paulis) = first(paulis)
    coeff = first(coeffs)
    k = 2
    @inbounds for j in 2:length(paulis)
        if paulis[j] != last_pauli
            last_pauli = paulis[k] = paulis[j]
            coeffs[k] = coeffs[j]
            k += 1
        else
            coeffs[k-1] += coeffs[j]
        end
    end
    resize!(paulis, k-1)
    resize!(coeffs, k-1)
    return nothing
end

## This is expensive. Most time is spent in sortperm.
## There is no ThreadsX.sortperm, only sort.
Base.sort!(psum::PauliSum) = (sort_sums!(psum.strings, psum.coeffs); psum)
# function sort_sums!(strings, coeffs)
#     p = sortperm(strings; alg=MergeSort) # alg=MergeSort is 50% faster for 1000x10 strings
#     permute!(strings, p)
#     permute!(coeffs, p)
#     return nothing
# end

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
#### IO
####

function Base.show(io::IO, psum::PauliSum)
    for i in eachindex(psum)
        show(io, psum[i])
        if i != lastindex(psum)
            print(io, "\n")
        end
    end
end

####
#### Container interface
####

## Fails for empty psum
Base.size(psum::PauliSum) = (length(psum), length(first(psum)))

Base.size(psum::PauliSum, i::Integer) = size(psum)[i]

Base.getindex(psum::PauliSum, j::Integer) = PauliTerm(psum.strings[j], psum.coeffs[j])

Base.getindex(psum::PauliSum, j::Integer, k::Integer) = psum.strings[j][k]
# TODO: Use already_sorted flag ?

Base.getindex(psum::PauliSum, inds) = PauliSum(psum.strings[inds], psum.coeffs[inds])

# Iterate uses getindex to return `PauliTerm`s.
function Base.iterate(psum::PauliSum, state=1)
    state > lastindex(psum) && return nothing
    return (psum[state], state + 1)
end

# Enables using `findall`, for instance.
# Fallback methods for `values` and `pairs` are OK.
Base.keys(psum::PauliSum) = eachindex(psum)

for func in (:length, :eachindex, :lastindex, :firstindex)
    @eval begin
        function Base.$func(ps::PauliSum, args...)
            return $func(ps.coeffs, args...)
        end
    end
end

"""
    reverse(ps::PauliSum)

Reverse qubit order in `ps` and sort terms.
"""
Base.reverse(ps::PauliSum) = reverse!(copy(ps))

"""
    reverse!(ps::PauliSum)

Reverse qubit order in `ps` in place and sort terms.
"""
function Base.reverse!(ps::PauliSum)
    strings = ps.strings
    @inbounds for i in eachindex(strings)
        reverse!(strings[i])
    end
    return Base.sort!(ps)
end

####
#### Updating / adding elements
####

function Base.deleteat!(ps::PauliSum, args...)
    deleteat!(ps.coeffs, args...)
    deleteat!(ps.strings, args...)
    return ps
end

"""
    insert!(ps::PauliSum, ind, p::PauliTerm)

Insert `p` into `ps` without sorting resulting `ps`.
"""
Base.insert!(ps::PauliSum, ind, p::PauliTerm) = insert!(ps, ind, (p.paulis, p.coeff))

@inline function Base.insert!(ps::PauliSum, ind, (paulis, coeff))
    insert!(ps.strings, ind, paulis)
    insert!(ps.coeffs, ind, coeff)
    return ps
end

# TODO: Should we check that the length of the PauliTerm is correct ?
"""
    push!(psum::PauliSum, ps::PauliTerm...)

Push `ps` to the end of `psum` without regard to order
or possible duplication.

See `sort_and_sum_duplicates!`.
"""
function Base.push!(psum::PauliSum, ps::PauliTerm...)
    for p in ps
        push!(psum.strings, p.paulis)
        push!(psum.coeffs, p.coeff)
    end
    return ps
end

function Base.push!(psum::PauliSum, (string, coeff))
    push!(psum.strings, string)
    push!(psum.coeffs, coeff)
end

####
#### Compare / predicates
####

# This will fail for an empty `psum`. Use type info instead.
# There no well-defined `one` for `PauliSum`. It depends on the
# width of the string.
function Base.one(psum::PauliSum)
    t = one(first(psum))
    already_sorted = true
    PauliSum([t.paulis], [t.coeff], already_sorted)
end

function Base.:(==)(psum1::PauliSum, psum2::PauliSum)
    if length(psum1) != length(psum2)
        return false
    end
    # This is 8x faster for large arrays, 10^4 or 5.
    return ThreadsX.all(i -> psum1[i] == psum2[i], eachindex(psum1))
end

####
#### Algebra / mathematical operations
####

"""
    add!(to::PauliSum, from::PauliSum)

Adds the terms in `from` to `to` in place. `to` is mutated. `from` is not.
"""
function add!(to::PauliSum, from::PauliSum)
    for p in from
        add!(to, p)
    end
    return to
end

"""
    add!(psum::PauliSum, pt::PauliTerm...)

Add `PauliTerm`s to `psum` in place, assuming `psum` is sorted and has no repeated
strings. Either a new term is inserted, or the coefficient is added to an existing
term. After adding the `pt`, `psum` will be left sorted, with no duplicates, and
no zero coefficients. Use `push!` to insert a term at the end of `psum` with no
simplification performed.
"""
function add!(psum::PauliSum, pt::PauliTerm...)
    for p in pt
        add!(psum, p)
    end
    return psum
end

add!(psum::PauliSum, pt::PauliTerm) = add!(psum, pt.paulis, pt.coeff)

function add!(psum::PauliSum, paulis, coeff)
    inds = searchsorted(psum.strings, paulis)
    if length(inds) == 0 # paulis not found, add a new term
        insert!(psum, first(inds), (paulis, coeff))
    elseif length(inds) == 1 # one element equal to paulis
        i = first(inds) # get the (single) index
        @inbounds psum.coeffs[i] += coeff # add p to existing term
        @inbounds if isapprox_zero(psum.coeffs[i])
            @inbounds deleteat!(psum, [i])
        end
    else
        throw(ErrorException("Duplicate terms found in `PauliSum`."))
    end
    return psum
end

# We use lmul! because that's how LinearAlgebra offers "scaling" of a Matrix (or rmul!)
"""
    lmul!(psum::PauliSum, n)

Left multiplies the coefficient of `psum` by `n` in place.
"""
function LinearAlgebra.lmul!(psum::PauliSum, n)
    @. psum.coeffs = n * psum.coeffs
    return psum
end

Base.:+(terms::T...) where {T <: PauliTerm} = PauliSum([terms...])

function Base.:+(ps0::PauliSum, pss::PauliSum...)
    ps_out = copy(ps0)
    for ps in pss
        add!(ps_out, ps)
    end
    return ps_out
end

## TODO: Do something more efficient here.
function Base.:-(pt1::PauliTerm, pt2::PauliTerm)
    return PauliSum([pt1, -one(pt2.coeff) * pt2])
end

function Base.:-(psum::PauliSum)
    already_sorted = true
    PauliSum(psum.strings, -one(eltype(psum.coeffs)) .* psum.coeffs, already_sorted)
end

function Base.:*(n::Number, psum::PauliSum)
    already_sorted = true
    PauliSum(psum.strings, n .* psum.coeffs, already_sorted)
end

function Base.:/(psum::PauliSum, n)
    already_sorted = true
    PauliSum(psum.strings, psum.coeffs ./ n, already_sorted)
end

## TODO: new_coeffs may be a Vector{<:Real}, but because of phase
## we need to set a coefficient to a complex type, which errors out.
## Currently, the user needs to make the coefficients of psum complex.
## We could also widen the type here somehow, say through promotion.
## A similar situation arises in other places in this library.
function Base.:*(pterm::PauliTerm, psum::PauliSum)
    new_coeffs = similar(psum.coeffs)
    new_strings = similar(psum.strings)
    ## Threading does not help; PauliSum() call is bottleneck
    ## And it is slower anyway for small sums
    @inbounds for j in eachindex(psum)
        new_term = pterm * psum[j]
        new_coeffs[j] = new_term.coeff
        new_strings[j] = new_term.paulis
    end
    return PauliSum(new_strings, new_coeffs)
end

"""
    *(ps1::PauliSum, ps2::PauliSum)

Multiply `ps1` and `ps2` returning another `PauliSum`

# Examples
```jldoctest
julia> a = rand_pauli_sum(5, 3);

julia> b = rand_pauli_sum(5, 3);

julia> a = rand_pauli_sum(5, 3); ma = Matrix(a);

julia> b = rand_pauli_sum(5, 3); mb = Matrix(b);

julia> Matrix(a * b * a) == ma * mb * ma
true
```
"""
function Base.:*(ps1::PauliSum, ps2::PauliSum)
    psum_out = ps1[1] * ps2
    for i in 2:length(ps1)
        add!(psum_out, ps1[i] * ps2)
    end
    return psum_out
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
