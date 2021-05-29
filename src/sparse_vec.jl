@inline neutral(x) = zero(x)
@inline isneutral(x) = iszero(x)

@inline neutral(T::Type{<:AbstractOp}) = one(T)
@inline neutral(x::AbstractOp) = one(x)
@inline isneutral(x::AbstractOp) = isone(x)

abstract type AbstractSparseVec{Tv, Ti} <:AbstractArray{Tv, 1} end

struct SparseVec{Tv,Ti<:Integer} <: AbstractSparseVec{Tv, Ti}
    n::Int            # Length of the sparse vector
    ind::Vector{Ti}   # Indices of stored values
    val::Vector{Tv}   # Stored values

    function SparseVec{Tv,Ti}(n::Integer, ind::Vector{Ti}, val::Vector{Tv}) where {Tv,Ti<:Integer}
        n >= 0 || throw(ArgumentError("The number of elements must be non-negative."))
        length(ind) == length(val) ||
            throw(ArgumentError("index and value vectors must be the same length"))
        new(convert(Int, n), ind, val)
    end
end

nvals(sp::AbstractSparseVec) = sp.n  # like nnz
inds(sp::AbstractSparseVec) = sp.ind  # like nonzeroinds
vals(sp::AbstractSparseVec) = sp.val  # like nonzerovals
#Base.length(sp::AbstractSparseVec) = sp.n

function SparseVec(ind::AbstractVector{Ti}, val::AbstractVector{Tv}; sorted=true) where {Ti, Tv}
    n = ind[end]
    return SparseVec{Tv, Ti}(n, ind, val)
end

SparseVec(n::Integer, ind::Vector{Ti}, val::Vector{Tv}) where {Tv,Ti} =
    SparseVec{Tv,Ti}(n, ind, val)

Base.size(x::SparseVec) = (getfield(x, :n),)
Base.count(f, x::SparseVec) = count(f, vals(x)) + f(neutral(eltype(x)))*(length(x) - nvals(x))

### getindex

function _spgetindex(m::Int, ind::AbstractVector{Ti}, val::AbstractVector{Tv}, i::Integer) where {Tv,Ti}
    ii = searchsortedfirst(ind, convert(Ti, i))
    (ii <= m && ind[ii] == i) ? val[ii] : neutral(Tv)
end

function Base.getindex(x::AbstractSparseVec, i::Integer)
    checkbounds(x, i)
    _spgetindex(nvals(x), inds(x), vals(x), i)
end

function _dense2indval!(nzind::Vector{Ti}, nzval::Vector{Tv}, s::AbstractArray{Tv}) where {Tv,Ti}
    Base.require_one_based_indexing(s)
    cap = length(nzind)
    @assert cap == length(nzval)
    n = length(s)
    c = 0
    @inbounds for (i, v) in enumerate(s)
        if !isneutral(v)
            if c >= cap
                cap = (cap == 0) ? 1 : 2*cap
                resize!(nzind, cap)
                resize!(nzval, cap)
            end
            c += 1
            nzind[c] = i
            nzval[c] = v
        end
    end
    if c < cap
        resize!(nzind, c)
        resize!(nzval, c)
    end
    return (nzind, nzval)
end

function _dense2sparsevec(s::AbstractArray{Tv}, initcap::Ti) where {Tv,Ti}
    ind, val = _dense2indval!(Vector{Ti}(undef, initcap), Vector{Tv}(undef, initcap), s)
    SparseVec(length(s), ind, val)
end

SparseVec{Tv,Ti}(s::AbstractVector{Tv}) where {Tv,Ti} =
    _dense2sparsevec(s, convert(Ti, max(8, div(length(s), 8))))

SparseVec{Tv}(s::AbstractVector{Tv}) where {Tv} = SparseVec{Tv,Int}(s)

SparseVec(s::AbstractVector{Tv}) where {Tv} = SparseVec{Tv,Int}(s)

function Base.show(io::IO, ::MIME"text/plain", x::AbstractSparseVec)
    xnnz = length(vals(x))
    print(io, length(x), "-element ", typeof(x), " with ", xnnz,
           " stored ", xnnz == 1 ? "entry" : "entries")
    if xnnz != 0
        println(io, ":")
        show(IOContext(io, :typeinfo => eltype(x)), x)
    end
end

show(io::IO, x::AbstractSparseVec) = show(convert(IOContext, io), x)
function show(io::IOContext, x::AbstractSparseVec)
    # TODO: make this a one-line form
    n = length(x)
    nzind = inds(x)
    nzval = vals(x)
    if isempty(nzind)
        return show(io, MIME("text/plain"), x)
    end
    limit = get(io, :limit, false)::Bool
    half_screen_rows = limit ? div(displaysize(io)[1] - 8, 2) : typemax(Int)
    pad = ndigits(n)
    if !haskey(io, :compact)
        io = IOContext(io, :compact => true)
    end
    for k = eachindex(nzind)
        if k < half_screen_rows || k > length(nzind) - half_screen_rows
            print(io, "  ", '[', rpad(nzind[k], pad), "]  =  ")
            if isassigned(nzval, Int(k))
                show(io, nzval[k])
            else
                print(io, Base.undef_ref_str)
            end
            k != length(nzind) && println(io)
        elseif k == half_screen_rows
            println(io, "   ", " "^pad, "   \u22ee")
        end
    end
end


