abstract type AbstractTerm end

## type params only to get correct dispatch. There must be a better way
function Base.show(io::IO, term::AbstractTerm) # where { T, V, CoeffT, AbstractTerm{T,V,CoeffT}}
    print(io, op_string(term))
    print(io, " * ")
    if term.coeff isa Real  # could use CoeffT here.
        print(io, term.coeff)
    else
        print(io, "(", term.coeff, ")")
    end
end

####
#### Container interface
####

# :popat!
for func in (:length, :size, :eltype, :eachindex, :axes, :splice!, :getindex,
             :setindex!, :iterate, :pop!, :popfirst!)
    @eval begin
        Base.$func(ps::AbstractTerm, args...) = $func(op_string(ps), args...)
    end
end

for func in (:push!, :pushfirst!, :insert!)
    @eval begin
        Base.$func(ps::AbstractTerm, args...) = ($func(op_string(ps), args...); ps)
    end
end
