import ..FermiOps: count_raise_lower

for f in (:count_raise_lower, )
    @eval $f(op::IndOp{FermiOp}) = $f(op.op)
end

do_phase_count(f1::IndOp{<:AbstractFermiOp}, f2::IndOp{<:AbstractFermiOp}) = do_phase_count(f1.op, f2.op)
