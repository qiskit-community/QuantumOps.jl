import ..Utils

Base.:*(z::PyCall.PyObject, ps::OpTerm) = typeof(ps)(ps.paulis, ps.coeff * z)
Base.:*(ps::OpTerm, z::PyCall.PyObject) = *(z, ps)

function Utils.isapprox_zero(x::PyCall.PyObject)
    return x == 0
end
