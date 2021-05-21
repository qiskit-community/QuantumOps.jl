####
#### FermiOp
####

struct FermiOp
    ind::UInt8
end

op_index(op::FermiOp) = op.ind

Base.isless(op1::FermiOp, op2::FermiOp) = isless(op_index(op1), op_index(op2))

const id_op = FermiOp(0)
const number_op = FermiOp(1)
const empty_op = FermiOp(2)
const raise_op = FermiOp(3)
const lower_op = FermiOp(4)
const no_op = FermiOp(5)

function Base.show(io::IO, op::FermiOp)
    _ferm_chars = ('I', 'N', 'I', '+', '-', 'X')
    print(io, _ferm_chars[op_index(op) + 1])
end
