struct IndOp{OpT} <: AbstractOp
    op::OpT
    ind::Int
end

#function OpTerm{IndOp{
