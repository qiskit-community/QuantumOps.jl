@testset "IndOp" begin
    op = IndOp(Paulis.X, 1)
    @test copy(op) == op
end
