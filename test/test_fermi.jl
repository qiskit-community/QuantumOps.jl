@testset "FermiOp" begin
    @test FermiOp(0) === QuantumOps.FermiOps.id_op
    @test FermiOp(1) === QuantumOps.FermiOps.number_op
end
