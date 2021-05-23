@testset "FermiOp" begin
    @test FermiOp(0) === QuantumOps.FermiOps.I_op
    @test FermiOp(1) === QuantumOps.FermiOps.number_op
end

@testset "FermiTerm" begin
    t = FermiTerm("++-")
    @test reverse(t) != t
    @test reverse(reverse(t)) == t
end
