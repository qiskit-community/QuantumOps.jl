@testset "FermiOp" begin
    @test FermiOp(0) === QuantumOps.FermiOps.I
    @test FermiOp(1) === QuantumOps.FermiOps.NumberOp
end

@testset "FermiTerm" begin
    t = FermiTerm("++-")
    @test reverse(t) != t
    @test reverse(reverse(t)) == t

    t1 = FermiTerm("+")
    t2 = FermiTerm("-")
    @test t1 * t2 == FermiTerm("N")

    t1 = FermiTerm("++")
    t2 = FermiTerm("--")
    @test t1 * t2 == FermiTerm("NN", -1)

    t1 = FermiTerm("+++")
    t2 = FermiTerm("---")
    @test t1 * t2 == FermiTerm("NNN", -1)

    t1 = FermiTerm("++++")
    t2 = FermiTerm("----")
    @test t1 * t2 == FermiTerm("NNNN", 1)

    t1 = FermiTerm("+E++")
    t2 = FermiTerm("-E--")
    @test t1 * t2 == FermiTerm("NENN", -1)
end
