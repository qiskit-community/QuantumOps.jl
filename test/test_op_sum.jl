@testset "OpSum" begin
    for _type = (Pauli, PauliI, FermiOp)
        m = rand(_type, 2, 3)
        @test size(OpSum(m)) == (2, 3)
    end

    for term in (PauliTerm("XYZ"), FermiTerm("--++"))
        @test reverse(term) != term
        @test reverse(reverse(term)) == term
        term2 = reverse(term)
        reverse!(term2)
        @test term2 == term
    end
end
