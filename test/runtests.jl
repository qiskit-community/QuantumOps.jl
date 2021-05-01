using PauliStrings
using Test

@testset "Pauli" begin
    p = Pauli(0)
    @test p == Pauli(:I)
end
