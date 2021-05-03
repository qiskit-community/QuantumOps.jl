using PauliStrings
using Test

@testset "Pauli" begin
    @test Pauli(0) == Pauli(:I)
    @test Pauli(1) == Pauli(:X)
    @test Pauli(2) == Pauli(:Y)
    @test Pauli(3) == Pauli(:Z)
    @test Pauli(0) == Pauli('I')
    @test Pauli(1) == Pauli('X')
    @test Pauli(2) == Pauli('Y')
    @test Pauli(3) == Pauli('Z')
    @test Pauli(1) == Pauli("X")

    for (a, b, c) in ((:I, :I, :I), (:I, :X, :X), (:I, :Y, :Y), (:I, :Z, :Z),
                      (:X, :I, :X), (:X, :X, :I), (:X, :Y, :Z), (:X, :Z, :Y),
                      (:Y, :I, :Y), (:Y, :X, :Z), (:Y, :Y, :I), (:Y, :Z, :X),
                      (:Z, :I, :Z), (:Z, :X, :Y), (:Z, :Y, :X), (:Z, :Z, :I))
        @test Pauli(a) * Pauli(b) == Pauli(c)
    end
end
