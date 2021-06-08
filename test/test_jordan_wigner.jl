import QuantumOps.FermiOps
import JLD2

@testset "Jordan-Wigner" begin
    s1 = PauliTerm("ZZZXIIII", complex(-1/2)) +  PauliTerm("ZZZYIIII", im/2)
    @test jordan_wigner(FermiOps.Raise, 4, 8) == s1

    s2 = PauliTerm("ZZZXIIII", complex(-1/2)) +  PauliTerm("ZZZYIIII", -im/2)
    @test jordan_wigner(FermiOps.Lower, 4, 8) == s2

    s3 = PauliTerm("IIIIIIII", 0.5 + 0.0im) +
        PauliTerm("IIIZIIII", -0.5 + 0.0im)
    @test jordan_wigner(FermiOps.NumberOp, 4, 8) == s3
    @test s3 == s1 * s2

    s4 = PauliTerm("IIIIIIII", 0.5 + 0.0im) +
        PauliTerm("IIIZIIII", 0.5 + 0.0im)
    @test jordan_wigner(FermiOps.Empty, 4, 8) == s4
    @test s4 == s2 * s1

    @test jordan_wigner_fermi(FermiOps.Raise, 4, 8) == FermiTerm("ZZZ+IIII", 1.0 + 0.0im)
    @test jordan_wigner_fermi(FermiOps.Lower, 4, 8) == FermiTerm("ZZZ-IIII", 1.0 + 0.0im)
    @test jordan_wigner_fermi(FermiOps.NumberOp, 4, 8) == FermiTerm("IIINIIII", 1.0 + 0.0im)

    h2_hamiltonian = FermiSum([FermiTerm("IIII", 0.7137539936876182),FermiTerm("IIIN", -0.47594871522096355),
                               FermiTerm("IINI", -0.47594871522096355),FermiTerm("IINN", 0.6973937674230275),
                               FermiTerm("INII", -1.2524635735648981),FermiTerm("ININ", 0.482179288212072),
                               FermiTerm("INNI", 0.663468096423568),FermiTerm("NIII", -1.2524635735648981),
                               FermiTerm("NIIN", 0.663468096423568),FermiTerm("NINI", 0.482179288212072),
                               FermiTerm("NNII", 0.6744887663568375),FermiTerm("++--", -0.181288808211496),
                               FermiTerm("+--+", 0.181288808211496),FermiTerm("-++-", 0.181288808211496),
                               FermiTerm("--++", -0.181288808211496)])

    data = JLD2.load("iop_h2.jld2")
    iop = data["iop"]
    @test OpSum{FermiOp}(iop) == h2_hamiltonian

    jw_h2_hamiltonian = PauliSum([PauliTerm("IIII", -0.09886396933545721 + 0.0im),PauliTerm("IIIZ", -0.2227859304041851 + 0.0im),PauliTerm("IIZI", -0.2227859304041851 + 0.0im),PauliTerm("IIZZ", 0.17434844185575687 + 0.0im),PauliTerm("IZII", 0.17119774903432972 + 0.0im),PauliTerm("IZIZ", 0.120544822053018 + 0.0im),PauliTerm("IZZI", 0.165867024105892 + 0.0im),PauliTerm("XXYY", -0.045322202052874 + 0.0im),PauliTerm("XYYX", 0.045322202052874 + 0.0im),PauliTerm("YXXY", 0.045322202052874 + 0.0im),PauliTerm("YYXX", -0.045322202052874 + 0.0im),PauliTerm("ZIII", 0.17119774903432966 + 0.0im),PauliTerm("ZIIZ", 0.165867024105892 + 0.0im),PauliTerm("ZIZI", 0.120544822053018 + 0.0im),PauliTerm("ZZII", 0.16862219158920938 + 0.0im)])

    @test jordan_wigner(h2_hamiltonian) == jw_h2_hamiltonian

    @test jordan_wigner(h2_hamiltonian, Pauli) == jordan_wigner(h2_hamiltonian, PauliI)
end
