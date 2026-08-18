@testset "gates" begin
    gs = [Id(), X(), Y(), Z(), H(), S(), Sdg(), T(), Tdg(),
          RX(0.7), RY(-1.3), RZ(2.1), PHASE(0.4), CNOT(), CZ(), SWAP()]
    for g in gs
        m = matrix(g)
        @test m * m' ≈ I(size(m, 1))
        @test size(m, 1) == 1 << nqubits(g)
    end

    # the identity the Gray-code multiplexor is built on
    Xm = matrix(X())
    for θ in (0.0, 0.3, -2.2, pi)
        @test Xm * matrix(RZ(θ)) * Xm ≈ matrix(RZ(-θ))
        @test Xm * matrix(RY(θ)) * Xm ≈ matrix(RY(-θ))
    end

    @test matrix(controlled(X())) ≈ matrix(CNOT())
    @test matrix(controlled(Z())) ≈ matrix(CZ())
    @test size(matrix(controlled(X(), 2))) == (8, 8)
    @test label(RZ(pi/2)) == "RZ(1.57)"
end
