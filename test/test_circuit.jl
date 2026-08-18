@testset "circuit" begin
    @testset "big-endian convention" begin
        c = Circuit(2); push!(c, X(), 1)
        @test statevector(c) ≈ ComplexF64[0, 0, 1, 0]   # |10⟩ has index 2
        c = Circuit(2); push!(c, X(), 2)
        @test statevector(c) ≈ ComplexF64[0, 1, 0, 0]   # |01⟩ has index 1
    end

    @testset "bell state" begin
        c = Circuit(2); push!(c, H(), 1); push!(c, CNOT(), 1, 2)
        @test statevector(c) ≈ ComplexF64[1, 0, 0, 1] ./ sqrt(2)
    end

    @testset "reversed CNOT wires" begin
        c = Circuit(2); push!(c, CNOT(), 2, 1)
        @test matrix(c) ≈ ComplexF64[1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]
    end

    @testset "matrix vs statevector" begin
        c = Circuit(3)
        push!(c, H(), 2); push!(c, CNOT(), 2, 3); push!(c, RY(0.8), 1)
        push!(c, SWAP(), 1, 3); push!(c, RZ(-0.4), 2)
        U = matrix(c)
        @test U * U' ≈ I(8)
        ψ0 = normalize!(randn(ComplexF64, 8))
        @test statevector(c, ψ0) ≈ U * ψ0
    end

    @testset "global phase" begin
        c = Circuit(1); c.global_phase = pi/3
        @test matrix(c) ≈ cis(pi/3) * I(2)
    end

    @testset "bookkeeping and errors" begin
        c = Circuit(2); push!(c, CNOT(), 1, 2); push!(c, CNOT(), 2, 1); push!(c, H(), 1)
        @test count_cnots(c) == 2
        @test count_gates(c, :H) == 1
        @test length(c) == 3
        @test_throws ArgumentError push!(c, CNOT(), 1)      # wrong arity
        @test_throws ArgumentError push!(c, CNOT(), 1, 1)   # repeated wire
        @test_throws ArgumentError push!(c, H(), 3)         # out of range
        @test sprint(draw, c) isa String
    end
end
