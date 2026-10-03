@testset "embeddings in real symmetric involutions" begin
    minusY = ComplexF64[1, -im] ./ sqrt(2)

    @testset "realify" begin
        for N in (1, 2, 3, 4)
            U, V = randu(N), randu(N)
            R = realify(U)
            @test eltype(R) <: Real && size(R) == (2N, 2N)
            @test R' * R ≈ I(2N)                              # orthogonal
            @test realify(U * V) ≈ R * realify(V)             # homomorphism
            @test realify(U') ≈ R'
            ψ = normalize!(randn(ComplexF64, N))
            @test R * kron(minusY, ψ) ≈ kron(minusY, U * ψ)  # acts as U on |−Y⟩
        end
        @test statevector(push!(Circuit(1), RX(π / 2), 1)) ≈ minusY
    end

    @testset "hermitian_dilation" begin
        for N in (1, 2, 4)
            R = randu(N)
            S = hermitian_dilation(R)
            @test S ≈ S' && S * S ≈ I(2N) && abs(tr(S)) < 1e-12
            ψ = normalize!(randn(ComplexF64, N))
            @test S * kron([0, 1], ψ) ≈ kron([1, 0], R * ψ)  # |1⟩ψ ↦ |0⟩Rψ
            Q = Matrix(qr(randn(N, N)).Q)
            @test eltype(hermitian_dilation(Q)) <: Real
            @test hermitian_dilation(Q) == transpose(hermitian_dilation(Q))
        end
        @test_throws ArgumentError hermitian_dilation(zeros(2, 3))
    end

    @testset "symmetric_embedding" begin
        for N in (1, 2, 3, 8)
            S = symmetric_embedding(randu(N))
            @test eltype(S) <: Real && size(S) == (4N, 4N)
            @test S == transpose(S)
            @test S * S ≈ I(4N)
            @test abs(tr(S)) < 1e-12
            λ = eigvals(Symmetric(S))                        # half +1, half -1
            @test count(>(0), λ) == count(<(0), λ) == 2N
            @test all(x -> isapprox(abs(x), 1; atol=1e-10), λ)
        end
    end

    @testset "embedding_circuit" begin
        for n in 0:3
            U = randu(1 << n)
            c = embedding_circuit(U)
            @test nqubits(c) == n + 2 && ancillas(c) == [n + 1, n + 2]
            @test implementation_error(c, U) < 1e-12
            @test is_clean(c)
            @test length(c) == 4                              # S plus three 1-qubit gates
            @test count(op -> length(op.qubits) == 1, c.ops) == 3
        end
        # splice in a synthesised S: anything that builds S builds U
        U = randu(4)
        sc = qsd(symmetric_embedding(U))
        sc.global_phase += 0.3                               # phase must carry through
        c = embedding_circuit(sc)
        @test implementation_error(c, cis(0.3) * U) < 1e-10
        @test count_cnots(c) == count_cnots(sc)
        # without the |−Y⟩ preparation the ancilla is not clean
        bad = embedding_circuit(randu(2)); deleteat!(bad.ops, 2); pop!(bad.ops)
        @test !is_clean(bad)
        @test_throws ArgumentError embedding_circuit(randu(3))
        @test_throws ArgumentError embedding_circuit(Circuit(1))
    end
end
