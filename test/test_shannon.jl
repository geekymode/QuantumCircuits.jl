@testset "cosine–sine and quantum Shannon decomposition" begin
    @testset "rand_unitary" begin
        for N in (2, 4, 8, 16)
            @test is_unitary(rand_unitary(N))
        end
    end

    @testset "cosine–sine, generic inputs" begin
        for N in (2, 4, 8, 16, 32), _ in 1:3
            U = rand_unitary(N)
            F = cosine_sine(U)
            @test Matrix(F) ≈ U
            @test F.c .^ 2 .+ F.s .^ 2 ≈ ones(N ÷ 2)
            @test all(is_unitary(M; atol=1e-8) for M in (F.L1, F.L2, F.R1, F.R2))
            @test all(0 .<= F.c .<= 1) && all(0 .<= F.s .<= 1)
        end
        @test_throws ArgumentError cosine_sine(randn(4, 4))
        @test_throws ArgumentError cosine_sine(rand_unitary(3))
    end

    @testset "cosine–sine, degenerate inputs" begin
        # every one of these has sines that are exactly zero — the case where
        # sqrt(1 - c^2) silently returns √eps instead of 0
        cases = ["identity" => Matrix{ComplexF64}(I, 8, 8),
                 "block diagonal" => cat(rand_unitary(4), rand_unitary(4); dims=(1, 2)),
                 "anti-diagonal" => kron(matrix(X()), Matrix{ComplexF64}(I, 4, 4)),
                 "CNOT" => matrix(CNOT()),
                 "CZ" => matrix(CZ()),
                 "SWAP" => matrix(SWAP()),
                 "diagonal" => Matrix(Diagonal(cis.(randn(8)))),
                 "permutation" => Matrix{ComplexF64}(I, 8, 8)[:, [3, 1, 2, 5, 4, 8, 6, 7]],
                 "repeated singular values" => kron(matrix(RY(0.7)), Matrix{ComplexF64}(I, 4, 4)),
                 "Hadamards" => kron_n(matrix(H()), 3)]
        for (name, U) in cases
            F = cosine_sine(U)
            @test Matrix(F) ≈ U
            @test all(is_unitary(M; atol=1e-8) for M in (F.L1, F.L2, F.R1, F.R2))
            @test F.c .^ 2 .+ F.s .^ 2 ≈ ones(size(U, 1) ÷ 2)
        end
        # a block-diagonal input must produce all-zero sines, not √eps
        F = cosine_sine(cat(rand_unitary(4), rand_unitary(4); dims=(1, 2)))
        @test all(F.s .< 1e-12)
    end

    @testset "the middle factor is a Gray-code multiplexor" begin
        for n in 2:4
            U = rand_unitary(1 << n)
            F = cosine_sine(U)
            C, S = Diagonal(F.c), Diagonal(F.s)
            mid = ComplexF64[C -S; S C]
            # ... exactly a uniformly controlled RY on the top wire
            uc = multiplexed_ry(csd_angles(F), 2:n, 1; n=n)
            @test matrix(uc) ≈ mid
            @test count_cnots(uc) == 1 << (n - 1)
        end
    end

    @testset "quantum Shannon decomposition" begin
        for n in 1:4, _ in 1:2
            U = rand_unitary(1 << n)
            c = qsd(U)
            @test matrix(c) ≈ U                                  # exact, phase included
            @test count_cnots(c) == qsd_cnot_count(n)            # (3/4)4ⁿ - (3/2)2ⁿ
            @test all(length(op.qubits) <= 2 for op in c.ops)    # nothing but 1q + CNOT
        end
        @test qsd_cnot_count.(1:5) == [0, 6, 36, 168, 720]

        for U in (Matrix{ComplexF64}(I, 8, 8), matrix(controlled(X(), 2)),
                  kron(matrix(CNOT()), Matrix{ComplexF64}(I, 2, 2)),
                  kron(matrix(SWAP()), Matrix{ComplexF64}(I, 2, 2)),
                  kron_n(matrix(H()), 3), Matrix(Diagonal(cis.(randn(8)))),
                  [cis(2π * i * j / 8) / sqrt(8) for i in 0:7, j in 0:7])
            @test matrix(qsd(U)) ≈ U
        end

        # beats the two-level route by a wide margin at n = 3
        U = rand_unitary(8)
        @test count_cnots(qsd(U)) < count_cnots(synthesize_unitary(U; lower=true))

        @test_throws ArgumentError qsd(rand_unitary(6))
        @test_throws ArgumentError qsd!(Circuit(3), rand_unitary(4), [1, 2, 3])
    end
end
