@testset "matrix decompositions" begin
    @testset "ZYZ (Euler)" begin
        recompose(e) = cis(e[1]) * matrix(RZ(e[2])) * matrix(RY(e[3])) * matrix(RZ(e[4]))
        for U in (randu(2), randu(2), randu(2), matrix(Id()), matrix(X()), matrix(Y()),
                  matrix(Z()), matrix(H()), matrix(S()), matrix(T()), matrix(RX(0.9)),
                  matrix(RY(-2.2)), matrix(RZ(1.4)), ComplexF64[0 1; -1 0])
            @test recompose(zyz(U)) ≈ U
        end
        @test_throws ArgumentError zyz(randn(3, 3))
        @test_throws ArgumentError zyz([1 1; 0 1])

        for U in (randu(2), matrix(H()), matrix(T()))
            @test matrix(decompose_1q(U)) ≈ U                 # global phase tracked
        end
        c = Circuit(3); decompose_1q!(c, matrix(H()), 2)
        @test matrix(c) ≈ embed(matrix(H()), [2], 3)
    end

    @testset "two-level (Givens) decomposition" begin
        for n in 1:3
            N = 1 << n
            U = randu(N)
            fs = two_level_decompose(U)
            @test length(fs) <= N * (N - 1) ÷ 2 + N
            @test all(is_unitary(t.V) for t in fs)
            @test all(t.a < t.b for t in fs)
            @test reduce(*, matrix(t, N) for t in fs) ≈ U     # factors multiply back
            @test matrix(synthesize_unitary(U)) ≈ U           # ... and synthesise back
        end
        # structured and degenerate inputs
        for U in (Matrix{ComplexF64}(I, 8, 8), kron(matrix(CZ()), matrix(Id())),
                  kron(matrix(SWAP()), matrix(Id())), kron_n(matrix(H()), 3),
                  matrix(diagonal([0.0, 0.3, 1.1, -0.7])))
            @test matrix(synthesize_unitary(U)) ≈ U
        end
        @test_throws ArgumentError two_level_decompose(randn(4, 4))
        @test_throws ArgumentError synthesize_unitary(randu(6))
    end

    @testset "two-level routing" begin
        # a two-level gate must be the identity outside span(|a⟩,|b⟩)
        V = randu(2)
        for (a, b) in ((0, 7), (1, 6), (2, 3), (0, 4), (5, 7))
            c = Circuit(3); two_level!(c, TwoLevel(a, b, V), 3)
            @test matrix(c) ≈ matrix(TwoLevel(a, b, V), 8)
        end
    end

    @testset "demultiplexing" begin
        for m in (2, 4, 8)
            U1, U2 = randu(m), randu(m)
            V, θ, W = demultiplex(U1, U2)
            @test is_unitary(V) && is_unitary(W)
            mid = Diagonal(vcat(cis.(-θ ./ 2), cis.(θ ./ 2)))     # D ⊕ D†
            @test kron(I(2), V) * mid * kron(I(2), W) ≈ cat(U1, U2; dims=(1, 2))
        end
        U = randu(4)                                              # degenerate: U₁ = U₂
        V, θ, W = demultiplex(U, U)
        @test kron(I(2), V) * Diagonal(vcat(cis.(-θ ./ 2), cis.(θ ./ 2))) * kron(I(2), W) ≈
              cat(U, U; dims=(1, 2))
    end

    @testset "uniformly controlled one-qubit gates" begin
        for (controls, target, n) in (([], 1, 1), ([1], 2, 2), ([2], 1, 2),
                                      ([1, 2], 3, 3), ([3, 1], 2, 3), ([1, 2, 3], 4, 4))
            k = length(controls)
            Us = [randu(2) for _ in 1:(1 << k)]
            c = multiplexed_1q(Us, controls, target; n=n)
            @test matrix(c) ≈ ref_multiplexed_gate(Us, controls, target, n)
            @test count_cnots(c) == (k == 0 ? 0 : 4 * (1 << k) - 2)
        end
        Us = [matrix(Id()), matrix(X())]
        @test matrix(multiplexed_1q(Us, [1], 2)) ≈ matrix(CNOT())
        @test_throws ArgumentError multiplexed_1q([randu(2)], [1], 2)
    end
end
