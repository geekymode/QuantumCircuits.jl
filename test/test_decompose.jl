# Reference implementation built straight from the *definition* of a
# uniformly controlled rotation, with no Gray code anywhere — this is what the
# synthesised circuit is checked against.
function ref_multiplexed(R, α, controls, target, n)
    N = 1 << n
    k = length(controls)
    U = zeros(ComplexF64, N, N)
    tshift = n - target
    for col in 0:N-1
        j = 0
        for (i, q) in enumerate(controls)
            j |= ((col >> (n - q)) & 1) << (k - i)
        end
        m = matrix(R(α[j+1]))
        b = (col >> tshift) & 1
        for bp in 0:1
            row = (col & ~(1 << tshift)) | (bp << tshift)
            U[row+1, col+1] = m[bp+1, b+1]
        end
    end
    U
end

@testset "decompositions" begin
    @testset "multiplex transform" begin
        for k in 0:5
            M = multiplex_matrix(k)
            @test M * M' ≈ (1 << k) * I(1 << k)          # orthogonal, Walsh-like
            α = randn(1 << k)
            @test M * multiplex_angles(α) ≈ α            # solves the system
        end
        @test_throws ArgumentError multiplex_angles(randn(3))
    end

    @testset "uniformly controlled rotations" begin
        for (kind, R) in ((:RY, RY), (:RZ, RZ))
            for (controls, target, n) in (([], 1, 1),
                                          ([1], 2, 2),
                                          ([2], 1, 2),
                                          ([1, 2], 3, 3),
                                          ([3, 1], 2, 3),
                                          ([4, 2, 1], 3, 4),
                                          ([1, 2, 3, 4], 5, 5))
                k = length(controls)
                α = randn(1 << k)
                c = Circuit(n)
                multiplexed_rotation!(c, kind, α, controls, target)
                @test matrix(c) ≈ ref_multiplexed(R, α, controls, target, n)
                @test count_cnots(c) == (k == 0 ? 0 : 1 << k)     # 2^k CNOTs, flat in k
                @test count_gates(c, kind) == 1 << k
            end
        end
        @test matrix(multiplexed_rz([0.3, -1.1], [1], 2)) ≈
              ref_multiplexed(RZ, [0.3, -1.1], [1], 2, 2)
        @test matrix(multiplexed_ry([0.3, -1.1], [2], 1)) ≈
              ref_multiplexed(RY, [0.3, -1.1], [2], 1, 2)
    end

    @testset "multiplexor errors" begin
        c = Circuit(3)
        @test_throws ArgumentError multiplexed_rotation!(c, :RX, [0.1, 0.2], [1], 2)
        @test_throws ArgumentError multiplexed_rotation!(c, :RZ, [0.1], [1], 2)
        @test_throws ArgumentError multiplexed_rotation!(c, :RZ, [0.1, 0.2], [1], 1)
    end

    @testset "diagonal unitaries" begin
        for n in 1:5
            φ = randn(1 << n)
            c = diagonal(φ)
            @test nqubits(c) == n
            @test matrix(c) ≈ Diagonal(cis.(φ))          # exact, incl. global phase
            @test count_cnots(c) == (1 << n) - 2         # 2^n - 2
        end
        @test matrix(diagonal([0.0, pi])) ≈ matrix(Z())
        @test matrix(diagonal([0.0, 0.0, 0.0, pi])) ≈ matrix(CZ())
        @test_throws ArgumentError diagonal(randn(6))
    end

    @testset "state preparation" begin
        for n in 1:5
            a = randn(ComplexF64, 1 << n)
            c = prepare_state(a)
            @test statevector(c) ≈ a ./ norm(a)          # exact, no phase slack
            @test count_cnots(c) == 2 * ((1 << n) - 2)
        end

        # sparse / degenerate inputs must not produce NaNs
        for a in (ComplexF64[1, 0, 0, 0], ComplexF64[0, 0, 0, 1],
                  ComplexF64[0, 1, 0, 0], ComplexF64[1, 0, 0, im],
                  ComplexF64[0, 0, 1, 1])
            c = prepare_state(a)
            ψ = statevector(c)
            @test all(isfinite, ψ)
            @test ψ ≈ a ./ norm(a)
        end

        # real, positive amplitudes need no phase-fixing beyond the diagonal
        a = abs.(randn(8))
        @test statevector(prepare_state(a)) ≈ a ./ norm(a)

        @test_throws ArgumentError prepare_state(zeros(ComplexF64, 4))
        @test_throws ArgumentError prepare_state(randn(ComplexF64, 3))
    end
end
