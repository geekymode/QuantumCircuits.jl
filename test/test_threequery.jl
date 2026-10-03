@testset "three-query synthesis" begin
    X2 = [0.0 1.0; 1.0 0.0]

    @testset "components" begin
        for K in (2, 4, 7, 8, 16)
            F = centered_dft(K)
            @test F' * F ≈ I(K) && F ≈ transpose(F)
            # Eq. 47: the inverse QFT between two diagonal phases, exactly
            QFTdg = [cis(-2π * u * w / K) for w in 0:K-1, u in 0:K-1] ./ sqrt(K)
            D1 = Diagonal([cis(2π / K * (K - 1) / 2 * (u - (K - 1) / 4)) for u in 0:K-1])
            @test F ≈ D1 * QFTdg * D1
        end
        E = qho_encoding(3, 8)
        @test size(E) == (8^3, 3) && E' * E ≈ I(3)
        # one-hot structure: column x is h0 ⊗ … ⊗ h1 (register x) ⊗ … ⊗ h0
        h = qho_encoding(1, 8)[:, 1]                       # h1 alone
        h0 = abs.(qho_encoding(2, 8)[1:8:end, 2])          # h0, read off register 1
        @test E[:, 1] ≈ kron(h, normalize(h0), normalize(h0))
        q = chirp_phases(X2, 8)
        @test length(q) == 64 && all(z -> abs(z) ≈ 1, q)
        @test chirp_phases(zeros(2, 2), 4) ≈ ones(16)
        @test_throws ArgumentError centered_dft(1)
        @test_throws ArgumentError qho_encoding(8, 16)    # 16^8 > maxdim
    end

    @testset "Theorem 4.1: exponential convergence" begin
        # Pauli X, and the 4×4 embedding of a random phase
        S4 = symmetric_embedding(randu(1))
        for S in (X2, S4)
            errs = [three_query_error(S, K) for K in (8, 16, 32)]
            @test errs[1] < 1e-2
            @test errs[3] < 1e-9
            # at least the paper's rate e^{-πK/8}; observed is about e^{-πK/4}
            @test errs[1] / errs[2] > exp(π * 8 / 8)
            @test errs[2] / errs[3] > exp(π * 16 / 8)
            @test errs[2] / errs[3] > exp(0.6 * 16)
        end
        # a 6×6 real involution
        Q = Matrix(qr(randn(3, 3)).Q)
        @test three_query_error(hermitian_dilation(Q), 10) < 2e-3
        # the result is real up to the discretisation error
        M = three_query(X2, 32)
        @test maximum(abs, imag(M)) < 1e-9
    end

    @testset "the trace condition is necessary" begin
        Q = Matrix(qr(randn(3, 3)).Q)
        St = Symmetric(Q * Diagonal([1.0, -1.0, 1.0]) * Q')   # trace 1
        @test_throws ArgumentError three_query(St, 8)
        e16 = three_query_error(St, 16; check=false)
        e24 = three_query_error(St, 24; check=false)
        @test e16 > 0.5 && abs(e16 - e24) < 1e-3              # does not converge
        @test_throws ArgumentError three_query([1.0 0.0; 0.0 -0.5], 8)    # not unitary
        @test_throws ArgumentError three_query([0.0 1.0; -1.0 0.0], 8)    # not symmetric
        @test_throws ArgumentError three_query(ComplexF64[0 im; im 0], 8) # not real
    end

    @testset "full pipeline" begin
        U = randu(2)
        # a one-qubit U already means 8 registers, so only coarse grids fit
        e5 = opnorm(three_query_synthesis(U, 5) - U)
        e6 = opnorm(three_query_synthesis(U, 6) - U)
        @test e6 < e5 < 0.1
        @test opnorm(three_query_synthesis(cis(0.7) * ones(1, 1), 24) - cis(0.7) * ones(1, 1)) < 1e-6
        @test_throws ArgumentError three_query_synthesis(U, 16)  # 16^8 registers
    end
end
