@testset "applications of Gray coding" begin
    @testset "matrix roots" begin
        for m in 0:3, U in (randu(2), randu(4), matrix(X()), matrix(H()), Matrix{ComplexF64}(I, 2, 2))
            V = matrix_root(U, m)
            @test is_unitary(V)
            @test V^(1 << m) ≈ U
        end
        @test_throws ArgumentError matrix_root(randu(2), -1)
    end

    @testset "multi-controlled gates" begin
        for k in 0:4
            U = randu(2)
            c = multicontrolled(U, 1:k, k + 1)
            @test matrix(c) ≈ (k == 0 ? U : matrix(controlled(Gate(:U, U), k)))
            if k >= 2                                     # the Gray-code construction
                @test count_cnots(c) == (1 << k) - 2
                @test length(c) - count_cnots(c) == (1 << k) - 1   # controlled-V gates
            end
            @test all(length(op.qubits) <= 2 for op in c.ops)      # elementary throughout
        end
        # named gates come out right
        @test matrix(multicontrolled(matrix(X()), [1, 2], 3)) ≈ matrix(controlled(X(), 2))
        @test matrix(multicontrolled(matrix(X()), [1], 2)) ≈ matrix(CNOT())
        @test matrix(multicontrolled(matrix(Z()), [1], 2)) ≈ matrix(CZ())
        # scattered / reversed wires
        c = Circuit(4); multicontrolled!(c, matrix(X()), [4, 1], 3)
        @test matrix(c) ≈ embed(matrix(controlled(X(), 2)), [4, 1, 3], 4)
        @test_throws ArgumentError multicontrolled(randu(4), [1], 2)
        @test_throws ArgumentError multicontrolled(randu(2), [1, 2], 2)
    end

    @testset "the code as a circuit" begin
        for n in 2:6
            E, D = gray_encoder(n), gray_decoder(n)
            @test count_cnots(E) == n - 1
            @test count_cnots(D) == n - 1
            @test length(E) == n - 1                       # nothing but CNOTs
            ME = matrix(E)
            for x in 0:(1 << n)-1                          # |x⟩ ↦ |gray(x)⟩
                @test ME[gray(x)+1, x+1] ≈ 1
            end
            @test matrix(D) * ME ≈ I(1 << n)               # exact inverses
        end
    end

    @testset "counters" begin
        for k in 2:4
            N = 1 << k
            M = matrix(increment(k))
            for x in 0:N-1
                @test M[(x + 1) % N + 1, x+1] ≈ 1
            end
            @test is_unitary(M)
        end
        for n in 2:4
            N = 1 << n
            M = matrix(gray_increment(n))
            for i in 0:N-1                                 # |gray(i)⟩ ↦ |gray(i+1)⟩
                @test M[gray((i + 1) % N)+1, gray(i)+1] ≈ 1
            end
            # every step of the counter changes exactly one bit
            for x in 0:N-1
                y = argmax(abs.(M[:, x+1])) - 1
                @test hamming(x, y) == 1
            end
        end
    end

    @testset "SELECT (LCU / QROM addressing)" begin
        for k in 1:3, m in 1:2
            Us = [randu(1 << m) for _ in 1:(1 << k)]
            ref = cat(Us...; dims=(1, 2))                  # controls are the leading wires
            cg = select(Us, 1:k, k+1:k+m; order=:gray)
            cn = select(Us, 1:k, k+1:k+m; order=:natural)
            @test matrix(cg) ≈ ref
            @test matrix(cn) ≈ ref
            @test count_gates(cg, :X) <= count_gates(cn, :X)
            @test count_gates(cg, Symbol("C"^k, "U")) == 1 << k   # one branch each way
        end
        # Gray addressing strictly wins once there is something to win
        for k in 3:6
            Us = [randu(2) for _ in 1:(1 << k)]
            @test count_gates(select(Us, 1:k, [k+1]; order=:gray), :X) <
                  count_gates(select(Us, 1:k, [k+1]; order=:natural), :X)
        end
        @test_throws ArgumentError select([randu(2)], [1, 2], [3])
        @test_throws ArgumentError select([randu(2), randu(2)], [1], [1])
        @test_throws ArgumentError select([randu(2), randu(2)], [1], [2]; order=:nope)
    end

    @testset "term ordering utilities" begin
        @test support_mask("III") == 0
        @test support_mask("ZII") == 0b100
        @test support_mask("IXY") == 0b011
        terms = ["IIZ" => 1.0, "ZZI" => 2.0, "IZZ" => 3.0, "ZII" => 4.0]
        ordered = gray_order(terms)
        @test sort(first.(ordered)) == sort(first.(terms))        # a permutation
        @test issorted([ungray(support_mask(first(t))) for t in ordered])
    end

    @testset "Walsh-series truncation" begin
        n = 6
        N = 1 << n
        φ = [2.0 * sin(2π * x / N) + 0.4 * (x / N)^2 for x in 0:N-1]
        pp = phase_polynomial(φ)
        prev = Inf
        for k in (2, 4, 8, 16, 32)
            tp = truncate_terms(pp, k)
            @test nterms(tp) <= k
            err = maximum(abs, phases(tp) .- φ)
            @test err <= prev + 1e-12                              # monotone improvement
            prev = err
            @test matrix(synthesize(tp)) ≈ Diagonal(cis.(phases(tp)))
        end
        @test phases(truncate_terms(pp, N)) ≈ φ                    # keeping everything is exact
        @test count_cnots(synthesize(truncate_terms(pp, 8))) < count_cnots(synthesize(pp))
    end

    @testset "lowering multi-controlled gates" begin
        for n in 1:3
            U = randu(1 << n)
            lo = synthesize_unitary(U; lower=true)
            @test matrix(lo) ≈ U
            @test all(length(op.qubits) <= 2 for op in lo.ops)
        end
        c2 = Circuit(3)
        V = randu(2)
        two_level!(c2, TwoLevel(1, 6, V), 3; lower=true)
        @test matrix(c2) ≈ matrix(TwoLevel(1, 6, V), 8)
        @test all(length(op.qubits) <= 2 for op in c2.ops)
    end
end
