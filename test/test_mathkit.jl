@testset "linear algebra toolkit" begin
    @testset "Walsh–Hadamard transform" begin
        for n in 0:6
            v = randn(1 << n)
            @test fwht(v) ≈ walsh_matrix(n) * v
            @test fwht(fwht(v)) ≈ (1 << n) .* v          # W² = N·I
        end
        a = [1.0, 2, 3, 4]
        @test fwht!(a) === a                             # mutates in place
        @test a == [10.0, -2, -4, 0]
        @test_throws ArgumentError fwht(randn(6))
    end

    @testset "Pauli basis" begin
        @test pauli("I") ≈ matrix(Id())
        @test pauli("X") ≈ matrix(X())
        @test pauli("XZ") ≈ kron(matrix(X()), matrix(Z()))
        @test length(pauli_strings(3)) == 64
        @test allunique(pauli_strings(3))
        for n in 1:3                                     # orthogonality: tr(P Q) = 2ⁿ δ
            ps = pauli_strings(n)
            @test all(abs(tr(pauli(p) * pauli(q))) ≈ (p == q ? 1 << n : 0)
                      for p in ps[1:4], q in ps[1:4])
        end
        for n in 1:3
            M = randn(ComplexF64, 1 << n, 1 << n)
            @test pauli_recompose(pauli_decompose(M)) ≈ M
        end
        H = randn(ComplexF64, 8, 8); H = H + H'
        @test all(isapprox(imag(c), 0; atol=1e-10) for (_, c) in pauli_decompose(H))
        @test pauli_decompose(matrix(Z())) == ["Z" => 1.0 + 0im]
        @test_throws ArgumentError pauli("Q")
    end

    @testset "embedding and predicates" begin
        @test embed(matrix(X()), [1], 2) ≈ kron(matrix(X()), I(2))
        @test embed(matrix(X()), [2], 2) ≈ kron(I(2), matrix(X()))
        @test embed(matrix(CNOT()), [1, 2], 3) ≈ kron(matrix(CNOT()), I(2))
        # reversed wires: control on qubit 2, target on qubit 1
        @test embed(matrix(CNOT()), [2, 1], 2) ≈ ComplexF64[1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]
        @test kron_n(matrix(H()), 3) ≈ kron(matrix(H()), matrix(H()), matrix(H()))

        @test is_unitary(matrix(H()))
        @test !is_unitary([1 1; 0 1])
        @test gate_fidelity(matrix(Z()), matrix(Z())) ≈ 1
        @test gate_fidelity(matrix(Z()), cis(0.7) .* matrix(Z())) ≈ 1   # phase-blind
        @test gate_fidelity(matrix(Z()), matrix(X())) ≈ 0 atol = 1e-12
        @test global_phase_between(cis(0.7) .* matrix(Z()), matrix(Z())) ≈ 0.7
    end

    @testset "Schmidt values and entropy" begin
        bell = ComplexF64[1, 0, 0, 1] ./ sqrt(2)
        @test schmidt_values(bell, 1) ≈ [1, 1] ./ sqrt(2)
        @test entanglement_entropy(bell, 1) ≈ 1

        prod2 = kron(ComplexF64[1, 0], normalize(randn(ComplexF64, 2)))
        @test entanglement_entropy(prod2, 1) ≈ 0 atol = 1e-12

        ghz = zeros(ComplexF64, 8); ghz[1] = ghz[8] = 1/sqrt(2)
        @test entanglement_entropy(ghz, 1) ≈ 1
        @test entanglement_entropy(ghz, 2) ≈ 1

        ψ = normalize(randn(ComplexF64, 16))
        @test sum(schmidt_values(ψ, 2) .^ 2) ≈ 1
        @test 0 <= entanglement_entropy(ψ, 2) <= 2
        @test_throws ArgumentError schmidt_values(ψ, 4)
    end
end
