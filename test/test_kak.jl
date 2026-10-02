@testset "two-qubit KAK decomposition" begin
    can(a, b, c) = exp(im * (a * pauli("XX") + b * pauli("YY") + c * pauli("ZZ")))
    loc() = kron(randu(2), randu(2))
    inchamber(k) = all(x -> -π / 4 - 1e-12 < x <= π / 4 + 1e-12, (k.a, k.b, k.c))

    @testset "canonical gate" begin
        for (a, b, c) in ((0.3, -0.2, 0.7), (0.0, 0.0, 0.0), (π / 4, π / 4, π / 4), (1.3, -2.0, 0.9))
            @test canonical_gate(a, b, c) ≈ can(a, b, c)
        end
        @test is_unitary(canonical_gate(0.1, 0.2, 0.3))
    end

    @testset "decomposition" begin
        for _ in 1:200
            U = randu(4)
            k = kak(U)
            @test Matrix(k) ≈ U                                 # exact, phase included
            @test inchamber(k)
            for A in (k.A1, k.A2, k.B1, k.B2)
                @test is_unitary(A)
            end
        end
        # degenerate and structured inputs, where eigenvalues coincide
        for U in (Matrix{ComplexF64}(I, 4, 4), matrix(CNOT()), matrix(CZ()), matrix(SWAP()),
                  loc(), matrix(controlled(RY(0.3))), can(0.3, 0.3, 0.3), can(π / 4, π / 4, π / 4),
                  can(0.2, 0.2, 0.0), kron(matrix(H()), matrix(H())), sqrt(matrix(SWAP())),
                  can(1.3, -2.0, 0.9), cis(0.4) * loc() * can(0.1, 0.0, -0.3) * loc())
            k = kak(U)
            @test Matrix(k) ≈ U
            @test inchamber(k)
        end
        # local gates have no entangling content
        k = kak(loc())
        @test (k.a, k.b, k.c) == (0.0, 0.0, 0.0)
        # the angles are local invariants: dressing U changes the factors, not
        # the multiset of |angles| (up to the chamber's symmetries)
        U = randu(4)
        @test sort(abs.([kak(U).a, kak(U).b, kak(U).c])) ≈
              sort(abs.(collect((k -> (k.a, k.b, k.c))(kak(loc() * U * loc())))))
        @test_throws ArgumentError kak(randu(3))
        @test_throws ArgumentError kak(ones(4, 4))
    end

    @testset "circuits use the fewest CNOTs" begin
        cases = [  # (unitary, CNOTs)
            (Matrix{ComplexF64}(I, 4, 4), 0), (loc(), 0),
            (matrix(CNOT()), 1), (matrix(push!(Circuit(2), CNOT(), 2, 1)), 1), (matrix(CZ()), 1),
            (loc() * matrix(CNOT()) * loc(), 1), (matrix(controlled(RY(π))), 1),
            (can(π / 4, 0, 0), 1), (can(-π / 4, 0, 0), 1), (can(0, π / 4, 0), 1), (can(0, 0, -π / 4), 1),
            (matrix(controlled(RY(0.3))), 2), (can(0.3, 0, 0), 2), (can(0.3, 0, 0.2), 2),
            (can(0, 0.3, 0.2), 2), (can(0.3, 0.2, 0), 2), (can(π / 4, π / 4, 0), 2),
            (matrix(SWAP()), 3), (sqrt(matrix(SWAP())), 3), (can(0.3, 0.2, 0.1), 3),
        ]
        for (U, n) in cases
            c = two_qubit(U)
            @test matrix(c) ≈ U
            @test count_cnots(c) == n
            @test all(length(op.qubits) == 1 || op.gate.name === :CNOT for op in c.ops)
        end
        for _ in 1:100
            U = randu(4)
            c = two_qubit(U)
            @test matrix(c) ≈ U && count_cnots(c) == 3
        end
        # on arbitrary wires of a larger register, either order
        U = randu(4)
        c = two_qubit!(Circuit(3), U, [3, 1])
        @test matrix(c) ≈ embed(U, [3, 1], 3)
        @test_throws ArgumentError two_qubit!(Circuit(2), U, [1, 1])
        @test_throws ArgumentError two_qubit!(Circuit(3), U, [1, 2, 3])
    end
end
