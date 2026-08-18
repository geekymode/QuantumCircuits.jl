@testset "phase polynomials" begin
    @testset "phase gadgets" begin
        for (qs, n) in (([1], 1), ([1, 2], 2), ([2, 1], 2), ([1, 2, 3], 3), ([1, 3], 3), ([3, 1, 2], 3))
            θ = 0.83
            c = phase_gadget(θ, qs; n=n)
            Zs = ["I"^(q - 1) * "Z" * "I"^(n - q) for q in qs]
            P = reduce(*, pauli.(Zs))
            @test matrix(c) ≈ exp(-im * θ / 2 * P)
            @test count_cnots(c) == 2 * (length(qs) - 1)
        end
        @test_throws ArgumentError phase_gadget(0.1, Int[]; n=2)
        @test_throws ArgumentError phase_gadget(0.1, [1, 1]; n=2)
    end

    @testset "Pauli rotations" begin
        for s in ("I", "X", "Y", "Z", "XX", "XY", "YZ", "ZI", "IZ", "XYZ", "IXI", "ZZZ", "YYYY")
            θ = 0.61
            c = Circuit(length(s)); pauli_rotation!(c, θ, s)
            @test matrix(c) ≈ exp(-im * θ / 2 * pauli(s))
        end
        c = Circuit(3); pauli_rotation!(c, 0.4, "XZ"; qubits=[3, 1])
        @test matrix(c) ≈ exp(-im * 0.4 / 2 * pauli("ZIX"))
        @test_throws ArgumentError pauli_rotation!(Circuit(2), 0.1, "XYZ")
    end

    @testset "Trotter steps" begin
        Hm = randn(ComplexF64, 4, 4); Hm = Hm + Hm'
        terms = pauli_decompose(Hm)
        dt = 0.03
        c = Circuit(2); trotter_step!(c, terms, dt; parity_network=false)
        # circuits apply terms left to right, matrices multiply right to left
        exact_product = reduce(*, (exp(-im * dt * real(cf) * pauli(s)) for (s, cf) in reverse(terms)))
        @test matrix(c) ≈ exact_product                       # exact, term by term

        # pooling the commuting diagonal terms is a different ordering, but
        # still a first-order step for the same Hamiltonian
        cp = Circuit(2); trotter_step!(cp, terms, dt; parity_network=true)
        @test gate_fidelity(matrix(cp), exp(-im * dt * Hm)) > 1 - 1e-3
        @test count_cnots(cp) <= count_cnots(c)
        @test gate_fidelity(matrix(c), exp(-im * dt * Hm)) > 1 - 1e-3   # ≈ exp(-iHt)
        @test_throws ArgumentError trotter_step!(Circuit(2), ["XY" => 1.0 + 1.0im], dt)
        @test_throws ArgumentError trotter_step!(Circuit(3), ["XY" => 1.0], dt; qubits=[1, 2, 3])
        # a short string defaults to the leading wires, which is fine
        @test trotter_step!(Circuit(3), ["XY" => 1.0], dt) isa Circuit
    end

    @testset "polynomial round trip" begin
        for n in 1:6
            φ = randn(1 << n)
            pp = phase_polynomial(φ)
            @test pp.n == n
            @test phases(pp) ≈ φ
            @test nterms(pp) == (1 << n) - 1                  # dense, generically
        end
        pp = phase_polynomial([0.0, 0.0, 0.0, 0.0])
        @test isempty(support(pp))
        @test_throws ArgumentError phase_polynomial(randn(3))
    end

    @testset "synthesis" begin
        for n in 1:5
            φ = randn(1 << n)
            pp = phase_polynomial(φ)
            D = Diagonal(cis.(φ))
            cg = synthesize(pp; order=:gray)
            cn = synthesize(pp; order=:gadgets)
            @test matrix(cg) ≈ D
            @test matrix(cn) ≈ D
            # the parity network matches the optimal recursive construction
            @test count_cnots(cg) == (1 << n) - 2
            @test count_cnots(cg) <= count_cnots(cn)
        end
        # sparse support
        n = 5
        co = zeros(1 << n)
        for S in (0b00011, 0b01110, 0b11100, 0b10001); co[S+1] = randn(); end
        pp = PhasePolynomial(n, co)
        for ord in (:gray, :gadgets)
            @test matrix(synthesize(pp; order=ord)) ≈ Diagonal(cis.(phases(pp)))
        end
        @test count_cnots(synthesize(pp; order=:gray)) <= count_cnots(synthesize(pp; order=:gadgets))
        @test_throws ArgumentError synthesize(phase_polynomial(randn(4)); order=:nope)
    end

    @testset "CNOT cancellation preserves the unitary" begin
        c = Circuit(3)
        push!(c, CNOT(), 1, 2); push!(c, H(), 3); push!(c, CNOT(), 1, 2)  # cancels around H(3)
        before = matrix(c)
        cancel_adjacent_cnots!(c)
        @test count_cnots(c) == 0
        @test matrix(c) ≈ before

        c = Circuit(2)
        push!(c, CNOT(), 1, 2); push!(c, H(), 2); push!(c, CNOT(), 1, 2)  # must NOT cancel
        before = matrix(c)
        cancel_adjacent_cnots!(c)
        @test count_cnots(c) == 2
        @test matrix(c) ≈ before

        for _ in 1:20                                          # random circuits
            c = Circuit(3)
            for _ in 1:12
                r = rand(1:3)
                if r == 1
                    a, b = randperm(3)[1:2]; push!(c, CNOT(), a, b)
                elseif r == 2
                    push!(c, RZ(randn()), rand(1:3))
                else
                    push!(c, H(), rand(1:3))
                end
            end
            before = matrix(c)
            n0 = count_cnots(c)
            cancel_adjacent_cnots!(c)
            @test matrix(c) ≈ before
            @test count_cnots(c) <= n0
        end
    end
end
