@testset "stabilizer simulation" begin
    samestate(a, b) = abs(abs(dot(a, b)) - 1) < 1e-9
    rng = Random.Xoshiro(2026)
    oneq = [X(), Y(), Z(), H(), S(), Sdg(), Id(), RZ(π / 2), RZ(π), RX(π / 2), RY(π / 2),
            RY(-π / 2), RX(3π / 2), PHASE(π / 2), PHASE(-π / 2)]
    twoq = [CNOT(), CZ(), SWAP(), controlled(X()), controlled(Y()), controlled(Z())]
    function randclifford(n, m)
        c = Circuit(n)
        for _ in 1:m
            if n > 1 && rand(rng) < 0.4
                a, b = randperm(rng, n)[1:2]
                push!(c, rand(rng, twoq), a, b)
            else
                push!(c, rand(rng, oneq), rand(rng, 1:n))
            end
        end
        c
    end

    @testset "agrees with the statevector simulator" begin
        for _ in 1:150
            n = rand(rng, 1:5)
            c = randclifford(n, rand(rng, 5:40))
            t = apply!(Tableau(n), c)
            ψ = statevector(c)
            @test samestate(statevector(t), ψ)
            for _ in 1:4
                P = PauliOp(join(rand(rng, ['I', 'X', 'Y', 'Z']) for _ in 1:n))
                P = rand(rng, Bool) ? -P : P
                @test expectation(t, P) == round(Int, real(dot(ψ, matrix(P) * ψ)))
            end
            # measuring collapses onto the projection for the reported outcome
            P = PauliOp(join(rand(rng, ['I', 'X', 'Y', 'Z']) for _ in 1:n))
            out = measure!(t, P; rng=rng)
            φ = (I + (out ? -1 : 1) * matrix(P)) / 2 * ψ
            @test norm(φ) > 1e-9 && samestate(statevector(t), φ / norm(φ))
        end
    end

    @testset "measurement" begin
        t = Tableau(2)
        @test measure!(t, 1) == false && measure!(t, PauliOp("ZZ")) == false   # deterministic
        apply!(t, X(), 2)
        @test measure!(t, 2) == true && expectation(t, PauliOp("ZZ")) == -1
        # a Bell pair: outcomes random but perfectly correlated
        for _ in 1:20
            t = Tableau(2); apply!(t, H(), 1); apply!(t, CNOT(), 1, 2)
            @test expectation(t, PauliOp("ZI")) == 0
            @test measure!(t, 1; rng=rng) == measure!(t, 2; rng=rng)
        end
        ones_ = count(_ -> measure!(Tableau(1), PauliOp("X"); rng=rng), 1:2000)
        @test 850 < ones_ < 1150                                           # fair coin
        # stabilizers after H⊗I, CNOT are XX and ZZ
        t = Tableau(2); apply!(t, H(), 1); apply!(t, CNOT(), 1, 2)
        @test Set(string.(stabilizers(t))) == Set(["+XX", "+ZZ"])
        # large: a 1000-qubit GHZ state, every qubit measures the same
        t = Tableau(1000); apply!(t, H(), 1)
        for q in 2:1000; apply!(t, CNOT(), q - 1, q); end
        outs = [measure!(t, q; rng=rng) for q in 1:1000]
        @test all(==(outs[1]), outs)
    end

    @testset "rejects non-Clifford input" begin
        @test_throws ArgumentError apply!(Tableau(1), T(), 1)
        @test_throws ArgumentError apply!(Tableau(1), RZ(0.3), 1)
        @test_throws ArgumentError apply!(Tableau(2), Circuit(3))
        @test_throws ArgumentError measure!(Tableau(1), PauliOp("iZ"))
        @test_throws ArgumentError measure!(Tableau(1), 2)
        @test_throws ArgumentError Tableau(0)
    end

    @testset "error correction at scale" begin
        # |0̄⟩ by measuring stabilizers, for CSS and non-CSS codes alike
        for code in (five_qubit_code(), steane_code(), shor_code(), quantum_reed_muller(4),
                     rotated_surface_code(5))
            t = prepare_logical_zero(code; rng=rng)
            @test all(expectation(t, s) == 1 for s in stabilizers(code))
            @test all(expectation(t, z) == 1 for z in logical_operators(code)[2])
        end
        # it is the encoder's state
        code = steane_code()
        c, inp = encoding_circuit(code)
        @test samestate(statevector(prepare_logical_zero(code; rng=rng)), statevector(c))
        # measured syndromes equal the computed ones, up to 97 qubits
        for d in (3, 5, 7)
            code = rotated_surface_code(d)
            n = length(code)
            for _ in 1:5
                E = PauliOp(rand(rng, n) .< 0.1, rand(rng, n) .< 0.1)
                @test sample_syndrome(code, E; rng=rng) == syndrome(code, E)
            end
        end
        # the five-qubit code too — no encoder needed
        code = five_qubit_code()
        for q in 1:5, P in ("X", "Y", "Z")
            E = PauliOp(join(j == q ? P : "I" for j in 1:5))
            @test sample_syndrome(code, E; rng=rng) == syndrome(code, E)
        end
    end

    @testset "logical error rates" begin
        # the decoder's guarantee, exactly: every error of weight ≤ 2 on d = 5
        code = rotated_surface_code(5)
        table = lookup_decoder(code; maxweight=2)
        L = vcat(logical_operators(code)...)
        for t in 1:2, pos in QuantumCircuits._combinations(25, t),
            kinds in Iterators.product(ntuple(_ -> 1:3, t)...)
            x = falses(25); z = falses(25)
            for (q, k) in zip(pos, kinds); x[q] = k != 3; z[q] = k != 1; end
            R = table[syndrome(code, PauliOp(x, z))] * PauliOp(x, z)
            @test all(l -> commutes(R, l), L)
        end
        @test logical_error_rate(code, 0.0; shots=100, decoder=table)[1] == 0.0
        # Monte Carlo, with wide margins: d = 5 beats d = 3 at p = 0.005
        r3 = logical_error_rate(rotated_surface_code(3), 0.005; shots=40_000, rng=Random.Xoshiro(1))[1]
        r5 = logical_error_rate(code, 0.005; shots=40_000, decoder=table, rng=Random.Xoshiro(1))[1]
        @test r5 < r3
        @test_throws ArgumentError logical_error_rate(code, 1.5; decoder=table)
    end
end
