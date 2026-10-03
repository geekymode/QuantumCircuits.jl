"""Basis state with the given wire => bit assignments, other wires 0."""
qec_basis(n, bits) = (v = zeros(ComplexF64, 1 << n); v[sum((b << (n - w) for (w, b) in bits); init=0)+1] = 1; v)

function qec_encoded(code, xs)
    c, inputs = encoding_circuit(code)
    statevector(c, qec_basis(length(code), zip(inputs, xs)))
end

@testset "error-correction circuits" begin
    @testset "encoders" begin
        for code in (steane_code(), shor_code(), rotated_surface_code(3), quantum_reed_muller(4),
                     repetition_code(4), css_code([1 1 1 1 0 0; 0 0 1 1 1 1], [1 1 1 1 0 0; 0 0 1 1 1 1]))
            k = dimension(code)
            c, inputs = encoding_circuit(code)
            @test all(op.gate.name in (:H, :CNOT, :SWAP) for op in c.ops)   # Clifford: H + CNOT
            X̄, Z̄ = logical_operators(code)
            for xs in Iterators.product(ntuple(_ -> 0:1, k)...)
                ψ = qec_encoded(code, collect(xs))
                @test all(s * ψ ≈ ψ for s in stabilizers(code))            # in the code space
                @test all(Z̄[j] * ψ ≈ (-1)^xs[j] * ψ for j in 1:k)          # reads back x
                for j in 1:k                                               # X̄ⱼ flips bit j
                    ys = collect(xs); ys[j] ⊻= 1
                    @test X̄[j] * ψ ≈ qec_encoded(code, ys)
                end
            end
        end
        # superpositions encode linearly
        code = steane_code()
        c, inputs = encoding_circuit(code)
        α, β = 0.6, 0.8im
        ψin = α * qec_basis(7, [inputs[1] => 0]) + β * qec_basis(7, [inputs[1] => 1])
        @test statevector(c, ψin) ≈ α * qec_encoded(code, [0]) + β * qec_encoded(code, [1])
        @test_throws ArgumentError encoding_circuit(five_qubit_code())
    end

    @testset "syndrome extraction" begin
        for code in (steane_code(), five_qubit_code(), rotated_surface_code(3))
            n, r = length(code), length(stabilizers(code))
            sc = syndrome_circuit(code)
            @test nqubits(sc) == n + r && ancillas(sc) == collect(n+1:n+r)
            # a generic encoded state: project a random state onto the code space
            ψ = normalize(randn(ComplexF64, 1 << n))
            for s in stabilizers(code); ψ = (ψ + s * ψ) / 2; end
            normalize!(ψ)
            for q in 1:n, P in ("X", "Y", "Z")
                E = PauliOp(join(j == q ? P : "I" for j in 1:n))
                s = syndrome(code, E)
                anc = sum((Int(s[i]) << (r - i) for i in 1:r); init=0)
                # data untouched, ancillas hold exactly the syndrome
                e_anc = zeros(ComplexF64, 1 << r); e_anc[anc+1] = 1
                @test statevector(sc, E * ψ) ≈ kron(E * ψ, e_anc)
            end
            # it is meant to write the ancillas (checked where the isometry is small)
            n <= 7 && @test !is_clean(sc)
            # no error: ancillas stay |0⟩, the code state is untouched
            @test statevector(sc, ψ) ≈ kron(ψ, qec_basis(r, Pair{Int,Int}[]))
        end
    end

    @testset "transversal T on the [[15,1,3]] quantum Reed–Muller code" begin
        Tn(ψ, n) = [cis(π / 4 * count_ones(b)) * ψ[b+1] for b in 0:(1 << n)-1]
        code = quantum_reed_muller(4)
        z, o = qec_encoded(code, [0]), qec_encoded(code, [1])
        # codeword weights: 0 mod 8 in |0̄⟩, 7 mod 8 in |1̄⟩ …
        @test all(count_ones(b) % 8 == 0 for b in 0:(1 << 15)-1 if abs(z[b+1]) > 1e-9)
        @test all(count_ones(b) % 8 == 7 for b in 0:(1 << 15)-1 if abs(o[b+1]) > 1e-9)
        # … so T on every qubit is exactly a logical T†
        @test Tn(z, 15) ≈ z
        @test Tn(o, 15) ≈ cis(-π / 4) * o
        # on the Steane code it is not a logical gate: it leaves the code space
        st = steane_code()
        z7 = qec_encoded(st, [0])
        T7 = Tn(z7, 7)
        @test norm(T7 - dot(z7, T7) * z7) > 0.1
    end
end
