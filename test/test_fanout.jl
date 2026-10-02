"""Reference matrix: XOR the bit on `src` wires' parity into each `dst` wire."""
function ref_xor(n, srcs, dsts)
    N = 1 << n
    U = zeros(ComplexF64, N, N)
    bit(x, q) = (x >> (n - q)) & 1
    for x in 0:N-1
        p = reduce(⊻, (bit(x, q) for q in srcs); init=0)
        y = x
        for q in dsts
            y ⊻= p << (n - q)
        end
        U[y+1, x+1] = 1
    end
    U
end

clog2(r) = r <= 1 ? 0 : 64 - leading_zeros(r - 1)

@testset "depth and fan-out" begin
    @testset "depth and layers" begin
        @test depth(Circuit(3)) == 0
        @test isempty(layers(Circuit(3)))

        c = Circuit(4)
        push!(c, H(), 1); push!(c, H(), 2); push!(c, H(), 3)   # one layer
        push!(c, CNOT(), 1, 2); push!(c, CNOT(), 3, 4)           # one layer
        push!(c, CNOT(), 2, 3)                                   # one layer
        push!(c, RZ(0.3), 1)                                     # fits in layer 3
        @test depth(c) == 3
        @test depth(c, :CNOT) == 2
        @test depth(c, :H) == 1
        @test depth(c, :SWAP) == 0
        L = layers(c)
        @test length(L) == depth(c)
        @test sum(length, L) == length(c)
        @test [length(l) for l in L] == [3, 2, 2]
        for l in L                                    # disjoint wires per layer
            @test allunique(reduce(vcat, (op.qubits for op in l)))
        end
        # layers in order reproduce the circuit
        c2 = Circuit(4)
        for l in L, op in l; push!(c2, op.gate, op.qubits...); end
        @test matrix(c2) ≈ matrix(c)

        # a ladder is fully sequential
        c = Circuit(5)
        for i in 1:4; push!(c, CNOT(), i, i + 1); end
        @test depth(c) == 4
    end

    @testset "FANOUT and PARITY gates" begin
        for r in 1:4
            @test matrix(FANOUT(r)) ≈ ref_xor(r + 1, [1], 2:r+1)
            @test matrix(PARITY(r)) ≈ ref_xor(r + 1, 1:r, [r + 1])
            @test is_unitary(matrix(FANOUT(r)))
            # Lemma 5.3 [Moo99, HŠ05]: parity is fan-out conjugated by Hadamards,
            # with the fan-out control becoming the parity target
            c = Circuit(r + 1)
            for q in 1:r+1; push!(c, H(), q); end
            push!(c, FANOUT(r), r + 1, 1:r...)
            for q in 1:r+1; push!(c, H(), q); end
            @test matrix(c) ≈ matrix(PARITY(r))
        end
        @test matrix(FANOUT(1)) ≈ matrix(CNOT())
        @test matrix(PARITY(1)) ≈ matrix(CNOT())
        @test_throws ArgumentError FANOUT(0)
        @test_throws ArgumentError PARITY(0)
        # dense storage caps the gate size; the CNOT styles have no cap
        @test_throws ArgumentError FANOUT(11)
        @test_throws ArgumentError fanout!(Circuit(13), 1, 2:13; style=:gate)
        c = fanout!(Circuit(17), 1, 2:17; style=:tree)
        @test (count_cnots(c), depth(c)) == (31, 9)
    end

    @testset "fan-out styles" begin
        n = 7
        for r in 1:6, style in (:ladder, :tree, :gate)
            ts = collect(2:r+1)
            c = fanout!(Circuit(n), 1, ts; style=style)
            @test matrix(c) ≈ ref_xor(n, [1], ts)     # every input, not just |0⟩ targets
            if style === :ladder
                @test count_cnots(c) == r
                @test depth(c) == r
            elseif style === :tree
                @test count_cnots(c) == 2r - 1
                @test depth(c) == 2clog2(r) + 1
            else
                @test count_cnots(c) == 0
                @test depth(c) == 1
            end
        end
        # arbitrary wire placement, control below its targets
        c = fanout!(Circuit(5), 4, [5, 1, 3]; style=:tree)
        @test matrix(c) ≈ ref_xor(5, [4], [5, 1, 3])
        @test_throws ArgumentError fanout!(Circuit(3), 1, Int[])
        @test_throws ArgumentError fanout!(Circuit(3), 1, [2, 2])
        @test_throws ArgumentError fanout!(Circuit(3), 1, [1, 2])
        @test_throws ArgumentError fanout!(Circuit(3), 1, [2]; style=:wide)
    end

    @testset "parity styles" begin
        n = 7
        for r in 1:6, style in (:ladder, :tree, :gate)
            ss = collect(1:r)
            c = parity!(Circuit(n), ss, n; style=style)
            @test matrix(c) ≈ ref_xor(n, ss, [n])
            if style === :ladder
                @test (count_cnots(c), depth(c)) == (r, r)
            elseif style === :tree
                @test (count_cnots(c), depth(c)) == (2r - 1, 2clog2(r) + 1)
            else
                @test (count_cnots(c), depth(c)) == (0, 1)
            end
        end
        c = parity!(Circuit(5), [5, 2, 4], 1)
        @test matrix(c) ≈ ref_xor(5, [5, 2, 4], [1])
        @test_throws ArgumentError parity!(Circuit(3), [1, 3], 3)
    end

    @testset "phase gadget styles" begin
        θ = 0.71
        for k in 1:6, style in (:ladder, :tree, :gate)
            qs = randperm(6)[1:k]
            c = phase_gadget(θ, qs; n=6, style=style)
            Zs = ["I"^(q - 1) * "Z" * "I"^(6 - q) for q in qs]
            @test matrix(c) ≈ exp(-im * θ / 2 * reduce(*, pauli.(Zs)))
            if k == 1
                @test (count_cnots(c), depth(c)) == (0, 1)
            elseif style === :ladder
                @test (count_cnots(c), depth(c)) == (2(k - 1), 2k - 1)
            elseif style === :tree
                @test (count_cnots(c), depth(c)) == (2(k - 1), 2clog2(k) + 1)
            else
                @test (count_cnots(c), depth(c)) == (0, 3)
            end
        end
        # style reaches through Pauli rotations
        for style in (:ladder, :tree, :gate)
            c = Circuit(5); pauli_rotation!(c, 0.4, "XYZZX"; style=style)
            @test matrix(c) ≈ exp(-im * 0.4 / 2 * pauli("XYZZX"))
        end
        @test_throws ArgumentError phase_gadget(0.1, [1, 2]; style=:wide)
    end

    @testset "drawing" begin
        c = Circuit(4); fanout!(c, 2, [1, 4]; style=:gate); parity!(c, [1, 3], 4; style=:gate)
        s = sprint(draw, c)
        @test occursin("●", s) && occursin("⊕", s)
    end
end
