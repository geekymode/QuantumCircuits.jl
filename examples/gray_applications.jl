# Fifteen applications of Gray coding in quantum circuit synthesis.
#
#   julia --project=. examples/gray_applications.jl
#
# Every number printed here is measured from a circuit that is checked against
# its own definition in the same breath — nothing is quoted from a paper.

using QuantumCircuits
using LinearAlgebra
using Random

Random.seed!(20260818)
randu(n) = (F = qr(randn(ComplexF64, n, n)); Matrix(F.Q) * Diagonal(cis.(2π .* rand(n))))
hdr(i, s) = println("\n", "═"^76, "\n ", i, ". ", s, "\n", "═"^76)
ok(x) = x ? "✓" : "✗ MISMATCH"

# ---------------------------------------------------------------------------
hdr(1, "The code: one bit changes per step, and the walk closes")

println(" i  binary   gray    flip")
for i in 0:7
    println(lpad(i, 2), "  ", bits(i, 3), "    ", bits(gray(i), 3),
            "    ", i == 0 ? "–" : string(gray_flip_position(i, 3)))
end
println("\n adjacent throughout: ", ok(all(gray_adjacent(gray(i), gray(i+1)) for i in 0:6)))
println(" bijection:           ", ok(sort(graycode(3)) == collect(0:7)))
println(" cycle closes:        ", ok(foldl(⊻, (1 << p for p in gray_flip_positions(3))) == 0))
println(" flip position is just trailing_zeros(i): ",
        ok(all(gray_flip_position(i) == trailing_zeros(i) for i in 1:100)))

# ---------------------------------------------------------------------------
hdr(2, "The code IS a circuit: |x⟩ ↦ |gray(x)⟩ in n-1 CNOTs")

for n in 2:6
    E = gray_encoder(n)
    M = matrix(E)
    correct = all(M[gray(x)+1, x+1] ≈ 1 for x in 0:(1 << n)-1)
    inverse = matrix(gray_decoder(n)) * M ≈ I(1 << n)
    println(" n=", n, "  CNOTs=", count_cnots(E), "   |x⟩↦|gray(x)⟩ ", ok(correct),
            "   decoder inverts ", ok(inverse))
end
println("\n gray(x) = x ⊻ (x>>1) is linear over GF(2), so it is pure CNOT —")
println(" n-1 of them, no ancillas, no multi-qubit controls.")
draw(gray_encoder(4))

# ---------------------------------------------------------------------------
hdr(3, "A Gray counter: |gray(i)⟩ ↦ |gray(i+1)⟩")

for n in 2:4
    N = 1 << n
    M = matrix(gray_increment(n))
    stepped = all(M[gray((i+1) % N)+1, gray(i)+1] ≈ 1 for i in 0:N-1)
    onebit = all(hamming(x, argmax(abs.(M[:, x+1])) - 1) == 1 for x in 0:N-1)
    println(" n=", n, "  ", ok(stepped), " counts in Gray order    ",
            ok(onebit), " every step changes exactly one bit")
end
println("\n Decode → increment → encode.  A binary counter can change every bit at")
println(" once (0111→1000); a Gray counter never changes more than one, which is")
println(" why hardware that samples it mid-update still reads a valid code word.")

# ---------------------------------------------------------------------------
hdr(4, "Uniformly controlled rotations: 2^k CNOTs, flat in k")

println(" k | branches | CNOTs | naive (one multi-controlled rotation each)")
for k in 1:6
    α = randn(1 << k)
    c = multiplexed_rz(α, 1:k, k+1)
    println(lpad(k, 2), " | ", lpad(1 << k, 8), " | ", lpad(count_cnots(c), 5),
            " | ", lpad(2 * k * (1 << k), 10), "   ", ok(count_cnots(c) == 1 << k))
end
draw(multiplexed_rz([0.1, 0.5, -1.2, 2.0], [1, 2], 3))
println(" CNOT controls read q2,q1,q2,q1 — the Gray flip sequence for k=2.")

# ---------------------------------------------------------------------------
hdr(5, "Uniformly controlled *arbitrary* one-qubit gates")

for k in 1:4
    Us = [randu(2) for _ in 1:(1 << k)]
    c = multiplexed_1q(Us, 1:k, k+1)
    ref = cat(Us...; dims=(1, 2))
    println(" k=", k, "  CNOTs=", lpad(count_cnots(c), 4), "  ", ok(matrix(c) ≈ ref),
            "   (three multiplexors RZ·RY·RZ + a diagonal on the controls)")
end

# ---------------------------------------------------------------------------
hdr(6, "Diagonal unitaries: 2^n - 2 CNOTs")

for n in 2:7
    φ = randn(1 << n)
    c = diagonal(φ)
    println(" n=", n, "  CNOTs=", lpad(count_cnots(c), 4), "  ",
            ok(count_cnots(c) == (1 << n) - 2 && matrix(c) ≈ Diagonal(cis.(φ))))
end

# ---------------------------------------------------------------------------
hdr(7, "State preparation: 2^(n+1) - 4 CNOTs")

for n in 2:7
    a = randn(ComplexF64, 1 << n)
    c = prepare_state(a)
    println(" n=", n, "  CNOTs=", lpad(count_cnots(c), 4), "  ",
            ok(statevector(c) ≈ a ./ norm(a)), "   exact, phases included")
end
draw(prepare_state([1, 2, 3im, 4]))

# ---------------------------------------------------------------------------
hdr(8, "Two-level unitaries: routing |a⟩ and |b⟩ to be Gray-adjacent")

V = randu(2)
for (a, b) in ((0, 7), (1, 6), (2, 3), (0, 4))
    c = Circuit(3)
    two_level!(c, TwoLevel(a, b, V), 3)
    println(" |", bits(a, 3), "⟩↔|", bits(b, 3), "⟩  hamming=", hamming(a, b),
            "  routing CNOTs=", count_cnots(c), "   ",
            ok(matrix(c) ≈ matrix(TwoLevel(a, b, V), 8)))
end
println("\n hamming(a,b) - 1 CNOTs each side bring the pair one bit apart;")
println(" then a single multi-controlled V does the work.")
c = Circuit(3); two_level!(c, TwoLevel(1, 6, V), 3)
draw(c)

# ---------------------------------------------------------------------------
hdr(9, "Arbitrary unitary synthesis: Givens factors + Gray routing")

for n in 1:3
    U = randu(1 << n)
    fs = two_level_decompose(U)
    hi = synthesize_unitary(U)
    lo = synthesize_unitary(U; lower=true)
    println(" n=", n, "  two-level factors=", lpad(length(fs), 3),
            "   high-level gates=", lpad(length(hi), 4),
            "   fully lowered: ", lpad(length(lo), 4), " gates / ",
            lpad(count_cnots(lo), 4), " CNOTs   ",
            ok(matrix(hi) ≈ U && matrix(lo) ≈ U))
end
println("\n every gate elementary after lowering: ",
        ok(all(length(op.qubits) <= 2 for op in synthesize_unitary(randu(8); lower=true).ops)))

# ---------------------------------------------------------------------------
hdr(10, "Multi-controlled C^n(U), ancilla-free, from a parity walk")

println(" AND(x) = 2^-(n-1) Σ_{S≠∅} (-1)^{|S|+1} parity_S(x), and the 2^n-1")
println(" parities are visited in Gray order — one CNOT per step.\n")
println(" controls | controlled-V gates | CNOTs | correct")
for k in 2:5
    local U = randu(2)
    local c = multicontrolled(U, 1:k, k+1)
    println(lpad(k, 9), " | ", lpad(length(c) - count_cnots(c), 18), " | ",
            lpad(count_cnots(c), 5), " | ",
            ok(matrix(c) ≈ matrix(controlled(Gate(:U, U), k))))
end
println("\n Toffoli recovered from U=X: ",
        ok(matrix(multicontrolled(matrix(X()), [1, 2], 3)) ≈ matrix(controlled(X(), 2))))
draw(multicontrolled(matrix(X()), [1, 2], 3))

# ---------------------------------------------------------------------------
hdr(11, "Phase polynomials: parity networks vs one gadget per term")

println(" n | terms | gadget each | parity network | optimal 2^n-2 | correct")
for n in 2:7
    φ = randn(1 << n)
    pp = phase_polynomial(φ)
    a = synthesize(pp; order=:gadgets)
    b = synthesize(pp; order=:gray)
    D = Diagonal(cis.(φ))
    println(lpad(n, 2), " | ", lpad(nterms(pp), 5), " | ", lpad(count_cnots(a), 11),
            " | ", lpad(count_cnots(b), 14), " | ", lpad((1 << n) - 2, 13),
            " | ", ok(matrix(a) ≈ D && matrix(b) ≈ D))
end

# ---------------------------------------------------------------------------
hdr(12, "Sparse phase polynomials")

n = 6
for (name, masks) in ("5 terms, weight 4" => [0b001111, 0b011110, 0b111100, 0b111001, 0b100111],
                      "8 scattered terms" => [0b000011, 0b010101, 0b111000, 0b001110,
                                              0b100001, 0b011011, 0b110110, 0b101101],
                      "nearest-neighbour chain" => [0b000011, 0b000110, 0b001100, 0b011000, 0b110000])
    co = zeros(1 << n)
    for S in masks; co[S+1] = randn(); end
    pp = PhasePolynomial(n, co)
    D = Diagonal(cis.(phases(pp)))
    a = synthesize(pp; order=:gadgets)
    b = synthesize(pp; order=:gray)
    println(" ", rpad(name, 24), " gadgets=", lpad(count_cnots(a), 3),
            "  parity network=", lpad(count_cnots(b), 3),
            "   ", ok(matrix(a) ≈ D && matrix(b) ≈ D))
end
println("\n The chain ties: weight-2 terms already cost two CNOTs apiece, so there")
println(" is nothing for the network to save. Gray coding is not magic.")

# ---------------------------------------------------------------------------
hdr(13, "Walsh-series truncation: trading CNOTs for accuracy")

n = 6
N = 1 << n
φ = [2.0 * sin(2π * x / N) + 0.4 * (x / N)^2 for x in 0:N-1]
pp = phase_polynomial(φ)
println(" a smooth phase f(x) = 2 sin(2πx/64) + 0.4 (x/64)²\n")
println(" terms kept | CNOTs | max phase error | fidelity vs exact")
for k in (2, 4, 8, 16, 32, N)
    local tp = truncate_terms(pp, k)
    local c = synthesize(tp)
    println(lpad(nterms(tp), 11), " | ", lpad(count_cnots(c), 5), " | ",
            lpad(string(round(maximum(abs, phases(tp) .- φ), sigdigits=3)), 15), " | ",
            round(gate_fidelity(matrix(c), Diagonal(cis.(φ))), digits=6))
end

# ---------------------------------------------------------------------------
hdr(14, "LCU / QROM: sweeping an address register in Gray order")

println(" A branch fires on |1⟩, so address j needs an X on every control where")
println(" j has a zero. Consecutive Gray addresses differ in one bit ⟹ one X.\n")
println(" address bits | Gray X gates | natural X gates | saving | correct")
for k in 2:7
    Us = [randu(2) for _ in 1:(1 << k)]
    g = select(Us, 1:k, [k+1]; order=:gray)
    nat = select(Us, 1:k, [k+1]; order=:natural)
    ref = cat(Us...; dims=(1, 2))
    println(lpad(k, 13), " | ", lpad(count_gates(g, :X), 12), " | ",
            lpad(count_gates(nat, :X), 15), " | ",
            lpad(string(round(Int, 100 * (1 - count_gates(g, :X) / count_gates(nat, :X)))) * "%", 6),
            " | ", ok(matrix(g) ≈ ref && matrix(nat) ≈ ref))
end

# ---------------------------------------------------------------------------
hdr(15, "Hamiltonian simulation: pooling the commuting diagonal terms")

function tfim(n; J=1.0, h=0.6)                     # transverse-field Ising
    t = Pair{String,Float64}[]
    for i in 1:n-1
        s = collect("I"^n); s[i] = 'Z'; s[i+1] = 'Z'; push!(t, String(s) => J)
    end
    for i in 1:n
        s = collect("I"^n); s[i] = 'X'; push!(t, String(s) => h)
    end
    t
end

function qaoa_cost(n)                              # a dense diagonal cost operator
    t = Pair{String,Float64}[]
    for S in 1:(1 << n)-1
        s = collect("I"^n)
        for b in 0:n-1
            ((S >> b) & 1) == 1 && (s[n-b] = 'Z')
        end
        push!(t, String(s) => randn())
    end
    t
end

println(" transverse-field Ising (low-weight Z couplings — nothing to pool):")
for n in 4:7
    local t = tfim(n)
    local a = Circuit(n); trotter_step!(a, t, 0.05; parity_network=false)
    local b = Circuit(n); trotter_step!(b, t, 0.05; parity_network=true)
    println("   n=", n, "  ", lpad(length(t), 3), " terms:  gadgets=", lpad(count_cnots(a), 4),
            "   pooled=", lpad(count_cnots(b), 4))
end
println("\n dense diagonal cost operator (QAOA phase separator — everything pools):")
for n in 3:6
    local t = qaoa_cost(n)
    local a = Circuit(n); trotter_step!(a, t, 0.05; parity_network=false)
    local b = Circuit(n); trotter_step!(b, t, 0.05; parity_network=true)
    println("   n=", n, "  ", lpad(length(t), 3), " terms:  gadgets=", lpad(count_cnots(a), 4),
            "   pooled=", lpad(count_cnots(b), 4), "   same operator ",
            ok(matrix(a) ≈ matrix(b)), " (they all commute)")
end

Hm = randn(ComplexF64, 8, 8); Hm = Hm + Hm'
terms = pauli_decompose(Hm)
c = Circuit(3); trotter_step!(c, terms, 0.02)
println("\n a dense random 3-qubit Hamiltonian: ", length(terms), " Pauli terms, ",
        count_cnots(c), " CNOTs")
println(" fidelity vs exp(-iH·dt): ",
        round(gate_fidelity(matrix(c), exp(-im * 0.02 * Hm)), digits=6))

println("\n", "═"^76)
println(" Illustrations of all of the above:  using CairoMakie; circuitfigure(c)")
println("═"^76)
