# Error correction

```@meta
CurrentModule = QuantumCircuits
```

A quantum error-correcting code spreads `k` logical qubits over `n` physical
ones so that errors on a few of them can be found and undone without ever
learning the encoded state. The codes here are *stabilizer codes*, and the
most useful of them are built from pairs of classical codes — so this page
starts classical, with the family that ties the whole package together:
**Reed–Muller codes**.

```@setup qec
using QuantumCircuits, LinearAlgebra, Random
Random.seed!(20261002)
```

## Classical codes over GF(2)

A binary linear code is a subspace of GF(2)ⁿ: codewords are the row span of
a generator matrix `G`, equivalently the null space of a parity-check matrix
`H`. Everything is row reduction mod 2 ([`gf2_rref`](@ref),
[`gf2_nullspace`](@ref)).

The `[7, 4, 3]` Hamming code has a parity-check matrix whose column `j` is `j`
in binary — so the syndrome of a single bit flip *is* the flipped position:

```@example qec
h = hamming_code(3)
w = encode(h, [1, 0, 1, 1])
y = copy(w); y[5] = !y[5]                     # flip bit 5
println("syndrome ", Int.(syndrome(h, y)), "  →  decoded correctly: ", syndrome_decode(h, y) == w)
```

## Reed–Muller codes

`RM(r, m)` is the set of evaluation tables, over all `2ᵐ` points of GF(2)ᵐ,
of Boolean polynomials of degree at most `r`. These are the same objects the
package already uses as **phase polynomials**: a diagonal unitary's parity
expansion is a Reed–Muller-style polynomial, and minimising the `T` count of
a `{CNOT, T}` circuit is a Reed–Muller decoding problem.

```@example qec
for (r, m) in ((1, 3), (1, 4), (2, 4), (1, 5), (2, 6), (3, 5))
    c = reed_muller(r, m)
    println("RM($r,$m) = [", length(c), ", ", dimension(c), ", ", minimum_distance(c), "]")
end
```

The distances come from an exhaustive walk over all `2ᵏ` codewords — `2²⁶`
of them for `RM(3, 5)` — taken **in Gray-code order**, so that each step is a
single XOR with one generator row.

Reed–Muller codes are *locally decodable*: a degree-`r` polynomial sums to
zero over any `(r+1)`-dimensional affine subspace, so one coordinate can be
recovered from `2^(r+1) - 1` others, with a majority vote over random
subspaces absorbing a few errors ([`rm_local_decode`](@ref)):

```@example qec
c = reed_muller(2, 6)
w = encode(c, rand(Bool, dimension(c)))
y = copy(w); for i in randperm(64)[1:2]; y[i] = !y[i]; end     # two errors
count(rm_local_decode(y, 2, 6, x; trials = 31) == w[x+1] for x in 0:63)
```

This is the classical shadow of the [three-query identity](@ref "Trading
depth for width"): over the reals, a degree-2 form is pinned by **three**
collinear points; over GF(2) the same job takes `2³ - 1 = 7`.

## Stabilizer codes

A Pauli operator ([`PauliOp`](@ref)) is a pair of bit vectors and a phase;
two Paulis commute exactly when their symplectic product is even. A
stabilizer code is the joint `+1` eigenspace of `r` commuting Paulis on `n`
qubits — `k = n - r` logical qubits — and a *CSS code* ([`css_code`](@ref))
takes its `X` checks from one classical code and its `Z` checks from another.

```@example qec
for code in (five_qubit_code(), steane_code(), shor_code(), quantum_reed_muller(4),
             rotated_surface_code(3), rotated_surface_code(5))
    println(rpad(string(code), 36), "d = ", code_distance(code), is_css(code) ? "   CSS" : "")
end
```

Logical operators are computed when not given — by symplectic Gram–Schmidt
on the normalizer, or for CSS codes from the two classical codes so that `X̄`
is pure `X` and `Z̄` pure `Z`:

```@example qec
X̄, Z̄ = logical_operators(rotated_surface_code(3))
string(X̄[1]), string(Z̄[1])
```

## Encoding and syndrome circuits

A CSS encoder is a network of Hadamards and CNOTs ([`encoding_circuit`](@ref))
— the same kind of object as the Gray-code encoder:

```@example qec
c, inputs = encoding_circuit(steane_code())
println("logical input on wire ", inputs[1], "; ", count_cnots(c), " CNOTs, depth ", depth(c))
draw(c)
```

Syndrome extraction ([`syndrome_circuit`](@ref)) gives every stabilizer its
own ancilla — the ancilla-aware circuits of the depth page. On an encoded
state hit by a Pauli error, the ancillas end in the syndrome and the data is
left alone; circuits stay unitary, so this is checked exactly rather than by
sampling measurements:

```@example qec
code = steane_code()
sc = syndrome_circuit(code)
enc, inp = encoding_circuit(code)
ψ = statevector(enc, (v = zeros(ComplexF64, 128); v[1] = 0.6; v[1 + (1 << (7 - inp[1]))] = 0.8im; v))
E = PauliOp("IIIIYII")                                  # a Y error on qubit 5
s = syndrome(code, E)
expected = kron(E * ψ, (a = zeros(ComplexF64, 64); a[1 + sum(Int(s[i]) << (6 - i) for i in 1:6)] = 1; a))
println("syndrome ", Int.(s), ";  circuit output exact: ", statevector(sc, E * ψ) ≈ expected)
```

A lookup table ([`lookup_decoder`](@ref)) maps each syndrome to the
lowest-weight error that explains it; applying that correction leaves at most
a stabilizer, which does nothing to the encoded state:

```@example qec
code = rotated_surface_code(3)
table = lookup_decoder(code)
fixed = all(q -> all(P -> begin
            E = PauliOp(join(j == q ? P : "I" for j in 1:9))
            C = table[syndrome(code, E)]
            !any(syndrome(code, C * E))
        end, ("X", "Y", "Z")), 1:9)
println(length(table), " syndromes; every single-qubit error corrected: ", fixed)
```

## Transversal T on the 15-qubit Reed–Muller code

Applying `T` to every qubit of a code is fault tolerant by construction — an
error on one qubit stays on one qubit — but it is rarely a logical gate. On
the `[[15, 1, 3]]` quantum Reed–Muller code it is, and the reason is a fact
about Reed–Muller weights: every computational-basis term of `|0̄⟩` has
Hamming weight `≡ 0 (mod 8)`, and of `|1̄⟩` weight `≡ 7`. `T⊗15` multiplies a
term of weight `w` by `e^{iπw/4}`:

```@example qec
code = quantum_reed_muller(4)
enc, inp = encoding_circuit(code)
basis1(x) = (v = zeros(ComplexF64, 1 << 15); v[1 + (x << (15 - inp[1]))] = 1; v)
z, o = statevector(enc, basis1(0)), statevector(enc, basis1(1))
weights(ψ) = sort(unique(count_ones(b) % 8 for b in 0:(1 << 15)-1 if abs(ψ[b+1]) > 1e-9))
Tn(ψ) = [cis(π / 4 * count_ones(b)) * ψ[b+1] for b in 0:(1 << 15)-1]
println("weights mod 8:  |0̄⟩ ", weights(z), "   |1̄⟩ ", weights(o))
println("T⊗15 |0̄⟩ = |0̄⟩: ", Tn(z) ≈ z, ";   T⊗15 |1̄⟩ = e^{-iπ/4} |1̄⟩: ", Tn(o) ≈ cis(-π / 4) * o)
```

So `T⊗15` is exactly a logical `T†` — the property that makes this code the
workhorse of magic-state distillation. The Steane code, its `m = 3` sibling,
has weights `{0, 4}` and `{3, 7}` and no such luck.

## Simulating at scale: the stabilizer tableau

Every circuit on this page is Clifford, and Clifford circuits on stabilizer
states can be simulated in polynomial time (Gottesman–Knill). A
[`Tableau`](@ref) stores `n` stabilizers and `n` destabilizers as packed bits —
`O(n²)` memory instead of `2ⁿ` amplitudes — and brings **measurement**, which
the unitary `Circuit` model leaves out. [`apply!`](@ref) runs any Clifford
circuit on it:

```@example qec
t = Tableau(1000)
apply!(t, H(), 1)
for q in 2:1000; apply!(t, CNOT(), q - 1, q); end          # a 1000-qubit GHZ state
outs = [measure!(t, q) for q in 1:1000]
println("first outcome ", outs[1], "; all 1000 agree: ", all(==(outs[1]), outs))
```

With measurement, a logical state can be prepared the way hardware does it —
measure every stabilizer, fix the `-1` outcomes — which works for any
stabilizer code, including the five-qubit code that has no encoder here
([`prepare_logical_zero`](@ref)). And syndromes can be *measured*: a
distance-7 surface code with its syndrome ancillas is 97 qubits, far past the
statevector, and [`sample_syndrome`](@ref) runs the real circuit:

```@example qec
code = rotated_surface_code(7)
E = PauliOp("I"^10 * "X" * "I"^20 * "Z" * "I"^17)
println(length(code) + length(stabilizers(code)), " qubits; measured syndrome = computed: ",
        sample_syndrome(code, E) == syndrome(code, E))
```

## Decoding by matching

A lookup table only knows errors up to weight `⌊(d-1)/2⌋`; anything heavier
is a guaranteed failure, so it wastes most of what a large code offers. The
surface code has more structure to use. Every qubit lies in at most two checks
of each type, so an `X` error is a set of *edges* in a graph whose nodes are
the `Z` checks, and its syndrome is the set of nodes where error chains end.
(A qubit in only one check gets an edge to a virtual boundary node.)

The most likely error under independent noise is the shortest set of chains
explaining the syndrome: pair up the flipped checks, or send them to the
boundary, at minimum total length. [`matching_decoder`](@ref) precomputes all
shortest paths, and [`decode`](@ref) finds that matching **exactly**, by
dynamic programming over subsets of the flipped checks — `O(2ᵐ m)` for `m`
defects, which covers the codes and noise rates here (past `maxdefects` it
falls back to greedy pairing). `Z` errors decode the same way on the `X`
checks. Whatever the error, the correction it returns explains the syndrome
exactly:

```@example qec
code = rotated_surface_code(7)
dec = matching_decoder(code)
E = PauliOp(rand(49) .< 0.1, rand(49) .< 0.1)                # a heavy random error
println("error weight ", weight(E), ";  correction explains the syndrome: ",
        syndrome(code, decode(dec, syndrome(code, E))) == syndrome(code, E))
```

## Logical error rates

[`logical_error_rate`](@ref) estimates how often decoding fails when every
qubit independently suffers `X`, `Y` or `Z` with probability `p/3` each
(code-capacity noise). A code is *working* when making it bigger makes it
better. With the lookup table that only holds below about 1%; with matching
it holds up to about 15%, matching the known threshold for decoding `X` and
`Z` independently (10.3% per type, which is about 15% depolarizing):

```@example qec
ps = (0.01, 0.03, 0.06, 0.1)
for d in (3, 5, 7)
    local code = rotated_surface_code(d)
    local dec = matching_decoder(code)
    local tab = lookup_decoder(code; maxweight = (d - 1) ÷ 2)
    m = [logical_error_rate(code, p; shots = 10_000, decoder = dec)[1] for p in ps]
    l = [logical_error_rate(code, p; shots = 10_000, decoder = tab)[1] for p in ps]
    println("d = ", d, "   matching ", round.(m; sigdigits = 2), "   lookup ", round.(l; sigdigits = 2))
end
```

(`p = 0.01, 0.03, 0.06, 0.1`.) Matching is not optimal: it treats a `Y`
error as an unrelated `X` and `Z`, and correlated decoding does better.

## What is not here yet

* Circuit-level noise and repeated syndrome rounds, where the matching graph
  gains a time direction.
* Encoders for non-CSS codes as circuits (the tableau prepares them by
  measurement).
* `T`-count optimisation as Reed–Muller decoding of phase polynomials.
