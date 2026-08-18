# Applications of Gray coding

```@meta
CurrentModule = QuantumCircuits
```

Fifteen places the one-bit-at-a-time property pays, from the code itself as a
circuit up to Hamiltonian simulation. Every number on this page is measured
from a circuit that is checked against its own definition in the same block —
nothing is quoted.

`examples/gray_applications.jl` in the repository runs all of it in one go.

```@setup apps
using QuantumCircuits, CairoMakie, LinearAlgebra, Random
CairoMakie.activate!(type = "png", px_per_unit = 2)
Random.seed!(20260818)
randu(n) = (F = qr(randn(ComplexF64, n, n)); Matrix(F.Q) * Diagonal(cis.(2π .* rand(n))))
```

## 1. The code itself

```@example apps
using QuantumCircuits
[(i, bits(i, 3), bits(gray(i), 3), i == 0 ? -1 : gray_flip_position(i, 3)) for i in 0:7]
```

Two properties carry everything downstream: which bit flips is
`trailing_zeros(i)` — free to compute — and the walk closes into a cycle, so
every bit flips an even number of times.

```@example apps
foldl(⊻, 1 << p for p in gray_flip_positions(3)) == 0     # the cycle closes
```

```@example apps
graycodefigure(4)
```

## 2. The code *is* a circuit

`gray(x) = x ⊻ (x >> 1)` is linear over GF(2), so the map `|x⟩ ↦ |gray(x)⟩` is
pure CNOT — `n-1` of them, no ancillas, no multi-qubit controls. Bit `b` of the
output is `x_b ⊻ x_{b+1}`, one CNOT each, swept upwards so every source bit is
still untouched when it is read.

```@example apps
circuitfigure(gray_encoder(4); title = "gray_encoder(4): |x⟩ ↦ |gray(x)⟩")
```

```@example apps
M = matrix(gray_encoder(4))
all(M[gray(x)+1, x+1] ≈ 1 for x in 0:15), matrix(gray_decoder(4)) * M ≈ I(16)
```

## 3. A Gray counter

[`gray_increment`](@ref) is decode → increment → encode, and maps
`|gray(i)⟩ ↦ |gray(i+1)⟩`. A binary counter can change every bit at once
(`0111 → 1000`); a Gray counter never changes more than one, which is why
hardware sampling it mid-update still reads a valid code word.

```@example apps
M = matrix(gray_increment(3))
all(hamming(x, argmax(abs.(M[:, x+1])) - 1) == 1 for x in 0:7)
```

## 4. Uniformly controlled rotations

The central primitive: `2ᵏ` CNOTs regardless of how the branches are spread —
see [Gray coding](@ref "Gray coding and efficient quantum circuit decomposition")
for the derivation.

```@example apps
[(k, count_cnots(multiplexed_rz(randn(1 << k), 1:k, k+1)), 2k * (1 << k)) for k in 1:6]
```

The three columns are `k`, the Gray-code CNOT count, and what one
multi-controlled rotation per branch would cost.

## 5. Uniformly controlled arbitrary one-qubit gates

Euler-decompose each branch and run the three angle families through three
multiplexors; the leftover per-branch phases are a diagonal on the controls.

```@example apps
Us = [randu(2) for _ in 1:4]
c = multiplexed_1q(Us, [1, 2], 3)
matrix(c) ≈ cat(Us...; dims = (1, 2)), count_cnots(c)
```

## 6, 7. Diagonal unitaries and state preparation

```@example apps
[(n, count_cnots(diagonal(randn(1 << n))), (1 << n) - 2) for n in 2:7]
```

```@example apps
[(n, count_cnots(prepare_state(randn(ComplexF64, 1 << n))), 2 * (1 << n) - 4) for n in 2:7]
```

## 8. Two-level unitaries: routing to Gray adjacency

A Givens factor acts on `|a⟩` and `|b⟩`, which generally differ in several bits.
CNOTs fanning out from the lowest differing bit bring them one bit apart, at
which point a single multi-controlled gate does the work; the routing is then
undone.

```@example apps
V = randu(2)
c = Circuit(3); two_level!(c, TwoLevel(1, 6, V), 3)
circuitfigure(c; title = "|001⟩ ↔ |110⟩ — three bits apart, routed to one")
```

```@example apps
[(bits(a, 3), bits(b, 3), hamming(a, b), count_cnots((c = Circuit(3);
    two_level!(c, TwoLevel(a, b, V), 3); c))) for (a, b) in ((0, 7), (2, 3), (0, 4))]
```

## 9. Arbitrary unitary synthesis

[`two_level_decompose`](@ref) gives the factors, [`two_level!`](@ref) routes each
one, and `lower = true` expands the multi-controlled gates so nothing bigger
than a CNOT survives.

```@example apps
U = randu(8)
lo = synthesize_unitary(U; lower = true)
matrix(lo) ≈ U, length(lo), count_cnots(lo), all(length(op.qubits) <= 2 for op in lo.ops)
```

## 10. Multi-controlled gates, ancilla-free

Write the AND of the control bits in the parity basis:

    x₁ x₂ ⋯ xₙ  =  2^{-(n-1)} · Σ_{S ≠ ∅} (-1)^{|S|+1} · (⊕_{i ∈ S} xᵢ)

so `C^n(U)` factors into one controlled-`V^{±1}` per non-empty subset, each
controlled on the *parity* of that subset rather than on individual wires. Those
`2ⁿ-1` parities are then walked in Gray order — one CNOT per step.

```@example apps
circuitfigure(multicontrolled(matrix(X()), [1, 2], 3); title = "Toffoli, from V² = X")
```

```@example apps
[(k, length(c) - count_cnots(c), count_cnots(c)) for k in 2:5
 for c in (multicontrolled(randu(2), 1:k, k+1),)]
```

Controlled-`V` gates and CNOTs, against `2ᵏ-1` and `2ᵏ-2`. Exponential in the
number of controls — that is the price of using no ancillas.

## 11, 12. Phase polynomials and parity networks

A diagonal unitary is `exp(i φ(x))`, and `φ` expands uniquely over parity
functions. Emitting one CNOT-ladder gadget per term rebuilds most of the same
parity each time; a **parity network** carries a running parity on an anchor
wire so each Gray-order step costs one CNOT.

```@example apps
[(n, count_cnots(synthesize(pp; order = :gadgets)),
     count_cnots(synthesize(pp; order = :gray)), (1 << n) - 2)
 for n in 2:7 for pp in (phase_polynomial(randn(1 << n)),)]
```

On dense polynomials the network lands exactly on the optimal `2ⁿ-2`. It is not
magic, though — a chain of weight-2 terms already costs two CNOTs apiece and
there is nothing to save:

```@example apps
n = 6
co = zeros(1 << n)
for S in (0b000011, 0b000110, 0b001100, 0b011000, 0b110000); co[S+1] = randn(); end
pp = PhasePolynomial(n, co)
count_cnots(synthesize(pp; order = :gadgets)), count_cnots(synthesize(pp; order = :gray))
```

## 13. Walsh-series truncation

Keeping only the largest coefficients of a smooth phase trades accuracy for
gates ([`truncate_terms`](@ref)).

```@example apps
N = 64
φ = [2.0 * sin(2π * x / N) + 0.4 * (x / N)^2 for x in 0:N-1]
pp = phase_polynomial(φ)
[(nterms(tp), count_cnots(c), round(maximum(abs, phases(tp) .- φ), sigdigits = 3),
  round(gate_fidelity(matrix(c), Diagonal(cis.(φ))), digits = 6))
 for k in (2, 4, 8, 16, 32, N) for tp in (truncate_terms(pp, k),) for c in (synthesize(tp),)]
```

Eight of 47 terms — 14 CNOTs instead of 52 — already reaches fidelity 0.998.

## 14. LCU and QROM: sweeping an address register

A branch of a SELECT fires on `|1⟩`, so address `j` must be reached by `X`-ing
every control where `j` has a zero. Consecutive **Gray** addresses differ in one
bit, so each step costs a single `X` instead of the Hamming distance between
neighbouring integers.

```@example apps
[(k, count_gates(select(Us, 1:k, [k+1]; order = :gray), :X),
     count_gates(select(Us, 1:k, [k+1]; order = :natural), :X))
 for k in 2:7 for Us in ([randu(2) for _ in 1:(1 << k)],)]
```

`2ᵏ` versus `2ᵏ⁺¹` — the saving approaches half as the address widens.

```@example apps
circuitfigure(select([matrix(X()), matrix(H()), matrix(Z()), matrix(T())], [1, 2], [3]);
              title = "SELECT over a 2-bit address, Gray-swept")
```

## 15. Hamiltonian simulation

Diagonal (`I`/`Z`-only) terms all commute, so they can be pooled and synthesised
together through the parity network. How much that buys depends entirely on how
dense the diagonal part is.

A transverse-field Ising Hamiltonian is all weight-2 couplings — nothing to
pool:

```@example apps
function tfim(n; J = 1.0, h = 0.6)
    t = Pair{String,Float64}[]
    for i in 1:n-1
        s = collect("I"^n); s[i] = 'Z'; s[i+1] = 'Z'; push!(t, String(s) => J)
    end
    for i in 1:n
        s = collect("I"^n); s[i] = 'X'; push!(t, String(s) => h)
    end
    t
end
[(n, count_cnots((c = Circuit(n); trotter_step!(c, tfim(n), 0.05; parity_network = false); c)),
     count_cnots((c = Circuit(n); trotter_step!(c, tfim(n), 0.05; parity_network = true); c)))
 for n in 4:7]
```

A dense diagonal cost operator — a QAOA phase separator — pools completely:

```@example apps
function qaoa_cost(n)
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
[(n, count_cnots((c = Circuit(n); trotter_step!(c, qaoa_cost(n), 0.05; parity_network = false); c)),
     count_cnots((c = Circuit(n); trotter_step!(c, qaoa_cost(n), 0.05; parity_network = true); c)))
 for n in 3:6]
```

258 CNOTs down to 62 at `n = 6`, for exactly the same operator — the terms
commute, so pooling them changes nothing but the gate count.

## What Gray coding does *not* buy

Reordering Trotter terms so consecutive supports are Gray-adjacent
([`gray_order`](@ref)) sounds like it should let neighbouring CNOT ladders
cancel. It does not: each term's gadget is bracketed by its own basis-change
gates, which sit on the same wires as the ladder ends, so
[`cancel_adjacent_cnots!`](@ref) never sees two CNOTs meet. Measured on dense
random Hamiltonians the reordering is neutral to slightly worse. The win comes
from pooling commuting terms, not from reordering non-commuting ones.

## Reference

```@docs
multicontrolled
multicontrolled!
matrix_root
gray_encoder
gray_decoder
increment
increment!
gray_increment
select
select!
support_mask
gray_order
truncate_terms
```
