# Gray coding and efficient quantum circuit decomposition

## 1. The code

The reflected binary Gray code is the map

```
gray(i) = i ⊻ (i >> 1)
```

It is a bijection on `0 … 2ⁿ-1`, and its defining property is that
`gray(i)` and `gray(i+1)` **differ in exactly one bit**. For `n = 3`:

| `i` | binary | `gray(i)` | bit flipped from the previous row |
|----:|:------:|:---------:|:---------------------------------:|
| 0 | 000 | 000 | — |
| 1 | 001 | 001 | 0 |
| 2 | 010 | 011 | 1 |
| 3 | 011 | 010 | 0 |
| 4 | 100 | 110 | 2 |
| 5 | 101 | 111 | 0 |
| 6 | 110 | 101 | 1 |
| 7 | 111 | 100 | 0 |

Two facts do all the work downstream:

1. **Which bit flips is free to compute.** Going from `gray(i-1)` to `gray(i)`
   flips bit `trailing_zeros(i)` — that is `gray_flip_position(i)`. No search,
   no table.
2. **The walk closes into a cycle.** `gray(2ⁿ-1) = 100…0` is one bit away from
   `gray(0) = 0`, so a full cyclic walk flips every bit an even number of
   times. `gray_flip_position(i, k)` returns `k-1` for the wrap-around step,
   and `gray_flip_positions(k)` gives the whole cycle.

Fact 1 makes the synthesis loop O(1) per gate. Fact 2 is what makes the
resulting circuit *unitary-correct* rather than merely close: everything the
construction scribbles on the control register is erased by the time the walk
returns to the start.

## 2. Why a quantum compiler cares

Almost every interesting unitary is specified as a *table of 2ⁿ things* — one
rotation angle per control branch, one phase per basis state, one amplitude per
basis state. A structurally naive compilation treats those 2ⁿ entries as 2ⁿ
independent multi-controlled gates. Each multi-controlled gate costs `O(n)`
CNOTs (or needs ancillas), so you pay `O(n · 2ⁿ)`.

But the entries are not independent — they are indexed by bitstrings, and you
get to choose the *order* in which you process them. Choose Gray-code order and
consecutive entries differ in one bit, so the transition between them costs a
**single CNOT** instead of an uncompute/recompute ladder. That is the entire
trick, and it recurs everywhere:

| Task | naive | Gray-code |
|---|---|---|
| Uniformly controlled rotation, `k` controls | `O(k·2ᵏ)` CNOTs | `2ᵏ` CNOTs |
| Diagonal unitary on `n` qubits | `O(n·2ⁿ)` | `2ⁿ - 2` |
| Arbitrary state preparation | `O(n·2ⁿ)` | `2ⁿ⁺¹ - 4` |
| Multi-controlled `U` (Barenco et al. Lemma 7.5) | — | `O(n²)`, ancilla-free |
| Trotterised Pauli-string exponentials | ladder per term | one CNOT per term |

On hardware where two-qubit gates dominate both the error budget and the
runtime, dropping the factor of `n` and every ancilla is the whole ballgame.

## 3. The uniformly controlled rotation

The central primitive, sometimes called a *multiplexed rotation* or
*multiplexor*:

```
UCR(α) = Σ_j |j⟩⟨j| ⊗ R(α_j),     j = 0 … 2ᵏ-1
```

with `k` control qubits and one target: "if the controls read `j`, rotate the
target by `α_j`". `R` is `RZ` or `RY`.

### The construction

Alternate rotations on the target with CNOTs into the target, `2ᵏ` of each:

```
R(θ₁) ─ CNOT(c₁) ─ R(θ₂) ─ CNOT(c₂) ─ … ─ R(θ_{2ᵏ}) ─ CNOT(c_{2ᵏ})
```

where `cᵢ` is the control qubit named by the Gray-code flip position at step
`i`. No multi-qubit control anywhere.

### Why it works

Every gate is diagonal on the control register, so fix the controls in a basis
state `|j⟩` and ask what happens to the target. Each CNOT is then just `X^b`,
with `b` the value of the relevant bit of `j`. The key identity is

```
X R(θ) X = R(-θ)          for R ∈ {RX, RY, RZ}
```

so a rotation sandwiched between an odd number of `X`s comes back with its
angle negated. Before step `i`, the target has been hit by the `X`s from
control bits `c₁ … c_{i-1}`, so the sign of `θᵢ` is `(-1)` to the parity of
those bits of `j`. The mask of bits toggled so far is exactly the accumulated
XOR of the flip masks, which is the definition of the Gray code:

```
mask after i-1 steps = gray(i-1)
```

All the rotations commute (same axis), so they simply add:

```
α_j = Σᵢ (-1)^(gray(i-1) · j) · θᵢ
```

where `·` is the bitwise dot product mod 2 (`parity` in this package). And by
fact 2, the mask after all `2ᵏ` steps is `0`: the CNOTs cancel and the control
register is left untouched. The circuit really is `UCR(α)`, not `UCR(α)` times
some leftover entanglement.

### Inverting the transform

Write the relation as `α = M θ` with `M[j+1, i] = (-1)^(gray(i-1) · j)`. The
rows of `M` are the characters of `(Z₂)ᵏ` evaluated at the points `gray(i-1)`,
which run over the whole group — so `M` is a column-permuted Walsh–Hadamard
matrix and

```
M Mᵀ = 2ᵏ I     ⟹     θ = Mᵀ α / 2ᵏ
```

Inverting a `2ᵏ × 2ᵏ` system costs nothing but a transpose. In this package
that is `multiplex_angles`:

```julia
julia> using QuantumCircuits

julia> c = multiplexed_rz([0.1, 0.5, -1.2, 2.0], [1, 2], 3)
Circuit(3 qubits, 8 gates, 4 CNOTs)
q1: ──────────────────────●──────────────────────●─
q2: ──────────●───────────│──────────●───────────│─
q3: ─RZ(0.35)─⊕─RZ(-0.90)─⊕─RZ(0.70)─⊕─RZ(-0.05)─⊕─
```

Four branches, four CNOTs — and it stays four per branch as `k` grows, versus
the `~2k` CNOTs *per branch* a controlled-rotation-per-branch compilation would
emit.

Note the CNOT control pattern reading left to right: q2, q1, q2, q1 — the Gray
flip sequence `0, 1, 0, 1` for `k = 2`, with the last step being the
wrap-around that closes the cycle.

## 4. Diagonal unitaries: `2ⁿ - 2` CNOTs

`diagonal(φ)` builds `diag(exp(i φ_j))`. Split the basis index `j = (j', b)` on
the last qubit and look at each pair of phases `φ(j',0)`, `φ(j',1)`:

* their **difference** is a relative phase between `|0⟩` and `|1⟩` of the last
  qubit, conditioned on `j'` — that is precisely a uniformly controlled `RZ`,
  one multiplexor, `2ⁿ⁻¹` CNOTs;
* their **mean** does not involve the last qubit at all — it is a diagonal
  unitary on `n-1` qubits. Recurse.

```
CNOTs(n) = CNOTs(n-1) + 2ⁿ⁻¹,  CNOTs(1) = 0    ⟹    CNOTs(n) = 2ⁿ - 2
```

The recursion bottoms out at a global phase, which `Circuit` tracks in
`global_phase` so that `matrix(c)` is exactly `diag(exp(i φ))` — not merely
equal up to phase. That matters as soon as the diagonal appears inside a larger
controlled block, where the "unobservable" global phase becomes an observable
relative phase.

## 5. State preparation

`prepare_state(a)` maps `|0…0⟩` to any normalised amplitude vector in two
passes, both Gray-coded:

1. **Magnitudes.** For each level `ℓ = 1 … n`, a uniformly controlled `RY` with
   controls `1 … ℓ-1` and target `ℓ`. The branch angle for prefix `p` is
   `2·atan(‖subtree p1‖, ‖subtree p0‖)`: split the amplitude mass of the
   subtree between the two children. This is a binary-tree walk over the
   amplitude vector, and level `ℓ` costs `2^(ℓ-1)` CNOTs — `2ⁿ - 2` in total.
2. **Phases.** One `diagonal(angle.(a))`, another `2ⁿ - 2` CNOTs.

Total `2ⁿ⁺¹ - 4` CNOTs, ancilla-free. (The theoretical lower bound is
`2ⁿ - n - 1`; the Gray-code construction is within a small constant factor of
it and is what most compilers actually ship.)

```julia
julia> c = prepare_state([1, 2, 3im, 4])
Circuit(2 qubits, 10 gates, 4 CNOTs)
q1: ─RY(2.30)──────────●──────────●─RZ(0.79)───────────●──────────●─
q2: ──────────RY(2.03)─⊕─RY(0.18)─⊕──────────RZ(-0.79)─⊕─RZ(0.79)─⊕─
global phase: 0.3927 rad

julia> round.(statevector(c), digits=4)
4-element Vector{ComplexF64}:
  0.1826 + 0.0im
  0.3651 - 0.0im
    -0.0 + 0.5477im
  0.7303 + 0.0im
```

## 6. Where else Gray coding shows up

Not yet implemented here, but the same lever:

* **Multi-controlled gates without ancillas.** Barenco et al. (1995), Lemma
  7.5: `Cⁿ(U)` from `2ⁿ` controlled-`V` gates where `V^(2ⁿ⁻¹) = U`, with the
  control patterns enumerated in Gray-code order so consecutive terms differ by
  one CNOT.
* **Trotterisation.** `exp(i θ P)` for a Pauli string `P` is a CNOT ladder onto
  a parity qubit, an `RZ`, and the ladder undone. Ordering the terms of a
  Hamiltonian so consecutive strings are Gray-adjacent lets the ladders
  partially cancel — the standard "Pauli-string ordering" optimisation.
* **Binary-to-unary / Gray-code qubit encodings** for simulating fermionic or
  bosonic modes, where the one-bit-difference property means a hop between
  adjacent occupation levels touches a single qubit.
* **Quantum arithmetic and QROM/LUT loading**, where an address register is
  swept in Gray order so the address decoder updates incrementally.

## References

* M. Möttönen, J. J. Vartiainen, V. Bergholm, M. M. Salomaa, *Transformation of
  quantum states using uniformly controlled rotations*, quant-ph/0407010 (2004).
* V. Bergholm, J. J. Vartiainen, M. Möttönen, M. M. Salomaa, *Quantum circuits
  with uniformly controlled one-qubit gates*, PRA 71, 052330 (2005).
* A. Barenco et al., *Elementary gates for quantum computation*, PRA 52, 3457
  (1995) — Lemma 7.5 for the Gray-code multi-controlled construction.
* A. Shende, S. Bullock, I. Markov, *Synthesis of quantum-logic circuits*, IEEE
  TCAD 25, 1000 (2006).
