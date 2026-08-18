# Linear algebra and decompositions

```@meta
CurrentModule = QuantumCircuits
```

Circuit synthesis is applied matrix factorisation. This page covers the
numerical layer: the three bases worth changing to, and the three
factorisations that turn a matrix into gates.

## Three bases

A quantum object is usually handed to you in the computational basis, and the
work happens in some other one.

| Basis | Change of basis | What it makes obvious |
|---|---|---|
| computational | — | amplitudes, block structure |
| parity / Walsh | [`fwht!`](@ref) | which parity functions a diagonal phase is built from |
| Pauli | [`pauli_decompose`](@ref) | which Hamiltonian terms to Trotterise |

The Walsh–Hadamard transform is the workhorse. `W[s+1,j+1] = (-1)^(s·j)` is its
own inverse up to a factor of `N`, and [`fwht!`](@ref) applies it in
`O(N log N)`:

```@example math
using QuantumCircuits
v = [1.0, 2, 3, 4]
fwht(v)
```

```@example math
fwht(fwht(v)) == 4 .* v      # W² = N·I
```

It is doing the real work inside [`multiplex_angles`](@ref): the branch angles
of a multiplexor are the Walsh coefficients of the applied rotation angles,
read out in Gray order.

The Pauli basis is orthogonal under the Hilbert–Schmidt inner product, so
expanding an operator in it is an exact change of basis:

```@example math
Hmat = [1.0 0.3; 0.3 -0.5]      # note: `H` itself is the Hadamard constructor
pauli_decompose(Hmat)
```

```@example math
pauli_recompose(pauli_decompose(Hmat)) ≈ Hmat
```

## Euler (ZYZ) decomposition

Every one-qubit unitary is three rotations and a phase, `U = e^{iα} R_z(β) R_y(γ) R_z(δ)`:

```@example math
using LinearAlgebra
α, β, γ, δ = zyz(matrix(H()))
cis(α) * matrix(RZ(β)) * matrix(RY(γ)) * matrix(RZ(δ)) ≈ matrix(H())
```

[`decompose_1q`](@ref) turns that straight into a circuit. The degenerate cases
(`γ ≈ 0`, `γ ≈ π`, where only `β ± δ` is determined) are pinned rather than left
to floating-point luck.

## Two-level (Givens) decomposition

Any `N × N` unitary is a product of rotations in two-dimensional coordinate
planes — QR by plane rotations, and the classical route to synthesising an
arbitrary unitary:

```@example math
U = matrix(prepare_state([1, 2, 3im, 4]))
factors = two_level_decompose(U)
length(factors), first(factors)
```

```@example math
reduce(*, matrix(t, 4) for t in factors) ≈ U
```

Each factor acts on the span of two basis states `|a⟩, |b⟩`. Those generally
differ in several bits, and this is exactly where the Gray code earns its keep:
[`two_level!`](@ref) fans CNOTs out from the lowest differing bit until `a` and
`b` are one bit apart, applies a single multi-controlled gate, and undoes the
routing.

```@example math
matrix(synthesize_unitary(U)) ≈ U
```

## Demultiplexing

The recursive step behind the quantum Shannon decomposition. A block-diagonal
(multiplexed) unitary splits into two smaller unitaries and one uniformly
controlled `RZ`, at the cost of a single eigendecomposition:

```math
U_1 \oplus U_2 \;=\; (I \otimes V)\,(D \oplus D^\dagger)\,(I \otimes W),
\qquad U_1 U_2^\dagger = V D^2 V^\dagger,\quad W = D V^\dagger U_2
```

```@example math
U1, U2 = matrix(RY(0.4)), matrix(RZ(1.1))
V, θ, W = demultiplex(U1, U2)
mid = Diagonal(vcat(cis.(-θ./2), cis.(θ./2)))
kron(I(2), V) * mid * kron(I(2), W) ≈ cat(U1, U2; dims=(1,2))
```

`D ⊕ D†` is a Gray-code multiplexor targeting the *control* wire — so this
factorisation bottoms out in the primitive from the
[Gray coding](@ref "Gray coding and efficient quantum circuit decomposition")
page. [`multiplexed_1q`](@ref) uses the same idea to realise a uniformly
controlled arbitrary one-qubit gate in `4·2ᵏ - 2` CNOTs.

## Phase polynomials and parity networks

A diagonal unitary is `exp(i φ(x))` for a real function on bitstrings, and every
such function expands uniquely over parity functions:

```math
\varphi(x) \;=\; \sum_{S} a_S \, (-1)^{S \cdot x}
```

Each term is one `exp(i a_S Z_S)` — a *phase gadget*: a CNOT ladder computing
the parity of `S` onto a wire, an `RZ`, and the ladder undone.

```@example math
φ = randn(16)
pp = phase_polynomial(φ)
phases(pp) ≈ φ
```

Emitting one independent gadget per term is wasteful, because consecutive
ladders re-derive most of the same parity. A **parity network** keeps the
running parity on an anchor wire, so moving from one term to the next costs one
CNOT per differing bit — and visiting the terms in Gray order keeps that
difference at one bit:

```@example math
naive = synthesize(pp; order=:gadgets)
gray  = synthesize(pp; order=:gray)
count_cnots(naive), count_cnots(gray)
```

For a dense polynomial the parity network lands on exactly `2ⁿ - 2`, matching
the optimal recursive construction in [`diagonal`](@ref). Both are correct to
the last phase:

```@example math
matrix(gray) ≈ Diagonal(cis.(φ))
```

## Hamiltonian simulation

Putting the Pauli basis and phase gadgets together gives a Trotter step
directly from a Hamiltonian matrix:

```@example math
Hm = randn(4, 4); Hm = Hm + Hm'
terms = pauli_decompose(Hm)
c = Circuit(2)
trotter_step!(c, terms, 0.05)
gate_fidelity(matrix(c), exp(-im * 0.05 * Hm))
```

## State measures

```@example math
bell = statevector((c = Circuit(2); push!(c, H(), 1); push!(c, CNOT(), 1, 2); c))
entanglement_entropy(bell, 1)
```
