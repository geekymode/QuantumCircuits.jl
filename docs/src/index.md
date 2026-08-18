# QuantumCircuits.jl

A small, dependency-free Julia package for building quantum circuits and
decomposing structured unitaries into elementary gates.

The theme of this first release is **Gray coding**: order the computational
basis so consecutive states differ in exactly one bit, and the CNOT ladders in
a decomposition collapse to a single CNOT per step. The
[Gray coding](@ref "Gray coding and efficient quantum circuit decomposition")
page has the derivation and the reason it matters.

## Installation

Only stdlib dependencies, so setup is instant:

```julia
julia> using Pkg; Pkg.add(url="https://github.com/geekymode/QuantumCircuits.jl")
julia> using QuantumCircuits
```

Or clone the repo and work inside it with `julia --project=.`.

## Quick start

```julia
using QuantumCircuits

c = Circuit(2)
push!(c, H(), 1)
push!(c, CNOT(), 1, 2)
statevector(c)        # (|00⟩ + |11⟩)/√2
matrix(c)             # the 4×4 unitary
```

## What it buys you

| Task | naive | Gray-code (this package) |
|---|---|---|
| Uniformly controlled rotation, `k` controls | ``O(k 2^k)`` CNOTs | ``2^k`` |
| Diagonal unitary on `n` qubits | ``O(n 2^n)`` | ``2^n - 2`` |
| Arbitrary state preparation | ``O(n 2^n)`` | ``2^{n+1} - 4`` |

### Uniformly controlled (multiplexed) rotations

"Rotate the target by `α[j+1]` when the controls read `j`", in ``2^k``
rotations and ``2^k`` CNOTs, with no multi-controlled gates and no ancillas:

```julia
julia> multiplexed_rz([0.1, 0.5, -1.2, 2.0], [1, 2], 3)
Circuit(3 qubits, 8 gates, 4 CNOTs)
q1: ──────────────────────●──────────────────────●─
q2: ──────────●───────────│──────────●───────────│─
q3: ─RZ(0.35)─⊕─RZ(-0.90)─⊕─RZ(0.70)─⊕─RZ(-0.05)─⊕─
```

The CNOT controls follow the Gray flip sequence; the rotation angles come from
a Walsh–Hadamard transform of the branch angles ([`multiplex_angles`](@ref)).

### Diagonal unitaries

```julia
c = diagonal([0.0, 0.3, 1.1, -0.7])                  # diag(exp(im*φ))
matrix(c) ≈ Diagonal(cis.([0.0, 0.3, 1.1, -0.7]))    # true, global phase included
```

### State preparation

```julia
julia> c = prepare_state([1, 2, 3im, 4])
Circuit(2 qubits, 10 gates, 4 CNOTs)
q1: ─RY(2.30)──────────●──────────●─RZ(0.79)───────────●──────────●─
q2: ──────────RY(2.03)─⊕─RY(0.18)─⊕──────────RZ(-0.79)─⊕─RZ(0.79)─⊕─
global phase: 0.3927 rad

julia> statevector(c) ≈ normalize([1, 2, 3im, 4])
true
```

## Conventions

* **Qubit 1 is the most significant bit** of a basis index: ``|q_1 q_2 \dots q_n\rangle``
  has index ``\sum_q b_q 2^{n-q}``. The top wire in a drawing is the leftmost
  bit in a ket.
* **`controls[1]` is the most significant bit** of a multiplexor's branch index
  `j`, matching the same reading order.
* Rotations are ``R(\theta) = \exp(-i\theta P/2)``.
* Circuits carry an explicit `global_phase`, so [`matrix`](@ref) and
  [`statevector`](@ref) are exact rather than correct-up-to-phase — which is
  what you need once a block is used as the target of a control.

## Running the examples and tests

```
julia --project=. examples/demo.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```
