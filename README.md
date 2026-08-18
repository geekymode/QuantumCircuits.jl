# QuantumCircuits.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://geekymode.github.io/QuantumCircuits.jl/)
[![CI](https://github.com/geekymode/QuantumCircuits.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/geekymode/QuantumCircuits.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

**Documentation: <https://geekymode.github.io/QuantumCircuits.jl/>**

A small, dependency-free Julia package for building quantum circuits and
decomposing structured unitaries into elementary gates.

The theme of this first release is **Gray coding**: order the computational
basis so consecutive states differ in exactly one bit, and the CNOT ladders in
a decomposition collapse to a single CNOT per step. See
[`docs/graycode.md`](docs/graycode.md) for the derivation and the reason it
matters.

## Getting started

Only stdlib dependencies, so setup is instant:

```bash
git clone https://github.com/geekymode/QuantumCircuits.jl.git
cd QuantumCircuits.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # resolves Manifest.toml
julia --project=. examples/demo.jl                    # guided tour
julia --project=. -e 'using Pkg; Pkg.test()'          # 750 tests
```

Interactive REPL, with the project environment active:

```bash
julia --project=.
```
```julia
julia> using QuantumCircuits
julia> multiplexed_rz([0.1, 0.5, -1.2, 2.0], [1, 2], 3)
```

Press `]` for pkg mode (`pkg> instantiate`, `pkg> test`, `pkg> st`) and
backspace to leave it. To use the package from anywhere instead of only inside
this folder, `dev` it into your default environment:

```julia
julia> using Pkg; Pkg.add(url="https://github.com/geekymode/QuantumCircuits.jl")
```

For an edit–run–edit loop without restarting Julia, `Pkg.add("Revise")` in your
default environment and `using Revise` **before** `using QuantumCircuits`.

## Quick start

```julia
using QuantumCircuits

c = Circuit(2)
push!(c, H(), 1)
push!(c, CNOT(), 1, 2)
statevector(c)        # (|00⟩ + |11⟩)/√2
matrix(c)             # the 4×4 unitary
```

### Gray code

```julia
gray.(0:7)                 # [0, 1, 3, 2, 6, 7, 5, 4]
gray_flip_position(6)      # 1  — bit that changes between gray(5) and gray(6)
gray_flip_positions(3)     # the full cyclic walk: [0,1,0,2,0,1,0,2]
gray_adjacent(0b011, 0b010) # true
```

### Uniformly controlled (multiplexed) rotations — `2ᵏ` CNOTs

"Rotate the target by `α[j+1]` when the controls read `j`", in `2ᵏ` rotations
and `2ᵏ` CNOTs, with no multi-controlled gates and no ancillas:

```julia
julia> multiplexed_rz([0.1, 0.5, -1.2, 2.0], [1, 2], 3)
Circuit(3 qubits, 8 gates, 4 CNOTs)
q1: ──────────────────────●──────────────────────●─
q2: ──────────●───────────│──────────●───────────│─
q3: ─RZ(0.35)─⊕─RZ(-0.90)─⊕─RZ(0.70)─⊕─RZ(-0.05)─⊕─
```

The CNOT controls follow the Gray flip sequence; the rotation angles come from
a Walsh–Hadamard transform of the branch angles (`multiplex_angles`).

### Diagonal unitaries — `2ⁿ - 2` CNOTs

```julia
c = diagonal([0.0, 0.3, 1.1, -0.7])   # diag(exp(im*φ))
matrix(c) ≈ Diagonal(cis.([0.0, 0.3, 1.1, -0.7]))   # true, global phase included
```

### State preparation — `2ⁿ⁺¹ - 4` CNOTs

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

* **Qubit 1 is the most significant bit** of a basis index: `|q1 q2 … qn⟩` has
  index `Σ bit_q · 2^(n-q)`. The top wire in a drawing is the leftmost bit in a
  ket.
* **`controls[1]` is the most significant bit** of a multiplexor's branch index
  `j`, matching the same reading order.
* Rotations are `R(θ) = exp(-i θ P / 2)`.
* Circuits carry an explicit `global_phase`, so `matrix` and `statevector` are
  exact rather than correct-up-to-phase — which is what you need once a block
  is used as the target of a control.

## API

| | |
|---|---|
| Gray code | `gray`, `ungray`, `graycode`, `gray_flip_position`, `gray_flip_positions`, `gray_adjacent`, `gray_walk`, `hamming`, `parity`, `bits` |
| Gates | `Id`, `X`, `Y`, `Z`, `H`, `S`, `Sdg`, `T`, `Tdg`, `RX`, `RY`, `RZ`, `PHASE`, `CNOT`, `CZ`, `SWAP`, `controlled` |
| Circuits | `Circuit`, `push!`, `append!`, `matrix`, `statevector`, `zero_state`, `apply!`, `draw`, `count_cnots`, `count_gates` |
| Decompositions | `multiplex_angles`, `multiplex_matrix`, `multiplexed_rotation!`, `multiplexed_ry`, `multiplexed_rz`, `diagonal`, `prepare_state` |

## Documentation

Published at **<https://geekymode.github.io/QuantumCircuits.jl/>**, rebuilt from
`main` on every push. To build the [Documenter](https://documenter.juliadocs.org)
site locally:

```bash
# one-time: install Documenter + LiveServer into the docs environment
julia --project=docs -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'

# build once, then open in a browser
julia --project=docs docs/make.jl
open docs/build/index.html

# or serve with live reload while editing (http://localhost:8000)
julia --project=docs -e 'using LiveServer; servedocs()'
```

`docs/build/` is gitignored; regenerate it whenever needed.

## Tests

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

Decompositions are checked against reference matrices built straight from the
definitions (no Gray code in the reference path), including exact global phase
and exact CNOT counts.

## Roadmap

* Multi-controlled gates via the Barenco Gray-code construction (ancilla-free)
* Uniformly controlled *general* `SU(2)` gates, and cosine–sine decomposition
  for arbitrary `n`-qubit unitaries
* Pauli-string exponentials with Gray-ordered term scheduling for Trotter steps
* CNOT-count optimisation passes (adjacent-cancellation, rotation merging)
