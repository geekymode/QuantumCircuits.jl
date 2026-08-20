# QuantumCircuits.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://geekymode.github.io/QuantumCircuits.jl/)
[![CI](https://github.com/geekymode/QuantumCircuits.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/geekymode/QuantumCircuits.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

**Documentation: <https://geekymode.github.io/QuantumCircuits.jl/>**

A dependency-free Julia package for building quantum circuits, decomposing
structured unitaries into elementary gates, and drawing the result.

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

### Linear algebra and matrix decompositions

Synthesis is applied matrix factorisation, so the factorisations are first-class:

```julia
zyz(matrix(H()))                    # Euler: U = e^{iα} RZ(β) RY(γ) RZ(δ)
two_level_decompose(U)              # QR by plane rotations (Givens)
synthesize_unitary(U)               # ... turned into gates by Gray-code routing
demultiplex(U1, U2)                 # U₁⊕U₂ = (I⊗V)(D⊕D†)(I⊗W), one eigendecomposition
multiplexed_1q(Us, [1,2], 3)        # uniformly controlled arbitrary SU(2)

fwht(v)                             # Walsh–Hadamard transform, O(N log N)
pauli_decompose(Hmat)               # exact expansion in the Pauli basis
trotter_step!(c, terms, dt)         # ... straight into a Trotter circuit
entanglement_entropy(ψ, 1)          # Schmidt values across a cut
```

Phase polynomials tie the algebra back to the Gray code. A diagonal unitary is
`exp(i φ(x))`, and `φ` expands uniquely over parity functions; each term is a
CNOT-ladder phase gadget. Visiting the terms in Gray order and carrying a
running parity turns each ladder into a single CNOT:

```julia
pp = phase_polynomial(randn(64))
count_cnots(synthesize(pp; order = :gadgets))   # 114 — one ladder per term
count_cnots(synthesize(pp; order = :gray))      #  62 — parity network = 2ⁿ-2, optimal
```

### The hardest case: an arbitrary unitary

Hand the compiler a `2ⁿ × 2ⁿ` unitary with no structure at all and you get the
hardest of the standard problems. The **quantum Shannon decomposition** is the
answer, and its whole skeleton is Gray-code multiplexors — one level peels the
top wire with a cosine–sine split plus two demultiplexings, emitting three
multiplexors and four `(n-1)`-qubit sub-problems:

```julia
U = rand_unitary(8)
c = qsd(U)                      # 36 CNOTs, matrix(c) ≈ U exactly
count_cnots(synthesize_unitary(U; lower = true))   # 98 — the two-level route

cosine_sine(U)                  # U = (L₁⊕L₂)·[C -S; S C]·(R₁⊕R₂)†
csdfigure(rand_unitary(16))     # ... and what that looks like
qsdfigure(3)                    # the recursion, with CNOT counts per level
```

`[C -S; S C]` *is* a uniformly controlled `RY` on the top wire, so every CNOT in
the output comes from a Gray-code walk. Costs `(3/4)·4ⁿ - (3/2)·2ⁿ`: 36 at
`n=3`, 168 at `n=4`, against `O(n·4ⁿ)` for the two-level method. See the
[Hardest case](https://geekymode.github.io/QuantumCircuits.jl/dev/shannon/) page
— including the numerical trap that makes a naive cosine–sine implementation
return non-unitary factors on block-diagonal input.

### More Gray-code applications

Fifteen worked applications, all measured and checked — see the
[Applications](https://geekymode.github.io/QuantumCircuits.jl/dev/applications/)
page, or run `julia --project=. examples/gray_applications.jl`:

```julia
gray_encoder(4)                     # |x⟩ ↦ |gray(x)⟩ — the code IS a CNOT circuit, n-1 gates
gray_increment(4)                   # a Gray counter: one bit changes per tick
multicontrolled(U, 1:4, 5)          # ancilla-free C^n(U) from a Gray walk over parities
select(Us, 1:k, [k+1])              # LCU/QROM SELECT, address swept in Gray order
truncate_terms(pp, 8)               # Walsh-series truncation: accuracy for gates
synthesize_unitary(U; lower = true) # arbitrary unitary, down to CNOTs and 1-qubit gates
```

A few of the measured results:

| Task | naive | Gray-code |
|---|---|---|
| Arbitrary 3-qubit unitary | 98 CNOTs (two-level) | **36** (Shannon) |
| `C⁵(U)`, no ancillas | multi-controlled decomposition | 31 controlled-`V` + 30 CNOTs |
| Dense diagonal, `n=7` | 240 CNOTs (gadget per term) | **126** = `2ⁿ-2`, optimal |
| QAOA cost operator, `n=6` | 258 CNOTs | **62**, same operator |
| SELECT over 7 address bits | 254 X gates | **140** |
| Smooth phase, 8 of 47 terms | 52 CNOTs exact | **14** at fidelity 0.998 |

And one thing it does *not* buy: reordering Trotter terms by Gray adjacency
does not reduce CNOTs, because each term's basis-change gates block ladder
cancellation. The win comes from pooling commuting terms, not reordering
non-commuting ones — the docs show the measurement.

### Illustrations

Plotting is a package extension — it loads when a Makie backend is present and
costs nothing otherwise:

```julia
using Pkg; Pkg.add("CairoMakie")
using CairoMakie, QuantumCircuits

circuitfigure(prepare_state([1, 2, 3im, 4]))    # publication-quality diagram
matrixfigure(U; part = :phase)                  # block structure of a unitary
graycodefigure(4)                               # the walk, with the flipped bit ringed
costfigure(2:9)                                 # CNOT count vs naive compilation
save("circuit.pdf", circuitfigure(c))           # vector output for a paper
```

Every figure takes `theme = :light | :dark`. See the
[Illustrations](https://geekymode.github.io/QuantumCircuits.jl/dev/plots/) page.

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
| Gray-code synthesis | `multiplex_angles`, `multiplex_matrix`, `multiplexed_rotation!`, `multiplexed_ry`, `multiplexed_rz`, `diagonal`, `prepare_state` |
| Linear algebra | `fwht`, `walsh_matrix`, `pauli`, `pauli_decompose`, `pauli_recompose`, `embed`, `kron_n`, `is_unitary`, `gate_fidelity`, `schmidt_values`, `entanglement_entropy` |
| Matrix decompositions | `zyz`, `decompose_1q`, `TwoLevel`, `two_level_decompose`, `two_level!`, `synthesize_unitary`, `demultiplex`, `multiplexed_1q` |
| Phase polynomials | `PhasePolynomial`, `phase_polynomial`, `phases`, `support`, `synthesize`, `phase_gadget!`, `pauli_rotation!`, `trotter_step!`, `cancel_adjacent_cnots!` |
| Shannon decomposition | `qsd`, `cosine_sine`, `CSD`, `csd_angles`, `qsd_cnot_count`, `rand_unitary` |
| Applications | `multicontrolled`, `matrix_root`, `gray_encoder`, `gray_decoder`, `increment`, `gray_increment`, `select`, `support_mask`, `gray_order`, `truncate_terms` |
| Illustrations (Makie ext) | `circuitfigure`, `circuitplot!`, `matrixfigure`, `graycodefigure`, `costfigure`, `csdfigure`, `qsdfigure` |

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

1498 tests. Decompositions are checked against reference matrices built straight
from the definitions (no Gray code in the reference path), including exact
global phase and exact CNOT counts. Plot tests need a Makie backend:

```
QC_TEST_PLOTS=true julia --project=docs -e 'using Pkg; Pkg.test("QuantumCircuits")'
```

## Roadmap

* Two-qubit KAK / Cartan decomposition and optimal 3-CNOT two-qubit synthesis —
  the largest remaining win, taking `qsd` from `(3/4)·4ⁿ` to the standard
  `(9/16)·4ⁿ`
* Full GraySynth term-ordering for sparse phase polynomials (the current parity
  network is greedy: it hits the optimum on dense inputs, not always on sparse)
* More optimisation passes (rotation merging, commutation-aware cancellation)
