# Illustrations

```@meta
CurrentModule = QuantumCircuits
```

Plotting is a **package extension**: it loads automatically the moment a Makie
backend is available, and costs nothing when it is not.

```julia
using Pkg; Pkg.add("CairoMakie")
using CairoMakie, QuantumCircuits      # the extension loads here
```

[CairoMakie](https://docs.makie.org/stable/explanations/backends/cairomakie)
gives vector output — `save("circuit.pdf", fig)` produces a figure you can drop
into a paper. GLMakie works too if you want an interactive window.

Every figure takes `theme = :light` (default) or `:dark`, and a `scale`
multiplier.

```@setup plots
using QuantumCircuits, CairoMakie, LinearAlgebra, Random
CairoMakie.activate!(type = "png", px_per_unit = 2)
Random.seed!(20260817)
randu(n) = (F = qr(randn(ComplexF64, n, n)); Matrix(F.Q) * Diagonal(cis.(2π .* rand(n))))
```

## Circuit diagrams

[`circuitfigure`](@ref) packs instructions into columns greedily — gates on
disjoint wires share a column, the way you would draw it by hand.

```@example plots
circuitfigure(multiplexed_rz([0.1, 0.5, -1.2, 2.0], [1, 2], 3))
```

Read the CNOT controls left to right: q2, q1, q2, q1 — the Gray flip sequence
for two controls, with the last step closing the cycle so the control register
comes out untouched.

State preparation, both passes visible — the `RY` cascade distributing
magnitude down the amplitude tree, then the diagonal stamping on phases:

```@example plots
circuitfigure(prepare_state([1, 2, 3im, 4]))
```

Multi-controlled gates and Gray-code routing, from the two-level decomposition:

```@example plots
c = Circuit(3)
two_level!(c, TwoLevel(1, 6, randu(2)), 3)
circuitfigure(c; title = "two-level |001⟩ ↔ |110⟩, routed to be Gray-adjacent")
```

The two CNOTs on each side are the routing: they bring `|001⟩` and `|110⟩` — three
bits apart — to a pair differing in one bit, so a single controlled `V` suffices.

Dark theme:

```@example plots
circuitfigure(diagonal(randn(8)); theme = :dark)
```

To compose several circuits into one figure, use
[`circuitplot!`](@ref) with your own axes.

## Matrix structure

[`matrixfigure`](@ref) shows what a decomposition is actually building. `part`
selects the quantity *and* the kind of colour scale it deserves: magnitude on a
single-hue sequential ramp, signed parts on a diverging ramp with a neutral
midpoint, and phase on a cyclic ramp — because phase wraps, so its colour must
wrap too.

A uniformly controlled one-qubit gate is block diagonal, one `2 × 2` block per
control state:

```@example plots
U = matrix(multiplexed_1q([randu(2) for _ in 1:4], [1, 2], 3))
matrixfigure(U; part = :abs, title = "uniformly controlled SU(2): one block per branch")
```

A diagonal unitary keeps all its information in the phase:

```@example plots
matrixfigure(matrix(diagonal(randn(8))); part = :phase, title = "diagonal unitary: arg U")
```

## The Gray code itself

[`graycodefigure`](@ref) draws the walk: the bit grid over all `2ⁿ` steps with
the changed bit ringed, and the flip-position sequence — `trailing_zeros(i)` —
beneath it.

```@example plots
graycodefigure(4)
```

The bottom panel is the CNOT schedule of a four-control multiplexor read
directly: bit 0 flips on every other step, bit 3 only twice, and step 16 closes
the cycle.

## Cost

```@example plots
costfigure(2:9)
```

## Reference

```@docs
circuitfigure
circuitplot!
matrixfigure
graycodefigure
costfigure
```
