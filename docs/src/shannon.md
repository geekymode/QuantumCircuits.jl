# The hardest case: an arbitrary unitary

```@meta
CurrentModule = QuantumCircuits
```

Every other decomposition in this package exploits structure. A multiplexor is
block diagonal. A diagonal unitary is a phase function. A state-preparation
target is a vector, not a matrix. Take all of that away — hand the compiler a
`2ⁿ × 2ⁿ` unitary with no structure at all — and you are left with the hardest
of the standard problems.

The **quantum Shannon decomposition** is the standard answer, and its entire
skeleton is Gray-code multiplexors.

```@setup shannon
using QuantumCircuits, CairoMakie, LinearAlgebra, Random
CairoMakie.activate!(type = "png", px_per_unit = 2)
Random.seed!(20260819)
```

## Why it is hard

The naive route is the one on the [Applications](@ref "Applications of Gray
coding") page: factor `U` into two-level (Givens) rotations and route each with
CNOTs. It is exact, and it is bad — there are `O(4ⁿ)` factors and each one
drags a multi-controlled gate behind it.

```@example shannon
U = rand_unitary(8)
tl = synthesize_unitary(U; lower = true)
matrix(tl) ≈ U, length(tl), count_cnots(tl)
```

98 CNOTs for three qubits, and it degrades fast. The parameter count says we
should do much better: a unitary in `U(2ⁿ)` has `4ⁿ` real parameters, and each
CNOT plus its surrounding one-qubit gates buys a constant number of them, so
`Θ(4ⁿ)` CNOTs is the floor. The two-level route pays an extra factor of `n`.

## One level of the recursion

Split the register at the top wire. Two classical factorisations then peel that
wire off entirely.

**Cosine–sine.** Any unitary of even dimension factors as

```math
U \;=\; (L_1 \oplus L_2)\;\begin{bmatrix} C & -S \\ S & C\end{bmatrix}\;(R_1 \oplus R_2)^\dagger
```

with `C` and `S` diagonal and `C² + S² = I`. Written on qubits, that middle
factor **is** a uniformly controlled `RY` on the top wire — a Gray-code
multiplexor, `2ⁿ⁻¹` CNOTs — and the outer factors are block-diagonal, i.e.
multiplexed unitaries.

```@example shannon
F = cosine_sine(U)
Matrix(F) ≈ U, round.(F.c, digits = 3), round.(F.s, digits = 3)
```

The claim about the middle factor is worth checking rather than believing:

```@example shannon
C, S = Diagonal(F.c), Diagonal(F.s)
matrix(multiplexed_ry(csd_angles(F), 2:3, 1; n = 3)) ≈ [C -S; S C]
```

```@example shannon
csdfigure(rand_unitary(16))
```

Left to right: a structureless `U`, then the two block-diagonal factors and the
four-diagonal cosine–sine pattern it resolves into. The dashed lines mark the
split at the top wire.

**Demultiplexing.** Each block-diagonal factor is then split by
[`demultiplex`](@ref), `L₁ ⊕ L₂ = (I ⊗ V)(D ⊕ D†)(I ⊗ W)`, and `D ⊕ D†` is a
uniformly controlled `RZ` on the same top wire — another Gray-code multiplexor.

So one level emits **three multiplexors** and four `(n-1)`-qubit unitaries:

```@example shannon
qsdfigure(3)
```

## The cost

```math
\mathrm{CNOTs}(n) = 4\,\mathrm{CNOTs}(n-1) + 3\cdot 2^{n-1},
\qquad \mathrm{CNOTs}(1) = 0
```

which solves to `(3/4)·4ⁿ - (3/2)·2ⁿ`.

```@example shannon
[(n, count_cnots(qsd(rand_unitary(1 << n))), qsd_cnot_count(n)) for n in 1:4]
```

Against the two-level route at `n = 3`: **36 CNOTs instead of 98**, and the gap
widens with `n` because the two-level method carries an extra factor of `n`
while this one does not.

```@example shannon
circuitfigure(qsd(rand_unitary(4)); title = "QSD of a Haar-random two-qubit unitary")
```

Every CNOT you see came out of a multiplexor, which is to say out of a Gray-code
walk. Strip the Gray coding out and each multiplexor reverts to one
multi-controlled rotation per branch — `O(k·2ᵏ)` instead of `2ᵏ` — and the
whole recursion picks up the factor of `n` again, landing back where the
two-level method was.

## Honest limits

* This implementation runs a constant factor above the literature's
  `(9/16)·4ⁿ`, which is reached by handling the two-qubit blocks with a KAK
  (Cartan) decomposition — three CNOTs each — rather than recursing into them.
  KAK is not implemented here; it is the single largest remaining win.
  Further optimisations reach `(23/48)·4ⁿ`.
* The decomposition is **structure-blind**. It spends the same 36 CNOTs on the
  identity as on a Haar-random unitary:

```@example shannon
count_cnots(qsd(Matrix{ComplexF64}(I, 8, 8))), count_cnots(qsd(rand_unitary(8)))
```

  When you know the structure, use the routine that knows it too —
  [`diagonal`](@ref) for a diagonal, [`prepare_state`](@ref) for a state,
  [`multiplexed_1q`](@ref) for a multiplexor.

## A numerical trap worth knowing about

The sines in the cosine–sine decomposition are *not* computed as
`sqrt(1 - c²)`. Near `c = 1` that subtraction loses every significant digit:
the true sine is `0`, the computed one is `√eps ≈ 1.5e-8`, and the
reconstruction formulas then divide by it. A block-diagonal input — where every
sine is exactly zero — comes back with non-unitary factors and an error of
order 1.

[`cosine_sine`](@ref) takes the sines as column norms of `U₂₁R₁` instead, which
is accurate to machine epsilon and vanishes exactly when it should.

```@example shannon
F = cosine_sine(cat(rand_unitary(4), rand_unitary(4); dims = (1, 2)))
maximum(F.s), all(is_unitary(M) for M in (F.L1, F.L2, F.R1, F.R2))
```

Columns whose sine is genuinely zero are a real gauge freedom — nothing outside
the top-left block constrains them — so those columns are completed to an
orthonormal basis and their partners read off where the cosine is near 1.

## Reference

```@docs
qsd
qsd!
qsd_cnot_count
cosine_sine
CSD
csd_angles
rand_unitary
csdfigure
qsdfigure
```
