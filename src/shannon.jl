# ---------------------------------------------------------------------------
# Cosine–sine decomposition and the quantum Shannon decomposition
#
# The hardest of the standard decompositions: an arbitrary n-qubit unitary,
# with no structure to exploit, into CNOTs and one-qubit gates.
#
# The skeleton is three Gray-code multiplexors per level of recursion:
#
#   U  =  (L₁ ⊕ L₂) · CS · (R₁ ⊕ R₂)†                       cosine–sine
#
#   CS = [ C  -S ]   is a uniformly controlled RY on the top wire       ← Gray
#        [ S   C ]
#
#   L₁ ⊕ L₂ = (I ⊗ V)(D ⊕ D†)(I ⊗ W)                        demultiplexing
#   D ⊕ D†  is a uniformly controlled RZ on the top wire                ← Gray
#
# and V, W, and the two halves of the right factor are (n-1)-qubit unitaries,
# so the whole thing recurses.  Every two-qubit gate in the output comes from a
# multiplexor, i.e. from a Gray-code walk.
# ---------------------------------------------------------------------------

"""
    CSD

The result of a [`cosine_sine`](@ref) decomposition of a `2m × 2m` unitary:

    U = (L1 ⊕ L2) · [C  -S; S  C] · (R1 ⊕ R2)'

with `C = Diagonal(c)`, `S = Diagonal(s)` and `c.^2 + s.^2 == 1`.
"""
struct CSD
    L1::Matrix{ComplexF64}
    L2::Matrix{ComplexF64}
    c::Vector{Float64}
    s::Vector{Float64}
    R1::Matrix{ComplexF64}
    R2::Matrix{ComplexF64}
end

Base.show(io::IO, F::CSD) = print(io, "CSD(", size(F.L1, 1), "×", size(F.L1, 1), " blocks)")

# orthonormal completion: k more columns orthogonal to those of V
function _complete_basis(V::AbstractMatrix, k::Integer)
    m = size(V, 1)
    k == 0 && return zeros(ComplexF64, m, 0)
    size(V, 2) == 0 && return Matrix{ComplexF64}(LinearAlgebra.I, m, m)[:, 1:k]
    W = LinearAlgebra.nullspace(Matrix(V'))
    size(W, 2) >= k || error("cannot complete basis: degenerate subspace is too small")
    ComplexF64.(W[:, 1:k])
end

"""
    cosine_sine(U; atol=1e-9) -> CSD

Cosine–sine decomposition of a unitary of even dimension, split into equal
halves:

    U = (L1 ⊕ L2) · [C -S; S C] · (R1 ⊕ R2)'

The middle factor is the interesting one: written on qubits, `[C -S; S C]` is
exactly a uniformly controlled `RY` with the *top* wire as target and the rest
as controls — a Gray-code multiplexor, `2ⁿ⁻¹` CNOTs.

Implementation: one SVD of the top-left block gives `L1`, `c` and `R1`; the
`s`-equations then give the columns of `L2` and `R2` wherever `s` is safely
non-zero. Columns with `s ≈ 0` are a genuine gauge freedom — nothing outside the
top-left block constrains them — so those columns of `R2` are completed to an
orthonormal basis and the matching columns of `L2` read off from the bottom-right
block, where `c ≈ 1` makes the division well conditioned.

The sines are taken as column norms of `U21·R1` rather than from
`sqrt(1 - c²)`: near `c = 1` that subtraction loses every significant digit and
returns `√eps ≈ 1.5e-8` where the true value is zero, which the `s`-equations
would then divide into. `atol` is the threshold below which a sine counts as
zero. Generic unitaries are nowhere near it; exactly-degenerate ones
(block-diagonal inputs, the identity, permutations) land on it squarely and take
the completion branch.
"""
function cosine_sine(U::AbstractMatrix; atol::Real=1e-10)
    N = size(U, 1)
    N == size(U, 2) && iseven(N) || throw(ArgumentError("cosine_sine needs a square matrix of even size"))
    is_unitary(U; atol=1e-8) || throw(ArgumentError("cosine_sine expects a unitary"))
    m = N ÷ 2
    A = ComplexF64.(U)
    U11 = A[1:m, 1:m];      U12 = A[1:m, m+1:N]
    U21 = A[m+1:N, 1:m];    U22 = A[m+1:N, m+1:N]

    F = LinearAlgebra.svd(U11)                 # U11 = L1 * Diagonal(c) * R1'
    L1 = Matrix{ComplexF64}(F.U)
    R1 = Matrix{ComplexF64}(F.V)
    c = clamp.(F.S, 0.0, 1.0)

    P = L1' * U12                              # = -Diagonal(s) * R2'
    Q = U21 * R1                               # =  L2 * Diagonal(s)
    # Take the sines from the data, not from sqrt(1 - c^2): near c = 1 that
    # subtraction loses every significant digit and returns √eps ≈ 1.5e-8 where
    # the true value is 0, which then gets divided into.  Column norms of
    # U21·R1 are accurate to machine epsilon and vanish exactly when they should.
    s = [LinearAlgebra.norm(view(Q, :, j)) for j in 1:m]
    L2 = zeros(ComplexF64, m, m)
    R2 = zeros(ComplexF64, m, m)
    solid = [j for j in 1:m if s[j] > atol]
    gauge = [j for j in 1:m if s[j] <= atol]
    for j in solid
        R2[:, j] = -conj.(P[j, :]) ./ s[j]
        # for L2 prefer whichever of s, c is larger — one of them is ≥ 1/√2
        L2[:, j] = s[j] >= c[j] ? Q[:, j] ./ s[j] : (U22 * R2[:, j]) ./ c[j]
    end
    if !isempty(gauge)
        R2[:, gauge] = _complete_basis(R2[:, solid], length(gauge))
        for j in gauge
            L2[:, j] = (U22 * R2[:, j]) ./ c[j]
        end
    end
    CSD(L1, L2, c, s, R1, R2)
end

"""
    Matrix(F::CSD) -> Matrix{ComplexF64}

Reassemble the unitary from its cosine–sine factors.
"""
function Base.Matrix(F::CSD)
    m = length(F.c)
    C = LinearAlgebra.Diagonal(F.c)
    S = LinearAlgebra.Diagonal(F.s)
    mid = [C -S; S C]
    cat(F.L1, F.L2; dims=(1, 2)) * mid * cat(F.R1, F.R2; dims=(1, 2))'
end

"""
    csd_angles(F::CSD) -> Vector{Float64}

The `RY` angles of the middle factor: `θⱼ = 2·atan(sⱼ, cⱼ)`, so that
`[C -S; S C]` is the uniformly controlled rotation `RY(θⱼ)` on the top wire.
"""
csd_angles(F::CSD) = [2 * atan(F.s[j], F.c[j]) for j in eachindex(F.c)]

# --- quantum Shannon decomposition -----------------------------------------

# append A ⊕ B (multiplexed by qubits[1]) using one demultiplexing step
function _demux!(c::Circuit, A::AbstractMatrix, B::AbstractMatrix, qubits::Vector{Int})
    V, θ, W = demultiplex(A, B)
    rest = qubits[2:end]
    qsd!(c, W, rest)                                     # right factor first
    multiplexed_rotation!(c, :RZ, θ, rest, qubits[1])    # D ⊕ D†  — Gray code
    qsd!(c, V, rest)
    c
end

"""
    qsd!(c, U, qubits) -> c

Append an arbitrary unitary on `qubits`, by the quantum Shannon decomposition.

One level peels the top wire with a [`cosine_sine`](@ref) split and two
[`demultiplex`](@ref) steps, emitting three Gray-code multiplexors — one `RY`
from the cosine–sine middle, one `RZ` from each side — and recursing on four
`(n-1)`-qubit unitaries. The recursion bottoms out at [`decompose_1q!`](@ref).

    CNOTs(n) = 4·CNOTs(n-1) + 3·2^(n-1),  CNOTs(1) = 0
             = (3/4)·4ⁿ - (3/2)·2ⁿ

which is 36 CNOTs at `n = 3` and 168 at `n = 4`, against roughly `O(n·4ⁿ)` for
the two-level route of [`synthesize_unitary`](@ref) — 98 CNOTs at `n = 3` and
far worse beyond. The literature's `(9/16)·4ⁿ` needs the two-qubit blocks
handled by a KAK decomposition instead of recursing; that is not implemented
here, so this runs a constant factor above the best known.
"""
function qsd!(c::Circuit, U::AbstractMatrix, qubits::AbstractVector{<:Integer})
    qs = collect(Int, qubits)
    n = length(qs)
    N = size(U, 1)
    N == 1 << n || throw(ArgumentError("matrix of size $N does not match $n qubits"))
    n == 1 && return decompose_1q!(c, U, qs[1])

    F = cosine_sine(U)
    _demux!(c, Matrix(F.R1'), Matrix(F.R2'), qs)                  # (R1 ⊕ R2)†
    multiplexed_rotation!(c, :RY, csd_angles(F), qs[2:end], qs[1]) # CS — Gray code
    _demux!(c, F.L1, F.L2, qs)                                    # L1 ⊕ L2
    c
end

"""
    qsd(U) -> Circuit

Quantum Shannon decomposition of an arbitrary `2ⁿ × 2ⁿ` unitary into CNOTs and
one-qubit gates. See [`qsd!`](@ref).

```julia
U = rand_unitary(8)
c = qsd(U)
matrix(c) ≈ U          # exact, global phase included
count_cnots(c)         # 36
```
"""
function qsd(U::AbstractMatrix)
    N = size(U, 1)
    ispow2(N) || throw(ArgumentError("need a 2ⁿ × 2ⁿ matrix"))
    qsd!(Circuit(trailing_zeros(N)), U, 1:trailing_zeros(N))
end

"""
    qsd_cnot_count(n) -> Int

The CNOT count this implementation's recursion produces for `n` qubits,
`(3/4)·4ⁿ - (3/2)·2ⁿ`. Useful for checking a synthesised circuit.
"""
qsd_cnot_count(n::Integer) = n <= 0 ? 0 : (3 * (1 << (2n)) ÷ 4) - 3 * (1 << n) ÷ 2

"""
    rand_unitary(N; rng=Random.default_rng()) -> Matrix{ComplexF64}

A Haar-random `N × N` unitary, from the QR of a complex Gaussian matrix with
the phase ambiguity fixed. Convenient for exercising the decompositions.
"""
function rand_unitary(N::Integer; rng=Random.default_rng())
    F = LinearAlgebra.qr(randn(rng, ComplexF64, N, N))
    Q = Matrix{ComplexF64}(F.Q)
    d = LinearAlgebra.diag(F.R)
    Q * LinearAlgebra.Diagonal(d ./ abs.(d))
end
