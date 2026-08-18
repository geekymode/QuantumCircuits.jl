# ---------------------------------------------------------------------------
# Matrix decompositions → circuits
#
# Three classical factorisations, each turning a matrix identity into gates:
#
#   ZYZ (Euler)      any U(2) = e^{iα} RZ(β) RY(γ) RZ(δ)
#   two-level/Givens any U(N) = a product of rotations in 2-dimensional
#                    coordinate planes  (QR by plane rotations)
#   demultiplexing   U₁ ⊕ U₂ = (I⊗V)(D ⊕ D†)(I⊗W), from one eigendecomposition
#
# The Gray code enters at the point where each factor becomes hardware gates:
# a two-level rotation acts on basis states |a⟩ and |b⟩, and the cheapest way
# to reach it is to route `a` and `b` to be Gray-adjacent — one bit apart —
# with CNOTs, act once, and undo the routing.
# ---------------------------------------------------------------------------

"""
    zyz(U) -> (α, β, γ, δ)

Euler decomposition of a `2 × 2` unitary:

    U = exp(im*α) * RZ(β) * RY(γ) * RZ(δ)

Handles the degenerate cases `γ ≈ 0` and `γ ≈ π` (where `β` and `δ` are
individually undetermined and only their sum or difference is fixed) by
pinning `δ = 0`.
"""
function zyz(U::AbstractMatrix)
    size(U) == (2, 2) || throw(ArgumentError("zyz expects a 2×2 matrix"))
    Uc = ComplexF64.(U)
    is_unitary(Uc; atol=1e-8) || throw(ArgumentError("zyz expects a unitary matrix"))
    α = angle(LinearAlgebra.det(Uc)) / 2
    V = Uc .* cis(-α)                       # V ∈ SU(2)
    c, s = abs(V[1, 1]), abs(V[2, 1])
    γ = 2 * atan(s, c)
    if c < 1e-9                             # γ ≈ π: only β - δ is determined
        β = 2 * angle(V[2, 1])
        δ = 0.0
    elseif s < 1e-9                         # γ ≈ 0: only β + δ is determined
        β = 2 * angle(V[2, 2])
        δ = 0.0
    else
        p = angle(V[2, 2])                  # (β + δ)/2
        m = angle(V[2, 1])                  # (β - δ)/2
        β = p + m
        δ = p - m
    end
    (α, β, γ, δ)
end

"""
    decompose_1q!(c, U, q) -> c

Append an arbitrary one-qubit unitary to `c` on wire `q` as `RZ·RY·RZ` plus a
global phase, via [`zyz`](@ref).
"""
function decompose_1q!(c::Circuit, U::AbstractMatrix, q::Integer)
    α, β, γ, δ = zyz(U)
    push!(c, RZ(δ), q)
    push!(c, RY(γ), q)
    push!(c, RZ(β), q)
    c.global_phase += α
    c
end

"""
    decompose_1q(U, q; n=q) -> Circuit

Standalone circuit for [`decompose_1q!`](@ref).
"""
decompose_1q(U::AbstractMatrix, q::Integer=1; n::Integer=q) = decompose_1q!(Circuit(n), U, q)

# --- two-level (Givens) decomposition --------------------------------------

"""
    TwoLevel(a, b, V)

A unitary that acts as the `2 × 2` matrix `V` on the span of basis states
`|a⟩, |b⟩` (0-based, `a < b`) and as the identity everywhere else.
"""
struct TwoLevel
    a::Int
    b::Int
    V::Matrix{ComplexF64}
end

Base.show(io::IO, t::TwoLevel) = print(io, "TwoLevel(|", t.a, "⟩,|", t.b, "⟩)")

"""
    matrix(t::TwoLevel, N) -> Matrix{ComplexF64}

Expand a two-level unitary to its full `N × N` matrix.
"""
function matrix(t::TwoLevel, N::Integer)
    M = Matrix{ComplexF64}(LinearAlgebra.I, N, N)
    M[t.a+1, t.a+1] = t.V[1, 1]; M[t.a+1, t.b+1] = t.V[1, 2]
    M[t.b+1, t.a+1] = t.V[2, 1]; M[t.b+1, t.b+1] = t.V[2, 2]
    M
end

"""
    two_level_decompose(U; atol=1e-12) -> Vector{TwoLevel}

Factor any `N × N` unitary into two-level unitaries, `U = T₁ T₂ ⋯ T_m` (in the
listed order).

This is QR by plane rotations: for each column, Givens rotations zero the
entries below the diagonal. A unitary upper-triangular matrix is diagonal, so
what survives is a diagonal of phases, emitted as a final run of two-level
phase gates. At most `N(N-1)/2 + N` factors.

The classical route to synthesising an *arbitrary* unitary — see
[`synthesize_unitary`](@ref), which turns each factor into gates using
Gray-code routing.
"""
function two_level_decompose(U::AbstractMatrix; atol::Real=1e-12)
    N = size(U, 1)
    size(U, 2) == N || throw(ArgumentError("matrix must be square"))
    is_unitary(U; atol=1e-8) || throw(ArgumentError("two_level_decompose expects a unitary"))
    A = ComplexF64.(U)
    undo = TwoLevel[]                       # the G's, as their inverses
    for col in 1:N-1, row in N:-1:col+1
        x, y = A[col, col], A[row, col]
        abs(y) <= atol && continue
        r = hypot(abs(x), abs(y))
        # g * (x, y)ᵀ = (r, 0)ᵀ, and g is unitary
        g = ComplexF64[conj(x)/r conj(y)/r; -y/r x/r]
        @views begin
            rows = A[[col, row], :]
            A[[col, row], :] = g * rows
        end
        push!(undo, TwoLevel(col - 1, row - 1, ComplexF64[g[1,1] g[1,2]; g[2,1] g[2,2]]'))
    end
    # A is now diag(d₁ … d_N) with |dᵢ| = 1, and Gₘ ⋯ G₁ U = A, so
    # U = G₁⁻¹ G₂⁻¹ ⋯ Gₘ⁻¹ A — the inverses in the order they were recorded.
    factors = undo
    for k in 1:N-1
        d = A[k, k]
        isapprox(d, 1; atol=atol) && continue
        push!(factors, TwoLevel(k - 1, k, ComplexF64[d 0; 0 1]))
    end
    dN = A[N, N]
    isapprox(dN, 1; atol=atol) || push!(factors, TwoLevel(N - 2, N - 1, ComplexF64[1 0; 0 dN]))
    factors
end

"""
    two_level!(c, t::TwoLevel, n) -> c

Append a two-level unitary on an `n`-qubit register, using Gray-code routing.

`|a⟩` and `|b⟩` generally differ in several bits. CNOTs fanning out from the
lowest differing bit `p` bring them to a Gray-adjacent pair — identical except
in bit `p` — at which point one multi-controlled `V` on wire `p` does the whole
job. The routing is then undone, so the net operator is the identity outside
`span(|a⟩, |b⟩)`.

`lower = true` expands that multi-controlled gate into CNOTs and one-qubit
gates with [`multicontrolled!`](@ref), giving a circuit made entirely of
elementary gates.
"""
function two_level!(c::Circuit, t::TwoLevel, n::Integer=c.nqubits; lower::Bool=false)
    a, b, V = t.a, t.b, t.V
    a == b && throw(ArgumentError("two-level gate needs distinct basis states"))
    d = a ⊻ b
    p = trailing_zeros(d)                   # bit position, LSB = 0
    qp = n - p                              # ... as a wire index
    if (a >> p) & 1 == 1                    # ensure |a⟩ is the bit-p = 0 partner
        a, b = b, a
        V = ComplexF64[V[2,2] V[2,1]; V[1,2] V[1,1]]
    end
    routing = Int[]
    for q in 0:n-1
        (q == p || ((d >> q) & 1) == 0) && continue
        push!(c, CNOT(), qp, n - q)         # copy bit p onto every other differing bit
        push!(routing, n - q)
    end
    ctrls = sort!([n - q for q in 0:n-1 if q != p])
    zeroctrl = [q for q in ctrls if ((a >> (n - q)) & 1) == 0]
    for q in zeroctrl                       # controls fire on |1⟩, so flip the |0⟩ ones
        push!(c, X(), q)
    end
    if lower
        multicontrolled!(c, V, ctrls, qp)      # ... down to CNOTs and one-qubit gates
    else
        push!(c, controlled(Gate(:V, V), length(ctrls)), ctrls..., qp)
    end
    for q in zeroctrl
        push!(c, X(), q)
    end
    for q in Iterators.reverse(routing)
        push!(c, CNOT(), qp, q)
    end
    c
end

"""
    synthesize_unitary(U) -> Circuit

Synthesise an arbitrary `2ⁿ × 2ⁿ` unitary: [`two_level_decompose`](@ref) to get
plane rotations, then [`two_level!`](@ref) to turn each one into gates via
Gray-code routing.

Exact but not cheap — `O(4ⁿ)` factors. By default the multi-controlled `V`
gates are left as single instructions; `lower = true` expands them into CNOTs
and one-qubit gates via [`multicontrolled!`](@ref). Use it as a reference
implementation and for small `n`; [`prepare_state`](@ref) and
[`diagonal`](@ref) are the efficient paths for the structured cases.
"""
function synthesize_unitary(U::AbstractMatrix; lower::Bool=false)
    N = size(U, 1)
    ispow2(N) || throw(ArgumentError("need a 2ⁿ × 2ⁿ matrix"))
    n = trailing_zeros(N)
    c = Circuit(n)
    # circuits apply left to right, matrices right to left
    for t in Iterators.reverse(two_level_decompose(U))
        two_level!(c, t, n; lower=lower)
    end
    c
end

# --- demultiplexing --------------------------------------------------------

"""
    demultiplex(U₁, U₂) -> (V, θ, W)

Split a two-branch multiplexed unitary into two unmultiplexed halves and one
uniformly controlled `RZ`:

    U₁ ⊕ U₂ = (I ⊗ V) · (D ⊕ D†) · (I ⊗ W)

with `D = diag(exp(im*θ/(-2)))` realised by a Gray-code `RZ` multiplexor
targeting the *control* wire. The construction is one eigendecomposition:
`U₁U₂† = V D² V†` and `W = D V† U₂`.

This is the recursive step of the quantum Shannon decomposition, and the
reason [`multiplexed_rotation!`](@ref) is the primitive everything else is
written in terms of.
"""
function demultiplex(U1::AbstractMatrix, U2::AbstractMatrix)
    size(U1) == size(U2) || throw(ArgumentError("blocks must have the same size"))
    A = ComplexF64.(U1) * ComplexF64.(U2)'
    # A is unitary, hence normal: its Schur form is diagonal and Z is unitary,
    # which `eigen` cannot promise for degenerate spectra.
    F = LinearAlgebra.schur(A)
    V = Matrix{ComplexF64}(F.Z)
    d2 = LinearAlgebra.diag(F.T)
    D = LinearAlgebra.Diagonal(cis.(angle.(d2) ./ 2))
    W = D * V' * ComplexF64.(U2)
    θ = -2 .* angle.(LinearAlgebra.diag(D))
    (V, θ, Matrix{ComplexF64}(W))
end

# --- uniformly controlled one-qubit gates ----------------------------------

"""
    multiplexed_1q!(c, Us, controls, target) -> c

Uniformly controlled *arbitrary* one-qubit gate: apply `Us[j+1]` to `target`
when the control register reads `j`.

Each branch is Euler-decomposed ([`zyz`](@ref)) and the three angle families
are realised as three Gray-code multiplexors, `RZ·RY·RZ`; the leftover per-branch
phases form a [`diagonal`](@ref) unitary on the controls. Cost is
`4·2ᵏ - 2` CNOTs, versus `O(k·2ᵏ)` for one multi-controlled gate per branch.
"""
function multiplexed_1q!(c::Circuit, Us::AbstractVector{<:AbstractMatrix},
                         controls::AbstractVector, target::Integer)
    ctrl = collect(Int, controls)
    k = length(ctrl)
    length(Us) == 1 << k || throw(ArgumentError("expected $(1 << k) matrices for $k controls, got $(length(Us))"))
    e = [zyz(U) for U in Us]
    multiplexed_rotation!(c, :RZ, [x[4] for x in e], ctrl, target)   # δ
    multiplexed_rotation!(c, :RY, [x[3] for x in e], ctrl, target)   # γ
    multiplexed_rotation!(c, :RZ, [x[2] for x in e], ctrl, target)   # β
    φ = [x[1] for x in e]
    c.global_phase += k == 0 ? φ[1] : diagonal!(c, φ, ctrl)
    c
end

"""
    multiplexed_1q(Us, controls, target; n) -> Circuit

Standalone circuit for [`multiplexed_1q!`](@ref).
"""
function multiplexed_1q(Us::AbstractVector{<:AbstractMatrix}, controls::AbstractVector, target::Integer;
                        n::Integer=max(Int(target), maximum(Int, controls; init=0)))
    multiplexed_1q!(Circuit(n), Us, controls, target)
end
