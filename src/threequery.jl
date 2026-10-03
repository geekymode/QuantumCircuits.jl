# ---------------------------------------------------------------------------
# Three-query unitary synthesis (Nehoran & Yuen, arXiv:2609.40351)
#
# For a real, symmetric, traceless unitary S ∈ O(N) — every unitary reduces to
# one, see `symmetric_embedding` —
#
#     S  ≈  i Ê† Q̂ F̂ Q̂ F̂ Q̂ Ê                                   (Theorem 4.1)
#
# on N registers of K grid points each:
#
#   Ê    encode |x⟩ in the single-excitation sector of N discretised harmonic
#        oscillators — a "one-hot" code with Gaussians h0, h1 for 0 and 1
#   Q̂    the chirp |v⟩ ↦ exp(i/2 · α(v)ᵀ S α(v)) |v⟩, α the grid positions
#   F̂    the centred discrete Fourier transform on every register
#
# The three queries evaluate the quadratic form vᵀSv at u, at v, and (through
# F̂ Q̂ F̂) at their midpoint — the polynomial identity of a Reed–Muller decoder,
# which recovers uᵀSv and puts S's entries into the amplitudes.  The error
# falls like exp(-cK), so K = O(log(N/ε)) grid points suffice.
#
# This is a reference implementation: it simulates the K^N-dimensional
# register space directly, so it is exact up to floating point and feasible
# only for small N.  Its job is to check the identity, not to be fast.
# ---------------------------------------------------------------------------

# Grid positions α(a) = (a - (K-1)/2)·√(2π/K), a ∈ 0:K-1 — spacing √(2π/K),
# symmetric about 0, self-dual under the centred DFT.
_qq_grid(K::Integer) = [(a - (K - 1) / 2) * sqrt(2π / K) for a in 0:K-1]

function _qq_check_size(N::Integer, K::Integer, maxdim::Integer)
    K >= 2 || throw(ArgumentError("need at least 2 grid points, got K = $K"))
    big(K)^N <= maxdim || throw(ArgumentError(
        "$N registers of $K points is a $(big(K)^N)-dimensional space, above maxdim = $maxdim; " *
        "lower K or raise maxdim"))
    nothing
end

"""
    centered_dft(K) -> Matrix{ComplexF64}

The discrete Fourier transform on the grid `α(a) = (a - (K-1)/2)·√(2π/K)`:
`F̂[w, u] = exp(-i α(u) α(w)) / √K` — Nehoran–Yuen Eq. 46.  Unitary and
symmetric.  It is the ordinary inverse QFT on `ℤ_K` between two copies of
the diagonal `D₁ = diag(exp(2πi/K · (K-1)/2 · (u - (K-1)/4)))` (Eq. 47) —
exactly, with no leftover phase, for every `K` — which is what lets it be
built from a standard QFT circuit.
"""
function centered_dft(K::Integer)
    K >= 2 || throw(ArgumentError("need at least 2 grid points, got K = $K"))
    g = _qq_grid(K)
    [cis(-g[w] * g[u]) for w in 1:K, u in 1:K] ./ sqrt(K)
end

# h0 and h1: the discretised ground and first excited oscillator states.
function _qq_oscillator(K::Integer)
    g = _qq_grid(K)
    h0 = exp.(-g .^ 2 ./ 2)
    h1 = g .* h0
    h0 ./ norm(h0), h1 ./ norm(h1)
end

# Register r (1 = most significant, as for qubits) is array dimension N - r + 1
# of a K×…×K array, so that `vec` gives the kron ordering.
_qq_shape(N, K, r) = ntuple(d -> d == N - r + 1 ? K : 1, N)

# The encoded basis state |ψ̂_x⟩ (x = 1…N) as a K×…×K array.
function _qq_encoded(N::Integer, K::Integer, x::Integer, h0, h1)
    A = ones(Float64, ntuple(_ -> K, N))
    for r in 1:N
        A .*= reshape(r == x ? h1 : h0, _qq_shape(N, K, r))
    end
    A
end

"""
    qho_encoding(N, K) -> Matrix{Float64}

The encoding isometry `Ê` of Nehoran–Yuen Eq. 44, as a `K^N × N` matrix.
Column `x` is `h₀ ⊗ … ⊗ h₁ ⊗ … ⊗ h₀`, with the first excited state `h₁` on
register `x`: a one-hot code in which `|0⟩` and `|1⟩` are replaced by the
discretised ground and first excited states of a harmonic oscillator,
`h_b(a) ∝ α(a)^b exp(-α(a)²/2)`.  The columns are exactly orthonormal.
"""
function qho_encoding(N::Integer, K::Integer; maxdim::Integer=1 << 22)
    _qq_check_size(N, K, maxdim)
    h0, h1 = _qq_oscillator(K)
    hcat((vec(_qq_encoded(N, K, x, h0, h1)) for x in 1:N)...)
end

"""
    chirp_phases(S, K) -> Vector{ComplexF64}

The diagonal of the quadratic-phase oracle `Q̂_S` (Nehoran–Yuen Eq. 45):
`exp(i/2 · α(v)ᵀ S α(v))` for every grid point `v ∈ [K]^N`, in the same
ordering as [`qho_encoding`](@ref) (register 1 most significant).
"""
function chirp_phases(S::AbstractMatrix, K::Integer; maxdim::Integer=1 << 22)
    N = size(S, 1)
    _qq_check_size(N, K, maxdim)
    g = _qq_grid(K)
    G = [reshape(g, _qq_shape(N, K, r)) for r in 1:N]
    φ = zeros(Float64, ntuple(_ -> K, N))
    Sr = real(S)
    for a in 1:N, b in 1:N
        iszero(Sr[a, b]) || (φ .+= Sr[a, b] .* G[a] .* G[b])
    end
    vec(cis.(φ ./ 2))
end

# Apply the same K×K matrix to every register of a K^N vector, in place.
function _qq_each_register!(ψ::Vector{ComplexF64}, F::Matrix{ComplexF64}, N::Integer, K::Integer)
    tmp = Matrix{ComplexF64}(undef, 0, 0)
    for d in 1:N
        L, R = K^(d - 1), K^(N - d)
        A = reshape(ψ, L, K, R)
        size(tmp) == (L, K) || (tmp = Matrix{ComplexF64}(undef, L, K))
        for r in 1:R
            @views mul!(tmp, A[:, :, r], transpose(F))
            @views A[:, :, r] .= tmp
        end
    end
    ψ
end

function _qq_check_target(S::AbstractMatrix, atol::Real)
    N = size(S, 1)
    size(S) == (N, N) || throw(ArgumentError("expected a square matrix, got $(size(S))"))
    all(x -> abs(imag(x)) <= atol, S) || throw(ArgumentError("S must be real"))
    isapprox(S, transpose(S); atol=atol) || throw(ArgumentError("S must be symmetric"))
    isapprox(S * S, I; atol=atol * N) || throw(ArgumentError("S must be unitary (S² = I)"))
    abs(tr(S)) <= atol * N || throw(ArgumentError(
        "S must be traceless (trace $(tr(S))); the identity needs it. " *
        "symmetric_embedding(U) gives a traceless S for any U"))
    nothing
end

"""
    three_query(S, K; check=true, maxdim=2^22) -> Matrix{ComplexF64}

Run Nehoran and Yuen's discretised three-query algorithm and return the
`N × N` matrix it implements,

    i Ê† Q̂_S F̂ Q̂_S F̂ Q̂_S Ê,

with `K` grid points per register — see [`qho_encoding`](@ref),
[`chirp_phases`](@ref) and [`centered_dft`](@ref).  For a real, symmetric,
traceless unitary `S` this approaches `S` exponentially fast in `K`
(their Theorem 4.1); [`three_query_error`](@ref) measures how fast.

`S` is checked for those properties unless `check=false`.  Turning the check
off is how you can watch the identity fail: drop `Tr S = 0` and the error
stays order one however fine the grid.

The `K^N`-dimensional register space is simulated directly, so this is for
small `N`; `maxdim` guards against accidental blow-up.
"""
function three_query(S::AbstractMatrix, K::Integer; check::Bool=true, atol::Real=1e-10,
                     maxdim::Integer=1 << 22)
    check && _qq_check_target(S, atol)
    N = size(S, 1)
    _qq_check_size(N, K, maxdim)
    h0, h1 = _qq_oscillator(K)
    q = chirp_phases(S, K; maxdim=maxdim)
    F = centered_dft(K)
    E = [vec(_qq_encoded(N, K, x, h0, h1)) for x in 1:N]
    M = Matrix{ComplexF64}(undef, N, N)
    for x in 1:N
        ψ = q .* E[x]                         # query
        _qq_each_register!(ψ, F, N, K)        # Fourier
        ψ .*= q                               # query
        _qq_each_register!(ψ, F, N, K)        # Fourier
        ψ .*= q                               # query
        for y in 1:N
            M[y, x] = im * dot(E[y], ψ)       # decode
        end
    end
    M
end

"""
    three_query_error(S, K; kwargs...) -> Float64

`‖three_query(S, K) - S‖` in operator norm.  Nehoran and Yuen bound it by
`O(N^{3/2} K^{3/4} e^{-πK/8})`; numerically it falls faster, at about
`e^{-πK/4}`.
"""
three_query_error(S::AbstractMatrix, K::Integer; kwargs...) =
    opnorm(three_query(S, K; kwargs...) - S)

"""
    three_query_synthesis(U, K; maxdim=2^22) -> Matrix{ComplexF64}

The whole pipeline for an arbitrary unitary `U`: embed it with
[`symmetric_embedding`](@ref), run [`three_query`](@ref) on the result, and
read `U` back out through the two Lemma 2.1 ancillas — prepare `|1⟩ ⊗ |−Y⟩`,
project onto `⟨0| ⊗ ⟨−Y|`.  The error is at most that of the embedded `S`.

The embedding quadruples the dimension, so even a one-qubit `U` means eight
registers: `K^8` grows fast, and only coarse grids fit under `maxdim`.
"""
function three_query_synthesis(U::AbstractMatrix, K::Integer; maxdim::Integer=1 << 22)
    n = size(U, 1)
    size(U) == (n, n) || throw(ArgumentError("expected a square matrix, got $(size(U))"))
    M = three_query(symmetric_embedding(U), K; maxdim=maxdim)
    # input |1⟩|−Y⟩|ψ⟩ selects columns 2n+1:4n as (cols₀ - i·cols₁)/√2
    C = (M[:, 2n+1:3n] .- im .* M[:, 3n+1:4n]) ./ sqrt(2)
    # output ⟨0|⟨−Y| keeps rows 1:2n as (rows₀ + i·rows₁)/√2
    (C[1:n, :] .+ im .* C[n+1:2n, :]) ./ sqrt(2)
end
