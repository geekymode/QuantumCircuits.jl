# ---------------------------------------------------------------------------
# Linear-algebra toolkit
#
# The numerical machinery that circuit synthesis rests on: transforms between
# the computational basis and the two other bases that matter (Walsh/parity and
# Pauli), plus the standard unitary-matrix predicates and state measures.
# ---------------------------------------------------------------------------

"""
    is_unitary(U; atol=1e-10) -> Bool

Whether `U'U` is the identity to within `atol`.
"""
is_unitary(U::AbstractMatrix; atol::Real=1e-10) = isapprox(U' * U, LinearAlgebra.I, atol=atol)

"""
    gate_fidelity(U, V) -> Float64

`|tr(U'V)| / d`, which is `1` exactly when `U` and `V` agree up to a global
phase and `0` when they are orthogonal.  The right way to compare two
synthesised circuits when you do not care about phase.
"""
gate_fidelity(U::AbstractMatrix, V::AbstractMatrix) = abs(LinearAlgebra.tr(U' * V)) / size(U, 1)

"""
    global_phase_between(U, V) -> Float64

The angle `ϕ` minimising `‖U - exp(im*ϕ) V‖`, i.e. the phase `V` must be
multiplied by to line it up with `U`.  Meaningless if the two are not
proportional — check [`gate_fidelity`](@ref) first.
"""
global_phase_between(U::AbstractMatrix, V::AbstractMatrix) = angle(LinearAlgebra.tr(V' * U))

"""
    embed(U, qubits, n) -> Matrix{ComplexF64}

Lift a `k`-qubit matrix acting on `qubits` to the full `2ⁿ × 2ⁿ` operator on an
`n`-qubit register.  `qubits[1]` is the most significant index bit of `U`.
"""
function embed(U::AbstractMatrix, qubits::AbstractVector{<:Integer}, n::Integer)
    N = 1 << n
    qs = collect(Int, qubits)
    Uc = ComplexF64.(U)
    M = Matrix{ComplexF64}(undef, N, N)
    ψ = Vector{ComplexF64}(undef, N)
    for j in 1:N
        fill!(ψ, 0)
        ψ[j] = 1
        apply!(ψ, Uc, qs, Int(n))
        @views M[:, j] .= ψ
    end
    M
end

"""
    kron_n(U, k) -> Matrix

`k`-fold Kronecker power of `U`.
"""
function kron_n(U::AbstractMatrix, k::Integer)
    k >= 1 || throw(ArgumentError("k must be ≥ 1"))
    M = ComplexF64.(U)
    for _ in 2:k
        M = kron(M, ComplexF64.(U))
    end
    M
end

# --- the Walsh / parity basis --------------------------------------------

"""
    fwht!(a) -> a

In-place fast Walsh–Hadamard transform (unnormalised, natural order):

    a[s+1] ← Σ_j (-1)^(s · j) a[j+1]

`O(N log N)` instead of the `O(N²)` matrix product.  This is the change of
basis between "a value per computational basis state" and "a coefficient per
parity function", and it is the transform sitting inside every multiplexed
rotation ([`multiplex_angles`](@ref)) and every phase polynomial
([`phase_polynomial`](@ref)).

Applying it twice multiplies by `N`.
"""
function fwht!(a::AbstractVector)
    N = length(a)
    ispow2(N) || throw(ArgumentError("length must be a power of two, got $N"))
    h = 1
    while h < N
        for i in 1:(2h):N
            @inbounds for j in i:(i+h-1)
                x, y = a[j], a[j+h]
                a[j] = x + y
                a[j+h] = x - y
            end
        end
        h *= 2
    end
    a
end

"""
    fwht(a) -> Vector

Out-of-place [`fwht!`](@ref).
"""
fwht(a::AbstractVector) = fwht!(collect(float.(a)))

"""
    walsh_matrix(n) -> Matrix{Float64}

The `2ⁿ × 2ⁿ` matrix `W[s+1, j+1] = (-1)^(s · j)` that [`fwht!`](@ref) applies.
Built explicitly for inspection and tests only.
"""
walsh_matrix(n::Integer) = [(-1.0)^parity(s, j) for s in 0:(1 << n)-1, j in 0:(1 << n)-1]

# --- the Pauli basis -------------------------------------------------------

const _PAULI1 = Dict('I' => ComplexF64[1 0; 0 1], 'X' => ComplexF64[0 1; 1 0],
                     'Y' => ComplexF64[0 -im; im 0], 'Z' => ComplexF64[1 0; 0 -1])

"""
    pauli(s::AbstractString) -> Matrix{ComplexF64}

The Pauli string operator, e.g. `pauli("XIZ") == kron(X, I, Z)`.  Character `i`
of `s` is the Pauli on qubit `i`, so the leftmost character is the most
significant factor — the same reading order as everything else here.
"""
function pauli(s::AbstractString)
    isempty(s) && throw(ArgumentError("empty Pauli string"))
    M = ComplexF64[1;;]
    for ch in s
        haskey(_PAULI1, ch) || throw(ArgumentError("bad Pauli character '$ch'; use I, X, Y or Z"))
        M = kron(M, _PAULI1[ch])
    end
    M
end

"""
    pauli_strings(n) -> Vector{String}

All `4ⁿ` Pauli strings on `n` qubits, in odometer order (`"II"`, `"IX"`, …).
"""
function pauli_strings(n::Integer)
    chars = ('I', 'X', 'Y', 'Z')
    [String([chars[(k ÷ 4^(n - i)) % 4 + 1] for i in 1:n]) for k in 0:4^n-1]
end

"""
    pauli_decompose(M; atol=1e-12) -> Vector{Pair{String,ComplexF64}}

Expand `M` in the Pauli basis: `M = Σ_P c_P P` with
`c_P = tr(P M) / 2ⁿ`, dropping terms below `atol`.

The Pauli basis is orthogonal under the Hilbert–Schmidt inner product, so this
is an exact change of basis, not a fit.  It is how a Hamiltonian becomes a list
of terms you can Trotterise with [`pauli_rotation!`](@ref).
"""
function pauli_decompose(M::AbstractMatrix; atol::Real=1e-12)
    N = size(M, 1)
    N == size(M, 2) && ispow2(N) || throw(ArgumentError("need a 2ⁿ × 2ⁿ matrix"))
    n = trailing_zeros(N)
    out = Pair{String,ComplexF64}[]
    for s in pauli_strings(n)
        c = LinearAlgebra.tr(pauli(s) * M) / N
        abs(c) > atol && push!(out, s => c)
    end
    out
end

"""
    pauli_recompose(terms) -> Matrix{ComplexF64}

Rebuild an operator from the output of [`pauli_decompose`](@ref).
"""
function pauli_recompose(terms::AbstractVector{<:Pair{<:AbstractString,<:Number}})
    isempty(terms) && throw(ArgumentError("no terms"))
    n = length(first(first(terms)))
    M = zeros(ComplexF64, 1 << n, 1 << n)
    for (s, c) in terms
        M .+= c .* pauli(s)
    end
    M
end

# --- state measures --------------------------------------------------------

"""
    schmidt_values(ψ, k) -> Vector{Float64}

Singular values of `ψ` reshaped across the cut after qubit `k`, i.e. the
Schmidt coefficients of the bipartition `(1:k, k+1:n)`.  They square to the
eigenvalues of the reduced density matrix.
"""
function schmidt_values(ψ::AbstractVector, k::Integer)
    N = length(ψ)
    ispow2(N) || throw(ArgumentError("state length must be a power of two"))
    n = trailing_zeros(N)
    1 <= k < n || throw(ArgumentError("cut must satisfy 1 ≤ k < n = $n"))
    # qubit 1 is the most significant bit, so the leading k bits index the left
    # factor; column-major reshape puts the trailing bits down the columns.
    LinearAlgebra.svdvals(reshape(collect(ComplexF64, ψ), 1 << (n - k), 1 << k))
end

"""
    entanglement_entropy(ψ, k) -> Float64

Von Neumann entropy (in bits) of the reduced state on qubits `1:k`.  `0` for a
product state across the cut, `k` for a maximally entangled one.
"""
function entanglement_entropy(ψ::AbstractVector, k::Integer)
    s = schmidt_values(ψ, k)
    p = s .^ 2
    p ./= sum(p)
    -sum(x > 0 ? x * log2(x) : 0.0 for x in p)
end
