# ---------------------------------------------------------------------------
# Pauli operators and stabilizer codes
#
# An n-qubit Pauli operator is i^phase · P₁ ⊗ … ⊗ Pₙ, stored as two bit
# vectors: qubit j carries X if x[j], Z if z[j], Y if both.  Products need
# only XORs and a phase count, and two Paulis commute exactly when their
# symplectic product  x·z′ + z·x′  is even.
#
# A stabilizer code is the joint +1 eigenspace of r commuting, independent
# Paulis on n qubits: it encodes k = n - r logical qubits.  Errors are found by
# which stabilizers they anticommute with (the syndrome); logical operators
# are the Paulis that commute with every stabilizer without being one.
# ---------------------------------------------------------------------------

"""
    PauliOp(s)
    PauliOp(x, z, phase=0)

An `n`-qubit Pauli operator `i^phase · P₁ ⊗ … ⊗ Pₙ`.  From a string,
`PauliOp("XIZY")`, with an optional sign prefix `+`, `-`, `i` or `-i`.
Qubit `j` carries `X` when `x[j]`, `Z` when `z[j]`, and `Y` when both.

Supports `*` (with the phase), `==`, [`commutes`](@ref), [`weight`](@ref),
`matrix`, and application to a statevector as `p * ψ`.
"""
struct PauliOp
    x::BitVector
    z::BitVector
    phase::Int            # power of i, mod 4
    function PauliOp(x::AbstractVector, z::AbstractVector, phase::Integer=0)
        length(x) == length(z) || throw(ArgumentError("x and z parts differ in length"))
        new(_gf2(x), _gf2(z), mod(phase, 4))
    end
end

function PauliOp(s::AbstractString)
    phase = 0
    for (pre, ph) in (("-i", 3), ("+i", 1), ("i", 1), ("-", 2), ("+", 0))
        if startswith(s, pre)
            phase = ph
            s = s[length(pre)+1:end]
            break
        end
    end
    all(in("IXYZ"), s) || throw(ArgumentError("bad Pauli string \"$s\""))
    PauliOp([ch in ('X', 'Y') for ch in s], [ch in ('Z', 'Y') for ch in s], phase)
end

Base.length(p::PauliOp) = length(p.x)
Base.:(==)(p::PauliOp, q::PauliOp) = p.x == q.x && p.z == q.z && p.phase == q.phase
Base.hash(p::PauliOp, h::UInt) = hash((p.x, p.z, p.phase), h)
Base.:-(p::PauliOp) = PauliOp(p.x, p.z, p.phase + 2)

_pauli_letters(p::PauliOp) =
    join((x && z ? 'Y' : x ? 'X' : z ? 'Z' : 'I') for (x, z) in zip(p.x, p.z))
Base.string(p::PauliOp) = ("+", "+i", "-", "-i")[p.phase+1] * _pauli_letters(p)
Base.show(io::IO, p::PauliOp) = print(io, "PauliOp(\"", string(p), "\")")

"""
    weight(p) -> Int

Number of qubits on which a Pauli acts non-trivially.
"""
weight(p::PauliOp) = count(p.x .| p.z)

"""
    commutes(p, q) -> Bool

Whether two Paulis commute: their symplectic product `x·z′ + z·x′` is even.
"""
function commutes(p::PauliOp, q::PauliOp)
    length(p) == length(q) || throw(ArgumentError("Paulis act on different numbers of qubits"))
    iseven(count(p.x .& q.z) + count(p.z .& q.x))
end

# Exponent of i in (single-qubit Pauli x1,z1)·(x2,z2) — Aaronson–Gottesman.
@inline function _pauli_g(x1::Bool, z1::Bool, x2::Bool, z2::Bool)
    !x1 && !z1 && return 0
    x1 && z1 && return Int(z2) - Int(x2)
    x1 && return Int(z2) * (2 * Int(x2) - 1)
    Int(x2) * (1 - 2 * Int(z2))
end

function Base.:*(p::PauliOp, q::PauliOp)
    length(p) == length(q) || throw(ArgumentError("Paulis act on different numbers of qubits"))
    ph = p.phase + q.phase
    for j in 1:length(p)
        ph += _pauli_g(p.x[j], p.z[j], q.x[j], q.z[j])
    end
    PauliOp(p.x .⊻ q.x, p.z .⊻ q.z, ph)
end

"""
    matrix(p::PauliOp) -> Matrix{ComplexF64}

The dense `2ⁿ × 2ⁿ` matrix, qubit 1 most significant.
"""
matrix(p::PauliOp) = im^p.phase .* pauli(_pauli_letters(p))

# P|b⟩ = i^(phase + #Y) (-1)^|b ∧ z| |b ⊕ x⟩, with Y = iXZ.
function Base.:*(p::PauliOp, ψ::AbstractVector)
    n = length(p)
    length(ψ) == 1 << n || throw(ArgumentError("state has length $(length(ψ)), expected 2^$n"))
    xm = sum((1 << (n - j) for j in 1:n if p.x[j]); init=0)
    zm = sum((1 << (n - j) for j in 1:n if p.z[j]); init=0)
    c = im^(p.phase + count(p.x .& p.z))
    out = similar(ψ, ComplexF64)
    for b in 0:(1 << n)-1
        out[(b ⊻ xm)+1] = c * (isodd(count_ones(b & zm)) ? -ψ[b+1] : ψ[b+1])
    end
    out
end

"""
    pauli!(c, p, qubits=1:length(p)) -> c

Append the Pauli `p` to a circuit as one-qubit `X`, `Y`, `Z` gates on
`qubits`, its phase going into the global phase.  For injecting errors and
applying corrections.
"""
function pauli!(c::Circuit, p::PauliOp, qubits::AbstractVector{<:Integer}=1:length(p))
    length(qubits) == length(p) || throw(ArgumentError("need $(length(p)) qubits"))
    for (j, q) in enumerate(qubits)
        x, z = p.x[j], p.z[j]
        x && z ? push!(c, Y(), q) : x ? push!(c, X(), q) : z ? push!(c, Z(), q) : nothing
    end
    c.global_phase += p.phase * π / 2
    c
end

# Symplectic row [x | z] and its commutation partner [z | x].
_symp(p::PauliOp) = vcat(p.x, p.z)
_symp_matrix(ps) = BitMatrix(reduce(vcat, permutedims(_symp(p)) for p in ps))
_symp_inner(u::AbstractVector, v::AbstractVector) =
    (n = length(u) ÷ 2; isodd(count(u[1:n] .& v[n+1:end]) + count(u[n+1:end] .& v[1:n])))
_pauli_from_symp(v::AbstractVector) = (n = length(v) ÷ 2; PauliOp(v[1:n], v[n+1:end]))

# --- stabilizer codes --------------------------------------------------------

"""
    StabilizerCode(stabilizers; logical_x, logical_z, name="", distance=nothing)

The code fixed by commuting, independent, Hermitian Pauli `stabilizers` on
`n` qubits: `k = n - r` logical qubits for `r` generators.

Logical operators may be given (they are checked: each commutes with every
stabilizer, `X̄ᵢ` anticommutes with `Z̄ᵢ` alone).  Otherwise they are computed
by symplectic Gram–Schmidt on the normalizer; CSS codes ([`css_code`](@ref))
get pure-`X` and pure-`Z` logicals.

`length(code)` is `n`; [`dimension`](@ref) is `k`; [`code_distance`](@ref)
is `d`.  A known `distance` can be recorded (the catalogue does) so that
simulations need not search for it; `code_distance` always searches.
"""
struct StabilizerCode
    stabilizers::Vector{PauliOp}
    logical_x::Vector{PauliOp}
    logical_z::Vector{PauliOp}
    name::String
    distance::Union{Int,Nothing}      # known distance, if recorded
end

function StabilizerCode(stabs::AbstractVector; logical_x=nothing, logical_z=nothing,
                        name::AbstractString="", distance::Union{Integer,Nothing}=nothing)
    S = [s isa PauliOp ? s : PauliOp(s) for s in stabs]
    isempty(S) && throw(ArgumentError("need at least one stabilizer"))
    n = length(S[1])
    all(length(s) == n for s in S) || throw(ArgumentError("stabilizers act on different numbers of qubits"))
    all(iseven(s.phase) for s in S) || throw(ArgumentError("stabilizers must be Hermitian (phase ±1)"))
    for i in eachindex(S), j in i+1:lastindex(S)
        commutes(S[i], S[j]) || throw(ArgumentError("stabilizers $(S[i]) and $(S[j]) anticommute"))
    end
    gf2_rank(_symp_matrix(S)) == length(S) || throw(ArgumentError("stabilizers are not independent"))

    if logical_x === nothing && logical_z === nothing
        LX, LZ = _logicals(S, n)
    else
        LX = [p isa PauliOp ? p : PauliOp(p) for p in logical_x]
        LZ = [p isa PauliOp ? p : PauliOp(p) for p in logical_z]
    end
    code = StabilizerCode(S, LX, LZ, String(name), distance === nothing ? nothing : Int(distance))
    _check_logicals(code)
    code
end

function _check_logicals(code::StabilizerCode)
    k = length(code) - length(code.stabilizers)
    length(code.logical_x) == length(code.logical_z) == k ||
        throw(ArgumentError("need $k logical X and $k logical Z operators"))
    L = vcat(code.logical_x, code.logical_z)
    for l in L, s in code.stabilizers
        commutes(l, s) || throw(ArgumentError("logical $l anticommutes with stabilizer $s"))
    end
    for i in 1:k, j in 1:k
        commutes(code.logical_x[i], code.logical_z[j]) == (i != j) ||
            throw(ArgumentError("logical X̄$i and Z̄$j must anticommute exactly when i == j"))
        i < j || continue
        commutes(code.logical_x[i], code.logical_x[j]) && commutes(code.logical_z[i], code.logical_z[j]) ||
            throw(ArgumentError("logical operators $i and $j must commute"))
    end
    nothing
end

# Logical operators by symplectic Gram–Schmidt on the normalizer.  CSS codes
# are handled separately so the logicals stay pure X and pure Z.
function _logicals(S::Vector{PauliOp}, n::Int)
    _iscss(S) && return _css_logicals(S, n)
    Sm = _symp_matrix(S)
    # normalizer: v with ⟨v, s⟩ = 0, i.e. [z_s | x_s] · v = 0
    N = gf2_nullspace(hcat(Sm[:, n+1:end], Sm[:, 1:n]))
    W = [N[i, :] for i in 1:size(N, 1)]
    span = copy(Sm)
    LX, LZ = PauliOp[], PauliOp[]
    while true
        i = findfirst(u -> !_gf2_inrowspace(span, u), W)
        i === nothing && break
        u = popat!(W, i)
        j = findfirst(w -> _symp_inner(u, w), W)
        j === nothing && throw(ArgumentError("normalizer element without a partner — not a valid code"))
        w = popat!(W, j)
        for t in eachindex(W)                      # make the rest orthogonal to u, w
            v = W[t]
            _symp_inner(v, w) && (v = v .⊻ u)
            _symp_inner(v, u) && (v = v .⊻ w)
            W[t] = v
        end
        push!(LX, _pauli_from_symp(u)); push!(LZ, _pauli_from_symp(w))
        span = vcat(span, permutedims(u), permutedims(w))
    end
    LX, LZ
end

_iscss(S) = all(s -> !any(s.x) || !any(s.z), S)

function _css_logicals(S::Vector{PauliOp}, n::Int)
    Hx = BitMatrix(reduce(vcat, [permutedims(s.x) for s in S if any(s.x)]; init=falses(0, n)))
    Hz = BitMatrix(reduce(vcat, [permutedims(s.z) for s in S if any(s.z)]; init=falses(0, n)))
    # X logicals: in ker(Hz), outside rowspace(Hx); Z logicals likewise
    pick(ker, rows) = begin
        out = BitVector[]
        span = rows
        for i in 1:size(ker, 1)
            v = ker[i, :]
            if !_gf2_inrowspace(span, v)
                push!(out, v)
                span = vcat(span, permutedims(v))
            end
        end
        out
    end
    Lx = pick(gf2_nullspace(Hz), Hx)
    Lz = pick(gf2_nullspace(Hx), Hz)
    k = length(Lx)
    k == 0 && return PauliOp[], PauliOp[]
    LxM = BitMatrix(reduce(vcat, permutedims.(Lx)))
    LzM = BitMatrix(reduce(vcat, permutedims.(Lz)))
    # pair them: Lz ← (M⁻¹)ᵀ Lz with M = Lx Lzᵀ, so that Lx Lzᵀ = I
    M = _gf2mul(LxM, permutedims(LzM))
    LzM = _gf2mul(permutedims(_gf2_inv(M)), LzM)
    [PauliOp(LxM[i, :], falses(n)) for i in 1:k], [PauliOp(falses(n), LzM[i, :]) for i in 1:k]
end

function _gf2_inv(M::AbstractMatrix)
    k = size(M, 1)
    R, piv = gf2_rref(hcat(_gf2(M), BitMatrix(LinearAlgebra.I(k))))
    piv[1:min(end, k)] == 1:k || throw(ArgumentError("matrix is singular over GF(2)"))
    R[:, k+1:end]
end

Base.length(c::StabilizerCode) = length(c.stabilizers[1])
dimension(c::StabilizerCode) = length(c) - length(c.stabilizers)

"""
    stabilizers(code) -> Vector{PauliOp}
"""
stabilizers(c::StabilizerCode) = copy(c.stabilizers)

"""
    logical_operators(code) -> (X̄, Z̄)

The logical `X` and `Z` operators, one pair per logical qubit.
"""
logical_operators(c::StabilizerCode) = (copy(c.logical_x), copy(c.logical_z))

"""
    is_css(code) -> Bool

Whether every stabilizer is purely `X`-type or purely `Z`-type — a
Calderbank–Shor–Steane code, built from two classical codes.
"""
is_css(c::StabilizerCode) = _iscss(c.stabilizers)

function Base.show(io::IO, c::StabilizerCode)
    print(io, isempty(c.name) ? "StabilizerCode" : c.name, " [[", length(c), ", ", dimension(c), "]]")
end

"""
    syndrome(code, E::PauliOp) -> BitVector

Bit `i` is set when the error `E` anticommutes with stabilizer `i` — the
pattern of `-1` outcomes a syndrome measurement would report.
"""
syndrome(c::StabilizerCode, E::PauliOp) = BitVector([!commutes(E, s) for s in c.stabilizers])

# Packed (x, z) words for fast commutation tests, n ≤ 64.
_pack(v::AbstractVector) = sum((UInt64(1) << (j - 1) for j in eachindex(v) if v[j]); init=UInt64(0))
_packed(p::PauliOp) = (_pack(p.x), _pack(p.z))

"""
    code_distance(code; maxweight=length(code)) -> Int

The smallest weight of a Pauli that commutes with every stabilizer without
being a stabilizer (up to phase) — an undetectable logical error.  A code of
distance `d` corrects any `⌊(d-1)/2⌋` errors.  Exhaustive over weights in
order, on packed bit words; fine for the small codes here.
"""
function code_distance(c::StabilizerCode; maxweight::Integer=length(c))
    n = length(c)
    n <= 64 || throw(ArgumentError("distance search supports n ≤ 64"))
    stabs = [_packed(s) for s in c.stabilizers]
    Sm = _symp_matrix(c.stabilizers)
    for t in 1:maxweight
        for pos in _combinations(n, t), kinds in Iterators.product(ntuple(_ -> 1:3, t)...)
            ex = UInt64(0); ez = UInt64(0)
            for (q, kd) in zip(pos, kinds)
                kd != 3 && (ex |= UInt64(1) << (q - 1))       # 1 = X, 2 = Y, 3 = Z
                kd != 1 && (ez |= UInt64(1) << (q - 1))
            end
            all(((count_ones(ex & sz) + count_ones(ez & sx)) & 1) == 0 for (sx, sz) in stabs) || continue
            v = vcat(BitVector([(ex >> (j - 1)) & 1 == 1 for j in 1:n]),
                     BitVector([(ez >> (j - 1)) & 1 == 1 for j in 1:n]))
            _gf2_inrowspace(Sm, v) || return t
        end
    end
    throw(ArgumentError("no logical operator of weight ≤ $maxweight"))
end

"""
    lookup_decoder(code; maxweight=1) -> Dict{BitVector, PauliOp}

A table from syndrome to the lowest-weight Pauli error producing it, over all
errors of weight up to `maxweight` — what a distance-3 code guarantees to fix
at the default.  Correct by applying the table's Pauli: if the real error had
weight ≤ `maxweight`, the product is a stabilizer.
"""
function lookup_decoder(c::StabilizerCode; maxweight::Integer=1)
    n, r = length(c), length(c.stabilizers)
    n <= 64 || throw(ArgumentError("lookup_decoder supports n ≤ 64"))
    stabs = [_packed(s) for s in c.stabilizers]
    table = Dict{BitVector,PauliOp}(falses(r) => PauliOp(falses(n), falses(n)))
    for t in 1:maxweight
        for pos in _combinations(n, t), kinds in Iterators.product(ntuple(_ -> 1:3, t)...)
            ex = UInt64(0); ez = UInt64(0)
            for (q, kd) in zip(pos, kinds)
                kd != 3 && (ex |= UInt64(1) << (q - 1))
                kd != 1 && (ez |= UInt64(1) << (q - 1))
            end
            s = BitVector([isodd(count_ones(ex & sz) + count_ones(ez & sx)) for (sx, sz) in stabs])
            haskey(table, s) && continue
            table[s] = PauliOp(BitVector([(ex >> (j - 1)) & 1 == 1 for j in 1:n]),
                               BitVector([(ez >> (j - 1)) & 1 == 1 for j in 1:n]))
        end
    end
    table
end

# --- the catalogue -----------------------------------------------------------

"""
    css_code(Hx, Hz; name="", distance=nothing) -> StabilizerCode

The CSS code with `X`-type stabilizers from the rows of `Hx` and `Z`-type
from the rows of `Hz` (which must satisfy `Hx Hzᵀ = 0`).  `Z` errors are
caught by `Hx`, `X` errors by `Hz`.
"""
function css_code(Hx::AbstractMatrix, Hz::AbstractMatrix; name::AbstractString="",
                  distance::Union{Integer,Nothing}=nothing)
    size(Hx, 2) == size(Hz, 2) || throw(ArgumentError("Hx and Hz have different lengths"))
    any(_gf2mul(Hx, permutedims(_gf2(Hz)))) && throw(ArgumentError("Hx Hzᵀ ≠ 0: the checks do not commute"))
    n = size(Hx, 2)
    S = vcat([PauliOp(_gf2(Hx[i, :]), falses(n)) for i in 1:size(Hx, 1)],
             [PauliOp(falses(n), _gf2(Hz[i, :])) for i in 1:size(Hz, 1)])
    StabilizerCode(S; name=name, distance=distance)
end

"""
    repetition_code(n) -> StabilizerCode

The `n`-qubit bit-flip code: stabilizers `Zᵢ Zᵢ₊₁`, logicals `X̄ = X⊗…⊗X`,
`Z̄ = Z₁`.  Corrects `⌊(n-1)/2⌋` bit flips but no phase flips, so its
quantum distance is 1.
"""
function repetition_code(n::Integer)
    n >= 2 || throw(ArgumentError("need at least 2 qubits"))
    S = [PauliOp(falses(n), BitVector(j in (i, i + 1) for j in 1:n)) for i in 1:n-1]
    StabilizerCode(S; logical_x=[PauliOp("X"^n)], logical_z=[PauliOp("Z" * "I"^(n - 1))],
                   name="repetition code", distance=1)
end

"""
    five_qubit_code() -> StabilizerCode

The `[[5, 1, 3]]` code, the smallest that corrects any single-qubit error:
stabilizers `XZZXI` and its cyclic shifts.  Not CSS.
"""
five_qubit_code() = StabilizerCode(["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"];
                                   logical_x=["XXXXX"], logical_z=["ZZZZZ"], name="five-qubit code", distance=3)

"""
    shor_code() -> StabilizerCode

Shor's `[[9, 1, 3]]` code: three bit-flip codes inside a phase-flip code.

Logicals are the CSS ones, `X̄` pure `X` and `Z̄` pure `Z`, like every CSS
code here.  Textbooks often swap them (`X̄ = Z⊗9`), which exchanges `|0̄⟩`
with `|+̄⟩`; it is the same code.
"""
shor_code() = StabilizerCode(["ZZIIIIIII", "IZZIIIIII", "IIIZZIIII", "IIIIZZIII",
                              "IIIIIIZZI", "IIIIIIIZZ", "XXXXXXIII", "IIIXXXXXX"]; name="Shor code", distance=3)

"""
    steane_code() -> StabilizerCode

Steane's `[[7, 1, 3]]` code: the CSS code of the `[7, 4, 3]` Hamming code
([`hamming_code`](@ref)`(3)`) against itself.  Equal to
[`quantum_reed_muller`](@ref)`(3)`.
"""
function steane_code()
    H = parity_check_matrix(hamming_code(3))
    css_code(H, H; name="Steane code", distance=3)
end

"""
    quantum_reed_muller(m) -> StabilizerCode

The `[[2ᵐ - 1, 1, 3]]` quantum Reed–Muller code, from punctured Reed–Muller
codes on the nonzero points of GF(2)ᵐ: `X` stabilizers are the `m` linear
monomials `xᵢ`, `Z` stabilizers the monomials of degree `1 … m-2`.  `m = 3` is
the Steane code; from `m = 4` (`[[15, 1, 3]]`) `T` applied to every qubit is
a logical gate — the reason it underlies magic-state distillation.
"""
function quantum_reed_muller(m::Integer)
    m >= 3 || throw(ArgumentError("quantum Reed–Muller codes need m ≥ 3"))
    pts = 1:(1 << m)-1                                 # punctured: drop x = 0
    row(S) = BitVector([(x & S) == S for x in pts])
    Hx = BitMatrix(reduce(vcat, permutedims(row(S)) for S in _monomials(1, m) if count_ones(S) == 1))
    Hz = BitMatrix(reduce(vcat, permutedims(row(S)) for S in _monomials(m - 2, m) if count_ones(S) >= 1))
    css_code(Hx, Hz; name="quantum Reed–Muller code", distance=3)
end

"""
    rotated_surface_code(d) -> StabilizerCode

The distance-`d` rotated surface code, `[[d², 1, d]]` for odd `d ≥ 3`: data
qubits on a `d × d` grid (qubit `(i, j)` is `(i-1)d + j`), weight-4 checks on
the faces in a checkerboard of `X` and `Z`, and weight-2 checks on the
boundary — `X` along the top and bottom, `Z` along the left and right.
"""
function rotated_surface_code(d::Integer)
    d >= 3 && isodd(d) || throw(ArgumentError("need odd d ≥ 3"))
    n = d^2
    q(i, j) = (i - 1) * d + j
    S = PauliOp[]
    for i in 0:d, j in 0:d                         # face with corners (i..i+1, j..j+1)
        isX = iseven(i + j)
        corners = [(a, b) for a in (i, i + 1), b in (j, j + 1) if 1 <= a <= d && 1 <= b <= d]
        length(corners) == 4 || length(corners) == 2 || continue
        if length(corners) == 2
            top_bottom = i == 0 || i == d
            (top_bottom && isX) || (!top_bottom && !isX) || continue
        end
        v = falses(n)
        for (a, b) in corners
            v[q(a, b)] = true
        end
        push!(S, isX ? PauliOp(v, falses(n)) : PauliOp(falses(n), v))
    end
    StabilizerCode(S; name="rotated surface code", distance=d)
end
