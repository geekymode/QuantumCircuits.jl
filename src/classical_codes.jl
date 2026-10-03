# ---------------------------------------------------------------------------
# Classical linear codes
#
# A binary linear [n, k, d] code is a k-dimensional subspace of GF(2)ⁿ,
# given by a generator matrix G (codewords are its row span) or a parity-check
# matrix H (codewords are its null space).  Quantum CSS codes are built from
# pairs of these, so they come first.
#
# Reed–Muller codes are the thread through this package: RM(r, m) is the set of
# evaluation tables of degree-≤r polynomials in m Boolean variables — the same
# objects as the parity expansions of `phase_polynomial` — and their local
# decodability is what Nehoran and Yuen's three-query identity exploits.
# ---------------------------------------------------------------------------

"""
    LinearCode(G)

The binary linear code spanned by the rows of `G` (read mod 2).  Dependent
rows are dropped; independent rows are kept as given, so structured
generators (like the monomials of [`reed_muller`](@ref)) survive.  The
parity-check matrix is computed.

`length(c)` is `n`, [`dimension`](@ref) is `k`.
"""
struct LinearCode
    G::BitMatrix          # k × n generator, full row rank
    H::BitMatrix          # (n-k) × n parity check
end

function LinearCode(G::AbstractMatrix)
    Gb = _gf2(G)
    # keep the given rows that are independent of the ones before them
    keep = Int[]
    for i in 1:size(Gb, 1)
        gf2_rank(Gb[vcat(keep, i), :]) > length(keep) && push!(keep, i)
    end
    Gk = Gb[keep, :]
    LinearCode(Gk, gf2_nullspace(Gk))
end

Base.length(c::LinearCode) = size(c.G, 2)

"""
    dimension(c) -> Int

Number of message bits `k` of a code (classical or quantum).
"""
dimension(c::LinearCode) = size(c.G, 1)

"""
    generator_matrix(c) -> BitMatrix

`k × n` generator matrix; codewords are `m' G` for messages `m`.
"""
generator_matrix(c::LinearCode) = copy(c.G)

"""
    parity_check_matrix(c) -> BitMatrix

`(n-k) × n` parity-check matrix `H`: `w` is a codeword iff `H w = 0`.
"""
parity_check_matrix(c::LinearCode) = copy(c.H)

Base.show(io::IO, c::LinearCode) = print(io, "LinearCode[", length(c), ", ", dimension(c), "]")

"""
    encode(c, m) -> BitVector

The codeword `mᵀ G` for message `m` (length `k`).
"""
function encode(c::LinearCode, m::AbstractVector)
    length(m) == dimension(c) || throw(ArgumentError("message must have $(dimension(c)) bits"))
    _gf2mul(permutedims(_gf2(m)), c.G)[1, :]
end

"""
    syndrome(c, w) -> BitVector

`H w` mod 2: zero exactly when `w` is a codeword, and for `w = codeword + e`
it depends only on the error `e`.
"""
function syndrome(c::LinearCode, w::AbstractVector)
    length(w) == length(c) || throw(ArgumentError("word must have $(length(c)) bits"))
    _gf2mul(c.H, _gf2(w))
end

"""
    iscodeword(c, w) -> Bool
"""
iscodeword(c::LinearCode, w::AbstractVector) = !any(syndrome(c, w))

"""
    codewords(c) -> Vector{BitVector}

All `2ᵏ` codewords.  Exponential; `k ≤ 20` only.
"""
function codewords(c::LinearCode)
    k = dimension(c)
    k <= 20 || throw(ArgumentError("2^$k codewords is too many to list"))
    [encode(c, digits(m; base=2, pad=k)) for m in 0:(1 << k)-1]
end

"""
    minimum_distance(c) -> Int

Smallest weight of a nonzero codeword.  A code of distance `d` detects
`d - 1` errors and corrects `⌊(d-1)/2⌋`.

Exhaustive, from whichever side is cheaper.  Up to `k = 26` it walks all
`2ᵏ` codewords **in Gray-code order**: consecutive messages differ in one bit,
so each step is a single XOR with one generator row
([`gray_flip_position`](@ref)), on words packed into machine integers.
Otherwise it finds the smallest set of columns of `H` summing to zero — a
weight-`d` codeword is exactly such a set.
"""
function minimum_distance(c::LinearCode)
    k, n = dimension(c), length(c)
    k == 0 && return typemax(Int)
    k <= 26 && n <= 128 && return _gray_min_weight(c)
    cols = [c.H[:, j] for j in 1:n]
    for t in 1:n
        for pos in _combinations(n, t)
            !any(reduce(.⊻, cols[pos])) && return t
        end
        binomial(n, t + 1) > 10^7 && throw(ArgumentError("distance search beyond weight $t is too large"))
    end
end

# Minimum nonzero weight over the code, by a Gray-code walk of the messages.
function _gray_min_weight(c::LinearCode)
    pack(row) = sum(UInt128(1) << (j - 1) for j in 1:length(row) if row[j]; init=UInt128(0))
    rows = [pack(c.G[i, :]) for i in 1:dimension(c)]
    w, best = UInt128(0), typemax(Int)
    for t in 1:(1 << dimension(c))-1
        w ⊻= rows[gray_flip_position(t)+1]       # gray(t) and gray(t-1) differ here
        best = min(best, count_ones(w))
    end
    best
end

"""
    dual(c) -> LinearCode

The dual code `C⊥ = {x : x·c = 0 for all c ∈ C}`: generator and parity-check
matrices swap roles.
"""
dual(c::LinearCode) = LinearCode(c.H, c.G)

"""
    puncture(c, positions) -> LinearCode

Delete the coordinates `positions` from every codeword.
"""
function puncture(c::LinearCode, positions)
    keep = setdiff(1:length(c), positions)
    LinearCode(c.G[:, keep])
end

"""
    hamming_code(r) -> LinearCode

The `[2ʳ - 1, 2ʳ - 1 - r, 3]` Hamming code: column `j` of the parity-check
matrix is `j` in binary, so the syndrome of a single bit flip *is* its
position.  (Not to be confused with [`hamming`](@ref), the Hamming distance.)
"""
function hamming_code(r::Integer)
    r >= 2 || throw(ArgumentError("Hamming codes need r ≥ 2"))
    n = (1 << r) - 1
    H = BitMatrix([(j >> (r - i)) & 1 == 1 for i in 1:r, j in 1:n])
    LinearCode(gf2_nullspace(H), H)
end

# Variable xᵢ is bit (m - i) of a point's index — x₁ most significant, as for
# qubits — and a monomial is the mask of its variables.  Ordered by degree,
# then by variable list: 1, x₁, …, x_m, x₁x₂, x₁x₃, …
_rm_vars(S::Integer, m::Integer) = [i for i in 1:m if (S >> (m - i)) & 1 == 1]
_monomials(r::Integer, m::Integer) =
    sort!([S for S in 0:(1 << m)-1 if count_ones(S) <= r]; by=S -> (count_ones(S), _rm_vars(S, m)))

# Evaluation table of a monomial: 1 at the points containing all its variables.
_rm_row(S::Integer, m::Integer) = BitVector([(x & S) == S for x in 0:(1 << m)-1])

"""
    reed_muller(r, m) -> LinearCode

The Reed–Muller code `RM(r, m)`: evaluation tables, over all `2ᵐ` points of
GF(2)ᵐ, of the Boolean polynomials of degree at most `r`.  Parameters
`[2ᵐ, Σ_{i≤r} C(m,i), 2^(m-r)]`, and `RM(r, m)⊥ = RM(m-r-1, m)`.

Rows of the generator matrix are the monomials, by degree — `1`, `x₁ … x_m`,
`x₁x₂ …` — so `encode(c, coeffs)` evaluates the polynomial with those
coefficients.  Coordinate `x + 1` is the point whose binary digits, most
significant first, are `x₁ … x_m`.
"""
function reed_muller(r::Integer, m::Integer)
    0 <= r <= m || throw(ArgumentError("need 0 ≤ r ≤ m, got r = $r, m = $m"))
    G = BitMatrix(reduce(vcat, permutedims(_rm_row(S, m)) for S in _monomials(r, m)))
    LinearCode(G, gf2_nullspace(G))
end

"""
    syndrome_decode(c, w; maxweight=…) -> Union{BitVector, Nothing}

Nearest-codeword decoding by syndrome lookup: find the lowest-weight error
with the syndrome of `w` and remove it.  Searches errors up to `maxweight`
(default `⌊(d-1)/2⌋`, what the code is guaranteed to correct) and returns
`nothing` if none matches.
"""
function syndrome_decode(c::LinearCode, w::AbstractVector;
                         maxweight::Integer=(minimum_distance(c) - 1) ÷ 2)
    s = syndrome(c, w)
    word = _gf2(w)
    !any(s) && return word
    n = length(c)
    for t in 1:maxweight
        for pos in _combinations(n, t)
            e = falses(n); e[pos] .= true
            syndrome(c, e) == s && return word .⊻ e
        end
    end
    nothing
end

# All t-subsets of 1:n, lexicographic.
function _combinations(n::Integer, t::Integer)
    out = Vector{Vector{Int}}()
    comb = collect(1:t)
    t == 0 && return [Int[]]
    t > n && return out
    while true
        push!(out, copy(comb))
        i = t
        while i >= 1 && comb[i] == n - t + i
            i -= 1
        end
        i == 0 && return out
        comb[i] += 1
        for j in i+1:t
            comb[j] = comb[j-1] + 1
        end
    end
end

"""
    rm_local_decode(w, r, m, x; trials=15, rng) -> Bool

Recover coordinate `x` (a point of GF(2)ᵐ given as an integer `0 … 2ᵐ-1`) of
a corrupted `RM(r, m)` codeword `w` while reading only a few other
coordinates — Reed's majority-logic decoding.

A degree-≤`r` polynomial sums to zero over every affine subspace of
dimension `r + 1`, so `f(x) = Σ f(x + v)` over the nonzero `v` of any
`(r+1)`-dimensional subspace: `2^(r+1) - 1` queries.  Each trial uses a random
subspace and the majority vote wins, so a small fraction of errors is
tolerated.

This is the classical shadow of Nehoran and Yuen's three-query identity: over
the reals a degree-2 form is pinned by **three** collinear points; over GF(2),
`2³ - 1 = 7` points of a 3-dimensional subspace are needed.
"""
function rm_local_decode(w::AbstractVector, r::Integer, m::Integer, x::Integer;
                         trials::Integer=15, rng=Random.default_rng())
    length(w) == 1 << m || throw(ArgumentError("an RM(r, $m) word has $(1 << m) bits"))
    0 <= x < 1 << m || throw(ArgumentError("point must lie in 0:$((1 << m) - 1)"))
    r + 1 <= m || throw(ArgumentError("local decoding needs r < m"))
    word = _gf2(w)
    votes = 0
    for _ in 1:trials
        # random (r+1)-dimensional subspace: draw independent vectors
        basis = Int[]
        while length(basis) < r + 1
            v = rand(rng, 1:(1 << m)-1)
            span = [reduce(⊻, (basis[i] for i in 1:length(basis) if (t >> (i - 1)) & 1 == 1); init=0)
                    for t in 0:(1 << length(basis))-1]
            v in span || push!(basis, v)
        end
        acc = false
        for t in 1:(1 << (r + 1))-1
            v = reduce(⊻, (basis[i] for i in 1:r+1 if (t >> (i - 1)) & 1 == 1); init=0)
            acc ⊻= word[(x ⊻ v) + 1]
        end
        votes += acc
    end
    2votes > trials
end
