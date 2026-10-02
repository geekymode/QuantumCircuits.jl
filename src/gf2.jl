# ---------------------------------------------------------------------------
# Linear algebra over GF(2)
#
# Codes are subspaces of GF(2)ⁿ, so everything in error correction — parity
# checks, syndromes, logical operators — is row reduction mod 2.  Matrices are
# `BitMatrix`; any integer or Bool matrix is accepted and read mod 2.
# ---------------------------------------------------------------------------

_gf2(A::AbstractMatrix) = BitMatrix(isodd.(A))
_gf2(v::AbstractVector) = BitVector(isodd.(v))

"""
    gf2_rref(A) -> (R, pivots)

Reduced row echelon form of `A` over GF(2): every pivot is the only `1` in its
column.  `R` keeps all rows (zero rows last); `pivots[i]` is the pivot column
of row `i` for the first `length(pivots)` rows.
"""
function gf2_rref(A::AbstractMatrix)
    R = _gf2(A)
    m, n = size(R)
    pivots = Int[]
    r = 1
    for col in 1:n
        r > m && break
        p = findfirst(i -> R[i, col], r:m)
        p === nothing && continue
        p += r - 1
        if p != r
            R[[r, p], :] = R[[p, r], :]
        end
        for i in 1:m
            if i != r && R[i, col]
                R[i, :] .⊻= R[r, :]
            end
        end
        push!(pivots, col)
        r += 1
    end
    R, pivots
end

"""
    gf2_rank(A) -> Int

Rank of `A` over GF(2).
"""
gf2_rank(A::AbstractMatrix) = length(gf2_rref(A)[2])

"""
    gf2_nullspace(A) -> BitMatrix

A basis of `{x : A x = 0 (mod 2)}`, one vector per row, so the result `N`
satisfies `A * N' .% 2 == 0`.  For a parity-check matrix this is a generator
matrix of the code, and vice versa.
"""
function gf2_nullspace(A::AbstractMatrix)
    R, pivots = gf2_rref(A)
    n = size(R, 2)
    free = setdiff(1:n, pivots)
    N = falses(length(free), n)
    for (j, f) in enumerate(free)
        N[j, f] = true
        for (i, p) in enumerate(pivots)
            N[j, p] = R[i, f]
        end
    end
    N
end

# Rows of A reduced to a basis: rref with the zero rows dropped.
_gf2_rowbasis(A::AbstractMatrix) = (R = gf2_rref(A); R[1][1:length(R[2]), :])

# Whether v lies in the row space of A.
_gf2_inrowspace(A::AbstractMatrix, v::AbstractVector) =
    gf2_rank(vcat(_gf2(A), permutedims(_gf2(v)))) == gf2_rank(A)

# A * B mod 2 for Bool/Bit matrices.
_gf2mul(A::AbstractMatrix, B::AbstractMatrix) = BitMatrix(isodd.(Int.(A) * Int.(B)))
_gf2mul(A::AbstractMatrix, v::AbstractVector) = BitVector(isodd.(Int.(A) * Int.(v)))
