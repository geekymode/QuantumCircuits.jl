# ---------------------------------------------------------------------------
# Gray code
#
# The reflected binary Gray code orders the integers 0 … 2^n-1 so that
# consecutive entries differ in exactly one bit.  In circuit synthesis that
# single-bit property is what turns an O(n) uncompute/recompute ladder into a
# single CNOT per step: the parity register you carry along only ever needs one
# incremental update.
# ---------------------------------------------------------------------------

"""
    gray(i) -> Int

Standard reflected binary Gray code of `i`: `i ⊻ (i >> 1)`.

`gray` is a bijection on `0:2^n-1` for every `n`, and `gray(i)` and `gray(i+1)`
always differ in exactly one bit.

```jldoctest
julia> gray.(0:7)
8-element Vector{Int64}:
 0
 1
 3
 2
 6
 7
 5
 4
```
"""
gray(i::Integer) = Int(i) ⊻ (Int(i) >> 1)

"""
    ungray(g) -> Int

Inverse of [`gray`](@ref): recovers `i` from `gray(i)` by prefix-XOR.
"""
function ungray(g::Integer)
    g >= 0 || throw(DomainError(g, "ungray expects a non-negative integer"))
    x = Int(g)
    sh = 1
    while sh < 8 * sizeof(x)
        x ⊻= (x >> sh)
        sh <<= 1
    end
    x
end

"""
    graycode(n) -> Vector{Int}

The full `n`-bit Gray code sequence, `gray.(0:2^n-1)`.
"""
graycode(n::Integer) = [gray(i) for i in 0:(1 << n)-1]

"""
    gray_flip_position(i) -> Int

Index (LSB = 0) of the single bit that changes between `gray(i-1)` and `gray(i)`.

This is just `trailing_zeros(i)`, which is the whole reason the Gray-code CNOT
ladder is cheap to generate: the control qubit for step `i` is a `ctz` away.
"""
gray_flip_position(i::Integer) = trailing_zeros(Int(i))

"""
    gray_flip_position(i, k) -> Int

Cyclic version over `k` bits.  For `i == 2^k` the sequence wraps from
`gray(2^k - 1)` back to `gray(0) == 0`, which flips the most significant bit,
so this returns `k-1`.  Closing the cycle is what makes a multiplexor's CNOTs
cancel out to the identity on the control register.
"""
function gray_flip_position(i::Integer, k::Integer)
    i == (1 << k) && return Int(k) - 1
    gray_flip_position(i)
end

"""
    gray_flip_positions(k) -> Vector{Int}

The `2^k` bit positions flipped by one full cyclic walk of the `k`-bit Gray
code.  `cumulative XOR` of `1 << p` over this list reproduces
`gray(1), gray(2), …, gray(2^k-1), 0`.
"""
gray_flip_positions(k::Integer) = [gray_flip_position(i, k) for i in 1:(1 << k)]

"""
    hamming(a, b) -> Int

Hamming distance between two non-negative integers.
"""
hamming(a::Integer, b::Integer) = count_ones(Int(a) ⊻ Int(b))

"""
    gray_adjacent(a, b) -> Bool

`true` when `a` and `b` differ in exactly one bit.
"""
gray_adjacent(a::Integer, b::Integer) = hamming(a, b) == 1

"""
    gray_walk(a, b) -> Vector{Int}

Bit positions (LSB = 0) to flip, one at a time, to walk from basis state `a` to
basis state `b`.  Each intermediate state is Gray-adjacent to the previous one.
Used by multi-controlled-gate constructions that need to visit a set of
computational basis states with one `X`/`CNOT` per hop.
"""
gray_walk(a::Integer, b::Integer) = [p for p in 0:(8 * sizeof(Int) - 1) if ((Int(a) ⊻ Int(b)) >> p) & 1 == 1]

"""
    bits(x, n) -> String

`n`-bit binary string of `x`, most significant bit first.
"""
bits(x::Integer, n::Integer) = string(Int(x); base=2, pad=n)

"""
    parity(a, b) -> Int

Bitwise dot product `a · b` mod 2, i.e. the parity of `a & b`.

The characters of ``Z_2^k`` are `(-1)^parity(a,b)`; the Walsh–Hadamard
transform behind every multiplexed rotation is built from exactly this.
"""
parity(a::Integer, b::Integer) = count_ones(Int(a) & Int(b)) & 1
