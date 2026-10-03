# ---------------------------------------------------------------------------
# Error correction on the tableau simulator
#
# The tableau handles codes far beyond the dense simulator — a distance-7
# surface code with its syndrome ancillas is 97 qubits — and brings real
# measurement: logical states are prepared by measuring stabilizers, and
# syndromes are read by measuring ancillas.
# ---------------------------------------------------------------------------

_pad(P::PauliOp, N::Int) = PauliOp(vcat(P.x, falses(N - length(P))), vcat(P.z, falses(N - length(P))), P.phase)

# Pure errors: F_i anticommutes with stabilizer i only, and commutes with every
# other stabilizer and every logical operator.
function _pure_errors(code::StabilizerCode)
    n = length(code)
    S = code.stabilizers
    L = vcat(code.logical_x, code.logical_z)
    # ⟨F, P⟩ = x_F·z_P + z_F·x_P: one equation per row [z_P | x_P]
    A = _symp_matrix(vcat(S, L))
    A = hcat(A[:, n+1:end], A[:, 1:n])
    map(eachindex(S)) do i
        b = BitVector(j == i for j in 1:size(A, 1))
        v = _gf2_solve(A, b)
        v === nothing && error("no pure error for stabilizer $i")
        _pauli_from_symp(v)
    end
end

"""
    prepare_logical_zero(code; ancillas=0, rng) -> Tableau

A tableau holding the encoded `|0̄…0̄⟩` of any stabilizer code — CSS or not —
prepared the way hardware would: start in `|0…0⟩`, measure every stabilizer,
and fix each `-1` outcome with a Pauli that flips that stabilizer alone; then
measure each `Z̄` and fix with `X̄`.  `ancillas` extra qubits in `|0⟩` are
appended after the data.
"""
function prepare_logical_zero(code::StabilizerCode; ancillas::Integer=0, rng=Random.default_rng())
    n = length(code)
    N = n + ancillas
    t = Tableau(N)
    F = _pure_errors(code)
    for (s, f) in zip(code.stabilizers, F)
        measure!(t, _pad(s, N); rng=rng) && pauli!(t, _pad(f, N))
    end
    for (z, x) in zip(code.logical_z, code.logical_x)
        measure!(t, _pad(z, N); rng=rng) && pauli!(t, _pad(x, N))
    end
    t
end

"""
    sample_syndrome(code, E; rng) -> BitVector

Run syndrome extraction for real on the tableau simulator: prepare `|0̄⟩`,
apply the error `E`, run [`syndrome_circuit`](@ref), and measure the
ancillas.  For a Pauli error the measured bits always equal
[`syndrome`](@ref)`(code, E)` — the point being that the circuit reports it.
"""
function sample_syndrome(code::StabilizerCode, E::PauliOp; rng=Random.default_rng())
    n, r = length(code), length(code.stabilizers)
    t = prepare_logical_zero(code; ancillas=r, rng=rng)
    pauli!(t, E, 1:n)
    apply!(t, syndrome_circuit(code))
    BitVector([measure!(t, n + i; rng=rng) for i in 1:r])
end

# iid depolarizing: each qubit gets X, Y or Z with probability p/3 each.
function _depolarizing(n::Int, p::Real, rng)
    x = falses(n); z = falses(n)
    for q in 1:n
        if rand(rng) < p
            k = rand(rng, 1:3)
            x[q] = k != 3; z[q] = k != 1
        end
    end
    PauliOp(x, z)
end

"""
    logical_error_rate(code, p; shots=10_000, decoder=lookup_decoder(code; maxweight), maxweight, rng)
        -> (rate, failures)

`decoder` is anything [`decode`](@ref) accepts: a lookup table, or a
[`MatchingDecoder`](@ref) for surface codes.

Monte Carlo estimate of the logical error rate under *code-capacity* noise:
each qubit independently suffers `X`, `Y` or `Z` with probability `p/3`
each; the syndrome is read perfectly, the lookup `decoder` proposes a
correction, and a shot fails if the residual error is a logical operator — or
if the syndrome is not in the table.

`maxweight` defaults to `⌊(d-1)/2⌋`, with `d` the code's recorded distance
(searched for when none is recorded): the decoder is bounded-distance, so it
gives up on heavier errors that a full decoder (minimum-weight matching, say)
might still fix.  The rate still falls as `p^(⌊(d-1)/2⌋+1)` for small `p`,
which is the point.
"""
function logical_error_rate(code::StabilizerCode, p::Real; shots::Integer=10_000,
                            maxweight::Integer=((code.distance === nothing ? code_distance(code) : code.distance) - 1) ÷ 2,
                            decoder=lookup_decoder(code; maxweight=maxweight),
                            rng=Random.default_rng())
    0 <= p <= 1 || throw(ArgumentError("p must lie in [0, 1]"))
    n = length(code)
    L = vcat(code.logical_x, code.logical_z)
    failures = 0
    for _ in 1:shots
        E = _depolarizing(n, p, rng)
        C = decode(decoder, syndrome(code, E))
        if C === nothing
            failures += 1
            continue
        end
        R = C * E
        any(l -> !commutes(R, l), L) && (failures += 1)
    end
    failures / shots, failures
end
