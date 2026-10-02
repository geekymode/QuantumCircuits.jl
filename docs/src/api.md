# API reference

```@meta
CurrentModule = QuantumCircuits
```

```@contents
Pages = ["api.md"]
Depth = 2
```

## Module

```@docs
QuantumCircuits
```

## Gray code

```@docs
gray
ungray
graycode
gray_flip_position
gray_flip_positions
gray_adjacent
gray_walk
hamming
parity
bits
```

## Gates

```@docs
Gate
nqubits
matrix
label
controlled
```

### Constructors

```@docs
Id
X
Y
Z
H
S
Sdg
T
Tdg
RX
RY
RZ
PHASE
CNOT
CZ
SWAP
FANOUT
PARITY
```

## Circuits

```@docs
Circuit
Instruction
Base.push!(::Circuit, ::Gate, ::Integer...)
Base.append!(::Circuit, ::Circuit)
statevector
zero_state
apply!
draw
count_gates
count_cnots
depth
layers
```

## Ancillas

```@docs
add_ancillas!
ancillas
data_qubits
isometry
logical_matrix
leakage
is_clean
implementation_error
```

## Two-qubit (KAK) decomposition

```@docs
KAK
kak
canonical_gate
two_qubit
two_qubit!
```

## GF(2) and classical codes

```@docs
gf2_rref
gf2_rank
gf2_nullspace
LinearCode
dimension
generator_matrix
parity_check_matrix
encode
syndrome
iscodeword
codewords
minimum_distance
dual
puncture
hamming_code
reed_muller
syndrome_decode
rm_local_decode
```

## Pauli operators and stabilizer codes

```@docs
PauliOp
weight
commutes
pauli!
StabilizerCode
stabilizers
logical_operators
is_css
code_distance
lookup_decoder
css_code
repetition_code
five_qubit_code
shor_code
steane_code
quantum_reed_muller
rotated_surface_code
encoding_circuit
syndrome_circuit
```

## Stabilizer simulation

```@docs
Tableau
apply!(::Tableau, ::Gate, ::Integer...)
measure!
expectation
stabilizers(::Tableau)
statevector(::Tableau)
pauli!(::Tableau, ::PauliOp, ::AbstractVector{<:Integer})
prepare_logical_zero
sample_syndrome
logical_error_rate
```

## Real symmetric embeddings

```@docs
realify
hermitian_dilation
symmetric_embedding
embedding_circuit
```

## Quantum Fourier transform

```@docs
qft
qft!
centered_dft!
```

## Three-query synthesis

```@docs
three_query
three_query_error
three_query_synthesis
qho_encoding
chirp_phases
centered_dft
```

## Fan-out, parity and AND

```@docs
fanout!
parity!
and!
```

## Linear algebra

```@docs
fwht
fwht!
walsh_matrix
pauli
pauli_strings
pauli_decompose
pauli_recompose
embed
kron_n
is_unitary
gate_fidelity
global_phase_between
schmidt_values
entanglement_entropy
```

## Matrix decompositions

```@docs
zyz
decompose_1q
decompose_1q!
TwoLevel
two_level_decompose
two_level!
synthesize_unitary
demultiplex
multiplexed_1q
multiplexed_1q!
```

## Phase polynomials

```@docs
PhasePolynomial
phase_polynomial
phases
support
nterms
synthesize
phase_gadget
phase_gadget!
pauli_rotation!
trotter_step!
cancel_adjacent_cnots!
```

## Gray-code decompositions

```@docs
multiplex_angles
multiplex_matrix
multiplexed_rotation!
multiplexed_ry
multiplexed_rz
diagonal
diagonal!
prepare_state
```
