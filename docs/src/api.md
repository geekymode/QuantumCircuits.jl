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
