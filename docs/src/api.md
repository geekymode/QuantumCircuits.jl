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
```

## Decompositions

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
