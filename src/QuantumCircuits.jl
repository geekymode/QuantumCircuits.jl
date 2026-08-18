"""
    QuantumCircuits

A small, dependency-free Julia package for building quantum circuits and
decomposing structured unitaries into elementary gates.

The organising idea of the first release is **Gray coding**: ordering the
computational basis so that consecutive states differ in exactly one bit turns
the CNOT ladders in a decomposition into a single CNOT per step.  See
`docs/graycode.md`.

Conventions
-----------
* Qubit `1` is the most significant bit of a basis index (`|q1 q2 … qn⟩`).
* `controls[1]` is the most significant bit of a multiplexor's branch index.
* Rotations are `R(θ) = exp(-i θ P / 2)`.
"""
module QuantumCircuits

using LinearAlgebra
using Printf: @sprintf

# Gray code
export gray, ungray, graycode, gray_flip_position, gray_flip_positions,
       gray_adjacent, gray_walk, hamming, bits, parity

# gates
export Gate, Id, X, Y, Z, H, S, Sdg, T, Tdg, RX, RY, RZ, PHASE,
       CNOT, CZ, SWAP, controlled, label

# circuits
export Circuit, Instruction, nqubits, matrix, statevector, zero_state,
       apply!, draw, count_gates, count_cnots

# decompositions
export multiplex_angles, multiplex_matrix, multiplexed_rotation!,
       multiplexed_ry, multiplexed_rz, diagonal, diagonal!, prepare_state

# linear algebra
export is_unitary, gate_fidelity, global_phase_between, embed, kron_n,
       fwht, fwht!, walsh_matrix, pauli, pauli_strings, pauli_decompose,
       pauli_recompose, schmidt_values, entanglement_entropy

# matrix decompositions
export zyz, decompose_1q, decompose_1q!, TwoLevel, two_level_decompose,
       two_level!, synthesize_unitary, demultiplex, multiplexed_1q,
       multiplexed_1q!

# applications of Gray coding
export matrix_root, multicontrolled, multicontrolled!, gray_encoder, gray_decoder,
       increment, increment!, gray_increment, select, select!, support_mask,
       gray_order, truncate_terms

# plotting (implemented by the Makie extension)
export circuitfigure, circuitplot!, matrixfigure, graycodefigure, costfigure

# phase polynomials
export phase_gadget, phase_gadget!, pauli_rotation!, trotter_step!,
       PhasePolynomial, phase_polynomial, phases, support, nterms,
       synthesize, cancel_adjacent_cnots!

include("graycode.jl")
include("gates.jl")
include("circuit.jl")
include("decompose.jl")
include("mathkit.jl")
include("matrixdecomp.jl")
include("phasepoly.jl")
include("applications.jl")
include("plots.jl")

end # module
