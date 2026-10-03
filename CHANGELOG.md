# Changelog

## Unreleased

### Error correction
- GF(2) linear algebra: `gf2_rref`, `gf2_rank`, `gf2_nullspace`.
- Classical codes: `LinearCode`, `hamming_code`, `reed_muller`, `encode`,
  `syndrome`, `syndrome_decode`, `minimum_distance` (Gray-code walk over all
  codewords), `dual`, `puncture`, and Reed's local decoder `rm_local_decode`.
- Pauli operators (`PauliOp`) and stabilizer codes (`StabilizerCode`,
  `css_code`), with logical operators computed when not given,
  `code_distance` and `lookup_decoder`.
- Catalogue: repetition, five-qubit, Shor, Steane, quantum Reed–Muller
  `[[2ᵐ-1, 1, 3]]`, rotated surface codes.
- Circuits: `encoding_circuit` (CSS, Hadamards and CNOTs) and
  `syndrome_circuit` (one ancilla per stabilizer), checked exactly.
- New docs page: *Error correction*, ending with transversal `T` on the
  `[[15, 1, 3]]` code.

### Stabilizer simulation
- `Tableau` — Aaronson–Gottesman Clifford simulation on packed bits, with
  `apply!` (any Clifford gate or `Circuit`), `measure!` (qubit or any Pauli),
  `expectation`, `stabilizers` and `statevector`. Circuits stay unitary;
  measurement lives here.
- `prepare_logical_zero` (any stabilizer code, by measuring stabilizers),
  `sample_syndrome` (runs the syndrome circuit and measures), and
  `logical_error_rate` (code-capacity Monte Carlo).
- Stabilizer codes can record a known `distance`; the catalogue does.

### Noisy syndrome extraction
- `NoiseModel`, `circuit_noise`, `phenomenological_noise`.
- `memory_experiment` — repeated noisy rounds of the syndrome circuit on the
  tableau, detectors as round-to-round changes, decoded by matching on the
  space-time graph. Phenomenological crossover near 3% (matching the known
  threshold); circuit-level about 0.3%, below tuned decoders, for documented
  reasons (no diagonal edges, uniform weights, sequential checks).
- `min_weight_perfect_matching` — Edmonds' blossom algorithm, `O(n³)`,
  checked against brute force; the matching decoder now uses it past a few
  defects instead of a greedy fallback, so decoding is always exact.

### Matching decoder
- `matching_decoder` / `decode` — minimum-weight matching for surface-code-like
  CSS codes, always exact: dynamic programming over defect subsets for a few
  defects, Edmonds' blossom algorithm (`min_weight_perfect_matching`,
  `O(n³)`) beyond. The correction always explains the syndrome; every error of
  weight `≤ ⌊(d-1)/2⌋` is corrected. Code-capacity crossover near 15%, against
  below 1% for the lookup table. `logical_error_rate` accepts either decoder.

### Two-qubit KAK decomposition
- `kak(U)` — `U = e^{iφ}(A₁⊗A₂)·exp(i(a XX + b YY + c ZZ))·(B₁⊗B₂)`, angles
  reduced to `(-π/4, π/4]`; `KAK` result type and `canonical_gate(a, b, c)`.
- `two_qubit!` / `two_qubit` — any two-qubit unitary in the fewest CNOTs its
  class allows: 0 (local), 1 (CNOT class), 2 (one canonical angle zero), 3.

### Quantum Shannon decomposition
- `qsd` stops at two-qubit KAK blocks: `(9/16)·4ⁿ - (3/2)·2ⁿ` CNOTs — 24 at
  `n = 3` (was 36), 120 at `n = 4` (was 168). Structured inputs now come in
  lower too (the identity on 3 qubits: 12). `kak=false` keeps the old
  recursion; `qsd_cnot_count` and `qsdfigure` take the same keyword.

## v0.2.0 — 2026-10-02

Depth, ancillas and parallel synthesis. Everything before this release
minimised CNOT *count* on exactly `n` wires; this release adds the other axis —
circuit **depth**, and what extra wires buy — built around Nehoran and Yuen,
[*All Unitaries Have Constant Depth Quantum Circuits*](https://arxiv.org/abs/2609.40351).
New docs page: *Depth and width*.

No breaking changes: existing circuits, defaults and printed output are
unchanged.

### Depth, fan-out and parity
- `depth(c)`, `depth(c, name)` (e.g. CNOT depth) and `layers(c)` — ASAP
  scheduling.
- `FANOUT(r)` and `PARITY(r)` gates, one step each (the unbounded fan-out
  model). Stored densely, so capped at 11 wires.
- `fanout!` and `parity!` with `style = :ladder | :tree | :gate`. The tree is
  exact on every input, not only on targets in `|0⟩`.
- `phase_gadget!` and `pauli_rotation!` take `style`; `:tree` keeps
  `2(k-1)` CNOTs at depth `2⌈log₂ k⌉ + 1` (8 qubits: depth 15 → 7).

### Ancillas
- `Circuit(n; ancillas)`, `add_ancillas!`, `ancillas`, `data_qubits`.
- `isometry` and `logical_matrix` — the circuit seen from the data wires.
- `leakage`, `is_clean`, `implementation_error` — `‖C(I⊗|0⟩) − U⊗|0⟩‖`.
- `and!` — `CᵏX` as a Toffoli tree through `k - 2` clean ancillas, depth
  `2⌈log₂ k⌉ - 1` (`k = 8`: depth 5, against 503 ancilla-free).
- `statevector` accepts a data-only input state; drawings label ancilla
  wires `a`.

### Real symmetric embeddings
- `realify`, `hermitian_dilation`, `symmetric_embedding` — any unitary as a
  real, symmetric, traceless involution of four times the dimension.
- `embedding_circuit` wraps any circuit for that involution into one that
  implements the original unitary with clean ancillas.

### Three-query synthesis
- `three_query`, `three_query_error`, `three_query_synthesis` — a reference
  implementation of the identity `S ≈ i Ê† Q̂ F̂ Q̂ F̂ Q̂ Ê`, with its components
  `qho_encoding`, `chirp_phases`, `centered_dft`.
- `threequeryfigure` — error against grid size (Makie extension).

### Quantum Fourier transform
- `qft!` and `qft`, with `inverse`, `swaps`, `cutoff` (approximate QFT) and
  `lower` (CNOTs only).
- `centered_dft!` — the three-query Fourier step as gates, exact.

### Fixes
- The README's plot-test command could not work (`Pkg.test` sandboxes the
  Makie backend away); it now runs `test/runtests.jl` in the docs environment.

## v0.1.0 — 2026-08-19

First release: Gray-code circuit decompositions.

- Gray code utilities; gates and a dense statevector simulator with exact
  global phase.
- Uniformly controlled rotations in `2ᵏ` CNOTs, diagonal unitaries in
  `2ⁿ - 2`, state preparation in `2ⁿ⁺¹ - 4`.
- Linear-algebra toolkit (Walsh–Hadamard and Pauli transforms), two-level and
  ZYZ decompositions, phase polynomials and parity networks.
- Fifteen worked applications of Gray coding.
- Cosine–sine and quantum Shannon decomposition, `(3/4)·4ⁿ - (3/2)·2ⁿ` CNOTs.
- Makie-based illustrations as a package extension.
