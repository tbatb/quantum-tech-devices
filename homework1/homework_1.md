# Homework 1: Quantum Wells and Perturbation Theory - Notebook Explanation

## Notebook Structure Overview
SciPy Docs: https://docs.scipy.org/doc/scipy/reference/
## Cell 1: imports - breakdown
- `import numpy as np`: Array operations, np.array converts lists to arrays.
- `scipy constant as const`: undefined numbers where truncation errors appear are brought into float form
- `scipy.sparse import diags`: physical and math units i.e pi, golden_ration etc
- `scipy.sparse.linalg import eigsh`: sparse linear algebra

---

## Cell 2: Markdown - Helper Functions Description
- **Type**: Markdown.
- **Content**:
  - Describes the `build_potential` function (detailed below).
  - Parameters: `layers` as list of tuples `(thickness, band_edge)` (thickness in nm, band_edge in eV).
  - Returns: Two arrays (`z` grid in space, potential `V` structure).
  - Example: `[(5, 1.0), (2, 0.3), (5, 1.0)]` for a barrier-well-barrier.
  - Tip: Suggests plotting a random structure for understanding.
- **Purpose**: Explains the upcoming code cell. References no external docs directly, but implies usage with NumPy arrays (see NumPy docs for `np.array`).

---

## Cell 3: Code - Imports and `build_potential` Function
- **Type**: Code (Python).
- **Execution Count**: `null` (not executed in the notebook).
- **Source Code**:
  ```python:disable-run
  import numpy as np
  import matplotlib.pyplot as plt
  import scipy.constants as const

  def build_potential(layers, dz_nm=0.05):
      """
      Build heterostructure potential from layers.
      layers: list of (thickness [nm], band edge [eV])
      dz_nm: grid spacing [nm]
      returns: z [m], V [J]
      """
      z = []
      V = []
      pos = 0.0
      for thickness_nm, band_edge_eV in layers:
          n_points = int(thickness_nm/dz_nm)
          for _ in range(n_points):
              z.append(pos * 1e-9)
              V.append(band_edge_eV * const.eV)
              pos += dz_nm
      return np.array(z), np.array(V)
  ```

### Breakdown
- **Imports**:
  - `import numpy as np`: Aliases NumPy for array operations. `np.array` converts lists to arrays (NumPy docs: [numpy.org/doc/stable/reference/generated/numpy.array.html](https://numpy.org/doc/stable/reference/generated/numpy.array.html)).
  - `import matplotlib.pyplot as plt`: Imports plotting submodule (not used here, but for later visualization; Matplotlib docs: [matplotlib.org/stable/api/pyplot_api.html](https://matplotlib.org/stable/api/pyplot_api.html)).
  - `import scipy.constants as const`: Aliases SciPy constants. `const.eV` is electronvolt in Joules (1.60217662e-19 J); `const.hbar` (Planck's constant / 2π) and `const.m_e` (electron mass) are used later (SciPy docs: [docs.scipy.org/doc/scipy/reference/constants.html#scipy.constants.physical_constants](https://docs.scipy.org/doc/scipy/reference/constants.html#scipy.constants.physical_constants)).

- **Function: `build_potential(layers, dz_nm=0.05)`**:
  - **Docstring**: Explains purpose, parameters, and returns. Units: input nm/eV, output m/J.
  - **Parameters**:
    - `layers`: List of tuples, e.g., `[(5.0, 1.0)]` (thickness in nm, band edge in eV). Represents layered structure (e.g., barriers/wells).
    - `dz_nm=0.05`: Default grid spacing in nm (scalar float).
  - **Local Variables**:
    - `z = []`: Empty list to store position grid (will become NumPy array in meters).
    - `V = []`: Empty list to store potential values (will become NumPy array in Joules).
    - `pos = 0.0`: Starting position (float, in nm; accumulates along z-direction).
    - `thickness_nm, band_edge_eV`: Unpacked from each tuple in `layers` loop (floats).
    - `n_points = int(thickness_nm/dz_nm)`: Integer number of grid points per layer (truncates via `int`; potential loss of precision for non-integer divisions).
    - `_`: Loop variable (unused, Python convention for dummy iterators).
  - **Logic**:
    - Outer loop: Iterates over each layer.
    - Inner loop: For each grid point in the layer, appends `pos * 1e-9` (nm to m conversion) to `z`, `band_edge_eV * const.eV` (eV to J) to `V`, then increments `pos` by `dz_nm`.
    - Creates a step-function potential on an equidistant grid (but note: total grid may not align perfectly across different `layers` calls due to `int` truncation).
  - **Returns**: Tuple `(np.array(z), np.array(V))` – 1D arrays of same length, suitable for finite-difference methods.
  - **Potential Issues**: Grid spacing is fixed per call; for multi-structure comparisons (e.g., Task 3), grids may misalign unless `dz_nm` is chosen carefully.
  - **Documentation Reference**: Relies on SciPy `const.eV` for unit conversion (see above). No external function calls beyond built-ins (`int`, `append`, `range` – Python stdlib docs: [docs.python.org/3/library/functions.html](https://docs.python.org/3/library/functions.html)).

---

## Cell 4: Markdown - Task 1: Finite Quantum Well
- **Type**: Markdown.
- **Content**:
  - Subtasks (a), (b), (c):
    - (a): Use `build_potential` for a 3 nm well, 1 eV offset.
    - (b): Implement `solve_schrodinger(z, V, n_eigs)`: Finite-difference Hamiltonian, electron mass `m_e`, return `n_eigs` energies (eV) and normalized wavefunctions (∫|ψ|² dx = 1). Boundary: ψ=0 at edges. Tip: Choose barrier thickness to avoid artifacts.
    - (c): Plot |ψ|² for first 3 states, compare energies to analytical (infinite well?). Suggest scaling ψ for plotting; use labels/legend.
- **Purpose**: Assignment instructions. References analytical finite well solutions implicitly (standard quantum mechanics texts, e.g., Griffiths' *Introduction to Quantum Mechanics* for transcendental equations).

---

## Cell 5: Code - Template for `solve_schrodinger`
- **Type**: Code.
- **Execution Count**: `null`.
- **Source Code**:
  ```python
  def solve_schrodinger(x, V, n_eigs=5):
      #Do some calculations
      
      return # put return values here
  ```

### Breakdown
- **Function: `solve_schrodinger(x, V, n_eigs=5)`**:
  - **Parameters**:
    - `x`: Likely position grid (array, e.g., from `build_potential`; note: docstring in Task says `z`, possible typo).
    - `V`: Potential array (in J, matching `x` length).
    - `n_eigs=5`: Number of lowest eigenvalues/wavefunctions to compute (int).
  - **Body**: Placeholder comment `#Do some calculations` and incomplete `return`.
  - **Expected Behavior** (per instructions): Build tridiagonal Hamiltonian matrix via finite differences (second derivative approximation). Use `const.m_e` for kinetic term. Solve via eigendecomposition (likely `scipy.linalg.eigh` for Hermitian matrix). Convert energies to eV (`/ const.e`). Normalize ψ via trapezoidal integration (∫|ψ|² dx ≈ sum(|ψ|² Δx) = 1). Boundary: Set ψ[0]=ψ[-1]=0, solve on interior points.
  - **Documentation Reference**: Would use SciPy linear algebra (e.g., `scipy.linalg.eigh` for symmetric matrices; docs: [docs.scipy.org/doc/scipy/reference/linalg.html](https://docs.scipy.org/doc/scipy/reference/linalg.html)). NumPy for arrays. Physical: `const.hbar`, `const.m_e`, `const.e` (SciPy constants).

---

## Cell 6: Markdown - Task 2: Coupled Quantum Wells
- **Type**: Markdown.
- **Content**:
  - Subtasks (a), (b), (c):
    - (a): Build/solve double well (same well widths); experiment to see subband splitting.
    - (b): Plot lowest 2 ψ; explain symmetric/antisymmetric.
    - (c): Plot splitting (ΔE between lowest 2 states) vs. barrier width.
- **Purpose**: Builds on Task 1 for tunneling effects. References quantum tunneling concepts (e.g., symmetric/antisymmetric via parity in coupled systems).

---

## Cell 7: Code - Empty for Task 2
- **Type**: Code.
- **Execution Count**: `null`.
- **Source**: Empty (`[]`).
- **Purpose**: Placeholder for student code (e.g., calling `build_potential`, `solve_schrodinger`, plotting with `plt`).

---

## Cell 8: Markdown - Task 3: Coupled States Using Perturbation Theory
- **Type**: Markdown.
- **Content**:
  - **Goal**: Use unperturbed (uncoupled) wells' ground states to approximate coupled states via 2x2 Hamiltonian in subspace: Ψ_coupled = c_L ψ_L + c_R ψ_R. Solve H_per c = E S c (generalized eigenproblem, S=overlap matrix).
  - Given `layers_left`, `layers_right`, `layers_coupled` (tuples: nm, eV).
  - **Steps**: 1-9 detailed (build potentials, solve uncoupled/coupled, compute H_per/S, solve eigensystem, reconstruct ψ, plot comparisons, extend Task 2(c) with perturbation).
  - References LaTeX equations for wavefunction and eigenproblem.
- **Purpose**: Introduces non-orthogonal basis perturbation theory. References quantum perturbation theory (e.g., time-independent, degenerate case for near-degenerate states).

---

## Cell 9: Code - Layer Definitions for Task 3
- **Type**: Code (ID: "9da3be02").
- **Execution Count**: `null`.
- **Source Code**:
  ```python
  # Define layer structures
  layers_left    = [(5.0, 0.0), (2.0, -0.5), (7.3, 0.0)]
  layers_right   = [(7.3, 0.0), (2.0, -0.5), (5.0, 0.0)]
  layers_coupled = [(5.0, 0.0), (2.0, -0.5), (0.3, 0.0), (2.0, -0.5), (5.0, 0.0)]

  dz_nm = 0.01
  ```

### Breakdown
- **Variables** (all global in this cell):
  - `layers_left`: List of 3 tuples – left barrier (5 nm, 0 eV), well (2 nm, -0.5 eV), right barrier (7.3 nm, 0 eV). Asymmetric to isolate left well.
  - `layers_right`: List of 3 tuples – left barrier (7.3 nm, 0 eV), well (2 nm, -0.5 eV), right barrier (5 nm, 0 eV). Mirror of left for right well isolation.
  - `layers_coupled`: List of 5 tuples – outer left barrier (5 nm, 0), left well (2 nm, -0.5), thin central barrier (0.3 nm, 0), right well (2 nm, -0.5), outer right barrier (5 nm, 0).
  - `dz_nm = 0.01`: Float, finer grid spacing (nm) for Task 3 to align grids across structures (addresses step 1 note on equidistant grids).
- **Comments**: `# Define layer structures` – descriptive header.
- **Purpose**: Predefines inputs for perturbation calculation. Note: Barrier thicknesses (5+7.3=12.3 nm total per side) ensure ψ≈0 at edges; central 0.3 nm for strong coupling.
- **Documentation Reference**: Tuples are Python built-ins (immutable lists; stdlib docs: [docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences](https://docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences)). Used as input to `build_potential`.

---

This document covers the entire notebook structure. The notebook emphasizes building from simple (single well) to advanced (perturbation) concepts, with code templates for implementation. For execution, cells would run sequentially in a Jupyter environment, building on prior definitions.
```