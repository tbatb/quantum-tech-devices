# Homework 1: Quantum Wells and Perturbation Theory - Notebook Explanation

## Notebook Structure Overview
SciPy Docs: https://docs.scipy.org/doc/scipy/reference/
## Cell: imports - breakdown
- `import numpy as np`: Array operations, np.array converts lists to arrays.
- `scipy constant as const`: undefined numbers where truncation errors appear are brought into float form
- `scipy.sparse import diags`: physical and math units i.e pi, golden_ration etc
- `scipy.sparse.linalg import eigsh`: sparse linear algebra

---

## Cell(s) Task1: Markdown - Helper Functions Description
- **Type**: Markdown.
- **Content**:
  - Describes the `build_potential` function (detailed below).
  - Parameters: `layers` as list of tuples `(thickness, band_edge)` (thickness in nm, band_edge in eV).
  - Returns: Two arrays (`z` grid in space, potential `V` structure).
  - Example: `[(5, 1.0), (2, 0.3), (5, 1.0)]` for a barrier-well-barrier.
  - Tip: Suggests plotting a random structure for understanding.
- **Purpose**: Explains the upcoming code cell. References no external docs directly, but implies usage with NumPy arrays (see NumPy docs for `np.array`).

Here is a detailed documentation for the code and concepts used in Task 1.

## Code Documentation for Task 1

This document provides a comprehensive explanation of the variables, functions, and scientific principles behind the numerical solution of a finite quantum well as implemented in Task 1.

---

### Library Imports

The script begins by importing necessary libraries.

* `numpy`: A fundamental package for scientific computing in Python. It's used here for creating and manipulating numerical arrays (`np.array`, `np.ones`, etc.), performing mathematical operations, and sorting.
* `matplotlib.pyplot`: A plotting library used to create visualizations of the potential well, energy levels, and wavefunctions.
* `scipy.constants`: A sub-module of the SciPy library that provides a comprehensive set of physical constants (e.g., `hbar`, electron mass `m_e`, `eV`). This ensures accuracy and avoids manual entry of these values.
* `scipy.sparse` and `scipy.sparse.linalg`: These modules are used to handle large, sparse matrices efficiently. A sparse matrix is one that is mostly filled with zeros. The Hamiltonian matrix we construct is tridiagonal, making it very sparse. Using sparse matrix tools saves significant memory and computational time compared to using dense NumPy arrays for the same problem.

---

### Function: `build_potential`

This function constructs the one-dimensional potential energy profile $V(z)$.

#### **Purpose**
To translate a high-level description of a semiconductor heterostructure—a series of layers with specified thicknesses and potential energies—into a discretized numerical array that can be used in the Schrödinger solver.

#### **Variables**
* **`layers`** (input): A Python `list` of `tuples`. Each tuple `(thickness_nm, band_edge_eV)` defines a single material layer.
    * `thickness_nm` (`float`): The physical thickness of the layer in nanometers.
    * `band_edge_eV` (`float`): The potential energy of that layer in electron-volts (eV).
* **`dz_nm`** (input): A `float` representing the grid spacing in nanometers. This determines the resolution of our simulation. A smaller `dz_nm` leads to a more accurate result at the cost of increased computation time.
* **`z`** (local list, returned as `np.array`): Stores the spatial coordinates of each grid point in **meters**.
* **`V`** (local list, returned as `np.array`): Stores the potential energy values at each grid point in **Joules**. The conversion from eV to Joules is essential for consistency with the SI units used in `scipy.constants`.
* **`pos`** (local `float`): A counter that keeps track of the current position along the z-axis as the structure is being built.
* **`n_points`** (local `int`): The number of discrete grid points that fit within a given layer's thickness.

#### **Mechanism**
The function iterates through each layer defined in the `layers` list. For each layer, it calculates how many grid points (`n_points`) are needed to represent it based on its thickness and the grid spacing `dz_nm`. It then loops `n_points` times, each time appending the current spatial position to the `z` list and the layer's potential energy to the `V` list. The position `pos` is incremented by `dz_nm` in each step. Finally, the Python lists are converted to more efficient NumPy arrays before being returned.

---

### Function: `solve_schrodinger`

This is the core of the exercise. It solves the 1D time-independent Schrödinger equation for a given potential `V(z)`.

[cite_start]$$-\frac{\hbar^2}{2m^*} \frac{\partial^2}{\partial z^2}\psi(z) + V(z)\psi(z) = E\psi(z) \quad [cite: 992]$$

This is an eigenvalue equation of the form $\hat{H}\psi = E\psi$, where $\hat{H}$ is the Hamiltonian operator. The goal is to find the allowed energies $E$ (eigenvalues) and their corresponding wavefunctions $\psi$ (eigenvectors).

#### **Mechanism: Finite Difference Method**
To solve this equation numerically, we discretize both space and the operators. The continuous function $\psi(z)$ becomes a vector `ψ[i]`, where `i` is the index of a grid point. [cite_start]The second derivative is approximated using a central finite difference formula[cite: 1017]:

$$\frac{\partial^2 \psi}{\partial z^2} \approx \frac{\psi(z+dz) - 2\psi(z) + \psi(z-dz)}{dz^2} \rightarrow \frac{\psi[i+1] - 2\psi[i] + \psi[i-1]}{dz^2}$$

[cite_start]Substituting this into the Schrödinger equation yields a system of linear equations that can be written as a single matrix eigenvalue problem[cite: 1059, 1088].

#### **Variables**
* **`z`, `V`** (input): The NumPy arrays for position (m) and potential (J) generated by `build_potential`.
* **`n_eigs`** (input): An `int` specifying how many of the lowest energy solutions to find.
* **`dz`**: The grid spacing in meters, calculated from the input `z` array.
* **`N`**: The total number of grid points in the simulation domain.
* **`T_const`**: A `float` constant that combines all the physical constants for the kinetic energy term: $-\frac{\hbar^2}{2m_e}$.
* **`T`** (`scipy.sparse.csc_matrix`): The **kinetic energy operator matrix**. It's a tridiagonal matrix where the main diagonal is `-2`, the diagonals above and below are `1`, and all other elements are zero. This structure comes directly from the finite difference formula for the second derivative. The entire matrix is then scaled by `T_const / dz**2`.
* **`U`** (`scipy.sparse.csc_matrix`): The **potential energy operator matrix**. In this basis, it's a simple diagonal matrix where the diagonal entries are the values from the input potential array `V`.
* **`H`** (`scipy.sparse.csc_matrix`): The **Hamiltonian matrix**, formed by adding the kinetic and potential matrices: `H = T + U`.
* **`eigenvalues`, `eigenvectors`**: The raw output from the `eigs` solver. Eigenvalues correspond to energies (in Joules), and eigenvectors are the (unnormalized) wavefunctions.
* **`energies`** (`np.array`): The final, sorted array of eigenenergies, converted to **eV**.
* **`wavefunctions`** (`np.array`): The final, sorted, and **normalized** wavefunctions.

#### **Execution Flow**
1.  **Construct Hamiltonian**: It builds the sparse matrices `T` and `U` as described above and sums them to get `H`.
2.  **Solve Eigenproblem**: It calls `scipy.sparse.linalg.eigs(H, k=n_eigs, which='SM')`. The argument `which='SM'` is crucial; it tells the solver to find the `k` eigenvalues with the **S**mallest **M**agnitude, which correspond to the lowest energy bound states.
3.  **Sort Results**: The solver doesn't guarantee the order of the results, so the code explicitly sorts the eigenvalues in ascending order and rearranges the eigenvectors to match.
4.  [cite_start]**Normalize Wavefunctions**: A physical wavefunction must be normalized such that the total probability of finding the particle is 1. The condition is $\int_{-\infty}^{\infty} |\psi(z)|^2 dz = 1$[cite: 800]. The discrete version is $\sum_i |\psi_i|^2 dz = 1$. The code enforces this by calculating the normalization constant for each wavefunction (`np.sqrt(np.sum(psi**2) * dz)`) and dividing the wavefunction by it.
5.  **Return Values**: The function returns the clean, sorted, and correctly unit-converted energies and wavefunctions.

---

### Main Script for Task 1

This part of the code uses the functions above to perform the specific simulation requested.

#### **(a) Build Potential**
* **`well_width_nm`, `barrier_width_nm`, `barrier_height_eV`**: These variables define the physical parameters of our specific quantum well: a 3 nm well surrounded by 10 nm barriers that are 1 eV high.
* **`layers_task1`**: This list organizes the parameters into the format required by the `build_potential` function.
* **`z, V = build_potential(...)`**: This line calls the function to generate the discrete potential profile.
* **`plt.plot(...)`**: This section plots the potential $V(z)$ to visually confirm that the structure was built correctly.

#### **(b) & (c) Solve and Analyze**
* **`energies, wavefunctions = solve_schrodinger(...)`**: This calls the solver to calculate the first 3 energy states (`n_eigs=3`).
* **Plotting Results**:
    * The potential is plotted again as a black line (`'k-'`).
    * The code then loops through the first three solutions.
    * For each solution, it plots the **probability density** $|\psi|^2$, not the wavefunction $\psi$ itself, as this density represents the probability of finding the particle at a given position.
    * The probability density is scaled by a factor `scaling` and shifted vertically by its corresponding energy `energies[i]` so that it appears centered on its energy level in the plot, which is a common and clear way to visualize quantum states.
    * A dashed horizontal line (`axhline`) is drawn at each energy level to guide the eye.
* **Analytical Comparison**:
    * This section serves as a sanity check for the numerical results.
    * [cite_start]It calculates the energy levels for a related but simpler problem: the **infinite potential well**, for which an exact analytical formula exists: $E_n = \frac{\hbar^2\pi^2 n^2}{2mL^2}$[cite: 803].
    * By comparing the numerical results of our *finite* well to the analytical results of an *infinite* well, we can check if they are reasonable. As observed, the finite well energies are lower, which is physically correct because the wavefunctions can penetrate the barriers, effectively widening the box.