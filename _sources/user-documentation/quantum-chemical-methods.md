# Quantum chemical methods


## Configuration interaction

CI methods find (one or a couple of the lowest) eigenvectors of the Hamiltonian expressed in a certain Fock subspace. GQCP offers the following ONV bases:


| Fock subspace | GQCP ONV basis | Use case |
| --- | --- | --- |
Full spin-resolved | `SpinResolvedONVBasis` | Full CI and Hubbard calculations |
Spin-resolved, selection | `SpinResolvedSelectedONVBasis` | Selected CI calculations |
Seniority-zero | `SeniorityZeroONVBasis` | DOCI, seniority-zero calculations

Let's explore a full CI calculation using the molecular Hamiltonian as an example, but note that you may replace the `SpinResolvedONVBasis` with a different one.

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
# Create the molecular Hamiltonian in the Löwdin spin-orbital basis.
molecule = gqcpy.Molecule.HChain(4, 0.742)
N_P = molecule.numberOfElectronPairs()

spinor_basis = gqcpy.RSpinOrbitalBasis(molecule, "STO-3G")
K = spinor_basis.numberOfSpatialOrbitals()

sq_hamiltonian = gqcpy.SQHamiltonian.Molecular(spinor_basis, molecule)

# Do a dense FCI calculation.
onv_basis = gqcpy.SpinResolvedONVBasis(K, N_P, N_P)

solver = gqcpy.EigenproblemSolver.Dense()
environment = gqcpy.CIEnvironment.Dense(sq_hamiltonian, onv_basis)

qc_structure = gqcpy.CI(onv_basis).optimize(solver, environment)
energy = qc_structure.groundStateEnergy()
linear_expansion = qc_structure.groundStateParameters()
```

<!--C++-->
```C++
// Create the molecular Hamiltonian in the Löwdin spin-orbital basis.
const auto molecule = GQCP::Molecule::HChain(4, 0.742);
const auto N_P = molecule.numberOfElectronPairs();

GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
const auto K = spinor_basis.numberOfSpatialOrbitals();

spinor_basis.lowdinOrthonormalize();
auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);

// Do a dense FCI calculation.
const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};

auto solver = GQCP::EigenproblemSolver::Dense();
auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);

const auto qc_structure = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment);
const auto energy = qc_structure.groundStateEnergy()
const auto linear_expansion = qc_structure.groundStateParameters();
```
<!--END_DOCUSAURUS_CODE_TABS-->