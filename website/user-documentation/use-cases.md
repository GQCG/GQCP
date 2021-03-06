# Common use cases

In this section, we'll provide a quick reference to some common use cases. Be sure to check out the [gqcpy example notebooks](https://github.com/GQCG/GQCP/tree/develop/gqcpy/examples) as well!

## Hartree-Fock calculations

### RHF

In this use case, we will calculate the restricted HF energy for H2 in an STO-3G basis set. 

#### Molecular setup

The first step consists of generating the [`Molecule`](#molecules.md) object from which we can obtain properties. The simplest way to do so is to provide an `xyz` file. 

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
molecule = gqcpy.Molecule.ReadXYZ("h2.xyz")
N = molecule.numberOfElectrons()
```

<!--C++-->
```C++
const auto molecule = GQCP::Molecule::ReadXYZ("data/h2.xyz");
const auto N = molecule.numberOfElectrons();
```
<!--END_DOCUSAURUS_CODE_TABS-->

Before we go over to the calculation of integrals, an [orbital basis](#orbital_bases.md) must be constructed. The RHF method uses a restricted spin-orbital basis. We can immediately use this object to find the overlap integrals.

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
spinor_basis = gqcpy.RSpinOrbitalBasis(molecule, "STO-3G")
S = spinor_basis.quantizeOverlapOperator().parameters()
```

<!--C++-->
```C++
const auto spinor_basis = GQCP::RSpinOrbitalBasis::ReadXYZ(molecule, "STO-3G";
const auto S = spinor_basis.overlap().parameters();
```
<!--END_DOCUSAURUS_CODE_TABS-->

This spinor basis is used to construct a Hamiltonian operator in second quantization (SQ). 

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
sq_hamiltonian = gqcpy.SQHamiltonian.Molecular(spinor_basis, molecule)  # in an AO basis
```

<!--C++-->
```C++
const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, h2);  // in an AO basis
```
<!--END_DOCUSAURUS_CODE_TABS-->

#### Plain RHF SCF

To solve the RHF SCF equations, we will need to set up a `Solver` and its associated `Environment`. The environment will provide all necessary intermediates for the solver to use. 

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
rhf_environment = gqcpy.RHFSCFEnvironment.WithCoreGuess(N, sq_hamiltonian, S)
plain_rhf_scf_solver = gqcpy.RHFSCFSolver.Plain()

gqcpy.RHF.optimize(solver, environment).groundStateParameters() # nog niet af!!!!
```

<!--C++-->
```C++
auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, spinor_basis.overlap().parameters());
auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
```
<!--END_DOCUSAURUS_CODE_TABS-->

In order to really confirm that the electronic structure model's parameter are 'optimal', in our case a diagonal Fock matrix, we must define an objective.

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
objective = gqcpy.DiagonalRHFFockMatrixObjective(sq_hamiltonian)  # use the default threshold of 1.0e-08
```

<!--C++-->
```C++
const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
```

<!--END_DOCUSAURUS_CODE_TABS-->

The objective, solver and environment are combined into the _optimize_ method of the RHF QCMethod, which returns a `QCStructure` containing the optimized RHF parameters that satisfy the objective. 

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
rhf_parameters = gqcpy.RHF.optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters()
```

<!--C++-->
```C++
const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();
```

<!--END_DOCUSAURUS_CODE_TABS-->

These optimized RHF parameters can be obtained through the `QCStructure` object. Examples are the coefficient matrix, with every spatial orbital as a column, and orbital energies.

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
C = rhf_parameters.expansion()
energies = rhf_parameters.orbitalEnergies()
```

<!--C++-->
```C++
const auto C = rhf_parameters.expansion();
const auto E = rhf_parameters.orbitalEnergies();
```

<!--END_DOCUSAURUS_CODE_TABS-->

### UHF

We can start off UHF SCF calculation in the same way as we did with RHF SCF calculations.

#### Molecular setup

The first step consists of generating the [`Molecule`](#molecules.md) object from which we can obtain properties. The simplest way to do so is to provide an `xyz` file. 

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
molecule = gqcpy.Molecule.ReadXYZ("h2.xyz")
N = molecule.numberOfElectrons()
N_alpha = N//2
N_beta = N//2
```

<!--C++-->
```C++
const auto molecule = GQCP::Molecule::ReadXYZ("data/h2.xyz");
const auto N = molecule.numberOfElectrons();
const auto N_alpha = N/2;
const auto N_beta = N/2;
```
<!--END_DOCUSAURUS_CODE_TABS-->

Before we go over to the calculation of integrals, an [orbital basis](#orbital_bases.md) must be constructed. The RHF method uses an AO basis. We can immediately use this object to find the overlap integrals.

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
spinor_basis = gqcpy.RSpinOrbitalBasis(molecule, "STO-3G")
S = spinor_basis.quantizeOverlapOperator().parameters()
```

<!--C++-->
```C++
const auto spinor_basis = GQCP::RSpinOrbitalBasis::ReadXYZ(molecule, "STO-3G";
const auto S = spinor_basis.overlap().parameters();
```
<!--END_DOCUSAURUS_CODE_TABS-->

This spinor basis is used to construct a Hamiltonian operator in second quantization (sq). 

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
sq_hamiltonian = gqcpy.SQHamiltonian.Molecular(spinor_basis, molecule)  # in an AO basis
```

<!--C++-->
```C++
const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, h2);  // in an AO basis
```
<!--END_DOCUSAURUS_CODE_TABS-->

#### Plain UHF SCF

To solve the UHF SCF equations, we will need to set up a `Solver` and its associated `Environment`. The environment will provide all necessary intermediates for the solver to use. 

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
uhf_environment = gqcpy.UHFSCFEnvironment.WithCoreGuess(N_alpha, N_beta, sq_hamiltonian, S)
plain_uhf_scf_solver = gqcpy.UHFSCFSolver.Plain()
```

<!--C++-->
```C++
auto uhf_environment = GQCP::UHFSCFEnvironment<double>::WithCoreGuess(N_alpha, N_beta,, sq_hamiltonian, spinor_basis.overlap().parameters());
auto plain_uhf_scf_solver = GQCP::UHFSCFSolver<double>::Plain();
```
<!--END_DOCUSAURUS_CODE_TABS-->
