# Introduction for developers

Welcome to the developer documentation! Here, we'll guide you through the design of the C++ library and the Python bindings.

GQCP has two main components: the C++ library `gqcp` and its associated Python bindings `gqcpy`, which is reflected in the folder structure of the GitHub: files that belong to `gqcp` can be found in the `./gqcp` folder and files that belong `gqcpy` can be found in the `./gqcpy` folder. Pretty straightforward.

> **Note**: Here `./` means the root folder of the repository, i.e. the one that is created after cloning this repository.

The `./gqcp`-folder has three main parts:
- `./gqcp/include`, which collects all C++ header files (`*.hpp`)
- `./gqcp/src`, which collects all C++ source files (`*.cpp`)
- `./gqcp/tests`, which collects all our C++ unit tests (`*_test.cpp`)

The Python bindings don't require any headers, so we only have two subfolders of `./gqcpy`:
- `./gqcpy/src`, which collects all C++ source files that are used to generate the bindings (`*_bindings.cpp`)
- `./gqcpy/examples`, which collects some illustrative calculations in Jupyter Notebooks.

Inside these folders, we have added the following structure:
- __Basis__: for everything related to spinors, spinor basis and transformations
- __DensityMatrix__: for everything related to density matrices
- __Mathematical__: for general mathematical utilities, not directly related to quantum chemistry
- __Molecule__: for everything related to molecules and nuclei
- __ONVBasis__: for everything related to ONVs and bases for Fock (sub)spaces
- __Operator__: for everything related to first- and second-quantized operators
- __Processing__: for collecting everything that happens _after_ the determination of the optimal values of the electronic structure model's parameters, like response calculations
- __QCMethod__: for the determination of the optimal parameters of an electronic structure model
- __Utilities__: for collecting general utilities that do not belong elsewhere


## General structure of the C++ library

In this section, we will go over the main aspects of the C++ library.
It is meant as an introduction to the source code for interested users, but mainly for C++ developers that are interested in how GQCP works under the hood.
By reading this section, we assume that you are familiar with C++.


### Second-quantized operators
Since the GQCP source code uses concepts related to second quantization, we should start by examining the core objects `SQOneElectronOperator` and `SQTwoElectronOperator`: the prefix 'SQ' means 'second-quantized', we have used it in order to somewhat abbreviate the class names. 
`SQOperator`, which is not a real class but is used to refer to either `SQOneElectronOperator` or `SQTwoElectronOperator`, encapsulates the integrals/parameters that are present in the corresponding second-quantized expression, so there must be a way to create an `SQOperator` from its integrals/parameters.
Let us go through a small snippet to create an `SQOneElectronOperator`. 


```cpp
#include <gqcp.hpp>


// Initialize a matrix and convert it into an operator
GQCP::SquareMatrix<double> M (2);  // the matrix representation
M << 1.0, 2.0,
     3.0, 4.0;

const GQCP::ScalarRSQOneElectronOperator<double> op {M};  // the operator itself
```


If the syntax for the matrix-related objects feels familiar for other C++ developers, it might be because we use the amazing linear algebra library [Eigen](http://eigen.tuxfamily.org) for representing vectors, matrices and tensors and performing BLAS/LAPACK operations.
In this example, we first create a matrix representation using Eigen's syntax and afterwards create an `SQOneElectronOperator` from it.
`SQOneElectronOperator` is actually a class template, which has two template arguments.
The first one is the underlying type of the matrix elements/integrals/parameters, and the second is the number of components the operator has.
We are most often working with `SQOperator`s with only one component, like the kinetic energy operator, or the Coulomb repulsion operator, but occasionally, we work with electronic dipole operators as well.

You might have guessed, but `ScalarRSQOperator` is just an alias, which for `SQOneElectronOperator` is written as follows:


```cpp
template <typename Scalar>
using ScalarRSQOneElectronOperator = RSQOneElectronOperator<Scalar, ScalarVectorizer>;
```


In order to access the underlying integrals/parameters/matrix elements, we use the API `.parameters(component_index)` or `.allParameters()`. Building on the previous example, we can now write the following.


```cpp
#include <gqcp.hpp>


// Initialize two matrices and convert them into the two components of an operator
GQCP::SquareMatrix<double> M1 (2);  // the first matrix representation
M1 << 1.0, 2.0,
      3.0, 4.0;

GQCP::SquareMatrix<double> M2 = GQCP::SquareMatrix<double>::Identity(2);  // the second first matrix representation


// Construct an operator with two components
const GQCP::SQOneElectronOperator<double, 2> op ({M1, 2*M1*M2});  // note that we can naturally work with multiplications of matrices thanks to SquareMatrix inheriting from Eigen


const auto all_components = op.allParameters();  // return all the components
std::cout << op.parameters(1)(0,1) << std::endl;  // access the element (0,1) of the first component: this should print '2'
```

We should note that this type of access both works in a read-only and a write way.


### Flexible solver algorithms

In GQCP, we have chosen for a flexible solver design instead of only providing our users with our implementations of certain optimization algorithms, like RHF SCF or Newton-step based minimizers.

The first class we will discuss is the `IterativeSolver`. At its core, it provides the implementation of `.perform()`, which iterates until the maximum number of allowed iterations is reached, or the convergence criterion is reached. In every iteration step, it will check if the convergence criterion is fulfilled, and if it is not, it will continue to execute all the steps in its `StepCollection`.

Convergence criteria can be implemented by deriving from `ConvergenceCriterion` and implementing its `isFulfilled()` method. An example for RHF SCF would be to check the norm on two subsequent density matrices. One important realization is that the iteration steps and the convergence criteria must be able to access the information that the algorithm in its entirety produces. For RHF SCF, this would be the coefficient matrices, the density matrices, the Fock matrices, etc. That is why every `IterativeSolver`, its `StepCollection`, its `Step`s and `ConvergenceCriterion` must all be defined with respect to an `Environment`, which is the template parameter that should be attached to each of these classes.

An example can make many things clear. Suppose we would like to do an RHF SCF calculation. The code that creates our suggested type of plain RHF SCF solver is the following:

```cpp
IterativeAlgorithm<Environment> Plain(const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {

    // Create the iteration cycle that effectively 'defines' a plain RHF SCF solver
    StepCollection<RHFSCFEnvironment<Scalar>> plain_rhf_scf_cycle {};
    plain_rhf_scf_cycle.add(RHFDensityMatrixCalculation<Scalar>())
                        .add(RHFFockMatrixCalculation<Scalar>())
                        .add(RHFFockMatrixDiagonalization<Scalar>())
                        .add(RHFElectronicEnergyCalculation<Scalar>());

    // Create a convergence criterion on the norm of subsequent density matrices
    const auto density_matrix_extractor = [] (const RHFSCFEnvironment<Scalar>& environment) { return environment.density_matrices; };
    const ConsecutiveIteratesNormConvergence<Orbital1DM<Scalar>, RHFSCFEnvironment<Scalar>> convergence_criterion (threshold, density_matrix_extractor);

    return IterativeAlgorithm<RHFSCFEnvironment<Scalar>>(plain_rhf_scf_cycle, convergence_criterion, maximum_number_of_iterations);
}
```

The reader might confirm that the code is clear: in every iteration, the algorithm will check the convergence on the density matrices (cfr `ConsecutiveIteratesNormConvergence` with a suitable iterate extractor). If the algorithm has not converged, an iteration cycle continues in which the following happens:
1. The RHF density matrix is calculated (from the most recent coefficient matrix)
1. The RHF Fock matrix is calculated (from the most recent density matrix);
1. The RHF Fock matrix is diagonalized (to yield a new coefficient matrix);
1. RHF energy is calculated (from the most recent coefficient matrix).

From this example, we can see that every `StepCollection` consists of `Step`s. Each of these `Step`s are instances of classes with an implemented `.execute(Environment)` method that usually 1) read from the environment, 2) calculate some value, 3) write to the environment, but the user is always free to implement his or her desires.

As a conclusion, we have achieved, in essence, a run-time specification of iterative algorithms and we therefore provide the highest amount of flexibility for the (knowing) user to experiment with different solver algorithms. We will continue to provide default implementations, but if a user requires a new kind of algorithm, i.e. one that requires a new kind of environment, he or she only has to implement a new type of `Environment` and the necessary `Step`s.



### Quantum chemical methods and models
We define a quantum chemical model as a set of related electronic structure models with the same classes of parameters. For example, the RHF wave function model can be mathematically written as the orbital rotation operator acting on the closed-shell reference determinant. Its parameters are all the (occupied-virtual) orbital rotation generators. When we optimize a quantum chemical method, such that a certain objective is fulfilled, we have found its 'optimal parameters'.

In the `C++` code, a `QCModel` is a hypothetical class (i.e. it does not really exist) that allows us to express all the variations for the electronic structure model. Coupled to `QCModel`, we inseparably have the class `QCObjective`. Let's take a look at the essential source code.

```cpp
template <typename _DerivedQCObjective>
class QCObjective :
    public CRTP<_DerivedQCObjective> {
public:
    using DerivedQCObjective = _DerivedQCObjective;

public:
    template <typename QCModel>
    bool isSatisfiedWith(const QCModel& model_parameters) const {
        this->derived().isSatisfiedWith(model_parameters);
    }
};
```

We see that we have made `QCObjective` a compile-time base class, whose derived classes should implement `isSatisfiedWith(const QCModel&)`. The optimized parameters are then an instance of `QCModel`, such that a corresponding objective `isSatisfiedWith` it.

Since we had already introduced the `namespace QCMethod`, we chose another name for the actual base class `QCMethodProtocol`, of which we'll inspect the source code.

```cpp
template <typename _QCModel, typename _DerivedQCMethod>
class QCMethodProtocol:
    public CRTP<_DerivedQCMethod> {

public:

    template <typename QCObjective, typename Solver, typename Environment>
    QCStructure<QCModel> optimize(const QCObjective& objective, Solver& solver, Environment& environment) {
        return this->derived().optimize(objective, solver, environment);
    }
};
```

This class is thus naturally coupled to a `QCModel`. Furthermore, classes that derive from `QCMethodProtocol`, or, better: conform to it, must implement an `optimize()` method, which takes an `Objective` to allow to check if the `solver`'s solution actually produces a solution that is consistent with the objective. It's the `QCMethod`'s responsibility to construct a `QCModel` from the `solver`'s solution, which is actually stored in the corresponding `environment` and put them in a `QCStructure`, which is a collection of ground- and possibly excited state model parameters.


## Usage in an external project

In order to use the C++ We have created a small [example](https://github.com/GQCG/GQCP-link) which showcases how to use `GQCP` in an external C++ project.
