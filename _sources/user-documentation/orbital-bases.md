# Orbital bases

## Creating an orbital basis

An intermediary step in-between the creation of a `Molecule` object and the calculation of the _integrals_ for this molecule, is the initialization of an orbital basis. In molecular electronic structure theory based on wave function models, we often express a set of orthonormal spinors or spin-orbitals as linear combinations of basis functions.

Here's where the spinor-basis objects come into play. They encapsulate the expansion coefficients, together with the atomic orbital basis. You can create a restricted spin-orbital basis (`RSpinOrbitalBasis`) by specifying the molecule and the name of the basisset that you would like to use.

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
spinor_basis = RSpinOrbitalBasis(molecule, "STO-3G")
```

<!--C++-->
```C++
const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};  // note that you have to specify the type of shell that underlies this spin-orbital basis
```
<!--END_DOCUSAURUS_CODE_TABS-->

The constructor then places scalar AOs on every nucleus in the molecule, according to the basisset specification.

> **Note**: Spinor bases constructed in this way represent non-orthogonal spin-orbital bases.

> **Additional information**: For more theoretical information, please visit our knowdes on [spinor bases](https://gqcg-res.github.io/knowdes/general-spinor-bases.html), [spin-orbital bases](https://gqcg-res.github.io/knowdes/spin-orbital-bases.html) and on [quantizing one- and two-electron operators](https://gqcg-res.github.io/knowdes/quantizing-one-and-two-electron-operators-in-general-spinor-bases.html).


## Calculating integrals

One of the most important ingredients for any quantum chemical method are its one- and two-electron integrals. Calculating the _integrals_, which is actually short for calculating the matrix representation of various first-quantized one- and two-electron operators in a spinor basis, is made fairly simple. If you are acquainted with the modern electronic structure framework of second quantization, you will find it very natural to use.

`SpinorBases` are responsible for expressing spinors in the underlying scalar orbitals. Subsequently, a `SpinorBasis` instance provides `.quantize()` methods that can create second-quantized representations of the corresponding first-quantized operators. More or less, this means that `SpinorBasis` instances provide the interface to calculate the matrix representations of one- and two-electron integrals.

Below, we can find examples how the overlap, kinetic energy, potential energy and interelectronic repulsion energy operators can be expressed in the recently initialized `spinor_basis`. The process of expressing these first-quantized one- and two-electron integrals in the spinor basis is what we call _quantizing_.

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
S_op = spinor_basis.quantizeOverlapOperator()
T_op = spinor_basis.quantizeKineticOperator()
V_op = spinor_basis.quantizeNuclearAttractionOperator(molecule)  # the nuclear attraction operator is defined with respect to the molecule's nuclear framework

g_op = spinor_basis.quantizeCoulombRepulsionOperator()
```

<!--C++-->
```C++
const auto S = spinor_basis.quantize(GQCP::Operator::Overlap());
const auto T = spinor_basis.quantize(GQCP::Operator::Kinetic());
const auto V = spinor_basis.quantize(GQCP::Operator::NuclearAttraction(molecule));  // the nuclear attraction operator is defined with respect to the molecule's nuclear framework

const auto g = spinor_basis.quantize(GQCP::Operator::Coulomb());
```
<!--END_DOCUSAURUS_CODE_TABS-->


The instances `S_op`, `T_op`, `V_op` and `g_op` represent the second-quantized representations of the overlap, kinetic energy, nuclear attraction and interelectronic repulsion operators and they encapsulate the corresponding matrix representations of the integrals. We can access the raw matrix representations using the `.parameters()` method.

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
S = S_op.parameters()
T = T_op.parameters()
V = V_op.parameters()

g = g_op.parameters()
```

<!--C++-->
```C++
const auto S = S_op.parameters();
const auto T = T_op.parameters();
const auto V = V_op.parameters();

const auto g = g_op.parameters();
```
<!--END_DOCUSAURUS_CODE_TABS-->

> **Note**: The matrices/tensors that we have calculated here are associated to the non-orthogonal AO basis.
