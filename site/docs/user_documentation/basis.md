---
id: user_basis
title: Orbital bases
sidebar_label: Orbital bases
---

## Generating a basis set

To represent an electronic wave function, a basis set containing a set of basis functions is used. The 

## Calculating integrals

One of the most important ingredients for any quantum chemical method are its one- and two-electron integrals. Calculating the `integrals', which is actually short for calculating the matrix representation of various first-quantized one- and two-electron operators in a spinor basis, is made fairly simple. If you are acquainted with the modern electronic structure framework of second quantization, you will find it very natural to use.

`SpinorBases` are responsible for expressing spinors in the underlying scalar orbitals. Subsequently, a `SpinorBasis` instance provides `.quantize()` methods that can create second-quantized representations of the corresponding first-quantized operators. More or less, this means that `SpinorBasis` instances provide the interface to calculate the matrix representations of one- and two-electron integrals.



Below, we can find examples how the overlap, kinetic energy, potential energy and interelectronic repulsion energy operators can be expressed in the recently initialized `spinor_basis`. The process of expressing these first-quantized one- and two-electron integrals in the spinor basis is what we call _quantizing_.

```python
spinor_basis = SpinorBasis(molecule = molecule, basisset_name = "STO-3G")



S_op = spinor_basis.quantizeOverlapOperator()
T_op = spinor_basis.quantizeKineticOperator()
V_op = spinor_basis.quantizeNuclearAttractionOperator(molecule)  # the nuclear attraction operator is defined with respect to the nuclear framework
g_op = spinor_basis.quantizeCoulombRepulsionOperator()
```

The instances `S_op`, `T_op`, `V_op` and `g_op` represent the second-quantized representations of the overlap, kinetic energy, nuclear attraction and interelectronic repulsion operators and they encapsulate the corresponding matrix representations of the integrals. We can access the raw matrix representations using the `.parameters()` method.

```python
S = S_op.parameters()
T = T_op.parameters()
V = V_op.parameters()
g = g_op.parameters()
```