# Quantum chemical models


In GQCP, a quantum chemical [(wave function) model](https://gqcg-res.github.io/knowdes/wave-function-models.html) is defined to be a parametrization of a certain kind of wave function. To encapsulate the optimizable parameters related to a certain quantum chemical method, we currently provide the following models:

- [CCD](#ccd)
- [CCSD](#ccsd)
- [CI](#ci)
- [AP1roG](#ap1rog)
- [vAP1roG](#vap1rog)
- [RHF](#rhf)
- [UHF](#uhf)


## CCD

## CCSD

## CI

All of the configuration interaction wave functions are expressed using the wave function model `LinearExpansion`. [Once obtained](user_quantum_chemical_methods), you can calculate the associated one- and two-electron density matrices straightforwardly. 

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
D = linear_expansion.calculate1DM()
d = linear_expansion.calculate2DM()
```

<!--C++-->
```C++
const auto D = linear_expansion.calculate1DM();
const auto d = linear_expansion.calculate2DM();
```

If the underlying ONV basis is spin-resolved, the C++ library also offers the alpha and beta components.
```C++
const auto D_spin_resolved = linear_expansion.calculateSpinResolved1DM();
const auto d_spin_resolved = linear_expansion.calculateSpinResolved2DM();
```
<!--END_DOCUSAURUS_CODE_TABS-->


## AP1roG

## vAP1roG

## RHF

## UHF