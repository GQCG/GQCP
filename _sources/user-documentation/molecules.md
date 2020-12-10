# Molecules


## Constructing molecules

In this section, we'll go over various ways to create molecules in gqcpy.


### XYZ-files

The simplest option is creating a geometric `.xyz`-file. This is a data file that specifies a molecular geometry, using coordinates in Ångstrom. The first line contains the number of atoms, while subsequent lines specify the coordinates of each element. If we would like to study a water molecule, here's what its `.xyz`-file would look like.

```
3

O          0.00000       -0.07579        0.00000
H          0.86681        0.60144        0.00000
H         -0.86681        0.60144        0.00000
```

> **Note:** The coordinates in an `.xyz`-file are specified in Ångstrom.

If we name the file `water.xyz`, we can read it in as follows.

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
water_molecule = gqcpy.Molecule.ReadXYZ("water.xyz")  # creates a neutral molecule
water_cation = gqcpy.Molecule.ReadXYZ("water.xyz", charge=+1)  # creates a charged ion
```

<!--C++-->
```C++
const auto water_molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");  // creates a neutral molecule
const auto water_cation = GQCP::Molecule::ReadXYZ("data/h2o.xyz", -1);  // creates a charged ion
```
<!--END_DOCUSAURUS_CODE_TABS-->



### Hydrogen toy systems

Often, when developing electronic structure models, we use toy systems that consist of only hydrogen atoms. We offer simplified methods for the creation of hydrogen chains and hydrogen rings.

For a chain of H-atoms, we specify the number of total hydrogens and the uniform separation. In order to create a H<sub>3</sub>-chain where adjacent hydrogens are 1 a.u. apart, we can use the following code.
<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
h3_chain = gqcpy.Molecule.HChain(3, 1.0)
```

<!--C++-->
```C++
const auto h3_chain = GQCP::Molecule::HChain(3, 1.0);
```
<!--END_DOCUSAURUS_CODE_TABS-->


If we would like to create a chain of H<sub>2</sub>-molecules, we can use the `H2Chain` method. Here, we'll have to specify (in order) the number of H<sub>2</sub>-molecules, the internuclear distance (i.e. the length of each H<sub>2</sub>-molecule) and the intermolecular distance (i.e. the distance between adjacent H<sub>2</sub>-molecules). In the following example, we'll create a molecule consisting of 4 hydrogens, where the internuclear distances are alternating between 1.0 and 1.5.

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
double_h2_chain = gqcpy.Molecule.H2Chain(2, 1.0, 1.5)
```

<!--C++-->
```C++
const auto double_h2_chain = GQCP::Molecule::H2Chain(2, 1.0, 1.5);
```
<!--END_DOCUSAURUS_CODE_TABS-->


We can create hydrogen in two ways. We can either specify the radius of the corresponding circle, or we can specify the internuclear distance between adjacent hydrogens. 

In order to create an H<sub>3</sub>-triangle, whose vertices are 1.0 a.u. away from the center, we can use `HRingFromRadius`.
<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
h3_ring = gqcpy.Molecule.HRingFromRadius(3, 1.0)
```

<!--C++-->
```C++
const auto h3_ring = GQCP::Molecule::HRingFromRadius(3, 1.0);
```
<!--END_DOCUSAURUS_CODE_TABS-->


To create an H<sub>4</sub>-square, whose edge lengths are 1.0 a.u., we can use `HRingFromDistance(4, 1.0)`.
<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
h4_ring = gqcpy.Molecule.HRingFromDistance(4, 1.0)
```

<!--C++-->
```C++
const auto h4_ring = GQCP::Molecule::HRingFromDistance(4, 1.0);
```
<!--END_DOCUSAURUS_CODE_TABS-->



### General construction
Molecules are a collection of nuclei and a number of electrons. In general, we can specify a molecule by supplying its nuclei and a charge. 

To generate a `Nucleus`, its coordinates and number of electrons are needed. For example, the nitrogen nucleus (with atom number 7) on position `(0.0, 0.0, 0.0)` is constructed as follows:

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
nucleus = Nucleus(7, 0.0,0.0,0.0)
```

<!--C++-->
```C++
const GQCP::Nucleus nucleus {7, 0.0,0.0,0.0};
```
<!--END_DOCUSAURUS_CODE_TABS-->

We can then combine nuclei in order to form a `Molecule`, given its charge. For example, we can create N<sub>2</sub> using two nitrogen nuclei. Here, we're using the explicit constructor with a charge equal to 0, but for neutral molecules, we may omit this argument.

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
# Here is the explicit code.
N_1 = gqcpy.Nucleus(7, 0.0,0.0,0.0)
N_2 = gqcpy.Nucleus(7, 0.0,0.0,1.0)

molecule = gqcpy.Molecule([N_1, N_2], charge=0)


# Here's a simple one-liner.
molecule = gqcpy.Molecule([Nucleus(7, 0.0,0.0,0.0), gqcpy.Nucleus(7, 0.0,0.0,1.0)], 0)
```

<!--C++-->
```C++
// Here is the explicit code.
const GQCP::Nucleus N_1 {7, 0.0,0.0,0.0};
const GQCP::Nucleus N_2 {7, 0.0,0.0,1.0};
const std::vector<GQCP::Nucleus> nuclei {N_1, N_2};

const GQCP::Molecule molecule {nuclei, 0};


// Here's a simple one-liner.
const GQCP::Molecule molecule {{{7, 0.0,0.0,0.0}, {7, 0.0,0.0,1.0}}, 0};
```
<!--END_DOCUSAURUS_CODE_TABS-->


## Inspecting molecules

After a `Molecule` object has been made, we can use this object to gather information on its properties. For a general output of a molecule's nuclear positions and number of electrons, we can print it to console.

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
print(molecule)
```

<!--C++-->
```C++
std::cout << molecule.description();
```
<!--END_DOCUSAURUS_CODE_TABS-->

> **Note:** To specify coordinates, GQCP internally uses Bohr (atomic units).


In order to obtain information on the number of electrons, or number of electron pairs, we can use following methods.

<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
N = molecule.numberOfElectrons()
N_P = molecule.numberOfElectronPairs()
```

<!--C++-->
```C++
const auto N = molecule.numberOfElectrons();
const auto N_P = molecule.numberOfElectronPairs();
```
<!--END_DOCUSAURUS_CODE_TABS-->

> **Note:** For a molecule with an odd number of electrons, `.numberOfElectronPairs()` returns the number of _completed_ electron pairs. For example, H<sub>2</sub><sup>-</sup> would still only have 1 electron pair.
