---
id: developer_ONV_path
title: ONV path manipulation
sidebar_label: ONV path manipulation
---

Bit string manipulations can be very hard to understand. This section should help aspiring developers how GQCP deals with the manipulation of ONV's (occupation number vectors) and ONV paths. It is necessary that you understand the theory first, since the code tries to follow the theory as closely as possible. To brush up on CI algorithms and ONV manipulation we highly recommend chapter 11 in [Molecular Electronic‐Structure Theory](https://onlinelibrary.wiley.com/doi/book/10.1002/9781119019572) by Helgaker, Jørgensen and Olsen.

## The key operators

OONV path manipulations are described using a second quantization approach. GQCP's API sticks as closely to this as possible. To manipulate a path, we use the elementary second quantization operators to annihilate or create arcs that represent occupied orbitals. In an ONV path, a diagonal arc represents an occupied orbital, while a vertical arc represents an unoccupied orbital.  We will use the path of ONV `10110` as an example.

![init_path](/GQCP/img/ONVPath_10110.png)

The function `annihilate()` removes a diagonal arc, starting at the given coordinate in the ONV path.

```C++
/**
*  Annihilate the diagonal arc that starts at the coordinate (q,n).
* 
*  @param q        the index of the orbital that should be annihilated
*  @param n        the number of electrons in the ONV/path up to the orbital index q
*/
void ONVPath::annihilate(const size_t q, const size_t n) {

    // During annihilation, we're removing an arc, so we'll have to update the current address by removing the corresponding arc weight.
    this->m_address -= this->onv_basis.arcWeight(q, n);

    // Update the next possible creation index. Since we're always constructing paths from the top-left to the bottom-right, we're only considering creation indices p > q.
    this->p = q + 1;
}
```

There are three diagonal arcs (starting at [0, 0], [2, 1] and [3, 2]). These are the coordinates where `annihilate()` can remove an electron. Let's say we `annihilate(0, 0)`. The diagonal arc starting at [0, 0] will then become a vertical arc, resulting in a new path:

![firs_annihilate](/GQCP/img/ONVPath_00110_1.png)

The ONV path is now 'open'. To 'close' it, the `create()` function can be used.

```C++
/**
*  Create the diagonal arc that starts at the coordinate (p, n).
* 
*  @param p        the index of the orbital that should be created
*  @param n        the number of electrons in the ONV/path up to the orbital index q, prior to the creation
*/
void ONVPath::create(const size_t p, const size_t n) {

    // During creation, we're adding an arc, so we'll have to update the current address by adding the corresponding arc weight.
    this->m_address += this->onv_basis.arcWeight(p, n);
}
```

 In this particular case, the path can be closed by creating at coordinate [1, 1], resulting in a new ONV: `01110`.

![firs_create](/GQCP/img/ONVPath_01110.png)

We will come back to this later, when looking at a complete example, side-by-side with the source code.  

This is not the only way to close the path. A more involved way consists of shifting the following diagonals (occupied orbitals) to the left until a vertical arc is encountered to close the gap. This is the use of `leftTranslate()`.

```C++
/**
  *  Translate the diagonal arc that starts at the coordinate (p, n) to the left.
  * 
  *  @param p        the index of the orbital that should be annihilated
  *  @param n        the number of electrons in the ONV/path up to the orbital index p
  */
 void ONVPath::leftTranslate(const size_t p, const size_t n) {

     // Translating a diagonal arc can be rewritten as an annihilation, followed by a creation.
     this->annihilate(p, n);
     this->create(p, n - 1);

     // Since a left-translation describes the process of 'encountering an electron/occupied orbital', the sign factor should be updated according to the fermionic anticommutation rules.
     this->m_sign *= -1;
}
```

In the initial example, the annihilation could be directly followed by a creation, since the next arc was vertical; i.e. an unoccupied orbital. This is not always the case however. If the next orbital is also an occupied orbital, the arc needs to be translated to the left, in order to check whether the subsequent arc represents an occupied or unoccupied orbital. Let us look back at the ONV `10110` after an `annihilate(0, 0)` operation, but assume we don't want to create at coordinate [1, 1], as previously shown. Since our second arc will remain vertical, it can be moved to the left without changing the address. The next index is the starting point of a diagonal arc, which means a creation is not possible there. But since we cannot translate this arc without changing the address, we need `leftTranslate()`. The resulting ONV paths we can be shown as:

![firs_annihilate_shift](/GQCP/img/ONVPath_00110_2_3_side.png)

with the left image showing the situation before `leftTranslate()` and the right image showing the situation after the operation. It is important to note that `leftTranslate()` multiplies the sign of the ONV path. This is done because of the fermionic anti-commutation relation which has to be taken into account each time an occupied orbital is encountered. 

## Example: ONV 10110

In this example, we'll show how the operators work on an ONV by visualizing each step, accompanied by the source code. The first step is to create an ONV basis. From this basis, an ONV, `10110`, can then be created. 

<!--DOCUSAURUS_CODE_TABS-->

<!--C++-->
```C++
// Set up a F(5,3) Fock space.
const size_t M = 5;
const size_t N = 3;

const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};

const auto I = GQCP::SpinUnresolvedONV::FromOccupiedIndices({0, 2, 3}, 5);  // |10110>
GQCP::ONVPath onv_path {onv_basis, I};
```

<!--Visual representation-->
```latex
ONV: 10110
address: 2
```

![init_path](/GQCP/img/ONVPath_10110.png)

<!--END_DOCUSAURUS_CODE_TABS-->

Let's annihilate the electron on the first orbital (index 0). This 'opens' the path, by turning a diagonal arc into a vertical one. 

<!--DOCUSAURUS_CODE_TABS-->

<!--C++-->
```C++
// Annihilate electron on the first orbital.
onv_path.annihilate(0, 0);
```

<!--Visual representation-->
```latex
ONV: 00110
address: 2
```

![firs_annihilate](/GQCP/img/ONVPath_00110_1.png)

<!--END_DOCUSAURUS_CODE_TABS-->

The next possible index to create an electron is [0, 1]. Creating an electron is turning a vertical arc into a diagonal one.

<!--DOCUSAURUS_CODE_TABS-->

<!--C++-->
```C++
// Create an electron on the second orbital.
onv_path.create(1, 0);
```

<!--Visual representation-->
```latex
ONV: 01110
address: 3
```

![first_create](/GQCP/img/ONVPath_01110.png)

<!--END_DOCUSAURUS_CODE_TABS-->

There is however another index where we can create an electron: index 4. In order to do this, let us undo the previous creation and start from ONV `00110`. The first arc is the unoccupied arc on which we do not want to create an electron. Moving this unoccupied orbital to the left does not change our ONV or it's address.

![init_path](/GQCP/img/ONVPath_00110_2.png)

This will be the starting point to find the other ONV coupled to `10110`. We are now in a situation where the next creation index corresponds to an occupied orbital. Since we can't create on these indices, we must translate the diagonal arc that starts at (2, 1) to (2, 0). A translation to the left means we've encountered an electron, and thus the sign should be updated.

<!--DOCUSAURUS_CODE_TABS-->

<!--C++-->
```C++
// Translate the diagonal arc that starts at (2, 1) to (2, 0).
onv_path.leftTranslate(2, 1);
```

<!--Visual representation-->
```latex
ONV: 00110
address: 3
```

![first_shift](/GQCP/img/ONVPath_00110_3.png)

<!--END_DOCUSAURUS_CODE_TABS-->

After the translation, the situation is still very similar. The next possible creation index is still an occupied index. We need to perform another translation. Once again, the sign is updated. Luckily, `leftTranslate()` does this automatically.

<!--DOCUSAURUS_CODE_TABS-->

<!--C++-->
```C++
// Translate the diagonal arc that starts at (3, 2) to (3, 1).
onv_path.leftTranslate(3, 2);
```

<!--Visual representation-->
```latex
ONV: 00110
address: 5
```

![first_shift](/GQCP/img/ONVPath_00110_4.png)

<!--END_DOCUSAURUS_CODE_TABS-->

We can now close the path by creating an electron, since the next index is now unoccupied.

<!--DOCUSAURUS_CODE_TABS-->

<!--C++-->
```C++
// Translate the diagonal arc that starts at (2, 1) to (2, 0).
onv_path.create(4, 2);
```

<!--Visual representation-->
```latex
ONV: 00111
address: 9
```

![first_shift](/GQCP/img/ONVPath_00111.png)

<!--END_DOCUSAURUS_CODE_TABS-->
