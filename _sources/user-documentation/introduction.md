# Introduction for users


## Introduction

GQCP naturally uses an object-oriented programming style in order to convey quantum chemical concepts into source code.
Since second quantization is a modern approach for wave function-based quantum chemistry, you will find many of the interfaces intuitive if you are familiar with it.

Using the Python bindings that we have written offers a couple of advantages as a user. First and foremost, they just _work_! We believe we've made it as simple as possible. If, however, you don't think this is the case and you are confused, please open a GitHub issue and we will help you to figure it out.

In all of the user documentation, we assume that you have imported the following libraries.


<!--DOCUSAURUS_CODE_TABS-->

<!--Python-->
```python
import gqcpy
import numpy as np
```


<!--C++-->
```C++
#include <gqcp.hpp>
```
<!--END_DOCUSAURUS_CODE_TABS-->
