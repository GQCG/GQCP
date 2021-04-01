# Tutorial

GQCP naturally uses an object-oriented programming style in order to convey quantum chemical concepts into source code. In the following tutorial, we assume that you have imported GQCP, either through its Python bindings

```python
import gqcpy
import numpy as np
```

or by directly linking to its C++ core

```C++
#include <gqcp.hpp>
```

We will first discuss [how to create molecules](molecules.md) and a corresponding [orbital basis](orbital-bases.md), after which we will introduce quantum chemical methods that can be used.
