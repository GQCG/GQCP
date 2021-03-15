/**
 *  Test driver
 */

#include <gqcp/gqcp.hpp>
#include <stdio.h>

int main() {
    GQCP::Nucleus nucleus (GQCP::elements::elementToAtomicNumber("N"), -1/2, 0, 0);
    std::cout << "Link successful." << std::endl;
    return 0;
}