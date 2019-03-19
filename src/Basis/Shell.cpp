#include "Basis/Shell.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param l                the angular momentum of the shell (x + y + z)
 *  @param atom             the atom on which the shell is centered
 *  @param exponents        the exponents, which are shared for every contraction
 *  @param coefficients     the contraction coefficients
 */
Shell::Shell(size_t l, const Atom& atom, const std::vector<double>& exponents, const std::vector<double>& coefficients) :
    l (l),
    atom (atom),
    exponents (exponents),
    coefficients (coefficients)
{
    if (exponents.size() != coefficients.size()) {
        throw std::invalid_argument("Shell(size_t, Atom, std::vector<double>, std::vector<double>): the exponents and contraction coefficients must match in size.");
    }
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the number of basis functions that are in this shell
 */
size_t Shell::numberOfBasisFunctions() const {
    return (this->l + 1) * (this->l + 2) / 2;  // Cartesian shell
}


}  // namespace GQCP
