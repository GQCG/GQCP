// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#ifndef Shell_hpp
#define Shell_hpp


#include "Atom.hpp"
#include "Basis/BasisFunction.hpp"


namespace GQCP {


/**
 *  A class that represents a shell. It is a specification of contraction coefficients and corresponding exponents for primitives with the same angular momentum, centered on an atom
 *
 *  Note that GQCP::Shell always represents a shell of Cartesian GTOs
 */
class Shell {
private:
    size_t l;  // angular momentum (x + y + z)
    Atom atom;  // atom on which the shell is centered
    std::vector<double> exponents;  // exponents, shared for every contraction
    std::vector<double> coefficients;  // contraction coefficients  // TODO: rename to contraction_coefficients?


public:
    // CONSTRUCTORS
    /**
     *  @param l                the angular momentum of the shell (x + y + z)
     *  @param atom             the atom on which the shell is centered
     *  @param exponents        the exponents, which are shared for every contraction
     *  @param coefficients     the contraction coefficients
     */
    Shell(size_t l, const Atom& atom, const std::vector<double>& exponents, const std::vector<double>& coefficients);


    // GETTERS
    size_t get_l() const { return this->l; }
    const Atom& get_atom() const { return this->atom; }
    const std::vector<double>& get_exponents() const { return this->exponents; }
    const std::vector<double>& get_coefficients() const { return this->coefficients; }


    // PUBLIC METHODS
    /**
     *  @return the number of basis functions that are in this shell
     */
    size_t numberOfBasisFunctions() const;

    /**
     *  @return the basis functions that are represented by this shell
     */
    std::vector<BasisFunction> basisFunctions() const;

    /**
     *  @return the length of the contraction in the shell, i.e. the number of primitives contracted in this shell
     */
    size_t contractionLength() const;
};


}  // namespace GQCP


#endif  /* Shell_hpp */
