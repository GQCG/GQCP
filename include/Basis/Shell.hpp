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


namespace GQCP {


/**
 *  A class that represents a shell of GTOs: it specifies in a condensed way which basis functions are on an atom
 *
 *  Note that an inclusion of any normalization factor inside the contraction coefficients is irrelevant if they do not change after the overlap matrix between the underlying basis functions has been calculated and all the integrals are calculated using the same contraction coefficients
 */
class Shell {
private:
    bool pure;  // true if spherical, false if Cartesian
    size_t l;  // the angular momentum of the shell
    Atom atom;  // atom on which the shell is centered
    std::vector<double> gaussian_exponents;  // Gaussian exponents (i.e. for the exponential), shared for every contraction
    std::vector<double> contraction_coefficients;


public:
    // CONSTRUCTORS
    /**
     *  @param l                            the angular momentum of the shell
     *  @param atom                         the atom on which the shell is centered
     *  @param gaussian_exponents           the Gaussian exponents, which are shared for every contraction
     *  @param contraction_coefficients     the contraction coefficients
     *  @param pure                         whether the shell is considered to be spherical or not
     */
    Shell(size_t l, const Atom& atom, const std::vector<double>& gaussian_exponents, const std::vector<double>& contraction_coefficients, bool pure=true);


    // GETTERS
    bool is_pure() const { return this->pure; }
    size_t get_l() const { return this->l; }
    const Atom& get_atom() const { return this->atom; }
    const std::vector<double>& get_gaussian_exponents() const { return this->gaussian_exponents; }
    const std::vector<double>& get_contraction_coefficients() const { return this->contraction_coefficients; }


    // OPERATORS
    /**
     *  @param rhs      the right-hand side of the operator ==
     *
     *  @return if this shell is considered equal to the other
     */
    bool operator==(const Shell& rhs) const;


    // PUBLIC METHODS
    /**
     *  @return the number of basis functions that are in this shell
     */
    size_t numberOfBasisFunctions() const;

    /**
     *  @return the size of the contraction in the shell, i.e. the number of primitives contracted in this shell
     */
    size_t contractionSize() const;
};


}  // namespace GQCP


#endif  /* Shell_hpp */
