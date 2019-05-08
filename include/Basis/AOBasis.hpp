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
#ifndef GQCP_AOBASIS_HPP
#define GQCP_AOBASIS_HPP


#include "Basis/ShellSet.hpp"
#include "Operator/OneElectronOperator.hpp"
#include "Operator/TwoElectronOperator.hpp"


namespace GQCP {


/**
 *  A class that represents an atomic orbital basis: it represents a collection of (scalar) atomic orbitals/basis functions
 */
class AOBasis {
private:
    ShellSet shell_set;  // the underlying collection of shells


public:
    // CONSTRUCTORS
    /**
     *  @param shell_set        the underlying collection of shells
     */
    AOBasis(const ShellSet& shell_set);

    /**
     *  Construct an AO basis by placing shells corresponding to the basisset specification on every atom of the molecule
     *
     *  @param molecule             the molecule containing the atoms on which the shells should be centered
     *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
     *
     *  Note that the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     */
    AOBasis(const Molecule& molecule, const std::string& basisset_name);


    // GETTERS
    const ShellSet& get_shell_set() const { return this->shell_set; }


    // PUBLIC METHODS
    /**
     *  @return the number of basis functions in this AO basis
     */
    size_t numberOfBasisFunctions() const;


    // PUBLIC METHODS - LIBINT2 INTEGRALS
    /**
     *  @return the matrix representation of the overlap operator in this AO basis, using the libint2 integral engine
     */
    OneElectronOperator<double> calculateLibintOverlapIntegrals() const;

    /**
     *  @return the matrix representation of the kinetic energy operator in this AO basis, using the libint2 integral engine
     */
    OneElectronOperator<double> calculateLibintKineticIntegrals() const;

    /**
     *  @return the matrix representation of the nuclear attraction operator in this AO basis, using the libint2 integral engine
     */
    OneElectronOperator<double> calculateLibintNuclearIntegrals() const;

    /**
     *  @param origin       the origin of the dipole
     *
     *  @return the matrix representation of the Cartesian components of the electrical dipole operator in this AO basis, using the libint2 integral engine
     */
    std::array<OneElectronOperator<double>, 3> calculateLibintDipoleIntegrals(const Vector<double, 3>& origin = Vector<double, 3>::Zero()) const;

    /**
     *  @return the matrix representation of the Coulomb repulsion operator in this AO basis, using the libint2 integral engine
     */
    TwoElectronOperator<double> calculateLibintCoulombRepulsionIntegrals() const;


    // PUBLIC METHODS - LIBCINT INTEGRALS
    /**
     *  @return the matrix representation of the overlap operator in this AO basis, using the libcint integral engine
     */
    OneElectronOperator<double> calculateLibcintOverlapIntegrals() const;

    /**
     *  @return the matrix representation of the kinetic energy operator in this AO basis, using the libcint integral engine
     */
    OneElectronOperator<double> calculateLibcintKineticIntegrals() const;

    /**
     *  @return the matrix representation of the nuclear attraction operator in this AO basis, using the libcint integral engine
     */
    OneElectronOperator<double> calculateLibcintNuclearIntegrals() const;

//    /**
//     *  @param origin       the origin of the dipole
//     *
//     *  @return the matrix representation of the Cartesian components of the electrical dipole operator in this AO basis, using the libcint integral engine
//     */
//    std::array<OneElectronOperator<double>, 3> calculateLibcintDipoleIntegrals(const Vector<double, 3>& origin = Vector<double, 3>::Zero()) const;
//
//    /**
//     *  @return the matrix representation of the Coulomb repulsion operator in this AO basis, using the libcint integral engine
//     */
//    TwoElectronOperator<double> calculateLibcintCoulombRepulsionIntegrals() const;
};


}  // namespace GQCP


#endif  // GQCP_AOBASIS_HPP
