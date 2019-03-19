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
 *  A class that represents an atomic orbital basis, i.e. the collection of (scalar) atomic orbitals/basis functions
 */
class AOBasis {
private:
    ShellSet shell_set;  // the underlying shell set


public:
    // CONSTRUCTORS
    /**
     *  Construct an AO basis by placing shells shells corresponding to the basisset information on every atom of the molecule
     *
     *  @param molecule             the molecule containing the atoms on which the shells should be centered
     *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
     */
    AOBasis(const Molecule& molecule, const std::string& basisset_name);


    // GETTERS
    const ShellSet& get_shell_set() const { return this->shell_set; }


    // PUBLIC METHODS
    /**
     *  @return the number of basis functions in this AO basis
     */
    size_t numberOfBasisFunctions() const;


    // PUBLIC METHODS - LIBINT INTEGRALS
    /**
     *  @return the matrix representation of the overlap operator in this AO basis
     */
    OneElectronOperator<double> calculateOverlapIntegrals() const;

    /**
     *  @return the matrix representation of the kinetic energy operator in this AO basis
     */
    OneElectronOperator<double> calculateKineticIntegrals() const;

    /**
     *  @return the matrix representation of the nuclear attraction operator in this AO basis
     */
    OneElectronOperator<double> calculateNuclearIntegrals() const;

    /**
     *  @param origin       the origin of the dipole
     *
     *  @return the matrix representation of the Cartesian components of the electrical dipole operator in this AO basis
     */
    std::array<OneElectronOperator<double>, 3> calculateDipoleIntegrals(const Vector<double, 3>& origin = Vector<double, 3>::Zero()) const;

    /**
     *  @return the matrix representation of the Coulomb repulsion operator in this AO basis
     */
    TwoElectronOperator<double> calculateCoulombRepulsionIntegrals() const;
};


}  // namespace GQCP


#endif  // GQCP_AOBASIS_HPP
