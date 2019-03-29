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
#include "Basis/AOBasis.hpp"
#include "Basis/LibintInterfacer.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param shell_set        the set of shells that are placed on the atoms
 */
AOBasis::AOBasis(const ShellSet& shell_set) :
    shell_set (shell_set)
{}


/**
 *  Construct an AO basis by placing shells corresponding to the basisset information on every atom of the molecule. The contraction coefficients in the underlying shells are modified such that the resulting spherical (or axis-aligned Cartesian) GTOs are normalized
 *
 *  @param molecule             the molecule containing the atoms on which the shells should be centered
 *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
 */
AOBasis::AOBasis(const Molecule& molecule, const std::string& basisset_name) :
    AOBasis(ShellSet(molecule, basisset_name))
{
    // At the moment, libint2 is used to read in the STO-3G basisset
    // This means that we should use libint2's normalization (i.e. only embed the normalization coefficients of the primitives, see LibintInterfacer's destructor)
    this->shell_set.embedNormalizationFactorsOfPrimitives();
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the number of basis functions in this AO basis
 */
size_t AOBasis::numberOfBasisFunctions() const {
    return this->shell_set.numberOfBasisFunctions();
}



/*
 *  PUBLIC METHODS - LIBINT INTEGRALS
 */

/**
 *  @return the matrix representation of the overlap operator in this AO basis
 */
OneElectronOperator<double> AOBasis::calculateOverlapIntegrals() const {

    auto libint_basisset = LibintInterfacer::get().interface(this->shell_set);
    return LibintInterfacer::get().calculateOneElectronIntegrals<1>(libint2::Operator::overlap, libint_basisset)[0];
}


/**
 *  @return the matrix representation of the kinetic energy operator in this AO basis
 */
OneElectronOperator<double> AOBasis::calculateKineticIntegrals() const {

    auto libint_basisset = LibintInterfacer::get().interface(this->shell_set);
    return LibintInterfacer::get().calculateOneElectronIntegrals<1>(libint2::Operator::kinetic, libint_basisset)[0];
}


/**
 *  @return the matrix representation of the nuclear attraction operator in this AO basis
 */
OneElectronOperator<double> AOBasis::calculateNuclearIntegrals() const {

    auto libint_basisset = LibintInterfacer::get().interface(this->shell_set);
    auto libint_atoms = LibintInterfacer::get().interface(this->shell_set.atoms());

    return LibintInterfacer::get().calculateOneElectronIntegrals<1>(libint2::Operator::nuclear, libint_basisset, make_point_charges(libint_atoms))[0];
}


/**
 *  @return the matrix representation of the Cartesian components of the electrical dipole operator
 */
std::array<OneElectronOperator<double>, 3> AOBasis::calculateDipoleIntegrals(const Vector<double, 3>& origin) const {

    std::array<double, 3> origin_array {origin.x(), origin.y(), origin.z()};
    auto libint_basisset = LibintInterfacer::get().interface(this->shell_set);

    auto all_integrals = LibintInterfacer::get().calculateOneElectronIntegrals<4>(libint2::Operator::emultipole1, libint_basisset, origin_array);  // overlap, x, y, z

    // Apply the minus sign which comes from the charge of the electrons -e
    return std::array<OneElectronOperator<double>, 3> {-all_integrals[1], -all_integrals[2], -all_integrals[3]};  // we don't need the overlap, so ignore [0]
}


/**
 *  @return the matrix representation of the Coulomb repulsion operator in this AO basis
 */
TwoElectronOperator<double> AOBasis::calculateCoulombRepulsionIntegrals() const {

    auto libint_basisset = LibintInterfacer::get().interface(this->shell_set);
    return LibintInterfacer::get().calculateTwoElectronIntegrals(libint2::Operator::coulomb, libint_basisset);
}


}  // namespace GQCP
