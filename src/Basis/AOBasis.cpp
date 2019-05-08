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
#include "Basis/LibcintInterfacer.hpp"


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
 *  Construct an AO basis by placing shells corresponding to the basisset specification on every atom of the molecule
 *
 *  @param molecule             the molecule containing the atoms on which the shells should be centered
 *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
 *
 *  Note that the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
 */
AOBasis::AOBasis(const Molecule& molecule, const std::string& basisset_name) :
    AOBasis(ShellSet(molecule, basisset_name))
{
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
 *  PUBLIC METHODS - LIBINT2 INTEGRALS
 */

/**
 *  @return the matrix representation of the overlap operator in this AO basis
 */
OneElectronOperator<double> AOBasis::calculateLibintOverlapIntegrals() const {

    auto libint_basisset = LibintInterfacer::get().interface(this->shell_set);
    return LibintInterfacer::get().calculateOneElectronIntegrals<1>(libint2::Operator::overlap, libint_basisset)[0];
}


/**
 *  @return the matrix representation of the kinetic energy operator in this AO basis
 */
OneElectronOperator<double> AOBasis::calculateLibintKineticIntegrals() const {

    auto libint_basisset = LibintInterfacer::get().interface(this->shell_set);
    return LibintInterfacer::get().calculateOneElectronIntegrals<1>(libint2::Operator::kinetic, libint_basisset)[0];
}


/**
 *  @return the matrix representation of the nuclear attraction operator in this AO basis
 */
OneElectronOperator<double> AOBasis::calculateLibintNuclearIntegrals() const {

    auto libint_basisset = LibintInterfacer::get().interface(this->shell_set);
    auto libint_atoms = LibintInterfacer::get().interface(this->shell_set.atoms());

    return LibintInterfacer::get().calculateOneElectronIntegrals<1>(libint2::Operator::nuclear, libint_basisset, make_point_charges(libint_atoms))[0];
}


/**
 *  @param origin       the origin of the dipole
 *
 *  @return the matrix representation of the Cartesian components of the electrical dipole operator in this AO basis
 */
std::array<OneElectronOperator<double>, 3> AOBasis::calculateLibintDipoleIntegrals(const Vector<double, 3>& origin) const {

    std::array<double, 3> origin_array {origin.x(), origin.y(), origin.z()};

    auto libint_basisset = LibintInterfacer::get().interface(this->shell_set);
    auto all_integrals = LibintInterfacer::get().calculateOneElectronIntegrals<4>(libint2::Operator::emultipole1, libint_basisset, origin_array);  // overlap, x, y, z

    // Apply the minus sign which comes from the charge of the electrons -e
    return std::array<OneElectronOperator<double>, 3> {-all_integrals[1], -all_integrals[2], -all_integrals[3]};  // we don't need the overlap, so ignore [0]
}


/**
 *  @return the matrix representation of the Coulomb repulsion operator in this AO basis
 */
TwoElectronOperator<double> AOBasis::calculateLibintCoulombRepulsionIntegrals() const {

    auto libint_basisset = LibintInterfacer::get().interface(this->shell_set);
    return LibintInterfacer::get().calculateTwoElectronIntegrals(libint2::Operator::coulomb, libint_basisset);
}



/*
 *  PUBLIC METHODS - LIBCINT INTEGRALS
 */

/**
 *  @return the matrix representation of the overlap operator in this AO basis, using the libcint integral engine
 */
OneElectronOperator<double> AOBasis::calculateLibcintOverlapIntegrals() const {

    const LibcintInterfacer libcint_interfacer;
    auto raw_container = libcint_interfacer.convert(this->shell_set);
    return libcint_interfacer.calculateOneElectronIntegrals<1>(cint1e_ovlp_cart, raw_container)[0];
}


/**
 *  @return the matrix representation of the kinetic energy operator in this AO basis, using the libcint integral engine
 */
OneElectronOperator<double> AOBasis::calculateLibcintKineticIntegrals() const {

    const LibcintInterfacer libcint_interfacer;
    auto raw_container = libcint_interfacer.convert(this->shell_set);
    return libcint_interfacer.calculateOneElectronIntegrals<1>(cint1e_kin_cart, raw_container)[0];
}


/**
 *  @return the matrix representation of the nuclear attraction operator in this AO basis, using the libcint integral engine
 */
OneElectronOperator<double> AOBasis::calculateLibcintNuclearIntegrals() const {

    const LibcintInterfacer libcint_interfacer;
    auto raw_container = libcint_interfacer.convert(this->shell_set);
    return libcint_interfacer.calculateOneElectronIntegrals<1>(cint1e_nuc_cart, raw_container)[0];
}


/**
 *  @param origin       the origin of the dipole
 *
 *  @return the matrix representation of the Cartesian components of the electrical dipole operator in this AO basis, using the libcint integral engine
 */
std::array<OneElectronOperator<double>, 3> AOBasis::calculateLibcintDipoleIntegrals(const Vector<double, 3>& origin) const {

    const LibcintInterfacer libcint_interfacer;
    auto raw_container = libcint_interfacer.convert(this->shell_set);
    const auto& all_integrals = libcint_interfacer.calculateOneElectronIntegrals<3>(cint1e_r_cart, raw_container);

    // Apply the minus sign which comes from the charge of the electrons -e
    return std::array<OneElectronOperator<double>, 3> {-all_integrals[0], -all_integrals[1], -all_integrals[2]};
}


/**
 *  @return the matrix representation of the Coulomb repulsion operator in this AO basis, using the libcint integral engine
 */
TwoElectronOperator<double> AOBasis::calculateLibcintCoulombRepulsionIntegrals() const {

    const LibcintInterfacer libcint_interfacer;
    auto raw_container = libcint_interfacer.convert(this->shell_set);
    return libcint_interfacer.calculateTwoElectronIntegrals(cint2e_cart, raw_container);
}


}  // namespace GQCP
