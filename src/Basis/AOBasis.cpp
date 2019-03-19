#include "Basis/AOBasis.hpp"
#include "LibintCommunicator.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Construct an AO basis by placing shells shells corresponding to the basisset information on every atom of the molecule
 *
 *  @param molecule             the molecule containing the atoms on which the shells should be centered
 *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
 */
AOBasis::AOBasis(const Molecule& molecule, const std::string& basisset_name) :
    shell_set (ShellSet(molecule, basisset_name))
{}


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

    auto libint_basisset = LibintCommunicator::get().interface(this->shell_set);
    return LibintCommunicator::get().calculateOneElectronIntegrals<1>(libint2::Operator::overlap, libint_basisset)[0];
}


/**
 *  @return the matrix representation of the kinetic energy operator in this AO basis
 */
OneElectronOperator<double> AOBasis::calculateKineticIntegrals() const {

    auto libint_basisset = LibintCommunicator::get().interface(this->shell_set);
    return LibintCommunicator::get().calculateOneElectronIntegrals<1>(libint2::Operator::kinetic, libint_basisset)[0];
}


/**
 *  @return the matrix representation of the nuclear attraction operator in this AO basis
 */
OneElectronOperator<double> AOBasis::calculateNuclearIntegrals() const {

    auto libint_basisset = LibintCommunicator::get().interface(this->shell_set);
    auto libint_atoms = LibintCommunicator::get().interface(this->shell_set.atoms());

    return LibintCommunicator::get().calculateOneElectronIntegrals<1>(libint2::Operator::nuclear, libint_basisset, make_point_charges(libint_atoms))[0];
}


/**
 *  @return the matrix representation of the Cartesian components of the electrical dipole operator
 */
std::array<OneElectronOperator<double>, 3> AOBasis::calculateDipoleIntegrals(const Vector<double, 3>& origin) const {

    std::array<double, 3> origin_array {origin.x(), origin.y(), origin.z()};
    auto libint_basisset = LibintCommunicator::get().interface(this->shell_set);

    auto all_integrals = LibintCommunicator::get().calculateOneElectronIntegrals<4>(libint2::Operator::emultipole1, libint_basisset, origin_array);  // overlap, x, y, z

    // Apply the minus sign which comes from the charge of the electrons -e
    return std::array<OneElectronOperator<double>, 3> {-all_integrals[1], -all_integrals[2], -all_integrals[3]};  // we don't need the overlap, so ignore [0]
}


/**
 *  @return the matrix representation of the Coulomb repulsion operator in this AO basis
 */
TwoElectronOperator<double> AOBasis::calculateCoulombRepulsionIntegrals() const {

    auto libint_basisset = LibintCommunicator::get().interface(this->shell_set);
    return LibintCommunicator::get().calculateTwoElectronIntegrals(libint2::Operator::coulomb, libint_basisset);
}


}  // namespace GQCP
