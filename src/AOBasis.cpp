#include "AOBasis.hpp"
#include "LibintCommunicator.hpp"

namespace GQCP {


AOBasis::AOBasis(const GQCP::Molecule& molecule, std::string basis_set) :
    atoms (molecule.atoms),
    basis_functions (libint2::BasisSet(std::move(basis_set), GQCP::LibintCommunicator::get().interface(this->atoms))),  // construct a libint2::BasisSet
    number_of_basis_functions (static_cast<size_t>(this->basis_functions.nbf()))
{}


}  // namespace GQCP
