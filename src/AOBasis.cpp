#include "AOBasis.hpp"
#include "LibintCommunicator.hpp"

namespace GQCG {


AOBasis::AOBasis(const GQCG::Molecule& molecule, std::string basis_set) :
    atoms (molecule.atoms),
    basis_functions (libint2::BasisSet(std::move(basis_set), GQCG::LibintCommunicator::get().interface(this->atoms)))
{
    //
}


}  // namespace GQCG
