#include "Basis/Basis.hpp"
#include "LibintCommunicator.hpp"


namespace GQCP {



    /*
     *  PUBLIC METHODS - LIBINT INTEGRALS
     */

    /**
     *  @return the matrix representation of the overlap operator in this AO basis
     */
    OneElectronOperator<double> Basis::calculateLibintOverlapIntegrals() const {

        auto libint_basisset = LibintCommunicator::get().interface(this->basisset);
        return LibintCommunicator::get().calculateOneElectronIntegrals<1>(libint2::Operator::overlap, libint_basisset)[0];
    }



    /*
     *  PUBLIC METHODS - LIBCINT INTEGRALS
     */


}  // namespace GQCP
