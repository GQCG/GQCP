#ifndef GQCP_BASEHAMILTONIANPARAMETERS_HPP
#define GQCP_BASEHAMILTONIANPARAMETERS_HPP


#include <memory>

#include "AOBasis.hpp"


namespace GQCP {


class BaseHamiltonianParameters {
protected:
    std::shared_ptr<GQCP::AOBasis> ao_basis;  // the initial atomic orbitals

public:
    // CONSTRUCTOR
    /**
     *  Constructor based on a given @param ao_basis
     */
    explicit BaseHamiltonianParameters(std::shared_ptr<GQCP::AOBasis> ao_basis);


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseHamiltonianParameters() = 0;
};


}  // namespace GQCP


#endif  // GQCP_BASEHAMILTONIANPARAMETERS_HPP
