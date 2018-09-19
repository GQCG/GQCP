#ifndef GQCG_BASEHAMILTONIANPARAMETERS_HPP
#define GQCG_BASEHAMILTONIANPARAMETERS_HPP


#include <memory>

#include "AOBasis.hpp"


namespace GQCG {


class BaseHamiltonianParameters {
protected:
    std::shared_ptr<GQCG::AOBasis> ao_basis_ptr;

public:
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseHamiltonianParameters() = 0;
};


}  // namespace GQCG


#endif  // GQCG_BASEHAMILTONIANPARAMETERS_HPP
