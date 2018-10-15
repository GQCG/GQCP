#ifndef GQCG_RDMBUILDER_HPP
#define GQCG_RDMBUILDER_HPP

#include "RDM/OneRDM.hpp"
#include "RDM/TwoRDM.hpp"

namespace GQCG {


/**
 *  RDMBuilder is an abstract base class for the calculation of a density matrix from a given wave function
 *  or coefficient expansion in a Fock space
 */
class RDMBuilder {
public:
    // CONSTRUCTOR
    RDMBuilder() = default;


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~RDMBuilder() = 0;


    // PURE VIRTUAL PUBLIC METHODS
    /**
     *  @return 1RDM from a coefficient vector @param x
     */
    virtual OneRDM construct1RDM(const Eigen::VectorXd& x) = 0;

    /**
     *  @return 2RDM from a coefficient vector @param x
     */
    virtual TwoRDM construct2RDM (const Eigen::VectorXd& x) = 0;

    /**
     *  @return the Fock space of the RDMBuilder
     */
    virtual BaseFockSpace* get_fock_space() = 0;
};


}  // namespace GQCG


#endif  // GQCG_RDMBUILDER_HPP
