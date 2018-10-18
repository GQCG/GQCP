#ifndef GQCG_BASERDMBUILDER_HPP
#define GQCG_BASERDMBUILDER_HPP


#include "RDM/OneRDM.hpp"
#include "RDM/TwoRDM.hpp"
#include "RDM/RDMs.hpp"
#include "FockSpace/BaseFockSpace.hpp"


namespace GQCG {


/**
 *  BaseRDMBuilder is an abstract base class for the calculation of a density matrix from a given wave function
 *  or coefficient expansion in a Fock space
 */
class BaseRDMBuilder {
public:
    // CONSTRUCTOR
    BaseRDMBuilder() = default;


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseRDMBuilder() = 0;


    // PURE VIRTUAL PUBLIC METHODS
    /**
     *  @return all 1-RDMs from a coefficient vector @param x
     */
    virtual OneRDMs calculate1RDMs(const Eigen::VectorXd& x) = 0;

    /**
     *  @return all 2-RDMs from a coefficient vector @param x
     */
    virtual TwoRDMs calculate2RDMs(const Eigen::VectorXd& x) = 0;

    /**
     *  @return the Fock space of the RDMBuilder
     */
    virtual BaseFockSpace* get_fock_space() = 0;
};


}  // namespace GQCG


#endif  // GQCG_BASERDMBUILDER_HPP
