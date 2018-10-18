#ifndef GQCG_RDMBUILDER_HPP
#define GQCG_RDMBUILDER_HPP

#include "RDM/BaseRDMBuilder.hpp"
#include "FockSpace/FockSpace.hpp"
#include "FockSpace/FockSpaceProduct.hpp"

#include <memory>


namespace GQCG {


/**
 *  RDMBuilder is a wrapper around the RDMBuilders that provides the functionality of the correct RDMBuilder
 *  for a given Fock space at compile- or runtime.
 */


class RDMBuilder {
private:
    std::shared_ptr<GQCG::BaseRDMBuilder> rdm_builder;
public:
    // CONSTRUCTOR
    /**
     *  Allocates a DOCIRDMBuilder based on @param fock_space
     */
    RDMBuilder(const FockSpace& fock_space);

    /**
     *  Allocates a FCIRDMBuilder based on @param fock_space
     */
    RDMBuilder(const FockSpaceProduct& fock_space);

    /**
     *  Allocates the correct derived BaseRDMBuilder based on @param fock_space
     */
    RDMBuilder(const BaseFockSpace& fock_space);


    // PUBLIC METHODS
    /**
     *  @return all 1-RDMs from a coefficient vector @param x
     */
    OneRDMs calculate1RDMs(const Eigen::VectorXd &x) { return rdm_builder->calculate1RDMs(x); }

    /**
     *  @return all 2-RDMs from a coefficient vector @param x
     */
    TwoRDMs calculate2RDMs(const Eigen::VectorXd &x) { return rdm_builder->calculate2RDMs(x); }
};


}  // namespace GQCG


#endif  // GQCG_RDMBUILDER_HPP
