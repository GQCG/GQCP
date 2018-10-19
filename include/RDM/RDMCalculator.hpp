#ifndef GQCP_RDMCALCULATOR_HPP
#define GQCP_RDMCALCULATOR_HPP

#include "RDM/BaseRDMBuilder.hpp"
#include "FockSpace/FockSpace.hpp"
#include "FockSpace/FockSpaceProduct.hpp"

#include <memory>


namespace GQCP {


/**
 *  RDMCalculator is a wrapper around the derived RDMBuilders that provides the functionality of the correct derived RDMBuilder
 *  for a given Fock space at compile- or runtime.
 */
class RDMCalculator {
private:
    std::shared_ptr<GQCP::BaseRDMBuilder> rdm_builder;

public:
    // CONSTRUCTOR
    /**
     *  Allocates a derived RDMBuilder based on the nature of the given @param fock_space
     */
    RDMCalculator(const FockSpace& fock_space);  // DOCIRDMBuilder
    RDMCalculator(const FockSpaceProduct& fock_space);  // FCIRDMBuilder
    RDMCalculator(const BaseFockSpace& fock_space);


    // PUBLIC METHODS
    /**
     *  @return all 1-RDMs from a coefficient vector @param x
     */
    OneRDMs calculate1RDMs(const Eigen::VectorXd& x);

    /**
     *  @return all 2-RDMs from a coefficient vector @param x
     */
    TwoRDMs calculate2RDMs(const Eigen::VectorXd& x);
};


}  // namespace GQCP


#endif  // GQCP_RDMCALCULATOR_HPP
