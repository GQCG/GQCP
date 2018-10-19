#include "RDM/RDMCalculator.hpp"

#include "RDM/DOCIRDMBuilder.hpp"
#include "RDM/FCIRDMBuilder.hpp"


namespace GQCG {

/*
 *  CONSTRUCTOR
 */

/**
 *  Allocates a DOCIRDMBuilder based on @param fock_space
 */
RDMCalculator::RDMCalculator(const FockSpace& fock_space) {
    rdm_builder = std::make_shared<GQCG::DOCIRDMBuilder>(fock_space);
}


/**
 *  Allocates a FCIRDMBuilder based on @param fock_space
 */
RDMCalculator::RDMCalculator(const FockSpaceProduct& fock_space) {
    rdm_builder = std::make_shared<GQCG::FCIRDMBuilder>(fock_space);
}


/**
 *  Allocates the correct derived BaseRDMBuilder based on @param fock_space
 */
RDMCalculator::RDMCalculator(const BaseFockSpace& fock_space) {

    switch (fock_space.get_type()){

        case FockSpaceType::FockSpace: {
            rdm_builder = std::make_shared<GQCG::DOCIRDMBuilder>(dynamic_cast<const GQCG::FockSpace&>(fock_space));

            break;
        }

        case FockSpaceType::FockSpaceProduct: {
            rdm_builder = std::make_shared<GQCG::FCIRDMBuilder>(dynamic_cast<const GQCG::FockSpaceProduct&>(fock_space));

            break;
        }
    }

}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return all 1-RDMs from a coefficient vector @param x
 */
OneRDMs RDMCalculator::calculate1RDMs(const Eigen::VectorXd& x) {
    return rdm_builder->calculate1RDMs(x);
}


/**
 *  @return all 2-RDMs from a coefficient vector @param x
 */
TwoRDMs RDMCalculator::calculate2RDMs(const Eigen::VectorXd& x) {
    return rdm_builder->calculate2RDMs(x);
}

}  // namespace GQCG