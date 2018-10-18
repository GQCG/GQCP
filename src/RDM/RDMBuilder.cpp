#include "RDM/RDMBuilder.hpp"

#include "RDM/DOCIRDMBuilder.hpp"
#include "RDM/FCIRDMBuilder.hpp"


namespace GQCG {

/*
 *  CONSTRUCTOR
 */

/**
 *  Allocates a DOCIRDMBuilder based on @param fock_space
 */
RDMBuilder::RDMBuilder(const FockSpace& fock_space) {
    rdm_builder = std::make_shared<GQCG::DOCIRDMBuilder>(fock_space);
}


/**
 *  Allocates a FCIRDMBuilder based on @param fock_space
 */
RDMBuilder::RDMBuilder(const FockSpaceProduct& fock_space) {
    rdm_builder = std::make_shared<GQCG::FCIRDMBuilder>(fock_space);
}


/**
 *  Allocates the correct derived BaseRDMBuilder based on @param fock_space
 */
RDMBuilder::RDMBuilder(const BaseFockSpace& fock_space) {

    FockSpaceType fock_space_type = fock_space.get_fock_space_type();

    switch (fock_space_type){

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

}  // namespace GQCG