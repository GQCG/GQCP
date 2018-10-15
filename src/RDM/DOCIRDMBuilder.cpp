#include "RDM/DOCIRDMBuilder.cpp"


namespace GQCG {


/*
 *  CONSTRUCTOR
 */
DOCIRDMBuilder::DOCIRDMBuilder(FockSpace fock_space) :
    fock_space(fock_space)
{}


/*
 *  OVERRIDEN PUBLIC METHODS
 */

/**
 *  @return 1RDM from a coefficient vector @param x
 */
OneRDM DOCIRDMBuilder::construct1RDM(const Eigen::VectorXd& x) {


}

/**
 *  @return 2RDM from a coefficient vector @param x
 */
TwoRDM DOCIRDMBuilder::construct2RDM(const Eigen::VectorXd& x) override;



};


}  // namespace GQCG
