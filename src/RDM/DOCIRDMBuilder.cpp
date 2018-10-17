#include "RDM/DOCIRDMBuilder.hpp"


namespace GQCG {


/*
 *  CONSTRUCTOR
 */
DOCIRDMBuilder::DOCIRDMBuilder(const FockSpace& fock_space) :
    fock_space (fock_space)
{}


/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @return 1RDM from a coefficient vector @param x
 */
OneRDMs DOCIRDMBuilder::calculate1RDMs(const Eigen::VectorXd& x) {
    // The formulas for the DOCI 1-RDMs can be found in (https://github.com/lelemmen/electronic_structure)

    size_t K = this->fock_space.get_K();
    size_t dim = this->fock_space.get_dimension();

    // For DOCI, one rdm covers both spins
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(K, K);

    // Create the first ONV (with address 0). In DOCI, the Fock space for alpha and beta is equal so we just use one
    ONV onv = this->fock_space.get_ONV(0);   


    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the spin strings

        for (size_t e1 = 0; e1 < this->fock_space.get_N(); e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupied_index(e1);  // retrieve the index of the orbital the electron occupies
            double c_I = x(I);  // coefficient of the I-th basis vector
            D(p,p) += 2*std::pow(c_I, 2);
        }
        
        if (I < dim-1) {
            this->fock_space.setNext(onv);
        }
    }

    OneRDM one_rdm (D);
    return OneRDMs (one_rdm);
}

/**
 *  @return 2RDM from a coefficient vector @param x
 */
TwoRDMs DOCIRDMBuilder::calculate2RDMs(const Eigen::VectorXd& x) {

    // The formulas for the DOCI 2-RDMs can be found in (https://github.com/lelemmen/electronic_structure)

    size_t K = this->fock_space.get_K();
    size_t dim = this->fock_space.get_dimension();
    
    Eigen::Tensor<double, 4> d_aaaa (K, K, K, K);
    d_aaaa.setZero();
    Eigen::Tensor<double, 4> d_aabb (K, K, K, K);
    d_aabb.setZero();

    // Create the first ONV (with address 0). In DOCI, the Fock space for alpha and beta is equal so we just use one
    ONV onv = this->fock_space.get_ONV(0);   


    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the spin strings
        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            if (onv.annihilate(p)) {  // if p is occupied in I

                double c_I = x(I);  // coefficient of the I-th basis vector
                double c_I_2 = std::pow(c_I, 2);  // square of c_I

                d_aabb(p,p,p,p) += c_I_2;

                for (size_t q = 0; q < p; q++) {  // q loops over SOs with an index smaller than p
                    if (onv.create(q)) {  // if q is not occupied in I
                        size_t J = this->fock_space.getAddress(onv);  // the address of the coupling string
                        double c_J = x(J);  // coefficient of the J-th basis vector

                        d_aabb(p,q,p,q) += c_I * c_J;
                        d_aabb(q,p,q,p) += c_I * c_J;  // since we're looping for q < p

                        onv.annihilate(q);  // reset the spin string after previous creation on q
                    }

                    else {  // if q is occupied in I
                        d_aaaa(p,p,q,q) += c_I_2;
                        d_aaaa(q,q,p,p) += c_I_2;  // since we're looping for q < p

                        d_aaaa(p,q,q,p) -= c_I_2;
                        d_aaaa(q,p,p,q) -= c_I_2;  // since we're looping for q < p

                        d_aabb(p,p,q,q) += c_I_2;
                        d_aabb(q,q,p,p) += c_I_2;  // since we're looping for q < p
                    }
                }
                onv.create(p);  // reset the spin string after previous annihilation on p
            }
        }

        if ( I < dim-1) {
            this->fock_space.setNext(onv);
        }
    }

    // For DOCI, we have additional symmetries (two_rdm_aaaa = two_rdm_bbbb, two_rdm_aabb = two_rdm_bbaa)
    TwoRDM two_rdm_aaaa (d_aaaa);
    TwoRDM two_rdm_aabb (d_aabb);
    TwoRDM two_rdm_bbaa (d_aabb);
    TwoRDM two_rdm_bbbb (d_aaaa);
    return TwoRDMs (two_rdm_aaaa, two_rdm_aabb, two_rdm_bbaa, two_rdm_bbbb);
}



}  // namespace GQCG
