#include "properties.hpp"


namespace GQCP {


/**
 *  Calculate the electronic energy as a result of the contraction of the 1- and 2-RDMs with the one- and two-electron integrals
 *
 *  @param ham_par      the Hamiltonian parameters containing the one- and two-electron integrals
 *  @param one_rdm      the 1-RDM
 *  @param two_rdm      the 2-RDM
 *
 *  @return the electronic energy
 */
double calculateElectronicEnergy(const GQCP::HamiltonianParameters& ham_par, const GQCP::OneRDM& one_rdm, const GQCP::TwoRDM& two_rdm) {

    auto h = ham_par.get_h().get_matrix_representation();
    auto g = ham_par.get_g().get_matrix_representation();

    auto D = one_rdm.get_matrix_representation();
    auto d = two_rdm.get_matrix_representation();


    double energy_by_contraction = (h * D).trace();

    // Specify the contractions for the relevant contraction of the two-electron integrals and the 2-RDM
    //      0.5 g(p q r s) d(p q r s)
    Eigen::array<Eigen::IndexPair<int>, 4> contractions = {Eigen::IndexPair<int>(0,0), Eigen::IndexPair<int>(1,1), Eigen::IndexPair<int>(2,2), Eigen::IndexPair<int>(3,3)};
    //      Perform the contraction
    Eigen::Tensor<double, 0> contraction = 0.5 * g.contract(d, contractions);

    // As the contraction is a scalar (a tensor of rank 0), we should access by (0).
    energy_by_contraction += contraction(0);

    return energy_by_contraction;
}


/**
 *  Calculate the electronic electric dipole moment
 *
 *  @param dipole_operator      the three components of the Cartesian dipole integrals in the orthonormal basis in which the 1-RDM is expressed
 *  @param one_rdm              the 1-RDM
 *
 *  @return the three Cartesian components of the electronic electric dipole moment
 */
Eigen::Vector3d calculateElectronicDipoleMoment(const std::array<GQCP::OneElectronOperator, 3>& dipole_operator, const GQCP::OneRDM& one_rdm) {

    Eigen::Vector3d electronic_dipole = Eigen::Vector3d::Zero();

    for (size_t i = 0; i < 3; i++) {

        Eigen::MatrixXd dipole_component = dipole_operator[i].get_matrix_representation();
        Eigen::MatrixXd D = one_rdm.get_matrix_representation();

        electronic_dipole(i) = (dipole_component * D).trace();
    }

    return electronic_dipole;
}


}  // namespace GQCP
