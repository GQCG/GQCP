// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#ifndef GQCP_LIBINTCOMMUNICATOR_HPP
#define GQCP_LIBINTCOMMUNICATOR_HPP


#include "AOBasis.hpp"
#include "Molecule.hpp"
#include "Operator/OneElectronOperator.hpp"
#include "Operator/TwoElectronOperator.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <libint2.hpp>


namespace GQCP {


/**
 *  A singleton class that takes care of interfacing with the Libint2 (version >2.2.0) C++ API
 *
 *  Singleton class template from (https://stackoverflow.com/a/1008289)
 */
class LibintCommunicator {
private:
    // SINGLETON METHODS
    /**
     *  Private constructor as required by the singleton class design
     */
    LibintCommunicator();

    /**
     *  Private destructor as required by the singleton class design
     */
    ~LibintCommunicator();


    // PRIVATE STRUCTS
    typedef struct {} empty;  // empty_pod is a private typedef for libint2::Engine, so we copy it over


    // PRIVATE METHODS
    /**
     *  @tparam N               the number of operator components
     *  @tparam Parameters      the type of the given parameters for the integral engine
     *
     *  @param operator_type    the name of the operator as specified by the enumeration
     *  @param basisset         the libint2 basis set representing the AO basis
     *
     *  @return an array of N OneElectronOperators corresponding to the matrix representations of the N components of the given operator type
     */
    template <size_t N, typename Parameters>
    std::array<GQCP::OneElectronOperator, N> calculateOneElectronIntegrals(libint2::Operator operator_type, const libint2::BasisSet& basisset, const Parameters& parameters = empty()) const;

    /**
     *  @return the TwoElectronOperator corresponding to the matrix representation of @param operator_type in the given
     *  @param ao_basis
     */
    GQCP::TwoElectronOperator calculateTwoElectronIntegrals(libint2::Operator operator_type, const GQCP::AOBasis& ao_basis) const;


public:
    // SINGLETON METHODS
    /**
     *  @return the static singleton instance
     */
    static LibintCommunicator& get();

    /**
     *  Remove the public copy constructor and a public assignment operator
     */
    LibintCommunicator(LibintCommunicator const& libint_communicator) = delete;
    void operator=(LibintCommunicator const& libint_communicator) = delete;


    // PUBLIC METHODS
    GQCP::OneElectronOperator calculateOverlapIntegrals(const GQCP::AOBasis& ao_basis) const;
    GQCP::OneElectronOperator calculateKineticIntegrals(const GQCP::AOBasis& ao_basis) const;
    GQCP::OneElectronOperator calculateNuclearIntegrals(const GQCP::AOBasis& ao_basis) const;

    GQCP::TwoElectronOperator calculateCoulombRepulsionIntegrals(const GQCP::AOBasis& ao_basis) const;

    /**
     *  @return a std::vector<libint2::Atom> based on a given std::vector<GQCP::Atom> @param atoms
     */
    std::vector<libint2::Atom> interface(const std::vector<GQCP::Atom>& atoms) const;
};



/*
 *  TEMPLATE IMPLEMENTATIONS
 */
/**
 *  @tparam N               the number of operator components
 *  @tparam Parameters      the type of the given parameters for the integral engine
 *
 *  @param operator_type    the name of the operator as specified by the enumeration
 *  @param basisset         the libint2 basis set representing the AO basis
 *  @param parameters       the parameters for the integral engine
 *
 *  @return an array of N OneElectronOperators corresponding to the matrix representations of the N components of the given operator type
 */
template <size_t N, typename Parameters>
std::array<GQCP::OneElectronOperator, N> LibintCommunicator::calculateOneElectronIntegrals(libint2::Operator operator_type, const libint2::BasisSet& basisset, const Parameters& parameters) const {

    // Use the basis_functions that is currently a libint2::BasisSet
    const auto nbf = static_cast<size_t>(basisset.nbf());  // nbf: number of basis functions in the basisset


    // Initialize the N components of the operator
    std::array<Eigen::MatrixXd, N> matrix_components;
    for (auto& matrix : matrix_components) {
        matrix = Eigen::MatrixXd::Zero(nbf, nbf);
    }


    // Construct the libint2 engine
    libint2::Engine engine (operator_type, basisset.max_nprim(), static_cast<int>(basisset.max_l()));

    // Set parameters for the nuclear attraction integrals
//    if (operator_type == libint2::Operator::nuclear) {
//        auto atoms = this->interface(ao_basis.get_atoms());  // convert from GQCP::Atoms to libint2::atoms
//        engine.set_params(make_point_charges(atoms));
//    }
    std::cout << "Got to before set_params." << std::endl;
    engine.set_params(parameters);



    std::array<const double*, N> calculated_integrals_components;


    const auto shell2bf = basisset.shell2bf();  // maps shell index to bf index

    const auto& buffer = engine.results();  // vector that holds pointers to computed shell sets
    // actually, buffer.size() is always 1, so buffer[0] is a pointer to
    //      the first calculated integral of these specific shells
    // the values that buffer[0] points to will change after every compute() call

    // ASSERT buffer.size() == N ??

    // One-electron integrals are between two basis functions, so we'll need two loops
    // Libint calculates integrals between libint2::Shells, so we will loop over the shells (sh) in the basisset
    const auto nsh = static_cast<size_t>(basisset.size());  // nsh: number of shells in the basisset
    for (auto sh1 = 0; sh1 != nsh; ++sh1) {  // sh1: shell 1
        for (auto sh2 = 0; sh2 != nsh; ++sh2) {  // sh2: shell 2
            // Calculate integrals between the two shells (basis_set is a decorated std::vector<libint2::Shell>)
            engine.compute(basisset[sh1], basisset[sh2]);


            for (size_t i = 0; i < N; i++) {
                calculated_integrals_components[i] = buffer[i];
            }
//            auto calculated_integrals = buffer[0];  // is actually a pointer: const double *


//            if (calculated_integrals == nullptr) {  // if the zeroth element is nullptr, then the whole shell has been exhausted
//                // or the libint engine predicts that the integrals are below a certain threshold
//                // in this case the value does not need to be filled in, and we are safe because we have properly initialized to zero
//                continue;
//            }

            // Extract the calculated integrals from calculated_integrals
            // In calculated_integrals, the integrals are stored in row major form
            auto bf1 = shell2bf[sh1];  // (index of) first bf in sh1
            auto bf2 = shell2bf[sh2];  // (index of) first bf in sh2

            auto nbf_sh1 = basisset[sh1].size();  // number of basis functions in first shell
            auto nbf_sh2 = basisset[sh2].size();  // number of basis functions in second shell

            for (auto f1 = 0; f1 != nbf_sh1; ++f1) {  // f1: index of basis function within shell 1
                for (auto f2 = 0; f2 != nbf_sh2; ++f2) { // f2: index of basis function within shell 2
                    for (size_t i = 0; i < N; i++) {
                        double computed_integral = calculated_integrals_components[i][f2 + f1 * nbf_sh2];  // integrals are packed in row-major form
                        matrix_components[i](bf1 + f1, bf2 + f2) = computed_integral;
                    }

//                    double computed_integral = calculated_integrals[f2 + f1 * nbf_sh2];  // integrals are packed in row-major form
//                    matrix(bf1 + f1, bf2 + f2) = computed_integral;
                }
            }  // data access loops

        }
    }  // shell loops


    std::array<GQCP::OneElectronOperator, N> operator_components;
    for (size_t i = 0; i < N; i++) {
        operator_components[i] = GQCP::OneElectronOperator(matrix_components[i]);  // convert to a GQCP::OneElectronOperator
    }

    return operator_components;
}



}  // namespace GQCP


#endif  // GQCP_LIBINTCOMMUNICATOR_HPP
