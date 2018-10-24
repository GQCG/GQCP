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
#include "LibintCommunicator.hpp"

#include <iostream>
#include <sstream>



namespace GQCP {


/*
 *  PRIVATE METHODS
 */

/**
 *  Private constructor as required by the singleton class design
 */
LibintCommunicator::LibintCommunicator() {
    libint2::initialize();
}


/**
 *  Private destructor as required by the singleton class design
 */
LibintCommunicator::~LibintCommunicator() {
    libint2::finalize();
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the static singleton instance
 */
LibintCommunicator& LibintCommunicator::get() {  // need to return by reference since we deleted the relevant constructor
    static LibintCommunicator singleton_instance;  // instantiated on first use and guaranteed to be destroyed
    return singleton_instance;
}


/**
 *  @return a std::vector<libint2::Atom> based on a given std::vector<GQCP::Atom> @param atoms
 */
std::vector<libint2::Atom> LibintCommunicator::interface(const std::vector<GQCP::Atom>& atoms) const {

    std::vector<libint2::Atom> libint_vector;  // start with an empty vector, we're doing push_backs later

    for (const auto& atom : atoms) {
        libint2::Atom libint_atom {static_cast<int>(atom.atomic_number), atom.x, atom.y, atom.z};
        libint_vector.push_back(libint_atom);
    }

    return libint_vector;
}


/**
 *  @return the OneElectronOperator corresponding to the matrix representation of @param operator_type in the given
 *  @param ao_basis
 */
GQCP::OneElectronOperator LibintCommunicator::calculateOneElectronIntegrals(libint2::Operator operator_type, const GQCP::AOBasis& ao_basis) const {

    // Use the basis_functions that is currently a libint2::BasisSet
    auto libint_basisset = ao_basis.get_basis_functions();
    const auto nbf = static_cast<size_t>(libint_basisset.nbf());  // nbf: number of basis functions in the basisset

    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(nbf, nbf);


    // Construct the libint2 engine
    libint2::Engine engine (operator_type, libint_basisset.max_nprim(), static_cast<int>(libint_basisset.max_l()));

    // Something extra for the nuclear attraction integrals
    if (operator_type == libint2::Operator::nuclear) {
        auto atoms = this->interface(ao_basis.get_atoms());  // convert from GQCP::Atoms to libint2::atoms
        engine.set_params(make_point_charges(atoms));
    }


    const auto shell2bf = libint_basisset.shell2bf();  // maps shell index to bf index

    const auto& buffer = engine.results();  // vector that holds pointers to computed shell sets
    // actually, buffer.size() is always 1, so buffer[0] is a pointer to
    //      the first calculated integral of these specific shells
    // the values that buffer[0] points to will change after every compute() call


    // One-electron integrals are between two basis functions, so we'll need two loops
    // Libint calculates integrals between libint2::Shells, so we will loop over the shells (sh) in the basisset
    const auto nsh = static_cast<size_t>(libint_basisset.size());  // nsh: number of shells in the basisset
    for (auto sh1 = 0; sh1 != nsh; ++sh1) {  // sh1: shell 1
        for (auto sh2 = 0; sh2 != nsh; ++sh2) {  // sh2: shell 2
            // Calculate integrals between the two shells (basis_set is a decorated std::vector<libint2::Shell>)
            engine.compute(libint_basisset[sh1], libint_basisset[sh2]);

            auto calculated_integrals = buffer[0];  // is actually a pointer: const double *

            if (calculated_integrals == nullptr) {  // if the zeroth element is nullptr, then the whole shell has been exhausted
                // or the libint engine predicts that the integrals are below a certain threshold
                // in this case the value does not need to be filled in, and we are safe because we have properly initialized to zero
                continue;
            }

            // Extract the calculated integrals from calculated_integrals
            // In calculated_integrals, the integrals are stored in row major form
            auto bf1 = shell2bf[sh1];  // (index of) first bf in sh1
            auto bf2 = shell2bf[sh2];  // (index of) first bf in sh2

            auto nbf_sh1 = libint_basisset[sh1].size();  // number of basis functions in first shell
            auto nbf_sh2 = libint_basisset[sh2].size();  // number of basis functions in second shell

            for (auto f1 = 0; f1 != nbf_sh1; ++f1) {  // f1: index of basis function within shell 1
                for (auto f2 = 0; f2 != nbf_sh2; ++f2) { // f2: index of basis function within shell 2
                    double computed_integral = calculated_integrals[f2 + f1 * nbf_sh2];  // integrals are packed in row-major form
                    matrix(bf1 + f1, bf2 + f2) = computed_integral;
                }
            }  // data access loops

        }
    }  // shell loops

    return GQCP::OneElectronOperator(matrix);
}


/**
 *  @return the TwoElectronOperator corresponding to the matrix representation of @param operator_type in the given
 *  @param ao_basis
 */
GQCP::TwoElectronOperator LibintCommunicator::calculateTwoElectronIntegrals(libint2::Operator operator_type, const GQCP::AOBasis& ao_basis) const {

    // Use the basis_functions that is currently a libint2::BasisSet
    auto libint_basisset = ao_basis.get_basis_functions();
    const auto nbf = static_cast<size_t>(libint_basisset.nbf());  // nbf: number of basis functions in the basisset


    // Initialize the rank-4 two-electron integrals tensor and set to zero
    Eigen::Tensor<double, 4> tensor (nbf, nbf, nbf, nbf);
    tensor.setZero();


    // Construct the libint2 engine
    libint2::Engine engine(libint2::Operator::coulomb, libint_basisset.max_nprim(), static_cast<int>(libint_basisset.max_l()));  // libint2 requires an int

    const auto shell2bf = libint_basisset.shell2bf();  // maps shell index to bf index

    const auto &buffer = engine.results();  // vector that holds pointers to computed shell sets
    // actually, buffer.size() is always 1, so buffer[0] is a pointer to
    //      the first calculated integral of these specific shells
    // the values that buffer[0] points to will change after every compute() call


    // Two-electron integrals are between four basis functions, so we'll need four loops
    // Libint calculates integrals between libint2::Shells, so we will loop over the shells (sh) in the basisset
    const auto nsh = static_cast<size_t>(libint_basisset.size());  // nsh: number of shells in the basisset
    for (auto sh1 = 0; sh1 != nsh; ++sh1) {  // sh1: shell 1
        for (auto sh2 = 0; sh2 != nsh; ++sh2) {  // sh2: shell 2
            for (auto sh3 = 0; sh3 != nsh; ++sh3) {  // sh3: shell 3
                for (auto sh4 = 0; sh4 != nsh; ++sh4) {  //sh4: shell 4
                    // Calculate integrals between the two shells (obs is a decorated std::vector<libint2::Shell>)
                    engine.compute(libint_basisset[sh1], libint_basisset[sh2], libint_basisset[sh3], libint_basisset[sh4]);

                    auto calculated_integrals = buffer[0];

                    if (calculated_integrals == nullptr) {  // if the zeroth element is nullptr, then the whole shell has been exhausted
                        // or the libint engine predicts that the integrals are below a certain threshold
                        // in this case the value does not need to be filled in, and we are safe because we have properly initialized to zero
                        continue;
                    }

                    // Extract the calculated integrals from calculated_integrals.
                    // In calculated_integrals, the integrals are stored in row major form.
                    auto bf1 = static_cast<long>(shell2bf[sh1]);  // (index of) first bf in sh1
                    auto bf2 = static_cast<long>(shell2bf[sh2]);  // (index of) first bf in sh2
                    auto bf3 = static_cast<long>(shell2bf[sh3]);  // (index of) first bf in sh3
                    auto bf4 = static_cast<long>(shell2bf[sh4]);  // (index of) first bf in sh4


                    auto nbf_sh1 = static_cast<long>(libint_basisset[sh1].size());  // number of basis functions in first shell
                    auto nbf_sh2 = static_cast<long>(libint_basisset[sh2].size());  // number of basis functions in second shell
                    auto nbf_sh3 = static_cast<long>(libint_basisset[sh3].size());  // number of basis functions in third shell
                    auto nbf_sh4 = static_cast<long>(libint_basisset[sh4].size());  // number of basis functions in fourth shell

                    for (auto f1 = 0L; f1 != nbf_sh1; ++f1) {
                        for (auto f2 = 0L; f2 != nbf_sh2; ++f2) {
                            for (auto f3 = 0L; f3 != nbf_sh3; ++f3) {
                                for (auto f4 = 0L; f4 != nbf_sh4; ++f4) {
                                    auto computed_integral = calculated_integrals[f4 + nbf_sh4 * (f3 + nbf_sh3 * (f2 + nbf_sh2 * (f1)))];  // integrals are packed in row-major form

                                    // Two-electron integrals are given in CHEMIST'S notation: (11|22)
                                    tensor(f1 + bf1, f2 + bf2, f3 + bf3, f4 + bf4) = computed_integral;
                                }
                            }
                        }
                    } // data access loops

                }
            }
        }
    } // shell loops

    return GQCP::TwoElectronOperator(tensor);
};


}  // namespace GQCP
