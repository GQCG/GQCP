// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
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
#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"

#include "Basis/CartesianGTO.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <iostream>
#include <sstream>



namespace GQCP {



/*
 *  PRIVATE METHODS - SINGLETON
 */

/**
 *  Private constructor as required by the singleton class design
 */
LibintInterfacer::LibintInterfacer() {

    // Set up the libint2 environment
    libint2::initialize();


    // By default, libint2 changes the contraction coefficients in a shell to correspond to a normalized basis function: undo this behavior
    // Note: libint2 uses normalization factors such that the axis-aligned (i.e. {2,0,0} or {0,3,0}) Cartesian functions are normalized
    libint2::Shell::do_enforce_unit_normalization(false);
}


/**
 *  Private destructor as required by the singleton class design
 */
LibintInterfacer::~LibintInterfacer() {

    // Finish the libint2 environment
    libint2::finalize();
}



/*
 *  PUBLIC METHODS - SINGLETON
 */

/**
 *  @return the static singleton instance
 */
LibintInterfacer& LibintInterfacer::get() {  // need to return by reference since we deleted the relevant constructor
    static LibintInterfacer singleton_instance;  // instantiated on first use and guaranteed to be destroyed
    return singleton_instance;
}



/*
 *  PUBLIC METHODS - INTERFACING (GQCP TO LIBINT)
 */

/**
 *  @param nucleus          the GQCP-nucleus that should be interfaced into a libint2::Atom
 *
 *  @return a libint2::Atom, interfaced from the given GQCP::Nucleus
 */
libint2::Atom LibintInterfacer::interface(const Nucleus& nucleus) const {

    libint2::Atom libint_atom {static_cast<int>(nucleus.charge()), nucleus.position().x(), nucleus.position().y(), nucleus.position().z()};

    return libint_atom;
}


/**
 *  @param nuclei           the GQCP-nuclei that should be interfaced
 *
 *  @return libint2::Atoms, interfaced from the given GQCP nuclei
 */
std::vector<libint2::Atom> LibintInterfacer::interface(const std::vector<Nucleus>& nuclei) const {

    std::vector<libint2::Atom> libint_vector;  // start with an empty vector, we're doing push_backs later
    libint_vector.reserve(nuclei.size());

    for (const auto& nucleus : nuclei) {
        libint_vector.push_back(this->interface(nucleus));
    }

    return libint_vector;
}


/**
 *  @param shell        the GQCP shell that should be interfaced
 *
 *  @return a libint2::Shell whose renorm()alization has been undone, interfaced from the GQCP Shell
 */
libint2::Shell LibintInterfacer::interface(const GTOShell& shell) const {

    // Part 1: exponents
    const std::vector<double>& libint_alpha = shell.get_gaussian_exponents();  // libint::Shell::real_t is double, so no need to use real_t


    // Part 2: contractions
    auto libint_l = static_cast<int>(shell.get_l());
    bool libint_pure = shell.is_pure();
    const std::vector<double>& libint_coeff = shell.get_contraction_coefficients();
    libint2::Shell::Contraction libint_contraction {libint_l, libint_pure, libint_coeff};


    // Part 3: origin
    const auto& position = shell.get_nucleus().position();
    std::array<double, 3> libint_O {position.x(), position.y(), position.z()};


    // Upon construction, libint2 renorm()alizes the contraction coefficients, so we want to undo this
    libint2::Shell libint_shell (libint_alpha, {libint_contraction}, libint_O);
    this->undo_renorm(libint_shell);
    return libint_shell;
}


/**
 *  @param shellset     the GQCP ShellSet that should be interfaced
 *
 *  @return a libint2::BasisSet (whose underlying libint2::Shells have been re-renorm()alized), interfaced from the GQCP ShellSet. Note that it is not possible to create libint2-sp-shells from a GQCP ShellSet
 */
libint2::BasisSet LibintInterfacer::interface(const ShellSet<GTOShell>& shellset) const {

    libint2::BasisSet libint_basisset;  // start with an empty vector, we're doing push_backs later
    libint_basisset.reserve(shellset.numberOfShells());

    for (const auto& shell : shellset.asVector()) {
        libint_basisset.push_back(this->interface(shell));
    }


    // At this point in the code, the libint2::BasisSet is 'uninitialized', i.e. its private member _nbf is -1, etc.
    // Therefore, we are using a hack to force the private libint2::BasisSet::init() being called

    std::vector<bool> pure_flags (libint_basisset.size());
    for (size_t i = 0; i < libint_basisset.size(); i++) {
        pure_flags[i] = libint_basisset[i].contr[0].pure;  // the libint2::BasisSet that was created can never have more than one libint2::Shell::Contraction (we split SP-shells)
    }

    // Call libint2::BasisSet::init() through set_pure, and re-set the 'pure' flags
    libint_basisset.set_pure(false);
    for (size_t i = 0; i < libint_basisset.size(); i++) {
        libint_basisset[i].contr[0].pure = pure_flags[i];
    }

    return libint_basisset;
}


/*
 *  PUBLIC METHODS - INTERFACING (LIBINT TO GQCP)
 */

/**
 *  Interface a libint2::Shell to the corresponding list of GQCP::Shells. Note that there is no one-to-one libint -> GQCP conversion, since GQCP does not support representing 'linked' sp-'shells'
 *
 *  @param libint_shell     the libint2 Shell that should be interfaced
 *  @param nuclei           the nuclei that can serve as centers of the Shells
 *  @param undo_renorm      if the libint2::Shell should be un-renorm()alized
 *
 *  @return a vector of GQCP::Shells
 */
std::vector<GTOShell> LibintInterfacer::interface(const libint2::Shell& libint_shell, const std::vector<Nucleus>& nuclei, bool undo_renorm) const {

    // If asked for, undo Libint2's default renorm()alization
    auto libint_shell_copy = libint_shell;
    if (undo_renorm) {
        this->undo_renorm(libint_shell_copy);
    }


    // Construct the corresponding GQCP::Shells
    std::vector<double> exponents = libint_shell_copy.alpha;

    std::vector<GTOShell> shells;
    shells.reserve(this->numberOfShells(libint_shell_copy));
    for (const auto& libint_contraction : libint_shell_copy.contr) {

        // Angular momentum and coefficients
        size_t l = libint_contraction.l;
        std::vector<double> coefficients = libint_contraction.coeff;

        // Libint2 only stores the origin of the shell, so we have to find the nucleus corresponding to the copied shell's origin
        Eigen::Map<const Eigen::Matrix<double, 3, 1>> libint_origin_map (libint_shell_copy.O.data());  // convert raw array data to Eigen
        Nucleus corresponding_nucleus;
        for (size_t i = 0; i < nuclei.size(); i++) {
            Nucleus nucleus = nuclei[i];

            if (nucleus.position().isApprox(libint_origin_map, 1.0e-06)) {  // tolerant comparison
                corresponding_nucleus = nucleus;
                break;
            }

            if (i == (nuclei.size() - 1)) {  // if we haven't broken out of the loop after exhausting the possible nuclei
                throw std::invalid_argument("LibintInterfacer::interface(libint2::Shell, std::vector<Nucleus>): No given nucleus matches the center of the libint2::Shell");
            }
        }

        // Other flags
        bool pure = libint_contraction.pure;
        shells.emplace_back(l, corresponding_nucleus, exponents, coefficients, pure, false, false);
    }

    return shells;
}


/**
 *  Interface a libint2::BasisSet to the corresponding GQCP::ShellSet and undo the libint2 renorm()alization
 *
 *  @param libint_basisset      the libint2 Shell that should be interfaced
 *  @param nuclei               the nuclei that can serve as centers of the Shells
 *
 *  @return a vector of GTOShells corresponding to the un-renorm()alized libint2::BasisSet
 */
std::vector<GTOShell> LibintInterfacer::interface(const libint2::BasisSet& libint_basisset, const std::vector<Nucleus>& nuclei) const {

    std::vector<GTOShell> shells;
    shells.reserve(this->numberOfShells(libint_basisset));
    for (const auto& libint_shell : libint_basisset) {
        for (const auto& shell : this->interface(libint_shell, nuclei)) {
            shells.push_back(shell);
        }
    }

    return shells;
}



/*
 *  PUBLIC METHODS - OTHER LIBINT2-RELATED FUNCTIONS
 */

/**
 *  @param libint_shell         the libint2::Shell
 *
 *  @return the number of true shells that are contained in the libint shell
 */
size_t LibintInterfacer::numberOfShells(const libint2::Shell& libint_shell) const {
    return libint_shell.ncontr();
}


/**
 *  @param libint_basisset      the libint2::BasisSet
 *
 *  @return the number of true shells that are contained in the libint2::BasisSet
 */
size_t LibintInterfacer::numberOfShells(const libint2::BasisSet& libint_basisset) const {

    size_t nsh {};  // number of shells

    for (const auto& libint_shell : libint_basisset) {
        nsh += this->numberOfShells(libint_shell);
    }

    return nsh;
}


/**
 *  Undo the libint2 default renormalization (see libint2::Shell::renorm())
 *
 *  @param libint_shell         the shell that should be un-renorm()alized
 */
void LibintInterfacer::undo_renorm(libint2::Shell& libint_shell) const {

    // Instead of multiplying (what libint2 does), divide each contraction coefficient by the normalization factor
    for (auto& contraction: libint_shell.contr) {
        for (size_t p = 0; p != libint_shell.nprim(); p++) {
            double alpha = libint_shell.alpha[p];
            size_t l = contraction.l;

            double N = CartesianGTO::calculateNormalizationFactor(alpha, CartesianExponents(l, 0, 0));

            contraction.coeff[p] /= N;
        }
    }
}


/*
 *  PUBLIC METHODS - ENGINES
 */


/**
 *  Construct a libint2 engine that corresponds to the given operator
 * 
 *  @param op               the overlap operator
 *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
 *  @param max_l            the maximum angular momentum of Gaussian shell
 * 
 *  @return the proper libint2 engine
 */
libint2::Engine LibintInterfacer::createEngine(const OverlapOperator& op, const size_t max_nprim, const size_t max_l) const {

    return libint2::Engine(libint2::Operator::overlap, max_nprim, static_cast<int>(max_l));
}


/**
 *  Construct a libint2 engine that corresponds to the given operator
 * 
 *  @param op               the kinetic operator
 *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
 *  @param max_l            the maximum angular momentum of Gaussian shell
 * 
 *  @return the proper libint2 engine
 */
libint2::Engine LibintInterfacer::createEngine(const KineticOperator& op, const size_t max_nprim, const size_t max_l) const {

    return libint2::Engine(libint2::Operator::kinetic, max_nprim, static_cast<int>(max_l));
}


/**
 *  Construct a libint2 engine that corresponds to the given operator
 * 
 *  @param op               the nuclear attraction operator
 *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
 *  @param max_l            the maximum angular momentum of Gaussian shell
 * 
 *  @return the proper libint2 engine
 */
libint2::Engine LibintInterfacer::createEngine(const NuclearAttractionOperator& op, const size_t max_nprim, const size_t max_l) const {

    return libint2::Engine(libint2::Operator::nuclear, max_nprim, static_cast<int>(max_l));
}


/**
 *  Construct a libint2 engine that corresponds to the given operator
 * 
 *  @param op               the electronic electric dipole operator
 *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
 *  @param max_l            the maximum angular momentum of Gaussian shell
 * 
 *  @return the proper libint2 engine
 */
libint2::Engine LibintInterfacer::createEngine(const ElectronicDipoleOperator& op, const size_t max_nprim, const size_t max_l) const {

    return libint2::Engine(libint2::Operator::emultipole1, max_nprim, static_cast<int>(max_l));
}


/**
 *  Construct a libint2 engine that corresponds to the given operator
 * 
 *  @param op               the Coulomb repulsion operator
 *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
 *  @param max_l            the maximum angular momentum of Gaussian shell
 * 
 *  @return the proper libint2 engine
 */
libint2::Engine LibintInterfacer::createEngine(const CoulombRepulsionOperator& op, const size_t max_nprim, const size_t max_l) const {

    return libint2::Engine(libint2::Operator::coulomb, max_nprim, static_cast<int>(max_l));
}



/*
 *  PUBLIC METHODS - INTEGRALS
 */

/**
 *  @param operator_type        the name of the operator as specified by the enumeration
 *  @param libint_basisset      the libint2 basis set representing the AO basis
 *
 *  @return the matrix representation of a two-electron operator in the given AO basis
 */
TwoElectronOperator<double> LibintInterfacer::calculateTwoElectronIntegrals(libint2::Operator operator_type, const libint2::BasisSet& libint_basisset) const {

    auto nbf = static_cast<size_t>(libint_basisset.nbf());  // nbf: number of basis functions in the basisset

    // Initialize the rank-4 two-electron integrals tensor and set to zero
    TwoElectronOperator<double> g (nbf);
    g.setZero();


    // Construct the libint2 engine
    libint2::Engine engine(libint2::Operator::coulomb, libint_basisset.max_nprim(), static_cast<int>(libint_basisset.max_l()));  // libint2 requires an int

    const auto& shell2bf = libint_basisset.shell2bf();  // maps shell index to bf index

    const auto& buffer = engine.results();  // vector that holds pointers to computed shell sets
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

                    const auto& calculated_integrals = buffer[0];

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
                                    const auto& computed_integral = calculated_integrals[f4 + nbf_sh4 * (f3 + nbf_sh3 * (f2 + nbf_sh2 * (f1)))];  // integrals are packed in row-major form

                                    // Two-electron integrals are given in CHEMIST'S notation: (11|22)
                                    g(f1 + bf1, f2 + bf2, f3 + bf3, f4 + bf4) = computed_integral;
                                }
                            }
                        }
                    } // data access loops

                }
            }
        }
    } // shell loops

    return g;
};


}  // namespace GQCP
