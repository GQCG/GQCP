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

#include "Basis/ScalarBasis/CartesianGTO.hpp"

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
LibintInterfacer& LibintInterfacer::get() {      // need to return by reference since we deleted the relevant constructor
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
    libint2::Shell libint_shell {libint_alpha, {libint_contraction}, libint_O};
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

    std::vector<bool> pure_flags(libint_basisset.size());  // create an array with the given size
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
        Eigen::Map<const Eigen::Matrix<double, 3, 1>> libint_origin_map(libint_shell_copy.O.data());  // convert raw array data to Eigen
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
    for (auto& contraction : libint_shell.contr) {
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


}  // namespace GQCP
