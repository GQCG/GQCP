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
#pragma once


#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"

#include <boost/preprocessor.hpp>  // include preprocessor before libint to fix libint-boost bug

#include <libint2.hpp>


namespace GQCP {


/**
 *  A singleton class that takes care of interfacing with the Libint2 (version 2.3.1) C++ API
 *
 *  Singleton class template from: https://stackoverflow.com/a/1008289
 */
class LibintInterfacer {
public:
    using libint_target_ptr_vec = libint2::Engine::target_ptr_vec;


private:
    // PRIVATE METHODS - SINGLETON
    /**
     *  Private constructor as required by the singleton class design
     */
    LibintInterfacer();

    /**
     *  Private destructor as required by the singleton class design
     */
    ~LibintInterfacer();


    // PRIVATE STRUCTS
    typedef struct {
    } empty;  // empty_pod is a private typedef for libint2::Engine, so we copy it over


public:
    // PUBLIC METHODS - SINGLETON

    /**
     *  @return the static singleton instance
     */
    static LibintInterfacer& get();

    /**
     *  Remove the public copy constructor and the public assignment operator
     */
    LibintInterfacer(LibintInterfacer const& libint_communicator) = delete;
    void operator=(LibintInterfacer const& libint_communicator) = delete;


    // PUBLIC METHODS - INTERFACING (GQCP TO LIBINT)

    /**
     *  @param nucleus          the GQCP-nucleus that should be interfaced into a libint2::Atom
     *
     *  @return a libint2::Atom, interfaced from the given GQCP::Nucleus
     */
    libint2::Atom interface(const Nucleus& nucleus) const;

    /**
     *  @param nuclei           the GQCP-nuclei that should be interfaced
     *
     *  @return libint2::Atoms, interfaced from the given GQCP nuclei
     */
    std::vector<libint2::Atom> interface(const std::vector<Nucleus>& nuclei) const;

    /**
     *  @param shell        the GQCP shell that should be interfaced
     *
     *  @return a libint2::Shell whose renorm()alization has been undone, interfaced from the GQCP GTOShell
     */
    libint2::Shell interface(const GTOShell& shell) const;

    /**
     *  @param shellset     the GQCP ShellSet that should be interfaced
     *
     *  @return a libint2::BasisSet (whose underlying libint2::Shells have been re-renorm()alized), interfaced from the GQCP ShellSet. Note that it is not possible to create libint2-sp-shells from a GQCP ShellSet
     */
    libint2::BasisSet interface(const ShellSet<GTOShell>& shellset) const;


    // PUBLIC METHODS - INTERFACING (LIBINT TO GQCP)
    /**
     *  Interface a libint2::Shell to the corresponding list of GQCP::Shells. Note that there is no one-to-one libint -> GQCP conversion, since GQCP does not support representing 'linked' sp-'shells'
     *
     *  @param libint_shell     the libint2 Shell that should be interfaced
     *  @param nuclei           the nuclei that can serve as centers of the Shells
     *  @param undo_renorm      if the libint2::Shell should be un-renorm()alized
     *
     *  @return a vector of GQCP::Shells corresponding to the given libint2::Shells
     */
    std::vector<GTOShell> interface(const libint2::Shell& libint_shell, const std::vector<Nucleus>& nuclei, bool undo_renorm = true) const;

    /**
     *  Interface a libint2::BasisSet to the corresponding GQCP::ShellSet and undo the libint2 renorm()alization
     *
     *  @param libint_basisset      the libint2 Shell that should be interfaced
     *  @param nuclei               the nuclei that can serve as centers of the Shells
     *
     *  @return a vector of GTOShells corresponding to the un-renorm()alized libint2::BasisSet
     */
    std::vector<GTOShell> interface(const libint2::BasisSet& libint_basisset, const std::vector<Nucleus>& nuclei) const;


    // PUBLIC METHODS - OTHER LIBINT2-RELATED FUNCTIONS
    /**
     *  @param libint_shell         the libint2::Shell
     *
     *  @return the number of true shells that are contained in the libint2::Shell
     */
    size_t numberOfShells(const libint2::Shell& libint_shell) const;

    /**
     *  @param libint_basisset      the libint2::BasisSet
     *
     *  @return the number of true shells that are contained in the libint2::BasisSet
     */
    size_t numberOfShells(const libint2::BasisSet& libint_basisset) const;

    /**
     *  Undo the libint2 default renormalization (see libint2::Shell::renorm())
     *
     *  @param libint_shell         the shell that should be un-renorm()alized
     */
    void undo_renorm(libint2::Shell& libint_shell) const;


    // PUBLIC METHODS - ENGINES

    /**
     *  Construct a libint2 engine that corresponds to the given operator
     * 
     *  @param op               the overlap operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return the proper libint2 engine
     */
    libint2::Engine createEngine(const OverlapOperator& op, const size_t max_nprim, const size_t max_l) const;

    /**
     *  Construct a libint2 engine that corresponds to the given operator
     * 
     *  @param op               the kinetic operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return the proper libint2 engine
     */
    libint2::Engine createEngine(const KineticOperator& op, const size_t max_nprim, const size_t max_l) const;

    /**
     *  Construct a libint2 engine that corresponds to the given operator
     * 
     *  @param op               the nuclear attraction operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return the proper libint2 engine
     */
    libint2::Engine createEngine(const NuclearAttractionOperator& op, const size_t max_nprim, const size_t max_l) const;

    /**
     *  Construct a libint2 engine that corresponds to the given operator
     * 
     *  @param op               the electronic electric dipole operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return the proper libint2 engine
     */
    libint2::Engine createEngine(const ElectronicDipoleOperator& op, const size_t max_nprim, const size_t max_l) const;

    /**
     *  Construct a libint2 engine that corresponds to the given operator
     * 
     *  @param op               the Coulomb repulsion operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return the proper libint2 engine
     */
    libint2::Engine createEngine(const CoulombRepulsionOperator& op, const size_t max_nprim, const size_t max_l) const;
};


}  // namespace GQCP
