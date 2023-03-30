// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/CoulombRepulsionOperator.hpp"
#include "Operator/FirstQuantized/ElectronicDipoleOperator.hpp"
#include "Operator/FirstQuantized/KineticOperator.hpp"
#include "Operator/FirstQuantized/NuclearAttractionOperator.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"

#include <functional>


extern "C" {


#include <cint.h>


/*
 *  The following functions are not inside the include header <cint.h>, so we should define them in order to be able to use them in our source code
 */
FINT cint1e_kin_cart(double* buf, const int* shls, const int* atm, int natm, const int* bas, int nbas, const double* env);   // kinetic energy
FINT cint1e_nuc_cart(double* buf, const int* shls, const int* atm, int natm, const int* bas, int nbas, const double* env);   // nuclear attraction energy
FINT cint1e_ovlp_cart(double* buf, const int* shls, const int* atm, int natm, const int* bas, int nbas, const double* env);  // overlap
FINT cint1e_r_cart(double* buf, const int* shls, const int* atm, int natm, const int* bas, int nbas, const double* env);     // dipole integrals


}  // extern "C"


namespace GQCP {


/*
 *  Aliases for functions that are used in Libcint
 */

using Libcint1eFunction = std::function<int(double*, const int*, const int*, int, const int*, int, const double*)>;
using Libcint2eFunction = std::function<int(double*, const int*, const int*, int, const int*, int, const double*, const CINTOpt*)>;
using Libcint2eOptimizerFunction = std::function<void(CINTOpt**, const int*, const int, const int*, const int, const double*)>;


/*
 *  Forward declarations
 */

namespace libcint {


class RawContainer;


}  // namespace libcint


/**
 *  A class that takes care of the interfacing with the libcint library
 */
class LibcintInterfacer {
public:
    // PUBLIC METHODS - INTERFACING

    /**
     *  @param shell_set        the GQCP::ShellSet whose information should be converted
     *
     *  @return the information in a GQCP::ShellSet as a libcint::RawContainer
     */
    libcint::RawContainer convert(const ShellSet<GTOShell>& shell_set) const;

    /**
     *  Set the origin for the calculation of all vector-related integrals
     *
     *  @param raw_container        the libcint::RawContainer that holds the data needed by libcint
     *  @param origin               the new origin for the calculation of all vector-related integrals
     */
    void setCommonOrigin(libcint::RawContainer& raw_container, const Vector<double, 3>& origin) const;


    //  PUBLIC METHODS - INTEGRAL FUNCTIONS

    /**
     *  @param op               the electronic electric dipole operator
     *
     *  @return the Libcint one-electron function that corresponds to the electronic electric dipole operator
     */
    Libcint1eFunction oneElectronFunction(const ElectronicDipoleOperator& op) const { return cint1e_r_cart; }

    /**
     *  @param op               the kinetic operator
     *
     *  @return the Libcint one-electron function that corresponds to the kinetic operator
     */
    Libcint1eFunction oneElectronFunction(const KineticOperator& op) const { return cint1e_kin_cart; }

    /**
     *  @param op               the nuclear attraction operator
     *
     *  @return the Libcint one-electron function that corresponds to the nuclear attraction operator
     */
    Libcint1eFunction oneElectronFunction(const NuclearAttractionOperator& op) const { return cint1e_nuc_cart; }

    /**
     *  @param op           the overlap operator
     *
     *  @return the Libcint one-electron function that corresponds to the overlap operator
     */
    Libcint1eFunction oneElectronFunction(const OverlapOperator& op) const { return cint1e_ovlp_cart; }

    /**
     *  @param op               the Coulomb repulsion operator
     *
     *  @return the Libcint two-electron function that corresponds to the Coulomb repulsion dipole operator
     */
    Libcint2eFunction twoElectronFunction(const CoulombRepulsionOperator& op) const { return cint2e_cart; }

    /**
     *  @param op               the Coulomb repulsion operator
     *
     *  @return the Libcint two-electron optimizer function that corresponds to the Coulomb repulsion dipole operator
     */
    Libcint2eOptimizerFunction twoElectronOptimizerFunction(const CoulombRepulsionOperator& op) const { return cint2e_cart_optimizer; }
};


namespace libcint {


/**
 *  C++ global variables instead of macro commands
 */
static constexpr int ptr_common_orig = PTR_COMMON_ORIG;  // an offset for the origin of vector operators
static constexpr int ptr_env_start = PTR_ENV_START;      // an offset such that libcint can retrieve the correct index inside the environment, starts at 20

static constexpr int atm_slots = ATM_SLOTS;  // the number of 'slots' (i.e. 'members') for an atom
static constexpr int charge_of = CHARGE_OF;  // slot offset for atomic charge
static constexpr int ptr_coord = PTR_COORD;  // slot offset for a 'pointer' of the atom inside the libcint environment


static constexpr int bas_slots = BAS_SLOTS;  // the number of 'slots' (i.e. 'members') for a basis function
static constexpr int atom_of = ATOM_OF;      // slot offset for the corresponding atom in the basis function
static constexpr int ang_of = ANG_OF;        // slot offset for the angular momentum
static constexpr int nprim_of = NPRIM_OF;    // slot offset for the number of primitives inside the basis function
static constexpr int nctr_of = NCTR_OF;      // slot offset for the number of contractions
static constexpr int ptr_exp = PTR_EXP;      // slot offset for a 'pointer' to the exponents of the shell inside the libcint environment
static constexpr int ptr_coeff = PTR_COEFF;  // slot offset for a 'pointer' to the contraction coefficients inside the libcint environment


/**
 *  A wrapper that owns raw libcint 'atm', 'bas' and 'env' arrays
 *
 *  @note The members are named after the libcint variables.
 */
class RawContainer {
private:
    int natm;  // number of atoms
    int nbf;   // number of basis functions
    int nsh;   // the number of shells

    int* libcint_atm;     // information about the atoms
    int* libcint_bas;     // information about the basis functions
    double* libcint_env;  // a raw block of doubles in which libcint (probably) places intermediary calculations


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Allocate memory for the raw libcint arrays
     *
     *  @param natm         the number of atoms
     *  @param nbf          the number of basis functions
     *  @param nsh          the number of shells
     */
    RawContainer(const size_t natm, const size_t nbf, const size_t nsh) :
        natm {static_cast<int>(natm)},
        nbf {static_cast<int>(nbf)},
        nsh {static_cast<int>(nsh)},
        libcint_atm {new int[this->natm * atm_slots]},
        libcint_bas {new int[this->nbf * bas_slots]},
        libcint_env {new double[10000]} {}


    /*
     *  DESTRUCTOR
     */

    /**
     *  Deallocate the memory used by the raw libcint arrays
     */
    ~RawContainer() {
        delete[] this->libcint_atm;
        delete[] this->libcint_bas;
        delete[] this->libcint_env;
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return (a pointer to) the atomic data in this libcint container
     */
    const int* atmData() const { return this->libcint_atm; }

    /**
     *  @return (a pointer to) the basis data in this libcint container
     */
    const int* basData() const { return this->libcint_bas; }

    /**
     *  @return (a pointer to) the environment data in this libcint container
     */
    const double* envData() const { return this->libcint_env; }

    /**
     *  @return the number of atoms in this libcint container
     */
    int numberOfAtoms() const { return this->natm; }

    /**
     *  @return the number of basis functions in this libcint container
     */
    int numberOfBasisFunctions() const { return this->nbf; }

    /**
     *  @return the number of shells in this libcint container
     */
    int numberOfShells() const { return this->nsh; }


    /*
     *  FRIENDS
     */
    friend class GQCP::LibcintInterfacer;
};


}  // namespace libcint
}  // namespace GQCP
