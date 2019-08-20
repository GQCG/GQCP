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
#ifndef LibcintInterfacer_hpp
#define LibcintInterfacer_hpp

#include "Basis/GTOShell.hpp"
#include "Basis/ShellSet.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Operator/OneElectronOperator.hpp"
#include "Operator/TwoElectronOperator.hpp"

#include <functional>
#include <unordered_map>



extern "C" {


#include <cint.h>


/*
 *  The following functions are not inside the include header <cint.h>, so we should define them in order to be able to use them in our source code
 */
FINT cint1e_ovlp_cart(double* buf, const int* shls, const int* atm, int natm, const int* bas, int nbas, const double* env);  // overlap
FINT cint1e_kin_cart(double* buf, const int* shls, const int* atm, int natm, const int* bas, int nbas, const double* env);  // kinetic energy
FINT cint1e_nuc_cart(double* buf, const int* shls, const int* atm, int natm, const int* bas, int nbas, const double* env);  // nuclear attraction energy
FINT cint1e_r_cart(double* buf, const int* shls, const int* atm, int natm, const int* bas, int nbas, const double* env);  // dipole integrals


}  // extern "C"



namespace GQCP {


class LibcintInterfacer;  // forward declaration before friending


namespace libcint {


/**
 *  C++ global variables instead of macro commands
 */
static constexpr int ptr_common_orig = PTR_COMMON_ORIG;  // an offset for the origin of vector operators
static constexpr int ptr_env_start = PTR_ENV_START;  // an offset such that libcint can retrieve the correct index inside the environment, starts at 20

static constexpr int atm_slots = ATM_SLOTS;  // the number of 'slots' (i.e. 'members') for an atom
static constexpr int charge_of = CHARGE_OF;  // slot offset for atomic charge
static constexpr int ptr_coord = PTR_COORD;  // slot offset for a 'pointer' of the atom inside the libcint environment


static constexpr int bas_slots = BAS_SLOTS;  // the number of 'slots' (i.e. 'members') for a basis function
static constexpr int atom_of = ATOM_OF;  // slot offset for the corresponding atom in the basis function
static constexpr int ang_of = ANG_OF;  // slot offset for the angular momentum
static constexpr int nprim_of = NPRIM_OF;  // slot offset for the number of primitives inside the basis function
static constexpr int nctr_of = NCTR_OF;  // slot offset for the number of contractions
static constexpr int ptr_exp = PTR_EXP;  // slot offset for a 'pointer' to the exponents of the shell inside the libcint environment
static constexpr int ptr_coeff = PTR_COEFF;  // slot offset for a 'pointer' to the contraction coefficients inside the libcint environment




/**
 *  A wrapper that owns raw libcint 'atm', 'bas' and 'env' arrays
 * 
 *  @note There is no easy way to ask a RawContainer the corresponding shell index given a GQCP::GTOShell. This is why we are explicitly holding such a map.
 */
class RawContainer {
private:
    int natm;  // number of atoms
    int nbf;  // number of basis functions
    int nsh;  // the number of shells

    int* libcint_atm;  // information about the atoms
    int* libcint_bas;  // information about the basis functions
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
        natm (static_cast<int>(natm)),
        nbf (static_cast<int>(nbf)),
        nsh (static_cast<int>(nsh)),
        libcint_atm (new int[this->natm * atm_slots]),
        libcint_bas (new int[this->nbf * bas_slots]),
        libcint_env (new double[10000])
    {}


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
     *  GETTERS
     */
    int numberOfAtoms() const { return this->natm; }
    int numberOfBasisFunctions() const { return this->nbf; }
    int numberOfShells() const { return this->nsh; }
    const int* atmData() const { return this->libcint_atm; }
    const int* basData() const { return this->libcint_bas; }
    const double* envData() const { return this->libcint_env; }

    /*
     *  FRIENDS
     */
    friend class GQCP::LibcintInterfacer;
};


}  // namespace libcint





using Libcint1eFunction = std::function<int (double*, const int*, const int*, int, const int*, int, const double*)>;
using Libcint2eFunction = std::function<int (double*, const int*, const int*, int, const int*, int, const double*, const CINTOpt*)>;




/**
 *  A class that takes care of the interfacing with the libcint library
 */
class LibcintInterfacer {
public:

    /*
     *  PUBLIC METHODS - INTERFACING
     */

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

    /**
     *  @param op           the overlap operator
     * 
     *  @return the Libcint one-electron function that corresponds to the overlap operator
     */
    Libcint1eFunction oneElectronFunction(const OverlapOperator& op) const;

    /**
     *  @param op               the kinetic operator
     * 
     *  @return the Libcint one-electron function that corresponds to the kinetic operator
     */
    Libcint1eFunction oneElectronFunction(const KineticOperator& op) const;

    /**
     *  @param op               the nuclear attraction operator
     * 
     *  @return the Libcint one-electron function that corresponds to the nuclear attraction operator
     */
    Libcint1eFunction oneElectronFunction(const NuclearAttractionOperator& op) const;

    /**
     *  @param op               the electronic electric dipole operator
     * 
     *  @return the Libcint one-electron function that corresponds to the electronic electric dipole operator
     */
    Libcint1eFunction oneElectronFunction(const ElectronicDipoleOperator& op) const;

    /**
     *  @param op               the Coulomb repulsion operator
     * 
     *  @return the Libcint two-electron function that corresponds to the Coulomb repulsion dipole operator
     */
    Libcint2eFunction twoElectronFunction(const CoulombRepulsionOperator& op) const;


    /*
     *  PUBLIC METHODS - INTEGRALS
     */

    /**
     *  @tparam N                   the number of libcint operator components
     *
     *  @param function             the libcint one-electron integral function
     *  @param raw_container        the data libcint needs to perform calculations
     *
     *  @return an array of N OneElectronOperators corresponding to the matrix representations of the N components of the given operator represented by the libcint function
     */
    template <size_t N>
    std::array<OneElectronOperator<double>, N> calculateOneElectronIntegrals(const Libcint1eFunction& function, libcint::RawContainer& raw_container) const {

        // Initialize the components to zero
        const auto& nbf = raw_container.nbf;  // number of basis functions
        std::array<OneElectronOperator<double>, N> operator_components;
        for (auto& op : operator_components) {
            op = OneElectronOperator<double>::Zero(nbf, nbf);
        }


        // Calculate the integrals over the shells
        size_t bf1 = 0;  // index of the first basis function in the first shell
        for (size_t sh1 = 0; sh1 < raw_container.nsh; sh1++) {

            int nbf_sh1 = CINTcgto_cart(static_cast<int>(sh1), raw_container.libcint_bas);  // number of basis functions in first shell

            size_t bf2 = 0;  // index of the first basis function in the second shell
            for (size_t sh2 = 0; sh2 < raw_container.nsh; sh2++) {

                int nbf_sh2 = CINTcgto_cart(static_cast<int>(sh2), raw_container.libcint_bas);  // number of basis functions in second shell


                int shell_indices[2];  // indices of the shells to be calculated over
                shell_indices[0] = static_cast<int>(sh1);
                shell_indices[1] = static_cast<int>(sh2);


                double buf[N * nbf_sh1 * nbf_sh2];  // buffer where the integrals are calculated to
                function(buf, shell_indices, raw_container.libcint_atm, raw_container.natm, raw_container.libcint_bas, raw_container.nbf, raw_container.libcint_env);  // TODO: what is result is zero: skip a shell?


                // Place the calculated integrals in the components of the Operator
                for (size_t f1 = 0; f1 < nbf_sh1; f1++) {
                    for (size_t f2 = 0; f2 < nbf_sh2; f2++) {
                        for (size_t i = 0; i < N; i++) {
                            const double& computed_integral = buf[f1 + nbf_sh1 * (f2 + nbf_sh2 * i)];  // integrals are packed in column-major form
                            operator_components[i](bf1 + f1, bf2 + f2) = computed_integral;
                        }
                    }
                }  // data access loops


                bf2 += nbf_sh2;  // update the 'first basis function' with the encountered number of basis functions
            }

            bf1 += nbf_sh1;  // update the 'first basis function' with the encountered number of basis functions
        }  // shell loops


        return operator_components;
    }


    /**
     *  @param function             the libcint two-electron integral function
     *  @param raw_container        the data libcint needs to perform calculations
     *
     *  @return an array of N OneElectronOperators corresponding to the matrix representations of the N components of the given operator represented by the libcint function
     */
    TwoElectronOperator<double> calculateTwoElectronIntegrals(const Libcint2eFunction& function, libcint::RawContainer& raw_container) const;
};


}  // namespace GQCP


#endif  /* LibcintInterfacer_hpp */
