#ifndef LibcintInterfacer_hpp
#define LibcintInterfacer_hpp

#include "Operator/OneElectronOperator.hpp"
#include "Operator/TwoElectronOperator.hpp"
#include "Molecule.hpp"
#include "Basis/ShellSet.hpp"

#include <functional>



extern "C" {

#include <cint.h>


/*
 *  FUNCTIONS THAT AREN'T INSIDE <cint.h>
 */
int cint1e_ovlp_cart(double* buf, int* shls, int* atm, int natm, int* bas, int nbas, double* env);  // overlap
int cint1e_kin_cart(double* buf, int* shls, int* atm, int natm, int* bas, int nbas, double* env);  // kinetic energy
int cint1e_nuc_cart(double* buf, int* shls, int* atm, int natm, int* bas, int nbas, double* env);  // nuclear attraction energy
int cint1e_r_cart(double* buf, int* shls, int* atm, int natm, int* bas, int nbas, double* env);  // dipole integrals, origin at zero


}  // extern "C"




namespace GQCP {



/*
 *  WRAPPERS AROUND C-STYLE LIBCINT ARRAYS
 */
namespace libcint {


    // https://github.com/sunqm/libcint/blob/master/examples/python_call.py

/**
 *  A struct wrapper around the libcint "atm" C array
 */
struct Atom {
public:
    int charge_of;  // nuclear charge/atomic number
    int ptr_coord;  // offset of the x-coordinate of the atom inside the libcint environment
};



/**
 *  A struct wrapper around the libcint "bas" C array
 */
struct BasisFunction {
public:
    int atom_of;  // index of the corresponding atom
    int ang_of;  // angular momentum
    int nprim_of;  // number of primitives
    int nctr_of;  // number of contractions
    int ptr_exp;  // offset of the exponents inside the libcint environment
    int ptr_coeff;  // offset to the contraction coefficients inside the libcint environment
};


/**
 *  A struct wrapper around the libcint "env" C array
 */
struct Environment {
public:

};


}  // namespace libcint





using Libint1eFunction = std::function<int (double*, int*, int*, int, int*, int, double*)>;
using Libint2eFunction = std::function<int (double*, int*, int*, int, int*, int, double*, CINTOpt*)>;




/**
 *  A class that takes care of the interfacing with the libcint library
 */
class LibcintInterfacer {
public:

    template <size_t N>
    std::array<OneElectronOperator<double>, N> calculateOneElectronIntegrals(const Libint1eFunction& function) const {

        /* general contracted DZ basis [3s1p/2s1p] for H2
         exponents    contract-coeff
         S   6.0          0.7               0.4
             2.0          0.6               0.3
             0.8          0.5               0.2
         P   0.9          1.
         */

        // Specify the example molecule
        Atom h1 (1,  0.0, 0.0, 0.8);  // coordinates in bohr
        Atom h2 (1,  0.0, 0.0, -0.8);


        // Put STO-3G on every H
        std::vector<double> gaussian_exponents {3.42525091, 0.623913730, 0.168855400};
        std::vector<double> contraction_coefficients {0.154328970, 0.535328140, 0.444634540};
        Shell shell1 {0, h1, gaussian_exponents, contraction_coefficients};
        Shell shell2 {0, h2, gaussian_exponents, contraction_coefficients};

        ShellSet shell_set {shell1, shell2};
        auto atoms = shell_set.atoms();


        int natm = static_cast<int>(atoms.size());
        int nbf = 2;  // number of basis functions


        // ATM_SLOTS = 6, BAS_SLOTS = 8 are declared inside <cint.h>

        // TODO: use std::vector and pass &[0] or .data()?
        int libcint_atm[natm * ATM_SLOTS];  // information about the atoms
        int libcint_bas[nbf * BAS_SLOTS];  // information about the basis functions
        double libcint_env[10000];  // a general container (env = environment) in which libcint (probably) places intermediary calculations


        /*
         *  ATM CONFIGURATION (ATOM)
         */

        // PTR_ENV_START = 20
        int offset = PTR_ENV_START;  // an offset such that libcint can retrieve the correct index inside the environment


        for (size_t i = 0; i < natm; i++) {
            libcint_atm[CHARGE_OF + ATM_SLOTS * i] = static_cast<int>(atoms[i].atomic_number);  // insert the charge/atomic number
            libcint_atm[PTR_COORD + ATM_SLOTS * i] = offset;  // pointer to the coordinates of the atom inside the libcint environment
            libcint_env[offset + 0] = atoms[i].position.x();  // insert the position of the atoms
            libcint_env[offset + 1] = atoms[i].position.y();
            libcint_env[offset + 2] = atoms[i].position.z();
            offset += 3;
        }



        /*
         *  BAS CONFIGURATION (BASIS)
         */


        int atom_index = 0;  // index of the atom the shell is centered on
        auto previous_atom = shell_set[0].get_atom();


        for (size_t n = 0; n < shell_set.numberOfShells(); n++) {

            const auto& current_shell = shell_set[n];
            const auto& contraction_coefficients = current_shell.get_contraction_coefficients();
            const auto& gaussian_exponents = current_shell.get_gaussian_exponents();


            // If there's a new atom, increment the index
            auto current_atom = current_shell.get_atom();
            if (current_atom != previous_atom) {
                atom_index++;
                previous_atom = current_atom;
            }
            libcint_bas[ATOM_OF  + BAS_SLOTS * n] = atom_index;


            libcint_bas[ANG_OF   + BAS_SLOTS * n] = static_cast<int>(current_shell.get_l());  // angular momentum
            libcint_bas[NPRIM_OF + BAS_SLOTS * n] = static_cast<int>(current_shell.contractionSize());  // number of primitives
            libcint_bas[NCTR_OF  + BAS_SLOTS * n] = 1;  // number of contractions  // FIXME ????

            libcint_bas[PTR_EXP  + BAS_SLOTS * n] = offset;  // pointer to the exponents of the shell inside the libcint environment
            for (size_t e = 0; e < gaussian_exponents.size(); e++, offset++) {  // also increment offset
                libcint_env[offset] = gaussian_exponents[e];
            }


            libcint_bas[PTR_COEFF + BAS_SLOTS * n] = offset;  // pointer to the contraction coefficients inside the libcint environment
            // input normalized coeff.
            for (size_t c = 0; c < contraction_coefficients.size(); c++, offset++) {  // also increment offset
                libcint_env[offset] = contraction_coefficients[c] * CINTgto_norm(libcint_bas[ANG_OF+BAS_SLOTS*n], libcint_env[libcint_bas[PTR_EXP+BAS_SLOTS*n]+c]);
            }

        }


        // CALCULATE ONE-ELECTRON INTEGRALS
        std::array<OneElectronOperator<double>, N> operator_components;
        for (auto& op : operator_components) {
            op = OneElectronOperator<double>::Zero(nbf, nbf);
        }

        for (size_t sh1 = 0; sh1 < shell_set.numberOfShells(); sh1++) {
            for (size_t sh2 = 0; sh2 < shell_set.numberOfShells(); sh2++) {

                int shls[2];
                shls[0] = static_cast<int>(sh1);
                shls[1] = static_cast<int>(sh2);
                int nbf_sh1 = CINTcgto_cart(static_cast<int>(sh1), libcint_bas);  // number of basis functions in first shell  // TODO: from our code?
                int nbf_sh2 = CINTcgto_cart(static_cast<int>(sh2), libcint_bas);  // number of basis functions in second shell  // TODO: from our code?
                //            auto nbf_sh1 = basisset[sh1].numberOfBasisFunctions();
                //            auto nbf_sh2 = basisset[sh2].numberOfBasisFunctions();


                double buf[N * nbf_sh1 * nbf_sh2];  // buffer where the integrals are calculated to
                //            double buf[n_components * nbf_sh1 * nbf_sh2];  // buffer where the integrals are calculated to

                function(buf, shls, libcint_atm, natm, libcint_bas, nbf, libcint_env);  // TODO: is zero?

                auto bf1 = shell_set.basisFunctionIndex(sh1);  // (index of) first bf in sh1
                auto bf2 = shell_set.basisFunctionIndex(sh2);  // (index of) first bf in sh2


                for (auto f1 = 0; f1 != nbf_sh1; ++f1) {  // f1: index of basis function within shell 1
                    for (auto f2 = 0; f2 != nbf_sh2; ++f2) { // f2: index of basis function within shell 2

                        double computed_integral = buf[f1 + f2 * nbf_sh1];  // integrals are packed in row-major form
                        operator_components[0](bf1 + f1, bf2 + f2) = computed_integral;

                    }
                }  // data access loops
            }
        }

        return operator_components;
    }





    TwoElectronOperator<double> calculateTwoElectronIntegrals(const Libint2eFunction& function) const;
};


}  // namespace GQCP


#endif  /* LibcintInterfacer_hpp */
