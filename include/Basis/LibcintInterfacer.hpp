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


/**
 *  Forward declaration
 */
class LibcintInterfacer;



/*
 *  WRAPPERS AROUND C-STYLE LIBCINT ARRAYS
 */
namespace libcint {


static constexpr int atm_slots = ATM_SLOTS;
static constexpr int bas_slots = BAS_SLOTS;


    // https://github.com/sunqm/libcint/blob/master/examples/python_call.py

/**
 *  A wrapper around the libcint "atm" C array
 */
struct Atom {
public:
    int charge_of;  // nuclear charge/atomic number
    int ptr_coord;  // offset of the x-coordinate of the atom inside the libcint environment
    int nuc_mod_of;
    int ptr_zeta;
    int reserve_atmlot1;
    int reserve_atmlot2;
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
    int kappa_of;
    int ptr_exp;  // offset of the exponents inside the libcint environment
    int ptr_coeff;  // offset to the contraction coefficients inside the libcint environment
    int reserve_baslot;
};


///**
// *  A struct wrapper around the libcint "env" C array
// */
//struct Environment {
//public:
//    double buffer [10000];
//};


///**
// *  A wrapper around libcint::Atoms, libcint::BasisFunctions and libcint::Environment
// */
//struct Container {
//public:
//    std::vector<Atom> atoms;
//    std::vector<BasisFunction> basis_functions;
//};


/**
 *  A wrapper that owns raw libcint 'atm', 'bas' and 'env' arrays
 */
class RawContainer {
private:
    // PRIVATE MEMBERS
    int natm;  // number of atoms
    int nbf;  // number of basis functions
    int nsh;  // the number of shells

    int* libcint_atm;  // information about the atoms
    int* libcint_bas;  // information about the basis functions
    double* libcint_env;  // a raw block of doubles in which libcint (probably) places intermediary calculations


public:
    // CONSTRUCTORS
    /**
     *  Allocate memory for the raw libcint arrays
     *
     *  @param natm         the number of atoms
     *  @param nbf          the number of basis functions
     *  @param nsh          the number of shells
     */
    RawContainer(size_t natm, size_t nbf, size_t nsh) :
        natm (static_cast<int>(natm)),
        nbf (static_cast<int>(nbf)),
        nsh (static_cast<int>(nsh)),
        libcint_atm (new int[this->natm * atm_slots]),
        libcint_bas (new int[this->nbf * bas_slots]),
        libcint_env (new double[10000])
    {}


    // DESTRUCTOR
    /**
     *  Deallocate the memory used by the raw libcint arrays
     */
    ~RawContainer() {
        delete[] this->libcint_atm;
        delete[] this->libcint_bas;
        delete[] this->libcint_env;
    }

    // FRIENDS
    friend class GQCP::LibcintInterfacer;
};





}  // namespace libcint





using Libint1eFunction = std::function<int (double*, int*, int*, int, int*, int, double*)>;
using Libint2eFunction = std::function<int (double*, int*, int*, int, int*, int, double*, CINTOpt*)>;




/**
 *  A class that takes care of the interfacing with the libcint library
 */
class LibcintInterfacer {
public:



    // PUBLIC METHODS

//    /**
//     *  @param libcint_container        the wrapper libcint container
//     *
//     *  @return the raw information about atoms that libcint can use
//     */
//    libcint::RawContainer interface(const libcint::Container& libcint_container) const;





    /**
     *  @param shell_set        the GQCP::ShellSet whose information should be converted
     *
     *  @return the information in a GQCP::ShellSet as a libcint::RawContainer
     */
    libcint::RawContainer convert(const ShellSet& shell_set) const;




    /**
     *  @tparam N                   the number of libcint operator components
     *
     *  @param function             the libcint integral function
     *  @param raw_container        the data libcint needs to perform calculations
     *
     *  @return an array of N OneElectronOperators corresponding to the matrix representations of the N components of the given operator represented by the libcint function
     */
    template <size_t N>
    std::array<OneElectronOperator<double>, N> calculateOneElectronIntegrals(const Libint1eFunction& function, libcint::RawContainer& raw_container) const {

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
                            double computed_integral = buf[f1 + nbf_sh1 * (f2 + nbf_sh2 * i)];  // integrals are packed in column-major form
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





    TwoElectronOperator<double> calculateTwoElectronIntegrals(const Libint2eFunction& function) const;
};


}  // namespace GQCP


#endif  /* LibcintInterfacer_hpp */
