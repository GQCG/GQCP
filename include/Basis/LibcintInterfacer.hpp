#ifndef LibcintInterfacer_hpp
#define LibcintInterfacer_hpp

#include "Operator/OneElectronOperator.hpp"

#include <functional>



extern "C" {

#include <cint.h>


/*
 *  FUNCTIONS THAT AREN'T INSIDE <cint.h>
 */
int cint1e_ovlp_cart(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env);  // ( \| \)
int cint1e_ipnuc_cart(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env);


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





using LibintFunction = std::function<int (double*, int*, int*, int, int*, int, double*)>;




/**
 *  A class that takes care of the interfacing with the libcint library
 */
class LibcintInterfacer {
public:
    OneElectronOperator<double> calculateOverlapIntegrals() const;


    OneElectronOperator<double> calculateOneElectronIntegrals(const LibintFunction& function) const;
};


}  // namespace GQCP


#endif  /* LibcintInterfacer_hpp */
