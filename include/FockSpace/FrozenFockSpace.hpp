#ifndef GQCP_FROZENFOCKSPACE_HPP
#define GQCP_FROZENFOCKSPACE_HPP

#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/FockPermutator.hpp"
#include "FockSpace/FockSpace.hpp"


namespace GQCP {


/**
 *  The frozen Fock space for a number of orbitals and number of electrons
 *
 *  Utilizes a full non frozen sub Fock space
 */
class FrozenFockSpace: public BaseFockSpace, public FockPermutator {
private:
    size_t N;  // number of electrons
    size_t X;  // number of frozen orbitals
    FockSpace fock_space; // non-frozen sub Fock space

public:
    // CONSTRUCTORS
    /**
     *  @param K        the number of orbitals
     *  @param N        the number of electrons
     *  @param X        the number of frozen orbitals
     */
    FrozenFockSpace(size_t K, size_t N, size_t X);

    /**
     *  @param fock_space       non-frozen sub Fock space
     *  @param X                the number of frozen orbitals
     */
    FrozenFockSpace(const FockSpace& fock_space, size_t X);


    // DESTRUCTORS
    ~FrozenFockSpace() override = default;


    // GETTERS
    size_t get_N() const { return this->N; }
    size_t get_X() const { return this->X; }
    FockSpaceType get_type() const override { return FockSpaceType::FrozenFockSpace; }


    // PUBLIC METHODS
    /**
     *  @param address      the address (i.e. the ordening number) of the ONV
     *
     *  @return the ONV with the corresponding address
     */
    ONV makeONV(size_t address) const override;

    /**
     *  Set the current ONV to the next ONV: performs ulongNextPermutation() and updates the corresponding occupation indices of the ONV occupation array
     *
     *  @param onv      the current ONV
     */
    void setNextONV(ONV& onv) const override;

    /**
     *  @param onv      the ONV
     *
     *  @return the address (i.e. the ordering number) of the given ONV
     */
    size_t getAddress(const ONV& onv) const override;

    /**
     *  Transform an ONV to one corresponding to the given address
     *
     *  @param onv          the ONV
     *  @param address      the address to which the ONV will be set
     */
    void transformONV(ONV& onv, size_t address) const override;

    /**
     *  @param onv       the ONV
     *
     *  @return the amount of ONVs (with a larger address) this ONV would couple with given a one electron operator
     */
    size_t countOneElectronCouplings(const ONV& onv) const;

    /**
     *  @param onv       the ONV
     *
     *  @return the amount of ONVs (with a larger address) this ONV would couple with given a two electron operator
     */
    size_t countTwoElectronCouplings(const ONV& onv) const;

    /**
     *  @return the amount non-zero (non-diagonal) couplings of a one electron coupling scheme in the Fock space
     */
    size_t countTotalOneElectronCouplings() const;

    /**
     *  @return the amount non-zero (non-diagonal) couplings of a two electron coupling scheme in the Fock space
     */
    size_t countTotalTwoElectronCouplings() const;
};


}  // namespace GQCP


#endif //GQCP_FROZENFOCKSPACE_HPP
