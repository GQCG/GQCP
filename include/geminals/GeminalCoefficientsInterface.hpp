#ifndef GeminalCoefficientsInterface_hpp
#define GeminalCoefficientsInterface_hpp


#include "WaveFunction/WaveFunction.hpp"
#include "FockSpace/FockSpace.hpp"


namespace GQCP {


class GeminalCoefficientsInterface {
public:
    // DESTRUCTOR
    virtual ~GeminalCoefficientsInterface();


    // PUBLIC METHODS
    /**
     *  @param onv      the ONV that is being projected on
     *
     *  @return the overlap of the APIG-like wave function with the given on, i.e. the projection of the APIG wave function onto that ONV
     */
    virtual double overlap(const ONV& onv) const = 0;

    /**
     *  @param fock_space       the seniority-zero Fock space the wave function should live in
     *
     *  @return the wave function expansion corresponding to the geminal coefficients
     */
    WaveFunction toWaveFunction(const FockSpace& fock_space) const;

};

}  // namespace GQCP



#endif /* GeminalCoefficientsInterface_hpp */
