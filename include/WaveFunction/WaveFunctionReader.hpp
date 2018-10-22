#ifndef GQCP_WAVEFUNCTIONREADER_HPP
#define GQCP_WAVEFUNCTIONREADER_HPP



#include <boost/dynamic_bitset.hpp>

#include "FockSpace/SelectedFockSpace.hpp"
#include "WaveFunction/WaveFunction.hpp"

#include <memory>


namespace GQCP {

/**
 *  Class that reads and stores a selected wavefunction expansion
 */
class WaveFunctionReader {
private:
    GQCP::SelectedFockSpace fock_space;
    Eigen::VectorXd coefficients;
    GQCP::WaveFunction wave_function;


public:
    /**
     *  Constructor based on a given @param GAMESS_filename
     */
    explicit WaveFunctionReader(const std::string& GAMESS_filename);


    // GETTERS
    SelectedFockSpace get_fock_space() const { return this->fock_space; }
    Eigen::VectorXd get_coefficients() const { return this->coefficients; }
    WaveFunction get_wave_function() const { return this->wave_function; }
};



}  // namespace GQCP



#endif //GQCP_WAVEFUNCTIONREADER_HPP
