#ifndef Matrix_hpp
#define Matrix_hpp


#include <fstream>

#include "typedefs.hpp"

#include <Eigen/Dense>

#include <boost/algorithm/string.hpp>



namespace GQCP {


constexpr auto Dynamic = Eigen::Dynamic;


/**
 *  An extension of the Eigen::Matrix class, with extra operations
 *
 *  We have decided to inherit from Eigen::Matrix, because we will use different hierarchies: see also: https://eigen.tuxfamily.org/dox-devel/TopicCustomizing_InheritingMatrix.html
 */
template <typename _Scalar = double, int _Rows = Dynamic, int _Cols = Dynamic>
class Matrix : public Eigen::Matrix<_Scalar, _Rows, _Cols> {
public:

    using Scalar = _Scalar;
    enum {
        Rows = _Rows,
        Cols = _Cols
    };

    using Self = Matrix<Scalar, Rows, Cols>;
    using Base = Eigen::Matrix<Scalar, Rows, Cols>;


public:

    /*
     *  CONSTRUCTORS
     */
    using Base::Base;  // inherit Base constructors


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Read a vector from a given file
     *
     *  @param filename     the name of the file to be read in
     *  @param rows         the number of expected rows
     */
    template <int Z = Cols>
    static enable_if_t<Z == 1, Self> FromFile(const std::string& filename, size_t rows) {

        // Initialize a zero vector
        Base result = Base::Zero(rows);

        std::ifstream file (filename);
        size_t index = 0;
        if (file.is_open()) {
            std::string line;

            while (std::getline(file, line)) {
                std::vector<std::string> splitted_line;  // create a container for the line to be split in

                // Split the line on any whitespace or tabs.
                boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

                if (splitted_line.size() != 1) {
                    throw std::runtime_error("Found a line that doesn't contain exactly 1 field delimited by whitespace.");
                }

                auto value = std::stod(splitted_line[0]);
                result(index) = value;

                ++index;
            }

            file.close();
        } else {
            throw std::runtime_error("Cannot open the given file. Maybe you specified a wrong path?");
        }

        return result;
    }
};



/*
 *  Convenience typedefs related to Matrix
 */

template <typename Scalar>
using MatrixX = Matrix<Scalar, Dynamic, Dynamic>;

template<typename Scalar, int Rows>
using Vector = Matrix<Scalar, Rows, 1>;

template <typename Scalar>
using VectorX = Vector<Scalar, Dynamic>;

using VectorXs = VectorX<size_t>;


}  // namespace GQCP



#endif  /* Matrix_hpp */
