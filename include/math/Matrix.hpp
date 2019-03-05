#ifndef Matrix_hpp
#define Matrix_hpp


#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "typedefs.hpp"

#include <Eigen/Dense>



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


private:

    static constexpr bool is_vector = (Cols == 1);
    static constexpr bool is_matrix = (Cols >= 2) || (Cols == Dynamic);  // if this is not a vector


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
    template <typename Z = Self>  // enable_if must have Z inside
    static enable_if_t<Self::is_vector, Z> FromFile(const std::string& filename, size_t rows) {

        // Initialize a zero vector
        Self result = Self::Zero(rows);

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


    /**
     *  Read a matrix from a given file
     *
     *  @param filename     the name of the file to be read in
     *  @param rows         the number of expected rows
     *  @param cols         the number of expected columns
     */
    template <typename Z = Self>  // enable_if must have Z inside
    static enable_if_t<Self::is_matrix, Z> FromFile(const std::string& filename, size_t rows, size_t cols) {

        Self result = Self::Zero(rows, cols);

        std::ifstream file (filename);
        if (file.is_open()) {
            std::string line;

            while (std::getline(file, line)) {
                std::vector<std::string> splitted_line;  // create a container for the line to be split in

                // Split the line on any whitespace or tabs.
                boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

                if (splitted_line.size() != 3) {
                    throw std::runtime_error("Found a line that doesn't contain exactly 3 fields delimited by whitespace.");
                }

                auto i = std::stoi(splitted_line[0]);
                auto j = std::stoi(splitted_line[1]);
                auto value = std::stod(splitted_line[2]);

                result(i,j) = value;
            }

            file.close();
        } else {
            throw std::runtime_error("Cannot open the given file. Maybe you specified a wrong path?");
        }

        return result;
    }



    /*
     *  PUBLIC METHODS
     */

    /**
     *  Print the contents of a matrix to an output filestream
     *
     *  @param output_stream        the stream used for outputting
     */
    template <typename Z = void>  // enable_if must have Z inside
    enable_if_t<Self::is_matrix, Z> print(std::ostream& output_stream = std::cout) const {

        for (size_t i = 0; i < this->rows(); i++) {
            for (size_t j = 0; j < this->cols(); j++) {
                output_stream << i << ' ' << j << "  " << this->operator()(i,j) << std::endl;
            }
        }
    }


    /**
     *  @param i    row index (starting from 0)
     *  @param j    column index (starting from 0)
     *
     *  @return the i-j minor (i.e. delete the i-th row and j-th column)
     */
    template <typename Z = Self>  // enable_if must have Z inside
    enable_if_t<Self::is_matrix, Z> minor(size_t i, size_t j) const {

        // Delete the i-th row
        Self A_i = Self::Zero(this->rows() - 1, this->cols());

        for (size_t i2 = 0; i2 < this->rows(); i2++) {
            if (i2 < i) {
                A_i.row(i2) = this->row(i2);
            } else if (i2 == i) {
                continue;
            } else if (i2 > i) {
                A_i.row(i2-1) = this->row(i2);
            }
        }

        // Delete the j-th column
        Self A_ij = Self::Zero(this->rows() - 1, this->cols() - 1);
        for (size_t j2 = 0; j2 < this->cols(); j2++) {
            if (j2 < j) {
                A_ij.col(j2) = A_i.col(j2);
            } else if (j2 == j) {
                continue;
            } else if (j2 > j) {
                A_ij.col(j2-1) = A_i.col(j2);
            }
        }

        return A_ij;
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
