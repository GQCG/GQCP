#ifndef GQCP_IO_HPP
#define GQCP_IO_HPP


#include <string>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>



namespace GQCP {


/**
  *  Read a vector from a given file and put the elements in the given vector
 *
 *  @param filename     the name of the file to be read in
 *  @param v            the vector that gets the contents of the file
 */
void readVectorFromFile(const std::string& filename, Eigen::VectorXd& v);

/**
 *  Read a matrix from a given file and put the elements in the given matrix
 *
 *  @param filename     the name of the file to be read in
 *  @param M            the matrix that gets the contents of the file
 */
void readArrayFromFile(const std::string& filename, Eigen::MatrixXd& M);

/**
 *  Read a tensor from a given file and put the elements in the given tensor
 *
 *  @param filename     the name of the file to be read in
 *  @param T            the tensor that gets the contents of the file
 */
void readArrayFromFile(const std::string& filename, Eigen::Tensor<double, 4>& T);

/**
 *  @param M        the matrix to be printed
 *
 *  Print the contents of a matrix in a fashionable way
 */
void print(const Eigen::MatrixXd& M);

/**
 *  Print the contents of a matrix to an output filestream
 *
 *  @param M                        the matrix to be printed
 *  @param output_filestream        the stream used for outputting
 */
void print(const Eigen::MatrixXd& M, std::ofstream& output_filestream);

/**
 *  Print the contents of a matrix to an output stream
 *
 *  @param T                        the tensor to be printed
 *  @param output_stream            the stream used for outputting
 */
void print(const Eigen::Tensor<double, 4>& T, std::ostream& output_stream);

/**
 *  Print the contents of a tensor in a fashionable way
 *
 *  @param T        the tensor to be printed
 */
void print(const Eigen::Tensor<double, 4>& T);


}  // namespace GQCP



#endif  // GQCP_IO_HPP
