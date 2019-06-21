#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace np = boost::python::numpy;

BOOST_PYTHON_MODULE(pygqcp) {
    Py_Initialize();
    np::initialize();
}