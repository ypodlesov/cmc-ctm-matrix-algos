#include "common.h"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/compute/algorithm.hpp>

using boost::numeric::ublas::vector;
using boost::numeric::ublas::inner_prod;
using boost::numeric::ublas::outer_prod;
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::triangular_matrix;
using boost::numeric::ublas::upper;
using boost::numeric::ublas::matrix_column;
using boost::numeric::ublas::matrix_vector_range;
using boost::numeric::ublas::matrix_vector_slice;
using boost::numeric::ublas::range;
using boost::numeric::ublas::slice;
using boost::numeric::ublas::norm_2;

bool QR_decomposition(matrix<double>& A, matrix<double>& Q, triangular_matrix<double, upper>& R);