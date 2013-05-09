/**
 * \file verify_accuracy.cpp
 * \brief unit test for SumByParts class
 * \version 1.0
 */

#include "../inner_prod_vector.hpp"
#include "../sum_by_parts.hpp"

#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;

using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char *argv[]) {

  const int order = 4;
  const int nodes_c = 41;
  const int nodes_f = 81;
  
  InnerProdVector func_c(nodes_c, 0.0), func_f(nodes_c, 0.0),
      dfunc_c(nodes_c, 0.0), dfunc_f(nodes_c, 0.0),
      dfdx_c(nodes_c, 0.0), dfdx_f(nodes_c, 0.0),
      ord(nodes_c, 0.0);

  // loop over the coarse grid nodes where we want to find the order
  for (int j = 0; j < nodes_c; j++) {      
    // for a given x_j, we find a function centered at x_j and
    // decrease the mesh spacing going to finer grid (number of nodes
    // remains fixed).  This allows us to find the order of accuracy
    // at node j.

    // define the function and derivative on the course grid
    double dx = 1.0/static_cast<double>(nodes_c-1);
    double x0 = dx*static_cast<double>(j);
    for (int i = 0; i < nodes_c; i++) {
      double x = dx*static_cast<double>(i-j);
      func_c(i) = exp(x);
      dfunc_c(i) = func_c(i);
    }
  
    // define the function on the fine grid
    dx *= 0.5;
    for (int i = 0; i < nodes_c; i++) {
      double x = dx*static_cast<double>(i-j);
      func_f(i) = exp(x);
      dfunc_f(i) = func_f(i);
    }
  
    // calculate the finite difference approximation
    SBP1stDerivative Dx(nodes_c, order);
    Dx.Apply(1, func_c, dfdx_c);
    Dx.Apply(1, func_f, dfdx_f);
    dfdx_c *= static_cast<double>(nodes_c-1);
    dfdx_f *= static_cast<double>(nodes_f-1);

    // calculate the error in the finite-difference approx
    double err_c(fabs(dfunc_c(j) - dfdx_c(j)));
    double err_f(fabs(dfunc_f(j) - dfdx_f(j)));
    ord(j) = (log(err_c/err_f)/log(2.0));
    cout << "node " << j << " err_c = " << err_c 
         << ": err_f " << err_f << endl;
  }

  for (int i = 0; i < nodes_c; i++) {
    cout << "node " << i << " order = " << ord(i) << endl;
  }
}
