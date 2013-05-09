/**
 * \file bspline.cpp
 * \brief function defintions for Bspline member functions
 * \author  Jason Hicken <jason.hicken@gmail.com>
 * \version 1.0
 */

#include "./bspline.hpp"

#include <math.h>

#include <ostream>
#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "./inner_prod_vector.hpp"

using std::cout;
using std::cerr;
using std::endl;

namespace ublas = boost::numeric::ublas;

// ======================================================================

Bspline::Bspline() {
  num_cpts_ = -1;
  order_ = -1;
}

// ======================================================================

Bspline::Bspline(const int & num_cpts, const int & order) {
  if ( (num_cpts <= 0) || (order <= 0) ) {
    cerr << "Error in Bspline(constructor): "
         << "number of control points or order invalid:" << endl;
    cerr << "num_cpts = " << num_cpts << ": order = " << order << endl;
    throw(-1);
  }
  num_cpts_ = num_cpts;
  order_ = order;
  cpts_.resize(num_cpts_);
  knot_.resize(num_cpts_ + order_);
  aj_.resize(order_);
  dl_.resize(order_-1);
  dr_.resize(order_-1);
  // set uniform knot values
  for (int j = 0; j < order_; j++)
    knot_(j) = 0.0;
  for (int j = order_; j < num_cpts_; j++)
    knot_(j) = static_cast<double>(j - order_ + 1)
        /static_cast<double>(num_cpts_ - order_ + 1);
  for (int j = 0; j < order_; j++)
    knot_(num_cpts_ + j) = 1.0;
}

// ======================================================================

Bspline::Bspline(const InnerProdVector & knots,
                 const InnerProdVector & cpts) {
  int size_knots = knots.size();
  int size_cpts = cpts.size();
  if (size_cpts <= 0) {
    cerr << "Error in Bspline(constructor): "
         << "number of control points invalid:" << endl;
    cerr << "cpts.size() = " << size_cpts << endl;
    throw(-1);
  }
  if (size_knots <= size_cpts) {
    cerr << "Error in Bspline(constructor): "
         << "number of knots invalid:" << endl;
    cerr << "knots.size() = " << size_knots << endl;
    throw(-1);
  }
  num_cpts_ = size_cpts;
  order_ = size_knots - num_cpts_;
  cpts_ = cpts;
  knot_ = knots;
  aj_.resize(order_);
  dl_.resize(order_-1);
  dr_.resize(order_-1);
  // need to check that knots are monotone increasing
}

// ======================================================================

void Bspline::SetCptsAndOrder(const InnerProdVector & cpts,
                              const int & order) {
  num_cpts_ = cpts.size();
  if ( (num_cpts_ <= 0) || (order <= 0) ) {
    cerr << "Error in Bspline(SetCptsAndOrder): "
         << "number of control points or order invalid:" << endl;
    cerr << "num_cpts_ = " << num_cpts_ << ": order = " << order << endl;
    throw(-1);
  }
  order_ = order;
  cpts_ = cpts;
  knot_.resize(num_cpts_ + order_);
  aj_.resize(order_);
  dl_.resize(order_-1);
  dr_.resize(order_-1);
  // set uniform knot values
  for (int j = 0; j < order_; j++)
    knot_(j) = 0.0;
  for (int j = order_; j < num_cpts_; j++)
    knot_(j) = static_cast<double>(j - order_ + 1)
        /static_cast<double>(num_cpts_ - order_ + 1);
  for (int j = 0; j < order_; j++)
    knot_(num_cpts_ + j) = 1.0;
}

// ======================================================================

void Bspline::FitCurve(const InnerProdVector & x,
                       const InnerProdVector & y) {
  if ( (x.size() != y.size()) || (x.size() < num_cpts_) ) {
    cerr << "Error in Bspline(FitCurve): "
         << "problem with number of data points." << endl;
    cerr << "x.size() = " << x.size() << endl;
    cerr << "y.size() = " << y.size() << endl;
    cerr << "num_cpts_ = " << num_cpts_ << endl;    
    throw(-1);
  }
  int num_pts = x.size();
  ublas::matrix<double> Ncoeff(num_pts, num_cpts_), A(num_cpts_, num_cpts_);
  ublas::vector<double> rhs(num_pts), b(num_cpts_);

  for (int i = 0; i < num_pts; i++) {
    rhs(i) = 0.0;
    for (int j = 0; j < num_cpts_; j++)
      Ncoeff(i, j) = 0.0;
  }

  for (int i = 0; i < num_pts; i++) {
    // evaluate the basis at x(i);
    int knot_interval = FindKnotInterval(x(i));
    InnerProdVector basis(order_, 0.0);
    BasisValues(x(i), knot_interval, basis);
    for (int j = 0; j < order_; j++)
      Ncoeff(i, j + knot_interval - order_) = basis(j);
    rhs(i) = y(i);
  }
  // form the least-squares system
  A = prod(trans(Ncoeff), Ncoeff);
  b = prod(trans(Ncoeff), rhs);
  int err = ublas::lu_factorize(A);
  if( err != 0 ) {
    cout << "Error in Bspline(FitCurve): "
         << "lu_factorize() returned with err = " << err << endl;
  }
  ublas::lu_substitute(static_cast<ublas::matrix<double> >(A), b);
  cpts_ = b;
}

// ======================================================================

int Bspline::FindKnotInterval(const double & x) const {
  // find the index such that knot_(index-1) <= x <= knot_(index)
  for (int j = order_; j <= num_cpts_; j++) {
    if ( (knot_(j-1) <= x) && (x <= knot_(j)) ) {
      return j;
    }
  }
}

// ======================================================================
void Bspline::BasisValues(const double & x, int & interval_index,
                          InnerProdVector & basis) {
  if ( (x < knot_(0)) || (x > knot_(knot_.size()-1)) ) {
    cerr << "Error in Bspline(BasisValues): "
         << "x evaluation point is outside knot range: x = " 
         << x << endl;
    throw(-1);
  }
  if (basis.size() != order_) {
    cerr << "Error in Bspline(BasisValues): "
         << "basis.size() is inconsistent with order_: "
         << "basis.size() = " << basis.size() 
         << ": order_ = " << order_ << endl;
    throw(-1);
  }
  interval_index = FindKnotInterval(x);

  basis(0) = 1.0;
  if (order_ > 1) {
    for (int k = 1; k <= order_ - 1; k++) {
      int km1 = k - 1;
      int kp1 = k + 1;
      dr_(km1) = knot_(interval_index + km1) - x;
      dl_(km1) = x - knot_(interval_index - k);
      double saved = 0.0;
      for (int i = 1; i <= k; i++) {
        double temp = basis(i-1)/(dr_(i-1) + dl_(k - i));
        basis(i-1) = saved + dr_(i-1)*temp;
        saved = dl_(k - i)*temp;
      }
      basis(k) = saved;
    }
  }

  // check on sum of basis = 1.0
  double temp = ublas::sum(basis);
  if (fabs(temp - 1.0) > 1.e-15) {
    cerr << "Error in Bspline(BasisValues): "
         << "basis values do not sum to 1.0: sum = "
         << temp << endl;
    throw(-1);
  }
}

// ======================================================================

double Bspline::EvalAtPoint(const double & x, const int & deriv_order) {
  // check that x is valid
  if ( (x < knot_(0)) || (x > knot_(knot_.size()-1)) ) {
    cerr << "Error in Bspline(EvalAtPoint): "
         << "x evaluation point is outside knot range: x = " 
         << x << endl;
    cerr << "knot_ = ";
    for (int j = 0; j < num_cpts_ + order_; j++)
      cerr << knot_(j) << " ";
    cerr << endl;
    throw(-1);
  }
  // check that deriv_order is valid
  if (deriv_order < 0) {
    cerr << "Error in Bspline(EvalAtPoint): "
         << "deriv_order deriviative order value is invalid: " << endl;
    cerr << "deriv_order = " << deriv_order << endl;
    throw(-1);
  }

  // find the left index such that knot_(left-1) <= x <= knot_(left)
  int left = FindKnotInterval(x);

  // we will store the (order) b-spline coefficients relevant
  // to the knot interval [knot_(left-1),knot_(left)]
  // in aj(0),...,aj(order-1)
  // and compute dl(j) = x - knot_(left-j)
  //             dr(j) = knot_(left+j-1) - x
  // for all j = 1,...,order_-1
  int km1 = order_ - 1;
  int jcmin = 1;
  int imk = left - order_;
  if (imk < 0) {
    // we are close to the left end of the knot interval, so some
    // of the aj will be set to zero later
    jcmin = 1 - imk;
    for (int j = 1; j <= left; j++)
      dl_(j-1) = x - knot_(left - j);
    for (int j = left; j <= km1; j++)
      dl_(j-1) = dl_(left-1);
  } else {
    for (int j = 1; j <= km1; j++)
      dl_(j-1) = x - knot_(left - j);
  }

  int jcmax = order_;
  int n = num_cpts_;
  int nmi = n - left;
  if (nmi < 0) {
    // we are close to the right end of the knot interval, so some
    // of the aj will be set to zero later
    jcmax = order_ + nmi;
    for (int j = 1; j <= jcmax; j++)
      dr_(j-1) = knot_(left + j - 1) - x;
    for (int j = jcmax; j <= km1; j++)
      dr_(j-1) = dr_(jcmax-1);
  } else {
    for (int j = 1; j <= km1; j++)
      dr_(j-1) = knot_(left + j - 1) - x;
  }

  // set all elements of aj(:,:,1) to zero, in case we are close to
  // the ends of the knot vector
  aj_ = 0.0;
  for (int jc = jcmin; jc <= jcmax; jc++)
    aj_(jc-1) = cpts_(imk + jc - 1);
  
  if (deriv_order != 0) {
    // derivative: apply the recursive formula X.12b from de Boor
    for (int j = 1; j <= deriv_order; j++) {
      int kmj = order_ - j;
      double dp_kmj = static_cast<double>(kmj);
      int ilo = kmj;
      for (int jj = 1; jj <= kmj; jj++) {
        aj_(jj-1) = (aj_(jj) - aj_(jj-1)) * dp_kmj 
            / (dl_(ilo-1) + dr_(jj-1));
        --ilo;
      }
    }
  }

  if (deriv_order != km1) {
    // if deriv_order /= order - 1, we need to apply the recursive
    // formula from de Boor
    for (int j = deriv_order+1; j <= km1; j++) {
      int kmj = order_ - j;
      int ilo = kmj;
      for (int jj = 1; jj <= kmj; jj++) {
        aj_(jj-1) = (aj_(jj)*dl_(ilo-1) + aj_(jj-1)*dr_(jj-1))
            / (dl_(ilo-1) + dr_(jj-1));
        --ilo;
      }
    }
  }
  return aj_(0);
}

// ======================================================================

void Bspline::EvalAtPoints(const InnerProdVector & x,
                           InnerProdVector & y) {
  if (x.size() != y.size()) {
    cerr << "Error in Bspline(EvalAtPoints): "
         << "output array y and input array x have different sizes." 
         << endl;
    throw(-1);
  }
  for (int i = 0; i < x.size(); i++)
    y(i) = EvalAtPoint(x(i));
}

// ======================================================================
