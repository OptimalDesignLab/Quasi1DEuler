/**
 * \file nozzle.cpp
 * \brief function defintions for Nozzle member functions
 * \author  Jason Hicken <jason.hicken@gmail.com>
 * \version 1.0
 */

#include "./nozzle.hpp"

#include <math.h>

#include <ostream>
#include <iostream>
#include <algorithm>

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "./inner_prod_vector.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::min;
namespace ublas = boost::numeric::ublas;

const double pi = boost::math::constants::pi<double>();

// ======================================================================

Nozzle::Nozzle() {
  num_coeff_ = -1;
}

// ======================================================================

void Nozzle::SetCoeff(const InnerProdVector & coeff) {
  if (coeff.size() <= 0) {
    cerr << "Nozzle(SetCoeff): size of coeff is zero" << endl;
    throw(-1);    
  }
  coeff_ = coeff;
  num_coeff_ = coeff_.size();
}

// ======================================================================

void Nozzle::GetCoeff(InnerProdVector & coeff) const {
  for (int i = 0; i < num_coeff_ ; i++)
    coeff(i) = coeff_(i);
}

// ======================================================================

void Nozzle::SetAreaAtEnds(const double & area_left,
                           const double & area_right) {
  if ( (area_left < 0.0) || (area_right < 0.0) ) {
    cerr << "Nozzle(SetAreaAtEnds): invalid areas" << endl;
    cerr << "area_left = " << area_left 
         << ": area_right = " << area_right << endl;
    throw(-1);
  }
  area_left_ = area_left;
  area_right_ = area_right;
}

// ======================================================================

InnerProdVector FourierNozzle::Area(const InnerProdVector & x_coord) {
  double x_min = 1.E+16;
  double x_max = -x_min;
  for (int i = 0; i < x_coord.size(); i++) {
    x_min = min(x_min, x_coord(i));
    x_max = max(x_max, x_coord(i));
  }
  if ( (x_max - x_min) < 1.E-16 ) {
    cerr << "FourierNozzle(Area): x_coord is invalid... "
         << "x_max - x_min < eps" << endl;
    throw(-1);
  }

  InnerProdVector area_at_x(x_coord.size(), 0.0);
  // add the constant and linear part of the nozzle function
  double fac = (area_right_ - area_left_)/(x_max - x_min);
  for (int i = 0; i < x_coord.size(); i++)
    area_at_x(i) += area_left_ + fac*(x_coord(i) - x_min);

  // add the sinusoidal part of the nozzle function
  for (int j = 0; j < num_coeff_; j++) {
    double scale = pi*static_cast<double>(j+1)/(x_max - x_min);
    for (int i = 0; i < x_coord.size(); i++) {
      double x_tmp = scale*(x_coord(i) - x_min);
      area_at_x(i) += coeff_(j)*sin(x_tmp);
    }
  }
  return area_at_x;
}

// ======================================================================

double FourierNozzle::Area(const double & x_min, const double & x_max,
                           const double & x) {
  if ( (x_max - x_min) < 1.E-16 ) {
    cerr << "FourierNozzle(Area): x_max - x_min < eps"
         << endl;
    throw(-1);
  }
  double area = 0.0;
  // add the constant and linear part of the nozzle function
  double x_hat = (x - x_min)/(x_max - x_min);
  area += area_left_ + (area_right_ - area_left_)*x_hat;

  // add the sinusoidal part of the nozzle function
  for (int j = 0; j < num_coeff_; j++) {
    double x_tmp = pi*static_cast<double>(j+1);
    x_tmp *= x_hat;
    area += coeff_(j)*sin(x_tmp);
  }
  return area;
}

// ======================================================================

InnerProdVector FourierNozzle::AreaForwardDerivative(
    const InnerProdVector & x_coord, const InnerProdVector & u) {
  // check for consistent sizes
  if (u.size() != num_coeff_) {
    cerr << "FourierNozzle(AreaForwardDerivative): "
         << "u.size() is incompatible with num_coeff_"
         << endl;
    throw(-1);
  }
  double x_min = 1.E+16;
  double x_max = -x_min;
  for (int i = 0; i < x_coord.size(); i++) {
    x_min = min(x_min, x_coord(i));
    x_max = max(x_max, x_coord(i));
  }
  if ( (x_max - x_min) < 1.E-16 ) {
    cerr << "FourierNozzle(AreaForwardDerivative): "
         << "x_coord is invalid... x_max - x_min < eps"
         << endl;
    throw(-1);
  }
  InnerProdVector d_area(x_coord.size(), 0.0);
  for (int j = 0; j < num_coeff_; j++) {
    double scale = pi*static_cast<double>(j+1)/(x_max - x_min);
    for (int i = 0; i < x_coord.size(); i++) {
      double x_tmp = scale*(x_coord(i) - x_min);
      d_area(i) += u(j)*sin(x_tmp);
    }
  }
  return d_area;
}

// ======================================================================

InnerProdVector FourierNozzle::AreaReverseDerivative(
    const InnerProdVector & x_coord, const InnerProdVector & u) {
  // check for consistent sizes
  if (u.size() != x_coord.size()) {
    cerr << "FourierNozzle(AreaReverseDerivative): "
         << "u.size() is incompatible with x_coord.size()"
         << endl;
    throw(-1);
  }
  double x_min = 1.E+16;
  double x_max = -x_min;
  for (int i = 0; i < x_coord.size(); i++) {
    x_min = min(x_min, x_coord(i));
    x_max = max(x_max, x_coord(i));
  }
  if ( (x_max - x_min) < 1.E-16 ) {
    cerr << "FourierNozzle(AreaReverseDerivative): "
         << "x_coord is invalid... x_max - x_min < eps"
         << endl;
    throw(-1);
  }
  InnerProdVector d_area(num_coeff_, 0.0);
  for (int j = 0; j < num_coeff_; j++) {
    double scale = pi*static_cast<double>(j+1)/(x_max - x_min);
    for (int i = 0; i < x_coord.size(); i++) {
      double x_tmp = scale*(x_coord(i) - x_min);
      d_area(j) += u(i)*sin(x_tmp);
    }
  }
  return d_area;
}

// ======================================================================

void BsplineNozzle::SetCoeff(const InnerProdVector & coeff) { //const int & order) {
  int order = 4;
  Nozzle::SetCoeff(coeff);
  InnerProdVector cpts(num_coeff_ + 2);
  cpts(0) = area_left_;
  for (int i = 1; i < num_coeff_ + 1; i++)
    cpts(i) = coeff(i-1);
  cpts(num_coeff_ + 1) = area_right_;
  spline_.SetCptsAndOrder(cpts, order);
}

// ======================================================================

void BsplineNozzle::FitNozzle(const InnerProdVector & x_coord,
                              const InnerProdVector & y_coord) {
  spline_.FitCurve(x_coord, y_coord);
  for (int i = 1; i < num_coeff_ + 1; i++)
    coeff_(i-1) = spline_.get_cpts(i);
}

// ======================================================================

InnerProdVector BsplineNozzle::Area(const InnerProdVector & x_coord) {
  double x_min = 1.E+16;
  double x_max = -x_min;
  for (int i = 0; i < x_coord.size(); i++) {
    x_min = min(x_min, x_coord(i));
    x_max = max(x_max, x_coord(i));
  }
  if ( (x_max - x_min) < 1.E-16 ) {
    cerr << "BsplineNozzle(Area): x_coord is invalid... "
         << "x_max - x_min < eps" << endl;
    throw(-1);
  }

  InnerProdVector area_at_x(x_coord.size(), 0.0);
  double fac = 1.0/(x_max - x_min);
  for (int i = 0; i < x_coord.size(); i++)
    area_at_x(i) = spline_.EvalAtPoint(fac*(x_coord(i) - x_min));
  return area_at_x;
}

// ======================================================================

double BsplineNozzle::Area(const double & x_min, const double & x_max,
                           const double & x) {
  if ( (x_max - x_min) < 1.E-16 ) {
    cerr << "BsplineNozzle(Area): x_max - x_min < eps"
         << endl;
    throw(-1);
  }
  double area = 0.0;
  double fac = 1.0/(x_max - x_min);
  area = spline_.EvalAtPoint(fac*(x - x_min));
  return area;
}

// ======================================================================

InnerProdVector BsplineNozzle::AreaForwardDerivative(
    const InnerProdVector & x_coord, const InnerProdVector & u) {
  // check for consistent sizes
  if (u.size() != num_coeff_) {
    cerr << "BsplineNozzle(AreaForwardDerivative): "
         << "u.size() is incompatible with num_coeff_"
         << endl;
    throw(-1);
  }
  double x_min = 1.E+16;
  double x_max = -x_min;
  for (int i = 0; i < x_coord.size(); i++) {
    x_min = min(x_min, x_coord(i));
    x_max = max(x_max, x_coord(i));
  }
  if ( (x_max - x_min) < 1.E-16 ) {
    cerr << "BsplineNozzle(AreaForwardDerivative): "
         << "x_coord is invalid... x_max - x_min < eps"
         << endl;
    throw(-1);
  }
  InnerProdVector d_area(x_coord.size(), 0.0);
  InnerProdVector cpts(num_coeff_ + 2);
  cpts(0) = 0.0;
  for (int i = 1; i < num_coeff_ + 1; i++)
    cpts(i) = u(i-1);
  cpts(num_coeff_ + 1) = 0.0;
  Bspline tmp_spline;
  int order = 4;
  tmp_spline.SetCptsAndOrder(cpts, order);

  double fac = 1.0/(x_max - x_min);
  for (int i = 0; i < x_coord.size(); i++)
    d_area(i) = tmp_spline.EvalAtPoint(fac*(x_coord(i) - x_min));
  return d_area;
}

// ======================================================================

InnerProdVector BsplineNozzle::AreaReverseDerivative(
    const InnerProdVector & x_coord, const InnerProdVector & u) {
  // check for consistent sizes
  if (u.size() != x_coord.size()) {
    cerr << "BsplineNozzle(AreaReverseDerivative): "
         << "u.size() is incompatible with x_coord.size()"
         << endl;
    throw(-1);
  }
  double x_min = 1.E+16;
  double x_max = -x_min;
  for (int i = 0; i < x_coord.size(); i++) {
    x_min = min(x_min, x_coord(i));
    x_max = max(x_max, x_coord(i));
  }
  if ( (x_max - x_min) < 1.E-16 ) {
    cerr << "BsplineNozzle(AreaReverseDerivative): "
         << "x_coord is invalid... x_max - x_min < eps"
         << endl;
    throw(-1); 
  }

  InnerProdVector d_area(num_coeff_, 0.0);
  InnerProdVector d_tmp(num_coeff_ + 2, 0.0);
  int order = spline_.get_order();
  InnerProdVector basis(order, 0.0);
  double fac = 1.0/(x_max - x_min);
  for (int i = 0; i < x_coord.size(); i++) {
    double u_tmp = u(i);
    int index;
    double x = fac*(x_coord(i) - x_min);
    spline_.BasisValues(x, index, basis);
    for (int j = 1; j <= order; j++)
      d_tmp(j + index - order - 1) += basis(j-1)*u_tmp;
  }
  for (int j = 0; j < num_coeff_; j++)
    d_area(j) = d_tmp(j+1);
  return d_area;

#if 0
  InnerProdVector d_area_2(num_coeff_, 0.0);
  InnerProdVector cpts(num_coeff_ + 2, 0.0);
  Bspline tmp_spline;
  order = spline_.get_order();
  fac = 1.0/(x_max - x_min);
  for (int j = 0; j < num_coeff_; j++) {
    cpts(j+1) = 1.0;
    tmp_spline.SetCptsAndOrder(cpts, order);    
    for (int i = 0; i < x_coord.size(); i++)
      d_area_2(j) += u(i)*tmp_spline.EvalAtPoint(fac*(x_coord(i) - x_min));
    cpts(j+1) = 0.0;
  }  
  for (int j = 0; j < num_coeff_; j++) {
    cout << "d_area = " << d_area(j) << ": d_area_2 = " << d_area_2(j)
         << endl;
  }
  throw(-1);
#endif
}

// ======================================================================
