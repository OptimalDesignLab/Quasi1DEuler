/**
 * \file bspline.hpp
 * \brief header file for Bspline class
 * \author  Jason Hicken <jason.hicken@gmail.com>
 * \version 1.0
 */

#pragma once

#include <math.h>

#include <ostream>
#include <iostream>

#include "./inner_prod_vector.hpp"

// ======================================================================

/*!
 * \class Bspline
 * \brief defines one dimensional Bsplines
 */
class Bspline {
 public:

  /*
   * \brief default constructor
   */
  Bspline();

  /*
   * \brief constructor for uniformly spaced knots and cpts_ = const
   * \param[in] num_cp - number of control points
   * \param[in] order - b-spline order (degree = order+1)
   */
  Bspline(const int & num_cp, const int & order);

  /*
   * \brief constructor based on given knots and control points
   * \param[in] knots - knot locations
   * \param[in] cpts - control point values
   */
  Bspline(const InnerProdVector & knots, const InnerProdVector & cpts);

  /*
   * \brief destructor
   */
  ~Bspline() {}

  /*
   * \brief set knot locations 
   * \param[in] knots - knot locations
   */
  void set_knots(const InnerProdVector & knots) { knot_ = knots; }

  /* \brief return the order of the b-spline
   * \returns the order_ member
   */
  int get_order() const { return order_; }

  /* \brief return an element of the control point vector
   * \param[in] i - desired index
   */
  const double & get_cpts(const int & i) const { return cpts_(i); }

  /*
   * \brief set the value of the control points and the order
   * \param[in] cpts - control point values
   * \param[in] order - the Bspline order
   */
  void SetCptsAndOrder(const InnerProdVector & cpts, const int & order);

  /*!
   * \brief set the b-spline coefficients by LS fitting to given data
   * \param[in] x - set of coordinates where the data are given
   * \param[in] y - set of data to fit
   */
  void FitCurve(const InnerProdVector & x, const InnerProdVector & y);

  /*
   * \brief find the index such that knot_(index-1) <= x <= knot_(index)
   * \param[in] x - coordinate whose knot interval we want
   */
  int FindKnotInterval(const double & x) const;

  /*
   * \brief returns the nonzero basis values at location x
   * \param[in] x - location where we want the basis
   * \param[out] interval_index - value returned by FindKnotInterval(x)
   * \param[in,out] basis - nonzero basis values
   */
  void BasisValues(const double & x, int & interval_index,
                   InnerProdVector & basis);

  /*
   * \brief evaluate the B-spline at a given x-coordinate
   * \param[in] x - location at which to evaluate B-spline
   * \param[in] deriv_order - derivative order >= 0
   */
  double EvalAtPoint(const double & x, const int & deriv_order = 0);

  /*
   * \brief evaluate the B-spline at a given x-coordinates
   * \param[in] x - locations at which to evaluate B-spline
   */
  void EvalAtPoints(const InnerProdVector & x, InnerProdVector & y);

  private:
  
  int num_cpts_; ///< number of control points
  int order_; ///< order of the B-spline
  InnerProdVector knot_; ///< knot locations
  InnerProdVector cpts_; ///< control point values
  InnerProdVector aj_; ///< work array
  InnerProdVector dl_; ///< work array
  InnerProdVector dr_; ///< work array
};
