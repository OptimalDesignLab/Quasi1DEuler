/**
 * \file nozzle.hpp
 * \brief header file for Nozzle class
 * \author  Jason Hicken <jason.hicken@gmail.com>
 * \version 1.0
 */

#pragma once

#include <math.h>

#include <ostream>
#include <iostream>

#include "./inner_prod_vector.hpp"
#include "./bspline.hpp"

// ======================================================================

/*!
 * \class Nozzle
 * \brief quasi-1d nozzle geometry definition and manipulation
 */
class Nozzle {
 public:

  /*!
   * \brief default cnstructor
   */
  Nozzle();

  /*!
   * \brief default destructor
   */
  virtual ~Nozzle() = 0;

  /*!
   * \brief set the area at the left and right ends of the domain
   * \param[in] area_left - area at the left of domain
   * \param[in] area_right - area at the right of domain
   */
  void SetAreaAtEnds(const double & area_left,
                     const double & area_right);

  /*!
   * \brief set the values of the coefficents that define the nozzle
   * \param[in] coeff - values of the design coefficents
   */
  virtual void SetCoeff(const InnerProdVector & coeff);

  /*!
   * \brief returns the (internal) b-spline coefficients
   * \param[out] coeff - coefficient values
   */
  void GetCoeff(InnerProdVector & coeff) const;

  /*!
   * \brief returns the nozzle area at a given set of x coordinates
   * \param[in] x_coord - x coordinates to evaluate nozzle area
   * \returns nozzle area at x_coord
   */
  virtual InnerProdVector Area(const InnerProdVector & x_coord) = 0;

  /*!
   * \brief returns the nozzle area at a given x coordinate
   * \param[in] x_min - minimum x coordinate
   * \param[in] x_max - maximum x coordinate
   * \param[in] x - x coordinate to evaluate nozzle area
   * \returns area at x
   */
  virtual double Area(const double & x_min, const double & x_max,
                      const double & x) = 0;

  /*!
   * \brief calculate the forward mode derivative of area w.r.t. coeff
   * \param[in] x_coord - x coordinates to evaluate area derivative
   * \param[in] u - the direction that the derivative is desired
   * \returns the resulting derivative components
   */
  virtual InnerProdVector AreaForwardDerivative(
      const InnerProdVector & x_coord, const InnerProdVector & u) = 0;

  /*!
   * \brief calculate the reverse mode derivative of area w.r.t. coeff
   * \param[in] x_coord - x coordinates to evaluate area derivative
   * \param[in] u - the direction that the derivative is desired
   * \returns the resulting derivative components
   */
  virtual InnerProdVector AreaReverseDerivative(
      const InnerProdVector & x_coord, const InnerProdVector & u) = 0;
  
 protected:

  int num_coeff_; ///< number of design coefficients
  double area_left_; ///< area at the left end of the domain
  double area_right_; ///< area at the right end of the domai
  InnerProdVector coeff_; ///< design coefficients
};
inline Nozzle::~Nozzle() {}

// ======================================================================

/*!
 * \class FourierNozzle
 * \brief quasi-1d nozzle geometry defined using a Fourier sine series
 */
class FourierNozzle : public Nozzle {
 public:

  /*!
   * \brief default cnstructor
   */
  FourierNozzle() : Nozzle() {}

  /*!
   * \brief default destructor
   */
  ~FourierNozzle() {}

  /*!
   * \brief returns the nozzle area at a given set of x coordinates
   * \param[in] x_coord - x coordinates to evaluate nozzle area
   * \returns nozzle area at x_coord
   */
  InnerProdVector Area(const InnerProdVector & x_coord);

  /*!
   * \brief returns the nozzle area at a given x coordinate
   * \param[in] x_min - minimum x coordinate
   * \param[in] x_max - maximum x coordinate
   * \param[in] x - x coordinate to evaluate nozzle area
   * \returns area at x
   */
  double Area(const double & x_min, const double & x_max,
                      const double & x);

  /*!
   * \brief calculate the forward mode derivative of area w.r.t. coeff
   * \param[in] x_coord - x coordinates to evaluate area derivative
   * \param[in] u - the direction that the derivative is desired
   * \returns the resulting derivative components
   */
  InnerProdVector AreaForwardDerivative(
      const InnerProdVector & x_coord, const InnerProdVector & u);

  /*!
   * \brief calculate the reverse mode derivative of area w.r.t. coeff
   * \param[in] x_coord - x coordinates to evaluate area derivative
   * \param[in] u - the direction that the derivative is desired
   * \returns the resulting derivative components
   */
  InnerProdVector AreaReverseDerivative(
      const InnerProdVector & x_coord, const InnerProdVector & u);
};

// ======================================================================
/*!
 * \class BsplineNozzle
 * \brief quasi-1d nozzle geometry defined using a Bspline
 */
class BsplineNozzle : public Nozzle {
 public:

  /*!
   * \brief default cnstructor
   */
  BsplineNozzle() : Nozzle() {}

  /*!
   * \brief default destructor
   */
  ~BsplineNozzle() {}

  /*!
   * \brief set the b-spline coefficients
   * \param[in] coeff - coefficient values
   * \param[in] order - order of the bspline
   */
  void SetCoeff(const InnerProdVector & coeff); //, const int & order = 4);

  /*!
   * \brief set the b-spline coefficients by LS fitting to given data
   * \param[in] x_coord - set of coordinates where the data are given
   * \param[in] y_coord - set of data to fit
   */
  void FitNozzle(const InnerProdVector & x_coord,
                 const InnerProdVector & y_coord);

  /*!
   * \brief returns the nozzle area at a given set of x coordinates
   * \param[in] x_coord - x coordinates to evaluate nozzle area
   * \returns nozzle area at x_coord
   */
  InnerProdVector Area(const InnerProdVector & x_coord);

  /*!
   * \brief returns the nozzle area at a given x coordinate
   * \param[in] x_min - minimum x coordinate
   * \param[in] x_max - maximum x coordinate
   * \param[in] x - x coordinate to evaluate nozzle area
   * \returns area at x
   */
  double Area(const double & x_min, const double & x_max,
                      const double & x);

  /*!
   * \brief calculate the forward mode derivative of area w.r.t. coeff
   * \param[in] x_coord - x coordinates to evaluate area derivative
   * \param[in] u - the direction that the derivative is desired
   * \returns the resulting derivative components
   */
  InnerProdVector AreaForwardDerivative(
      const InnerProdVector & x_coord, const InnerProdVector & u);

  /*!
   * \brief calculate the reverse mode derivative of area w.r.t. coeff
   * \param[in] x_coord - x coordinates to evaluate area derivative
   * \param[in] u - the direction that the derivative is desired
   * \returns the resulting derivative components
   */
  InnerProdVector AreaReverseDerivative(
      const InnerProdVector & x_coord, const InnerProdVector & u);

 protected:

  Bspline spline_;
};
