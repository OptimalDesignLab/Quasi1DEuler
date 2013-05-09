/**
 * \file exact_solution.hpp
 * \brief functions used to find the exact solution of the nozzle flow
 * \author  Jason Hicken <jason.hicken@gmail.com>
 * \version 1.0
 */

#pragma once

#include <math.h>

#include <ostream>
#include <iostream>
#include <limits>

#include <boost/math/tools/roots.hpp>

template <class T>
class MachRelation {
 public:

  /*!
   * \brief constructor for class
   * \param[in] gam - the heat capacity ratio
   * \parma[in] area_star - the critical area for the nozzle
   * \param[in] area - the area at a particular location
   */
  MachRelation(const T & gam, const T & area_star, const T & area) :
      gamma(gam), Astar(area_star), Ax(area) {}
    
  /*!
   * \brief returns a tuple containing the function and its derivative
   * \param[in] M - the current estimate for the Mach number
   * \result the value of the Mach relation and its derivative at M
   */
  boost::math::tuple<T, T> operator()(const T & M) {
    T fun = (1.0/M)*pow((2.0/(gamma + 1.0))
                        *(1.0 + ((gamma - 1.0)/2.0)*M*M), 
                        (gamma + 1.0)/(2.0*(gamma - 1.0))) - Ax/Astar;
    T dfun = -(fun+Ax/Astar)/M 
        + pow((2.0/(gamma + 1.0))*(1.0 + ((gamma - 1.0)/2.0)*M*M), 
              (3.0 - gamma)/(2.0*(gamma - 1.0)));
    return boost::math::make_tuple(fun, dfun);
  }

 private:
  T gamma; ///< value of the heat capcity ratio
  T Astar; ///< critical area for the nozzle
  T Ax; ///< area at the location of interest
};

/*!
 * \brief returns the exact Mach number solution for a nozzle flow
 * \param[in] gamma - the heat capacity ratio
 * \parma[in] area_star - the critical area for the nozzle
 * \param[in] area - the area at a particular location
 * \param[in] subsonic - chooses the subsonic and supersonic branches
 */
template <class T>
T CalcMachExact(const T & gamma, const T & area_star, const T & area,
                bool subsonic) { 
  T min_range, max_range, guess;
  if (subsonic) { // subsonic branch
    min_range = 0.001;
    max_range = 1.0;
    guess = 0.5;
  } else { // supersonic branch
    min_range = 1.0;
    max_range = 4.0;
    guess = 2.0;
  }
  // Maximum possible binary digits accuracy for type T.
  int digits = std::numeric_limits<T>::digits; 
  return boost::math::tools::
      newton_raphson_iterate(MachRelation<T>(gamma, area_star, area),
                             guess, min_range, max_range, digits);
}

/*!
 * \brief computes the exact density, momentum, and energy
 * \param[in] gamma - the heat capacity ratio
 * \param[in] R_gas - gas constant
 * \parma[in] area_star - the critical area for the nozzle
 * \param[in] area - the area at a particular location
 * \param[in] subsonic - chooses the subsonic and supersonic branches
 * \param[in] temp_stag - stagnation temperature
 * \param[in] press_stag - stagnation pressure
 * \param[out] rho - exact density
 * \param[out] rho_u - exact momentum
 * \param[out] e - exact energy
 */
template <class T>
void CalcFlowExact(const T & gamma, const T & R_gas, 
                   const T & area_star, const T & area,
                   bool subsonic, const T & temp_stag,
                   const T & press_stag, T & rho, T & rho_u, T & e) {
  T Mach = CalcMachExact(gamma, area_star, area, subsonic);
  T temp = 1.0/(1.0 + ((gamma - 1.0)/2.0)*(Mach*Mach));
  T press = press_stag*pow(temp, gamma/(gamma-1.0));
  temp *= temp_stag;
  rho = press/(R_gas*temp);
  T u = Mach*sqrt(gamma*press/rho);
  rho_u = rho*u;
  e = 0.5*rho*u*u + press/(gamma-1.0);
}
