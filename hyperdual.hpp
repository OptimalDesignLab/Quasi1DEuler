/*
 * Class: HyperDual
 * 
 * Implementation of a class where i^2 = j^2 = (ij)^2 = 0
 * instead of quaternions where i^2 = j^2 = k^2 = -1
 *
 * Only functions required by rangeboom.cpp and JOE have been implemented so far
 * Notes on some functions:
 *    pow(x,a) uses a small tolerance value, tol, to avoid division by zero
 *      Warning: This is not yet implemented in other functions such as asin
 *    fabs(x) uses simply -x if x<0, ignoring derivative discontinuities at x=0
 *
 * Written by: Jeffrey Fike
 * Stanford University, Department of Aeronautics and Astronautics
 * 
 * Last modified:  2/26/09  1:45 pm
 * Added "<<" and unary "+" operators in order to compile JOE
 * Added "const" in the appropriate functions in order to compile JOE
 * Broke into .h and .cpp in order to compile JOE
 * 
 * Last modified:  10/22/10  11:15 am
 * Added acos function
 */

#pragma once

#include <iostream>
//#include <math.h>
//using namespace std;
using std::ostream;

class HyperDual{
  double f0,f1,f2,f12;
 public:
  //create new and set values
  HyperDual();
  HyperDual(double x1,double x2,double x3,double x4);
  HyperDual(double x1);
  void setvalues(double x1,double x2,double x3,double x4);

  //examine values
  void view(void);
  const double & realpart(void) const;
  const double & ipart(void) const;
  const double & jpart(void) const;
  const double & ijpart(void) const;

  // set values
  double & realpart(void);
  double & ipart(void);
  double & jpart(void);
  double & ijpart(void);

  //basic manipulation
  //do not need to overload = operator?  fq=fq and fq=double seem to work
	
  HyperDual operator+ () const;
  HyperDual operator+ (const HyperDual rhs) const;
  friend HyperDual operator+ (const double lhs, const HyperDual rhs);
	
  HyperDual operator- () const;
  HyperDual operator- (const HyperDual rhs) const;
  friend HyperDual operator- (const double lhs, const HyperDual rhs);

  HyperDual operator* (const HyperDual rhs)const;
  friend HyperDual operator* (const double lhs, const HyperDual rhs);

  friend HyperDual operator/ (const HyperDual lhs, const HyperDual rhs);
  friend HyperDual operator/ (const double lhs, const HyperDual rhs);
  friend HyperDual operator/ (const HyperDual lhs, const double rhs);
	
  HyperDual& operator+= (HyperDual rhs);
  HyperDual& operator-= (HyperDual rhs);
  HyperDual& operator*= (HyperDual rhs);
  HyperDual& operator*= (double rhs);
  HyperDual& operator/= (double rhs);

  //math.h functions
  friend HyperDual pow (HyperDual x, double a);
  friend HyperDual pow (HyperDual x, HyperDual a);
  friend HyperDual exp(HyperDual x);
  friend HyperDual log(HyperDual x);
  friend HyperDual sin(HyperDual x);
  friend HyperDual cos(HyperDual x);
  friend HyperDual tan(HyperDual x);
  friend HyperDual asin(HyperDual x);
  friend HyperDual acos(HyperDual x);
  friend HyperDual atan(HyperDual x);
  friend HyperDual sqrt(HyperDual x);
  friend HyperDual fabs(HyperDual x);
  friend HyperDual max(HyperDual x1, HyperDual x2);
  friend HyperDual max(HyperDual x1, double x2);
  friend HyperDual max(double x1, HyperDual x2);
  friend HyperDual min(HyperDual x1, HyperDual x2);
  friend HyperDual min(HyperDual x1, double x2);
  friend HyperDual min(double x1, HyperDual x2);

  //comparisons
  friend bool operator> (HyperDual lhs, HyperDual rhs);
  friend bool operator> (double lhs, HyperDual rhs);
  friend bool operator> (HyperDual lhs, double rhs);
  friend bool operator>= (HyperDual lhs, HyperDual rhs);
  friend bool operator>= (double lhs, HyperDual rhs);
  friend bool operator>= (HyperDual lhs, double rhs);
  friend bool operator< (HyperDual lhs, HyperDual rhs);
  friend bool operator< (double lhs, HyperDual rhs);
  friend bool operator< (HyperDual lhs, double rhs);
  friend bool operator<= (HyperDual lhs, HyperDual rhs);
  friend bool operator<= (double lhs, HyperDual rhs);
  friend bool operator<= (HyperDual lhs, double rhs);
  friend bool operator== (HyperDual lhs, HyperDual rhs);
  friend bool operator== (double lhs, HyperDual rhs);
  friend bool operator== (HyperDual lhs, double rhs);
  friend bool operator!= (HyperDual lhs, HyperDual rhs);
  friend bool operator!= (double lhs, HyperDual rhs);
  friend bool operator!= (HyperDual lhs, double rhs);

  friend ostream& operator<<(ostream& output, const HyperDual& rhs);


};



