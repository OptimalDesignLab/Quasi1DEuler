/**
 * \file inner_prod_vector.hpp
 * \brief header file for InnerProdVector class
 * \author  Jason Hicken <jason.hicken@gmail.com>
 * \version 1.0
 */

#pragma once

#include <ostream>
#include <iostream>
#include <fstream>
#include <string>

#include <boost/numeric/ublas/vector.hpp>

using std::fstream;
using std::ofstream;
using std::ifstream;
using std::string;
using std::ios;
using std::cerr;
using std::endl;
namespace ublas = boost::numeric::ublas;

/*!
 * \class InnerProdVector
 * \brief defines a vector from a finite-dim inner-product space
 *
 * This is simply a front end for ublas::vector<double>, with some
 * additional functionality that we need for the Krylov iterative
 * solvers.  Also provides a simplified notation.
 */
class InnerProdVector : public ublas::vector<double> {
 public:

  /*!
   * \brief default constructor, creates an empty vector
   */
  InnerProdVector() : ublas::vector<double>() {}

  /*!
   * \brief constructor, creates and initializes a vector of given size 
   * \param[in] size - number of elements in vector
   * \param[in] val - all elements are given value val
   */
  InnerProdVector(const int & size, const double & val = 0.0) : 
      ublas::vector<double>(size, val) {}

  /*!
   * \brief constructor, creates a vector from a ublas::vector<double>
   * \param[in] u - ublas vector to initialize vector to
   */
  InnerProdVector(const ublas::vector<double> & u) : 
      ublas::vector<double>(u) {}

  /*!
   * \brief copy constructor
   * \param[in] u - existing vector that we want to copy
   */
  InnerProdVector(const InnerProdVector & u) :
      ublas::vector<double>(static_cast<ublas::vector<double> >(u)) {}

  /*!
   * \brief default destructor
   */
  ~InnerProdVector() {}

  /*!
   * \brief inner product between two InnerProdVectors
   * \param[in] u - first InnerProdVector in dot product
   * \param[in] v - second InnerProdVector in dot product
   */
  friend double InnerProd(const InnerProdVector & u, 
                   const InnerProdVector & v) {
    return ublas::inner_prod(u, v); 
  }

  /*!
   * \brief the L2 norm of the InnerProdVector
   */
  double Norm2() const {
    return norm_2(*this);
  }

  /*!
   * \brief general linear combination of two InnerProdVectors
   * \param[in] a - scalar factor for x
   * \param[in] x - first InnerProdVector in linear combination
   * \param[in] b - scalar factor for y
   * \param[in] y - second InnerProdVector in linear combination
   */
  void EqualsAXPlusBY(const double & a, const InnerProdVector & x,
                      const double & b, const InnerProdVector & y) {
    for (int i = 0; i < size(); i++)
      this->operator()(i) = a*x(i) + b*y(i);
  }
  
  /*!
   * \brief assign all elements of a vector the same value
   * \param[in] val - scalar value that we want the elements to be set to
   */
  void operator=(const double & val) {
    ublas::vector<double>::operator=(
        ublas::scalar_vector<double>(size(), val));
  }

  void BinaryWrite(ofstream & fout) const {
    if (!fout.is_open()) {
      cerr << "InnerProdVector::BinaryWrite(): fout is not open!" << endl;
      throw(-1);
    }
    // using &data() causes problems for BinaryRead
    fout.write(reinterpret_cast<const char*>(&operator[](0)),
               size()*sizeof(double));
  }

  void BinaryWrite(ofstream & fout, const unsigned long & ptr) const {
    if (!fout.is_open()) {
      cerr << "InnerProdVector::BinaryWrite(): fout is not open!" << endl;
      throw(-1);
    }
    // using &data() causes problems for BinaryRead
    fout.seekp(ptr);
    fout.write(reinterpret_cast<const char*>(&operator[](0)),
               size()*sizeof(double));
  }

  void BinaryWrite(fstream & fout) const {
    if (!fout.is_open()) {
      cerr << "InnerProdVector::BinaryWrite(): fout is not open!" << endl;
      throw(-1);
    }
    // using &data() causes problems for BinaryRead
    fout.write(reinterpret_cast<const char*>(&operator[](0)),
               size()*sizeof(double));
  }

  void BinaryWrite(fstream & fout, const unsigned long & ptr) const {
    if (!fout.is_open()) {
      cerr << "InnerProdVector::BinaryWrite(): fout is not open!" << endl;
      throw(-1);
    }
    // using &data() causes problems for BinaryRead
    fout.seekp(ptr);
    fout.write(reinterpret_cast<const char*>(&operator[](0)),
               size()*sizeof(double));
  }

  void BinaryRead(ifstream & fin) {
    if (!fin.is_open()) {
      cerr << "InnerProdVector::BinaryRead(): fin is not open!" << endl;
      throw(-1);
    }
    // using &data() causes problems for BinaryRead
    fin.read(reinterpret_cast<char*>(&operator[](0)), size()*sizeof(double));
  }

  void BinaryRead(ifstream & fin, const unsigned long & ptr) {
    if (!fin.is_open()) {
      cerr << "InnerProdVector::BinaryRead(): fin is not open!" << endl;
      throw(-1);
    }
    // using &data() causes problems for BinaryRead
    fin.seekg(ptr);
    fin.read(reinterpret_cast<char*>(&operator[](0)), size()*sizeof(double));
  }

};
