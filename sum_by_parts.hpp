/**
 * \file sum_by_parts.hpp
 * \brief header file for SumByParts class
 * \author  Jason Hicken <jason.hicken@gmail.com>
 * \version 1.0
 */

#pragma once

#include <math.h>

#include <ostream>
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include "./inner_prod_vector.hpp"

namespace ublas = boost::numeric::ublas;

/*!
 * \class SumByParts
 * \brief defines summation by parts operators
 */
class SumByParts {
 public:

  /*!
   * \brief default constructor
   * \param[in] nodes - the number of nodes in the grid
   * \param[in] order - the order of the diagonal norm
   */
  SumByParts(int nodes, int order);

  /*!
   * \brief default destructor
   */
  virtual ~SumByParts() = 0;

  /*!
   * \brief sets the order and size (nodes) of the operator
   * \param[in] nodes - the number of nodes in the grid
   * \param[in] order - the order of the diagonal norm
   *
   * This is useful when changing the order or size of the grid
   */
  void Define(int nodes, int order);

  /*!
   * \brief returns the order of the operator
   * \returns the value of the order_ member
   */
  int order() const { return order_; }

  /*!
   * \brief returns the value of the diagonal norm at the first/last node
   * \result the diagonal norm value at the first/last node
   */
  double Hinv() const { return hinv_(0); }

  /*!
   * \brief returns the H inner product of two vectors
   * \param[in] num_var - number of variables at each node
   * \param[in] u - first vector in the inner product
   * \param[in] v - second vector in the inner product
   * \result the H inner product of u and v
   */
  double InnerProductSBP(const int & num_var,
                         const InnerProdVector & u,
                         const InnerProdVector & v) const;

  /*!
   * \brief returns the L2-norm of a vector using the H weight matrix
   * \param[in] num_var - number of variables at each node
   * \param[in] u - vector whose L2-norm is being calculated
   * \result the L2-norm of u estimated using the H quadrature rule
   */
  double NormSBP(const int & num_var,
                 const InnerProdVector & u) const {
    return sqrt(InnerProductSBP(num_var, u, u));
  }

  /*!
   * \brief multipy a vector by the H norm
   * \param[in] num_var - number of variables at each node
   * \parma[in, out] u - vector that is being multiplied
   * \param[in, out] v - result of the multiplication
   *
   * The routine is in-place so v = u is a possibility
   */
  void HTimesVector(const int & num_var, InnerProdVector & u,
                    InnerProdVector & v) const;

  /*!
   * \brief multipy a vector by the Hinv norm
   * \param[in] num_var - number of variables at each node
   * \parma[in, out] u - vector that is being multiplied
   * \param[in, out] v - result of the multiplication
   *
   * The routine is in-place so v = u is a possibility
   */
  void HinvTimesVector(const int & num_var, InnerProdVector & u,
                    InnerProdVector & v) const;

 protected:

  int num_nodes_; ///< the number of nodes
  int order_; ///< the order of accuracy of the operator
  int numh_; ///< the number of non-unity values in the norm
  ublas::vector<int> ja_; ///< the column offsets, relative to diagonal
  ublas::vector<int> ibeg_; ///< pointer array to beginning of row info
  ublas::vector<int> iend_; ///< pointer array to end of row info
  InnerProdVector as_; ///< the SBP operator matrix entries
  InnerProdVector hinv_; ///< the inverse of the relevant diagonal norm

};
inline SumByParts::~SumByParts() {}

/*!
 * \class SBP1stDerivative
 * \brief defines first-derivative SBP operators
 */
class SBP1stDerivative : public SumByParts {
 public:

  /*!
   * \brief default constructor
   * \param[in] nodes - the number of nodes in the grid
   * \param[in] order - the order of the operator
   */
  SBP1stDerivative(int nodes, int order);

  /*!
   * \brief default destructor
   */
  ~SBP1stDerivative() {}

  /*!
   * \brief sets the order and size (nodes) of the operator
   * \param[in] nodes - the number of nodes in the grid
   * \param[in] order - the order of the diagonal norm
   *
   * This is useful when changing the order or size of the grid
   */
  void Define(int nodes, int order);

  /*!
   * \brief applies the SBP matrix operator to a given vector
   * \param[in] num_var - number of variables at each node
   * \param[in] u - vector which is being multiplied by SBP operator
   * \param[out] v - result of matrix vector product
   */
  void Apply(const int & num_var, const InnerProdVector & u,
                  InnerProdVector & v) const;

  /*!
   * \brief applies D^{T} = Q^{T} H^{-1} to a given vector
   * \param[in] num_var - number of variables at each node
   * \param[in] u - vector which is being multiplied by transponsed SBP
   * \param[out] v - result of matrix vector product
   */
  void ApplyTranspose(const int & num_var, const InnerProdVector & u,
                  InnerProdVector & v) const;
};

/*!
 * \class SBPDissipation
 * \brief defines numerical dissipation SBP operators
 */
class SBPDissipation : public SumByParts {
 public:

  /*!
   * \brief default constructor
   * \param[in] nodes - the number of nodes in the grid
   * \param[in] norm_order - the order of the diagonal norm
   * \param[in] order - the order of the operator
   */
  SBPDissipation(int nodes, int norm_order, int order);

  /*!
   * \brief default destructor
   */
  ~SBPDissipation() {}

  /*!
   * \brief sets the order and size (nodes) of the operator
   * \param[in] nodes - the number of nodes in the grid
   * \param[in] norm_order - the order of the diagonal norm
   * \param[in] order - the order of the diagonal norm
   *
   * This is useful when changing the order or size of the grid
   */
  void Define(int nodes, int norm_order, int order);

  /*!
   * \brief applies the SBP dissipation to a given vector
   * \param[in] num_var - number of variables at each node
   * \param[in] u - vector which is being multiplied by SBP operator
   * \param[in] b - diagonal of positive scaling values (e.g spec. rad.)
   * \param[out] v - result of matrix vector product
   */
  void Apply(const int & num_var, const InnerProdVector & u,
                  const InnerProdVector & b,
                  InnerProdVector & v) const;

  /*!
   * \brief applies transpose of SBP dissipation to a given vector
   * \param[in] num_var - number of variables at each node
   * \param[in] u - vector which is being multiplied by SBP operator
   * \param[in] b - diagonal of positive scaling values (e.g spec. rad.)
   * \param[out] v - result of matrix vector product
   */
  void ApplyTranspose(const int & num_var, const InnerProdVector & u,
                  const InnerProdVector & b,
                  InnerProdVector & v) const;

};
