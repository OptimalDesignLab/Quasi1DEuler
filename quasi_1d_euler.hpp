/**
 * \file quasi_1d_euler.hpp
 * \brief header file for Quasi1DEuler
 * \author  Jason Hicken <jason.hicken@gmail.com>
 * \version 1.0
 */

#pragma once

#include <math.h>

#include <ostream>
#include <iostream>
#include <complex>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <krylov.hpp>

#include "./inner_prod_vector.hpp"
#include "./sum_by_parts.hpp"
#include "./hyperdual.hpp"

namespace ublas = boost::numeric::ublas;

typedef std::complex<double> complex;

const double kGamma = 1.4; ///< specific heat ratio

enum objective {
  inverse = 0,       ///< inverse design objective
  total_kinetic = 1  ///< total kinetic energy objective
};

enum norm_type {
  inf = 0,           ///< infinity (max) norm
  L2 = 2             ///< L2 norm
};

/*! 
 * \brief complex version of fabs for complexified functions
 * \param[in] z - complex variable whose absolute value is being taken
 */
complex fabs(const complex & z);

// ======================================================================

/*!
 * \class JacobianVectorProduct
 * \brief specialization of matrix-vector product for InnerProdVectors
 */
class Quasi1DEuler;
class JacobianVectorProduct :
    public kona::MatrixVectorProduct<InnerProdVector> {
 public:

  /*!
   * \brief default constructor
   * \param[in] euler_solver - a Quasi1DEuler solver (defines product)
   */
  JacobianVectorProduct(Quasi1DEuler * euler_solver) {
    solver = euler_solver; 
  } 

  ~JacobianVectorProduct() {} ///< class destructor

  /*!
   * \brief operator that defines the Jacobian-Vector product
   * \param[in] u - vector that is being multiplied by the Jacobian
   * \param[out] v - vector that is the result of the product
   */
  void operator()(const InnerProdVector & u, InnerProdVector & v);

 private:
  Quasi1DEuler * solver; ///< used to access the Jacobian-Vector routine
};

// ======================================================================

/*!
 * \class ApproxJacobian
 * \brief specialization of preconditioner for InnerProdVectors
 */
class ApproxJacobian : 
    public kona::Preconditioner<InnerProdVector> {
 public:

  /*!
   * \brief default constructor
   * \param[in] euler_solver - a Quasi1DEuler solver to access precond.
   */
  ApproxJacobian(Quasi1DEuler * euler_solver) {
    solver = euler_solver; 
  }

  ~ApproxJacobian() {} ///< class destructor

  /*!
   * \brief operator that applies the approximate Jacobian LU decomp.
   * \param[in] u - vector that is being preconditioned
   * \param[out] v - vector that is the result of the preconditioning
   */
  void operator()(InnerProdVector & u, InnerProdVector & v);

 private:
  Quasi1DEuler * solver; ///< used to access the preconditioner
};

// ======================================================================

/*!
 * \class JacobianTransposedVectorProduct
 * \brief specialization of matrix-vector product for InnerProdVectors
 */
class Quasi1DEuler;
class JacobianTransposedVectorProduct :
    public kona::MatrixVectorProduct<InnerProdVector> {
 public:

  /*!
   * \brief default constructor
   * \param[in] euler_solver - a Quasi1DEuler solver (defines product)
   */
  JacobianTransposedVectorProduct(Quasi1DEuler * euler_solver) {
    solver = euler_solver; 
  } 

  ~JacobianTransposedVectorProduct() {} ///< class destructor

  /*!
   * \brief operator that defines the Jacobian-transposed-Vector product
   * \param[in] u - vector that is being multiplied by the Jacobian
   * \param[out] v - vector that is the result of the product
   */
  void operator()(const InnerProdVector & u, InnerProdVector & v);

 private:
  Quasi1DEuler * solver; ///< used to access Jacobian-transposed product
};

// ======================================================================

/*!
 * \class ApproxJacobianTransposed
 * \brief specialization of preconditioner for InnerProdVectors
 */
class ApproxJacobianTransposed : 
    public kona::Preconditioner<InnerProdVector> {
 public:

  /*!
   * \brief default constructor
   * \param[in] euler_solver - a Quasi1DEuler solver to access precond.
   */
  ApproxJacobianTransposed(Quasi1DEuler * euler_solver) {
    solver = euler_solver;
  }

  ~ApproxJacobianTransposed() {} ///< class destructor

  /*!
   * \brief operator that solves the approximate Jacobian transposed
   * \param[in] u - vector that is being preconditioned
   * \param[out] v - vector that is the result of the preconditioning
   */
  void operator()(InnerProdVector & u, InnerProdVector & v);

 private:
  Quasi1DEuler * solver; ///< used to access the preconditioner
};

// ======================================================================

/*!
 * \class UnsteadyJacobianVectorProduct
 * \brief specialization of matrix-vector product for InnerProdVectors
 */
class UnsteadyJacobianVectorProduct :
    public kona::MatrixVectorProduct<InnerProdVector> {
 public:

  /*!
   * \brief default constructor
   * \param[in] euler_solver - a Quasi1DEuler solver (defines product)
   */
  UnsteadyJacobianVectorProduct(Quasi1DEuler * euler_solver) {
    solver = euler_solver; 
  } 

  ~UnsteadyJacobianVectorProduct() {} ///< class destructor

  /*!
   * \brief operator that defines the unsteady Jacobian-Vector product
   * \param[in] u - vector that is being multiplied by the Jacobian
   * \param[out] v - vector that is the result of the product
   */
  void operator()(const InnerProdVector & u, InnerProdVector & v);

 private:
  Quasi1DEuler * solver; ///< used to access the Jacobian-Vector routine
};

// ======================================================================

/*!
 * \class UnsteadyJacTransVectorProduct
 * \brief specialization of matrix-vector product for InnerProdVectors
 */
class UnsteadyJacTransVectorProduct :
    public kona::MatrixVectorProduct<InnerProdVector> {
 public:

  /*!
   * \brief default constructor
   * \param[in] euler_solver - a Quasi1DEuler solver (defines product)
   */
  UnsteadyJacTransVectorProduct(Quasi1DEuler * euler_solver) {
    solver = euler_solver; 
  }

  ~UnsteadyJacTransVectorProduct() {} ///< class destructor

  /*!
   * \brief operator that defines the unsteady transposed Jacobian-Vector
   * \param[in] u - vector that is being multiplied by the Jacobian
   * \param[out] v - vector that is the result of the product
   */
  void operator()(const InnerProdVector & u, InnerProdVector & v);

 private:
  Quasi1DEuler * solver; ///< used to access the transposed product routine
};

// ======================================================================

/*!
 * \class Quasi1DEuler
 * \brief defines and solves a steady quasi-1d Euler problem
 */
class Nozzle;
class Quasi1DEuler {
 public:

  /*! 
   * \brief default constructor
   */
  //Quasi1DEuler() {
  //  num_nodes_ = -1;    
  //}

  /*!
   * \brief class constructor
   * \param[in] num_nodes - number of nodes in the domain
   * \param[in] order - order of accuracy of SBP operator
   */
  Quasi1DEuler(int num_nodes, int order) : 
      sbp_deriv_(num_nodes, order),
      sbp_diss_(num_nodes, order, order),      
      prec_(3*num_nodes, 3*num_nodes, 5, 5),
      x_coord_(num_nodes, 0.0),
      met_jac_(num_nodes, 0.0),
      area_(num_nodes, 0.0),
      press_(num_nodes, 0.0),
      press_targ_(num_nodes, 0.0),
      sndsp_(num_nodes, 0.0),
      spect_(num_nodes, 0.0),
      q_(3*num_nodes, 1.0),
      q_old_(3*num_nodes, 1.0),
      res_(3*num_nodes, 0.0),
      psi_(3*num_nodes, 0.0),
      bc_left_(),
      bc_right_() {
    num_nodes_ = num_nodes;
    dxi_ = 1.0/static_cast<double>(num_nodes_-1);
    ClearPreconditionerCallCounts();
  }
  
  /*!
   * \brief default destructor
   */
  ~Quasi1DEuler() {}

  /*!
   * \brief returns the number of nodes in the mesh
   * \returns num_nodes_ member value
   */
  const int & get_num_nodes() const { return num_nodes_;}

  /*!
   * \brief returns a const reference to the mesh spacing
   */
  const double & dxi() const {
    return dxi_;
  }

  /*!
   * \brief returns a const reference to the time step
   */
  const double & dt() const {
    return dt_;
  }
  
  /*!
   * \brief returns the number of calls to the primal preconditioner
   * \returns num_primal_precond_ member value
   */
  int get_num_primal_precond() { return num_primal_precond_; }

  /*!
   * \brief returns the number of calls to the adjoint preconditioner
   * \returns num_adjoint_precond_ member value
   */
  int get_num_adjoint_precond() { return num_adjoint_precond_; }

  /*!
   * \brief get the sum of primal and adjoint preconditioner calls
   * \returns sum of num_adjoint_precond_ and num_adjoint_adjoint_
   */
  int TotalPreconditionerCalls() {
    return num_primal_precond_ + num_adjoint_precond_;
  }

  /*!
   * \breif sets the number of calls to the preconditioners to zero
   */
  void ClearPreconditionerCallCounts() {
    num_primal_precond_ = 0;
    num_adjoint_precond_ = 0;
  }

  /*!
   * \brief set the coordinates of the nodes
   * \param[in] coord - values used to set coordinates
   */
  void set_x_coord(InnerProdVector & coord) { 
    x_coord_ = coord;
    sbp_deriv_.Apply(1, x_coord_, met_jac_);
  }

  /*!
   * \brief returns a const ref to the coordinates of the nodes
   * \returns a vector of nodal x coordinates
   */
  const InnerProdVector & get_x_coord() const { return x_coord_; }

  /*!
   * \brief set the area of the nozzle at each node
   * \param[in] A - values used to set areas
   */
  void set_area(const InnerProdVector & A) { area_ = A; }

  /*!
   * \brief set the momentum source for each time
   * \param[in] x - x location (origin) of source
   * \param[in] sigma - width (std) of source
   * \param[in] source - values used to set source (between time steps)
   */
  void set_source(const double & x, const double & sigma, 
                  const InnerProdVector & source) { 
    src_x_ = x;
    src_sig2_ = sigma*sigma;
    src_ = source; 
  }

  /*!
   * \brief set the value of the numerical dissipation coefficient
   * \param[in] val - value of the coefficient
   */
  void set_diss_coeff(const double & val) { diss_coeff_ = val; }
  
  /*!
   * \brief set the boundary conditions at the left boundary
   * \param[in] rho - value of the density at the left boundary
   * \param[in] rho_u - value of the momentum at the left boundary
   * \param[in] e - value of the total energy per unit volume 
   */
  void set_bc_left(const double & rho, const double & rho_u, 
                   const double & e) { 
    bc_left_(0) = rho;
    bc_left_(1) = rho_u;
    bc_left_(2) = e;
  }

  /*!
   * \brief set the boundary conditions at the right boundary
   * \param[in] rho - value of the density at the right boundary
   * \param[in] rho_u - value of the momentum at the right boundary
   * \param[in] e - value of the total energy per unit volume 
   */
  void set_bc_right(const double & rho, const double & rho_u, 
                   const double & e) { 
    bc_right_(0) = rho;
    bc_right_(1) = rho_u;
    bc_right_(2) = e;
  }

  /*!
   * \brief sets the flow to a given state
   * \param[in] q_new - desired state that the flow is set to
   *
   * This member function is useful for evaluating Jacobian-vector
   * products at the point of linearization given by q_new
   */
  void set_q(const InnerProdVector & q_new) { q_ = q_new; }

  /*!
   * \brief sets the previous flow state to a given state
   * \param[in] q_old_new - desired state that the flow is set to
   *
   * This member function is useful for evaluating Jacobian-vector
   * products at the point of linearization given by q_old_new
   */
  void set_q_old(const InnerProdVector & q_old_new) { q_old_ = q_old_new; }

  /*!
   * \brief sets the target pressure for inverse design problems
   * \param[in] press_targ - target pressure
   */
  void set_press_targ(const InnerProdVector & press_targ) {
    press_targ_ = press_targ;
  }

  /*!
   * \brief return a const reference to the solution vector
   * \returns the solution vector
   */
  const InnerProdVector & get_q() const { return q_; }

  /*!
   * \brief return a const reference to the previous solution vector
   * \returns the previous solution vector
   */
  const InnerProdVector & get_q_old() const { return q_old_; }

  /*!
   * \brief return a const reference to the adjoint solution vector
   * \returns the adjoint solution vector
   */
  const InnerProdVector & get_psi() const { return psi_; }
  
  /*!
   * \brief return a const reference to the residual vector
   * \returns the residual vector
   */
  const InnerProdVector & get_res() const { return res_; }

  /*!
   * \brief return a const reference to the pressure
   * \returns the pressure at each node
   */
  const InnerProdVector & get_press() const { return press_; }
  
  /*!
   * \brief set a uniform flow based on the input values
   * \param[in] rho - value of the density 
   * \param[in] rho_u - value of the momentum
   * \param[in] e - value of the total energy per unit volume 
   */
  void InitialCondition(const double & rho, const double & rho_u,
                        const double & e) {
    for (int i = 0; i < num_nodes_; i++) {
      q(i,0) = rho; 
      q(i,1) = rho_u;
      q(i,2) = e;
    }
  }

  void InitialCondition(const InnerProdVector & q_init) {
    q_old_ = q_init;
    q_ = q_init;
  }

  /*!
   * \brief sets the adjoint variable (useful for partitions in time)
   * \param[in] psi_init - initial value for the adjoint
   */
  void AdjointInitialCondition(const InnerProdVector & psi_init) {
    psi_ = psi_init;
  }
  
  /*!
   * \brief resizes the number of nodes on the grid
   * \param[in] coord - the new grid coordinates
   */
  void ResizeGrid(const InnerProdVector & coord);

  /*!
   * \brief calculate the residual vector
   * \result residual based on q_ is calculated and stored in res_
   */
  void CalcResidual();

  /*!
   * \brief calculate the unsteady residual vector
   * \result residual based on q_ and q_old is calculated and stored in res_
   */
  void CalcUnsteadyResidual();

  /*!
   * \brief adds the momentum source to the residual vector
   * \param[in] n - indicates the time iteration that is begin added [n,n+1]
   */
  void AddUnsteadySource(const int & n);

  /*!
   * \brief calculates the product (dpress/dq)*u
   * \param[in] u - vector that is being multiplied (size = 3*num_nodes)
   * \param[out] v - vector that holds resulting product (size = num_nodes)
   */
  void CalcDPressDQProduct(const InnerProdVector & u,
                           InnerProdVector & v);

  /*!
   * \brief calculates the product (dpress/dq)^{T}*u
   * \param[in] u - vector that is being multiplied (size = num_nodes)
   * \param[out] v - vector that holds resulting product (size = 3*num_nodes)
   */
  void CalcDPressDQTransposedProduct(const InnerProdVector & u,
                                     InnerProdVector & v);

  /*!
   * \brief checks the products involving dPress/dQ, using v^T*(dPress/dQ)*u
   */
  void TestDPressDQProducts();
  
  /*!
   * \brief calculate the Jacobian-vector product
   * \param[in] u - vector that is being multiplied
   * \param[out] v - vector that holds the resulting product
   * \pre the state held in q_ determines the Jacobian
   */
  void JacobianStateProduct(const InnerProdVector & u, 
                            InnerProdVector & v);
  
  /*! 
   * \brief tests the JacobianStateProduct() routine using FD
   */
  void TestJacobianStateProduct();
  
  /*!
   * \brief calculate the Jacobian-transposed-vector product
   * \param[in] u - vector that is being multiplied
   * \param[out] v - vector that holds the resulting product
   * \pre the state held in q_ determines the Jacobian
   */
  void JacobianTransposedStateProduct(const InnerProdVector & u,
                                      InnerProdVector & v);

  /*! 
   * \brief tests the JacobianTransposedStateProduct() routine
   */
  void TestJacobianTransposedStateProduct();

  /*!
   * \brief calculate the (unsteady) Jacobian-vector product
   * \param[in] u - vector that is being multiplied
   * \param[out] v - vector that holds the resulting product
   * \param[in] plus_time - if true, I*u is added, otherwise, I*u is subtracted
   * \pre the states held in q_ and q_old_ determines the Jacobian
   */
  void UnsteadyJacobianStateProduct(const InnerProdVector & u,
                                    InnerProdVector & v,
                                    const bool & plus_time = true);

  /*!
   * \brief calculate the (unsteady) approximate Jacobian-vector product
   * \param[in] u - vector that is being multiplied
   * \param[out] v - vector that holds the resulting product
   * \param[in] plus_time - if true, I*u is added, otherwise, I*u is subtracted
   * \pre the states held in q_ and q_old_ determines the approximate Jacobian
   */
  void UnsteadyApproxJacStateProduct(const InnerProdVector & u,
                                     InnerProdVector & v,
                                     const bool & plus_time = true);
  
  /*!
   * \brief calculate the (unsteady) transposed Jacobian-vector product
   * \param[in] u - vector that is being multiplied
   * \param[out] v - vector that holds the resulting product
   * \param[in] plus_time - if true, I*u is added, otherwise, I*u is subtracted
   * \pre the states held in q_ and q_old_ determines the Jacobian
   */
  void UnsteadyJacTransStateProduct(const InnerProdVector & u,
                                    InnerProdVector & v, 
                                    const bool & plus_time = true);

  /*!
   * \brief calculate the (unsteady) transposed approximate Jacobian-vector prod
   * \param[in] u - vector that is being multiplied
   * \param[out] v - vector that holds the resulting product
   * \param[in] plus_time - if true, I*u is added, otherwise, I*u is subtracted
   * \pre the states held in q_ and q_old_ determines the approximate Jacobian
   */
  void UnsteadyApproxJacTransStateProduct(const InnerProdVector & u,
                                          InnerProdVector & v,
                                          const bool & plus_time = true);

  /*! 
   * \brief tests the UnsteadyJacTransStateProduct() routine
   */
  void TestUnsteadyJacTransStateProduct();

  /*!
   * \brief calculate the residual Hessian-vector product
   * \param[in] psi - vector that left-multiplies 
   * \param[in] w - vector that right-multiplies 
   * \param[out] v - vector that holds the resulting product
   * \pre the state held in q_ determines the Hessian
   *
   * Note that w is not a const reference so that we can use a ublas
   * vector_range on it inside the routine.
   */
  void ResidualHessianProduct(const InnerProdVector & psi,
                              InnerProdVector & w,
                              InnerProdVector & v);
  
  /*! 
   * \brief tests the ResidualHessianProduct() routine
   */
  void TestResidualHessianProduct();

  /*! 
   * \brief constructs an LU-preconditioner based on a first-order approx
  */
  void BuildAndFactorPreconditioner();

  /*! 
   * \brief constructs an LU-preconditioner based on a first-order approx
   * \param[in] factor - if true factors the matrix, otherwise does not factor
  */
  void BuildAndFactorUnsteadyPreconditioner(const bool & factor = true);

  /*!
   * \brief apply the first-order LU-preconditioner to a vector
   * \param[in] u - vector that is being preconditioned
   * \param[out] v - preconditioned vector
   */
  void Precondition(const InnerProdVector & u, 
                    InnerProdVector & v);

  /*!
   * \brief apply the transposed first-order LU-preconditioner
   * \param[in] u - vector that is being preconditioned
   * \param[out] v - preconditioned vector
   */
  void PreconditionTransposed(const InnerProdVector & u,
                              InnerProdVector & v);

  /*!
   * \brief multiply a vector by the first-order LU-preconditioner
   * \param[in] u - vector that is being multiplied
   * \param[out] v - product
   * 
   * For testing
   */
  void PreconditionerMultiply(const InnerProdVector & u, 
                              InnerProdVector & v);

  /*!
   * \brief calculate the Jacobian (w.r.t the area) product
   * \param[in] u - vector that is being multiplied
   * \param[out] v - vector that holds the resulting product
   * \pre q_ and area_ determine the Jacobian
   */
  void JacobianAreaProduct(const InnerProdVector & u, 
                           InnerProdVector & v);

  /*!
   * \brief calculate the Jacobian (w.r.t the area) transposed product
   * \param[in] u - vector that is being multiplied
   * \param[out] v - vector that holds the resulting product
   * \pre q_ and area_ determine the Jacobian
   */
  void JacobianTransposedAreaProduct(const InnerProdVector & u,
                                     InnerProdVector & v);

  /*! 
   * \brief tests JacobianAreaProduct() and JacobianTransposedProduct()
   */
  void TestJacobianAreaProducts();

  /*!
   * \brief solves for the flow using explicit Euler time marching
   * \param[in] max_iter - maximum number of iterations permitted
   * \param[in] target_cfl - a target CFL number to use
   * \param[in] tol - tolerance with which to solve the system
   */
  void ExplicitEuler(const int & max_iter, const double & target_cfl, 
                     const double & tol);

  /*!
   * \brief solves for the flow using a Newton-Krylov algorithm
   * \param[in] max_iter - maximum number of iterations permitted
   * \param[in] tol - tolerance with which to solve the system
   * \returns - total number of preconditioner calls
   */
  int NewtonKrylov(const int & max_iter, const double & tol);

  /*!
   * \brief solves for the adjoint variables using a Krylov solver
   * \param[in] max_iter - maximum number of iterations permitted
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] dJdQ - the rhs of the adjoint linear system
   * \returns total number of preconditioner calls
   */
  int SolveAdjoint(const int & max_iter, const double & tol,
                    const InnerProdVector & dJdQ,
                    InnerProdVector & psi);
  
  /*!
   * \brief solves the linearized state equation using a Krylov solver
   * \param[in] max_iter - maximum number of iterations permitted
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] rhs - the rhs of the linearized system
   * \returns total number of preconditioner calls
   */
  int SolveLinearized(const int & max_iter, const double & tol,
                      const InnerProdVector & rhs,
                      InnerProdVector & dq);

  /*!
   * \brief solves an unsteady problem using the midpoint rule
   * \param[in] iter - number of time steps
   * \param[in] Time - period of time to solve
   * \param[in] tol - tolerance with which to solve at each iteration
   * \param[in] store - if true, the solution is stored to file (for adjoint)
   * \param[in] flow_file - name of the file used to store the flow
   * \param[in] tec_write - if true, a Tecplot solution file is written
   */
  int SolveUnsteady(const int & iter, const double & Time, const double & tol,
                    const bool & store = false, 
                    const string & flow_file = string("save_flow.bin"),
                    const bool & tec_write = false);

  /*!
   * \brief solve an unsteady adjoint problem using the midpoint rule
   * \param[in] flow_file - name of the binary file storing the flow
   * \param[in] max_krylov - maximum number of Krylov iterations allowed per step
   * \param[in] tol - tolerance with which to solve at each iteration
   * \param[in] dJdQ_file - name of the binary file storing the dJ/dQ derivative
   * \param[in] psi_file - name of the binary file where the adjoint is stored
   * \param[in] init_adjoint - if true, the value in psi_ is used as I.C.
   */
  int SolveUnsteadyAdjoint(const string & flow_file, const int & max_krylov, 
                           const double & tol, const string & dJdQ_file,
                           const string & psi_file,
                           const bool & init_adjoint = false);

  /*!
   * \brief solves one iteration of the unsteady linearized state equation
   * \param[in] max_iter - maximum number of iterations permitted
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] rhs - the rhs of the linearized system
   * \returns total number of preconditioner calls
   * \pre dt_, q_old_ and q_ are defined
   *
   * This solves the linearized system at ONE particular iteration, not the
   * complete unsteady system.
   */  
  int SolveUnsteadyIterLinearized(const int & max_iter,
                                  const double & tol,
                                  const InnerProdVector & rhs,
                                  InnerProdVector & dq);
  
  /*!
   * \brief solves one iteration of the unsteady adjoint equation
   * \param[in] max_iter - maximum number of iterations permitted
   * \param[in] tol - tolerance with which to solve the system
   * \param[in] rhs - the rhs of the adjoint system
   * \returns total number of preconditioner calls
   * \pre dt_, q_old_ and q_ are defined
   *
   * This solves the unsteady adjont system at ONE particular iteration, not the
   * complete unsteady adjoint system.
   */
  int SolveUnsteadyAdjointIter(const int & max_iter,
                               const double & tol,
                               const InnerProdVector & rhs,
                               InnerProdVector & psi);

  /*!
   * \brief write the solution to a Tecplot file
   * \param[in] rhoL - the density at the inlet (for dimensionalization)
   * \param[in] aL - sound speed at the inlet (for dimensionalization)
   */
  void WriteTecplot(const double & rhoL, const double & aL,
                    const string & filename = "quasi1d.dat");

  /*!
   * \brief writes the solution for the current time step (q_) to file
   * \param[in] iter - current iteration
   * \param[in] dt - time step
   * \param[in] fout - output stream
   */
  void WriteUnsteadyTecplot(const int & iter, const double & dt,
                            ostream & fout);

  /*!
   * \brief writes the solution stored in flow_file to a Tecplot data file
   * \param[in] flow_file - a binary file name storing the unsteady solution
   * \param[in] tec_file - the output file name to write the solution to
   */
  void WriteUnsteadyTecplot(const string & flow_file,
                            const string & tec_file);

  /*!
   * \brief writes the adjoint solution stored in psi_file to a Tecplot data file
   * \param[in] psi_file - a binary file name storing the unsteady adjoint
   * \param[in] tec_file - the output file name to write the solution to
   */
  void WriteUnsteadyAdjointTecplot(const string & psi_file,
                                   const string & tec_file);

  /*!
   * \brief writes the given vector to fout in binar
   * \param[in] fout - file to write to
   * \param[in] q_var - solution vector to write
   */
  void SaveSolution(ofstream & fout, const InnerProdVector & q_var) const; 

  /*!
   * \brief computes the error in the Mach number using exact solution
   * \param[in] area_star - the critical area for the nozzle
   * \param[in] subsonic - is the flow entirely subsonic
   * \param[out] L2_error - the L2 error in the Mach number
   * \param[out] max_error - the infinity norm error in the Mach number
   */
  void CalcMachError(const double & area_star, const bool & subsonic,
                     double & L2_error, double & max_error);

  /*!
   * \brief computes the integral of the energy over the whole domain
   * \param[in] sbp_quad - if true, uses SBP quadrature, else Simpsons
   * \returns the energy integral computed using SBP-norm quadrature
   */
  double CalcTotalEnergy(const bool & sbp_quad);

  /*!
   * \brief calculates gradient w.r.t. Q of the total energy functional 
   * \param[out] dJdQ - gradient of the objective
   */
  void CalcTotalEnergydJdQ(InnerProdVector & dJdQ);

  /*!
   * \brief calculates the Hessian of J w.r.t. Q and multiplies by w
   * \param[in] w - the vector multiplying d2JdQ2
   * \param[out] dJ2dQ2 - the product d2J/dQ2 * w
   */
  void CalcTotalEnergyd2JdQ2(const InnerProdVector & w,
                             InnerProdVector & dJdQ);

  /*!
   * \brief computes the pressure inverse design objective
   */
  double CalcInverseDesign();

  /*!
   * \brief calculates gradient of the pressure inverse design objective
   * \param[out] dJdQ - gradient of the objective
   */
  void CalcInverseDesigndJdQ(InnerProdVector & dJdQ);

  /*!
   * \brief calculates the Hessian of J w.r.t. Q and multiplies by w
   * \param[in] w - the vector multiplying d2JdQ2
   * \param[out] dJ2dQ2 - the product d2J/dQ2 * w
   */
  void CalcInverseDesignd2JdQ2(const InnerProdVector & w,
                               InnerProdVector & dJdQ);

  /*!
   * \brief Tests the accuracy of the d2J/dQ2*w product using FD
   * \param[in] obj - the desired objective function
   */
  void Testd2JdQ2(const objective & obj);

  /*!
   * \brief calculates the pressure sensor (unsteady) objective
   * \param[in] x_center - x coordinate where the sensor is centered
   * \param[in] sigma - standard deviation of the sensor kernel
   * \param[in] press_targ - the constant target pressure
   * \param[in] reg_param - regularization parameter
   * \param[in] flow_file - name of the file that stores solution
   */
  double CalcSensor(const double & x_center, const double & sigma,
                    const double & press_targ, const double & reg_param,
                    const string & flow_file);

  /*!
   * \brief calculates the (partial) derivative of the sensor objective
   * \param[in] x_center - x coordinate where the sensor is centered
   * \param[in] sigma - standard deviation of the sensor kernel
   * \param[in] press_targ - the constant target pressure
   * \param[in] reg_param - regularization parameter
   * \param[in] flow_file - name of the file that stores solution
   * \param[in] dJdQ_file - name of the file to store the derivative in
   */
  void CalcSensordJdQ(const double & x_center, const double & sigma,
                      const double & press_targ, const double & reg_param,
                      const string & flow_file,
                      const string & dJdQ_file);

  /*!
   * \brief tests the routine CalcSensordJdQ
   * \param[in] x_center - x coordinate where the sensor is centered
   * \param[in] sigma - standard deviation of the sensor kernel
   * \param[in] press_targ - the constant target pressure
   * \param[in] reg_param - regularization parameter
   * \param[in] flow_file - name of the file that stores solution
   */
  void TestSensordJdQ(const double & x_center, const double & sigma,
                      const double & press_targ, const double & reg_param,
                      const string & flow_file);

  /*!
   * \brief calculates the gradient of a specified objective function
   * \param[in] obj - the desired objective function
   * \param[in] nozzle_shape - a class defining the nozzle shape
   * \param[out] dJdX - the gradient of the objective
   */
  void CalcGradient(const objective & obj, Nozzle & nozzle_shape,
                    InnerProdVector & dJdX);

  /*!
   * \brief calculates the gradient of the sensor objective w.r.t. the source
   * \param[in] x_center - x coordinate where the sensor is centered
   * \param[in] sigma - standard deviation of the sensor kernel
   * \param[in] press_targ - the constant target pressure
   * \param[in] reg_param - regularization parameter
   * \param[in] max_krylov - maximum number of Krylov iterations allowed per step
   * \param[in] flow_file - file name where flow solution is stored
   * \param[out] dJdX - gradient of the sensor objective
   */
  void CalcSensorGradient(const double & x_center, const double & sigma,
                          const double & press_targ, const double & reg_param,
                          const int & max_krylov, const string & flow_file,
                          InnerProdVector & dJdX);

  /*!
   * \brief adds the partial deriv of the sensor w.r.t. the source to given vec
   * \param[in] reg_param - regularization parameter
   * \param[in] dJdX - vector that partial derivative is added to
   */
  void AddSensordJdSource(const double & reg_param, InnerProdVector & dJdX);
  
  /*!
   * \brief product of the Unsteady Jacobian, dRes/dx, with design variables
   */
  void UnsteadyJacobianSourceProduct(const InnerProdVector & u,
                                     const string & prod_file);
  
  /*!
   * \brief product of the transposed Unsteady Jacobian, dRes/dx, with adjoint
   * \param[in] psi_file - file name where adjoint solution is stored
   * \param[in] dJdX - gradient of the sensor objective
   */
  void UnsteadyJacTransSourceProduct(const string & psi_file,
                                     InnerProdVector & dJdX);

  /*!
   * \brief tests CalcSensorGradient() using a finite-difference approximation
   */
  void TestSensorGradient();

  /*!
   * \brief calculates an estimate of the error in the gradient norm
   * \param[in] norm - type of norm used on the gradient (L2, inf)
   * \param[in] obj - objective function corresponding to the gradient
   * \param[in] nozzle_shape - a class defining the nozzle shape
   * \param[in] dJdX - the gradient of the objective
   * \returns the estimate of the error using adjoint-weighted residual
   */
  double EstimateGradientError(const norm_type & norm,
                               const objective & obj,
                               Nozzle & nozzle_shape,
                               const InnerProdVector & dJdX);

  /*!
   * \brief calculates an estimate of the error in the gradient norm
   * \param[in] norm - type of norm used on the gradient (L2, inf)
   * \param[in] obj - objective function corresponding to the gradient
   * \param[in] nozzle_shape - a class defining the nozzle shape
   * \param[in] dJdX - the gradient of the objective
   * \returns the estimate of the error using adjoint-weighted residual
   *
   * This version uses finite-difference approximations to compute the
   * right-hand-side of one of the adjoint problems
   */
  double EstimateGradientErrorFD(const norm_type & norm,
                                 const objective & obj,
                                 Nozzle & nozzle_shape,
                                 const InnerProdVector & dJdX);

  /*!
   * \brief evaluates the adjoint residual for objective obj
   * \param[in] obj - dual functional of interest
   * \param[out] adj_res - adjoint residual evaluated at q_ and psi_
   */
  void EvaluateAdjointResidual(const objective & obj,
                               InnerProdVector & adj_res);

 private:

  /*!
   * \brief returns a reference to a flow variable at a particular node
   * \param[in] node - node at which variable is wanted
   * \param[in] var - index of desired variable
   */
  double & q(const int & node, const int & var) {
    return q_(3*node+var);
  }

  /*!
   * \brief returns a const reference to a flow variable at a node
   * \param[in] node - node at which variable is wanted
   * \param[in] var - index of desired variable
   */
  const double & q(const int & node, const int & var) const {
    return q_(3*node+var);
  }

  /*!
   * \brief returns a reference to all variables at a particular node
   * \param[in] node - node at which variable is wanted
   */
  ublas::vector_range<ublas::vector<double> > q(const int & node) {
    ublas::vector_range<ublas::vector<double> > 
        q_at_node(q_, ublas::range(3*node, 3*node + 3));
    return q_at_node;
  }

  /*!
   * \brief returns a reference to a old flow variable at a particular node
   * \param[in] node - node at which variable is wanted
   * \param[in] var - index of desired variable
   */
  double & q_old(const int & node, const int & var) {
    return q_old_(3*node+var);
  }

  /*!
   * \brief returns a const reference to an old flow variable at a node
   * \param[in] node - node at which variable is wanted
   * \param[in] var - index of desired variable
   */
  const double & q_old(const int & node, const int & var) const {
    return q_old_(3*node+var);
  }

  /*!
   * \brief returns a reference to all old variables at a particular node
   * \param[in] node - node at which variable is wanted
   */
  ublas::vector_range<ublas::vector<double> > q_old(const int & node) {
    ublas::vector_range<ublas::vector<double> > 
        q_at_node(q_old_, ublas::range(3*node, 3*node + 3));
    return q_at_node;
  }

  /*!
   * \brief returns a reference to a residual at a particular node
   * \param[in] node - node at which residual is wanted
   * \param[in] var - index of desired residual
   */
  double & res(const int & node, const int & var) {
    return res_(3*node+var); 
  }

  /*!
   * \brief returns a reference to all residuals at a particular node
   * \param[in] node - node at which residual is wanted
   */
  ublas::vector_range<ublas::vector<double> > res(const int & node) {
    ublas::vector_range<ublas::vector<double> > 
        res_at_node(res_, ublas::range(3*node, 3*node + 3));
    return res_at_node;
  }

  /*!
   * \brief calculates and stores the Euler flux based on value in q_
   * \param[in] q_var - the vector of the flow state used to evaluate flux
   * \param[out] flux - vector of fluxes at each node
   * 
   * Note: press_ is used inside this member function, so it must be consistent
   * with q_var
   */
  void CalcEulerFlux(const InnerProdVector & q_var, InnerProdVector & flux);

  /*!
   * \brief calculates the pressure, sound speed, and spectral radius
   * \param[in] q_var - the state used to compute the auxilliary variables
   */
  void CalcAuxiliaryVariables(const InnerProdVector & q_var);

  /*!
   * \brief calculates the flux Jacobian at a given state
   * \param[in] area - nozzle area at particular location
   * \param[in] q - flow state
   * \param[out] flux_jac - the computed flux Jacobian
   */
  void CalcFluxJacobian(
      const double & area,
      const ublas::vector_range<ublas::vector<double> > q,
      ublas::matrix<double> & flux_jac) const;

  /*!
   * \brief calculates the product of the flux Hessian and a vector
   * \param[in] area - nozzle area at particular location
   * \param[in] q - flow state
   * \param[in] w - a vector to multiply the flux Hessian
   * \param[out] flux_hess - the computed flux Hessian product
   */
  void CalcFluxHessianProduct(
      const double & area,
      const ublas::vector_range<ublas::vector<double> > q,
      const ublas::vector_range<ublas::vector<double> > w,
      ublas::matrix<double> & flux_hess) const;

  /*!
   * \brief calculates the product of the flux Hessian and a vector
   * \param[in] area - nozzle area at particular location
   * \param[in] q - flow state
   * \param[in] w - a vector to multiply the flux Hessian
   * \param[out] flux_hess - the computed flux Hessian product
   *
   * This version uses hyperdual numbers for testing
   */
  void CalcFluxHessianProductHD(
      const double & area,
      const ublas::vector_range<ublas::vector<double> > q,
      const ublas::vector_range<ublas::vector<double> > w,
      ublas::matrix<double> & flux_hess) const;

  /*!
   * \breif templated routine to calculate the Euler flux (for debugging)
   * \param[in] area - nozzle area at particular location
   * \param[in] q - flow state
   * \param[out] flux - inviscid flux vector
   */
  template <class T>
  void CalcFlux(const T & area, const ublas::bounded_vector<T,3> & q,
                ublas::bounded_vector<T,3> & flux) const {
    T rho = q(0);
    T u = q(1)/rho;
    T e = q(2);
    T press = static_cast<T>(kGamma-1.0)*(e - 0.5*rho*u*u);
    flux(0) = rho * u * area;
    flux(1) = (rho * u * u + press) * area;
    flux(2) = u * ( e + press) * area;
  }

  /*!
   * \brief calculates a simultaneous approximation term (SAT)
   * \param[in] bc - boundary condition state
   * \param[in] q - flow state
   * \param[out] sat - the SAT penalty term
   */
  template <class T>
  void CalcSAT(const ublas::bounded_vector<T,3> & bc,
               const T & area, const T & sgn,
               const ublas::bounded_vector<T,3> & q,
               ublas::bounded_vector<T,3> & sat) const {
    // calculate primative variables
    T rho = q(0);
    T u = q(1)/rho;
    T E = q(2);
  
    const T Gamma = static_cast<T>(kGamma);
    T p = (Gamma - 1.0)*(E - 0.5*rho*u*u);
    T a = sqrt(Gamma*p/rho);
    T H = (E + p)/rho;
    T phi = 0.5*u*u;
  
    // calculate the wave speeds
    T lam1 = 0.5*area*(fabs(u + a) + sgn*(u + a));
    T lam2 = 0.5*area*(fabs(u - a) + sgn*(u - a));
    T lam3 = 0.5*area*(fabs(u) + sgn*(u));

#if 0
    T spec = fabs(u) + a;
    T lam1 = 0.5*area*(max(fabs(u + a),kVn*spec) + sgn*(u + a));
    T lam2 = 0.5*area*(max(fabs(u - a),kVn*spec) + sgn*(u - a));
    T lam3 = 0.5*area*(max(fabs(u),kVl*spec) + sgn*(u));
#endif 

    // calculate the differences
    T dq1 = q(0) - bc(0);
    T dq2 = q(1) - bc(1);
    T dq3 = q(2) - bc(2);

    sat(0) = lam3*dq1;
    sat(1) = lam3*dq2;
    sat(2) = lam3*dq3;

    // temporary vectors
    ublas::bounded_vector<T, 3> E1dq, E2dq;
  
    // get E1 times dq
    E1dq(0) = phi*dq1 - u*dq2 + dq3;
    E1dq(1) = E1dq(0)*u;
    E1dq(2) = E1dq(0)*H;

    // get E2 times dq
    E2dq(0) = 0.0;
    E2dq(1) = -u*dq1 + dq2;
    E2dq(2) = E2dq(1)*u;

    // add to pen
    T tmp1 = 0.5*(lam1 + lam2) - lam3;
    T tmp2 = (Gamma-1.0)/(a*a);
    sat += tmp1*(tmp2*E1dq + E2dq);

    // get E3 times dq
    E1dq(0) = -u*dq1 + dq2;
    E1dq(1) = E1dq(0)*u;
    E1dq(2) = E1dq(0)*H;
  
    // get E4 times dq
    E2dq(0) = 0.0;
    E2dq(1) = phi*dq1 - u*dq2 + dq3;
    E2dq(2) = E2dq(1)*u;

    // add to sat
    tmp1 = 0.5*(lam1 - lam2)/a; 
    sat += tmp1*(E1dq + (Gamma-1.0)*E2dq);
  }

  template <class T>
  void CalcSAT(const ublas::bounded_vector<T,3> & bc, 
               const T & area, const T & sgn,
               ublas::vector_range<ublas::vector<T> > q,
               ublas::bounded_vector<T,3> & sat) const {
    ublas::bounded_vector<T,3> q_temp;
    q_temp(0) = q(0); q_temp(1) = q(1); q_temp(2) = q(2);
    CalcSAT(bc, area, sgn, q_temp, sat);
  }

#if 0  
  template <class T>
  void CalcSAT(ublas::bounded_vector<T,3> & bc,
               const T & area, const T & sgn,
               ublas::bounded_vector<T,3> & q,
               ublas::bounded_vector<T,3> & sat) const {
    ublas::vector_range<ublas::vector<T> >
        q_range(q, ublas::range(0, 3)),
        bc_range(bc, ublas::range(0, 3));
    CalcSAT(bc_range, area, sgn, q_range, sat);
  }
#endif

#if 0
  /*!
   * \brief complex version of CalcSAT
   * \param[in] bc - boundary condition state
   * \param[in] q - flow state
   * \param[out] sat - the SAT penalty term
   */
  void CalcSAT_Complex(const ublas::bounded_vector<complex,3> & bc,
                       const complex & area, const complex & sgn,
                       const ublas::bounded_vector<complex,3> & q,
                       ublas::bounded_vector<complex,3> & sat) const;
#endif

  /*!
   * \brief returns the L2 norm of the residual, based on the SBP norm
   * \result the L2 norm of res_ calculated using the SBP quadrature
   */
  double ResidualNorm() {
    InnerProdVector prod(num_nodes_, 0.0);
    for (int i = 0; i < num_nodes_; i++) {
      for (int j = 0; j < 3; j++)
        prod(i) += (res(i,j)*res(i,j));
    }
    return sqrt(sbp_deriv_.InnerProductSBP(1, met_jac_, prod));
  }

  /*!
   * \brief returns a global time step based on a target CFL number
   * \param[in] CFL - target CFL number
   * \result a time step that keeps the CFL below the target CFL
   */
  double CalcDtFromCFL(const double & CFL) {
    double dt = 1E+15;
    //double dx = 1.0/(static_cast<double>(num_nodes_-1));
    for (int i = 0; i < num_nodes_; i++) {
      double dx = met_jac_(i); 
      if (dt > CFL*dx/spect_(i)) 
        dt = CFL*dx/spect_(i);
    }
  }

  int num_nodes_; ///< number of nodes
  int num_primal_precond_; ///< num of calls to primal preconditioner
  int num_adjoint_precond_; ///< num of calls to adjoint preconditioner
  double diss_coeff_; ///< numerical dissipation coefficient 
  double dxi_; ///< mesh spacing for domain = [0,1] (or in computational space)
  double dt_; ///< time step size
  ublas::bounded_vector<double,3> bc_left_; ///< left-end boundary conditions
  ublas::bounded_vector<double,3> bc_right_; ///< right-end boundary conditions
  InnerProdVector x_coord_; ///< nodal coordinate locations
  InnerProdVector met_jac_; ///< determinant of the metric Jacobian
  InnerProdVector area_; ///< nozzle area at each node
  InnerProdVector q_; ///< solution vector of conservative variables
  InnerProdVector q_old_; ///< previous solution vector of conservative variables
  InnerProdVector press_; ///< pressure at each node
  InnerProdVector press_targ_; ///< target pressure at nodes, for inverse
  InnerProdVector sndsp_; ///< sound speed at each node
  InnerProdVector spect_; ///< flux Jacobian spectral radius
  InnerProdVector res_; ///< residual vector
  InnerProdVector psi_; ///< adjoint vector
  InnerProdVector src_; ///< momentum source vector (between time steps)
  double src_x_; ///< origin of momentum source
  double src_sig2_; ///< square of standard deviation of momentum source
  int src_indx_; ///< node index where momentum source is applied
  SBP1stDerivative sbp_deriv_; ///< operator for first derivative
  SBPDissipation sbp_diss_; ///< numerical dissipation operator
  ublas::banded_matrix<double> prec_; ///< holds the preconditioner
};
