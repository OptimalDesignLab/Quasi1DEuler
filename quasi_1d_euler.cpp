/**
 * \file quasi_1d_euler.cpp
 * \brief function defintions for Quasi1DEuler member functions
 * \author  Jason Hicken <jason.hicken@gmail.com>
 * \version 1.0
 */

#include "./quasi_1d_euler.hpp"

#include <math.h>

#include <ostream>
#include <iostream>
#include <fstream>
#include <string>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "../krylov.hpp"

#include "./exact_solution.hpp"
#include "./nozzle.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ofstream;
namespace ublas = boost::numeric::ublas;

// ======================================================================

complex fabs(const complex & z) {
  if (z.real() < 0.0) {
    return -z;
  } else {
    return z;
  }
}

// ======================================================================

void Quasi1DEuler::ResizeGrid(const InnerProdVector & coord) {
  num_nodes_ = coord.size();
  if (num_nodes_ == 0) {
    cerr << "Quasi1DEuler::ResizeGrid(): "
         << "coord is empty";
    throw(-1);
  }
  int order = sbp_deriv_.order();
  sbp_deriv_.Define(num_nodes_, order);
  sbp_diss_.Define(num_nodes_, order, order);
  prec_.clear();
  prec_.resize(3*num_nodes_, 3*num_nodes_, 5, 5, false);
  x_coord_.resize(num_nodes_);
  x_coord_ = coord;
  met_jac_.resize(num_nodes_);
  sbp_deriv_.Apply(1, x_coord_, met_jac_);
  area_.resize(num_nodes_);
  press_.resize(num_nodes_);
  press_targ_.resize(num_nodes_);
  sndsp_.resize(num_nodes_);
  spect_.resize(num_nodes_);
  q_.resize(3*num_nodes_);
  q_old_.resize(3*num_nodes_);
  res_.resize(3*num_nodes_);
  psi_.resize(3*num_nodes_);
}

// ======================================================================

void Quasi1DEuler::CalcResidual() {
  res_ = 0.0;

  // evaluate the pressure, soundspeed, and spectral radius
  CalcAuxiliaryVariables(q_);

  // evaluate the flux at each node and apply SBP first-derivative
  InnerProdVector flux(3*num_nodes_, 0.0);
  CalcEulerFlux(q_, flux);
  sbp_deriv_.Apply(3, flux, res_);

  // add the source term
  InnerProdVector work(num_nodes_, 0.0);
  sbp_deriv_.Apply(1, area_, work);
  for (int i = 0; i < num_nodes_; i++)
    res(i,1) -= work(i)*press_(i);

  // add numerical dissipation
  for (int i = 0; i < num_nodes_; i++)
    work(i) = diss_coeff_;
  sbp_diss_.Apply(3, q_, work, flux);
  res_ += flux;

  // add the SAT boundary penalty terms
  double dx = 1.0/static_cast<double>(num_nodes_-1);
  ublas::bounded_vector<double, 3> sat;
  CalcSAT(bc_left_, area_(0), 1.0, q(0), sat);
  res(0) += sbp_deriv_.Hinv()*sat;
  int nm1 = num_nodes_-1;
  CalcSAT(bc_right_, area_(nm1), -1.0, q(nm1), sat);
  res(nm1) += sbp_deriv_.Hinv()*sat;
}

// ======================================================================

void Quasi1DEuler::CalcUnsteadyResidual() {
  res_ = 0.0;

  // evaluate the pressure, soundspeed, and spectral radius
  InnerProdVector q_mid(3*num_nodes_,0.0);
  q_mid.EqualsAXPlusBY(0.5, q_, 0.5, q_old_);
  CalcAuxiliaryVariables(q_mid);

  // evaluate the flux at each node and apply SBP first-derivative
  InnerProdVector flux(3*num_nodes_, 0.0);
  CalcEulerFlux(q_mid, flux);
  sbp_deriv_.Apply(3, flux, res_);

  // add numerical dissipation
  InnerProdVector work(num_nodes_, 0.0);
  for (int i = 0; i < num_nodes_; i++)
    work(i) = diss_coeff_;
  sbp_diss_.Apply(3, q_mid, work, flux);
  res_ += flux;

  // add the SAT boundary penalty terms
  ublas::bounded_vector<double, 3> sat;
  ublas::vector_range<ublas::vector<double> >
      q_mid_at_node(q_mid, ublas::range(3*0, 3*0 + 3));
  CalcSAT(bc_left_, area_(0), 1.0, q_mid_at_node, sat);
  res(0) += sbp_deriv_.Hinv()*sat;
  int nm1 = num_nodes_-1;
  q_mid_at_node = ublas::vector_range<ublas::vector<double> >
      (q_mid, ublas::range(3*nm1, 3*nm1 + 3));
  CalcSAT(bc_right_, area_(nm1), -1.0, q_mid_at_node, sat);
  res(nm1) += sbp_deriv_.Hinv()*sat;

  // finally, scale residual by dt/dxi and add solution difference
  double dx = 1.0/static_cast<double>(num_nodes_-1);
  res_ *= dt()/dxi();
  for (int i = 0; i < num_nodes_; i++)
    res(i) += q(i) - q_old(i);
}

// ======================================================================

void Quasi1DEuler::AddUnsteadySource(const int & n) {
  double src_max = dt()*src_(n);
  for (int i = 0; i < num_nodes_; i++) {
    double dx = x_coord_(i) - src_x_;
    res_(3*i+1) -= src_max*exp(-dx*dx/src_sig2_);
  }
}

// ======================================================================

void Quasi1DEuler::CalcDPressDQProduct(const InnerProdVector & u,
                                       InnerProdVector & v) {
  // check for consistent sizes
  if ( (u.size() != 3*num_nodes_) || (v.size() != num_nodes_) ) {
    cerr << "Quasi1DEuler::CalcDPressDQProduct(): "
         << "inconsistent sizes.";
    throw(-1);
  }
  for (int i = 0; i < num_nodes_; i++) {
    double vel = q(i,1)/q(i,0);
    v(i) = (kGamma-1.0)*(0.5*vel*vel*u(3*i) - vel*u(3*i+1) + u(3*i+2));
  }
}

// ======================================================================

void Quasi1DEuler::CalcDPressDQTransposedProduct(
    const InnerProdVector & u, InnerProdVector & v) {
  // check for consistent sizes
  if ( (u.size() != num_nodes_) || (v.size() != 3*num_nodes_) ) {
    cerr << "Quasi1DEuler::CalcDPressDQTransposedProduct(): "
         << "inconsistent sizes.";
    throw(-1);
  }
  for (int i = 0; i < num_nodes_; i++) {
    double vel = q(i,1)/q(i,0);
    v(3*i) = (kGamma-1.0)*0.5*vel*vel*u(i);
    v(3*i+1) = -(kGamma-1.0)*vel*u(i);
    v(3*i+2) = (kGamma-1.0)*u(i);
  }
}

// ======================================================================

void Quasi1DEuler::TestDPressDQProducts() {
  // create random vectors to apply dPress/dQ to from either side
  InnerProdVector u(3*num_nodes_, 0.0), v(num_nodes_, 0.0),
      w(3*num_nodes_, 0.0), z(num_nodes_, 0.0);
  boost::random::mt19937 gen;
  boost::random::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (int i = 0; i < 3*num_nodes_; i++)
    u(i) = dist(gen);
  for (int i = 0; i < num_nodes_; i++)
    v(i) = dist(gen);
  // evaluate dPress/dQ-vector product and contract with v
  CalcDPressDQProduct(u, z);
  double forward = InnerProd(z, v);
  // evaluate the transposed-dPress/dQ-vector product and contract with u
  CalcDPressDQTransposedProduct(v, w);
  double backward = InnerProd(w, u);
  cout << "Quasi1DEuler::TestDPressDQProducts:" << endl;
  cout << "\tDifference between forward and backward (relative error) = "
       << fabs((backward - forward)/forward) << endl;
}

// ======================================================================

void Quasi1DEuler::JacobianStateProduct(const InnerProdVector & u,
                                        InnerProdVector & v) {
  // check for consistent sizes
  if ( (u.size() != 3*num_nodes_) || (v.size() != 3*num_nodes_) ) {
    cerr << "Quasi1DEuler::JacobianStateProduct(): "
         << "inconsistent sizes.";
    throw(-1);
  }
  v = 0.0;

  // compute the product of the flux Jacobian matrices with u
  InnerProdVector Au(3*num_nodes_, 0.0);
  ublas::matrix<double> flux_jac(3, 3, 0.0);
  for (int i = 0; i < num_nodes_; i++) {
    CalcFluxJacobian(area_(i), q(i), flux_jac);
    ublas::range irange(3*i, 3*(i+1));
    Au(irange) = ublas::prod(flux_jac, u(irange));
  }
  // apply the SBP first derivative to the vector diag(A)*u
  sbp_deriv_.Apply(3, Au, v);

  // add terms corresponding to source term
  InnerProdVector work(num_nodes_, 0.0);
  sbp_deriv_.Apply(1, area_, work);
  for (int i = 0; i < num_nodes_; i++) {
    int indx = 3*i+1;
    double vel = q(i,1)/q(i,0);
    v(indx) -= work(i)*(kGamma-1.0)*(
        0.5*vel*vel*u(3*i) - vel*u(3*i+1) + u(3*i+2));
  }

  // add terms corresponding to numerical dissipation
  for (int i = 0; i < num_nodes_; i++)
    work(i) = diss_coeff_;
  sbp_diss_.Apply(3, u, work, Au);
  v += Au;

  // add terms corresponding to the SAT boundary penalties
  // Here we use the complex-step method
  double ceps = 1.E-30;
  ublas::bounded_vector<complex, 3> sat_c, bc_c, q_c;
  complex area_c, sgn_c;
  for (int i = 0; i < 3; i++) {
    q_c(i) = complex(q(0,i), ceps*u(i));
    bc_c(i) = complex(bc_left_(i), 0.0);
  }
  area_c = complex(area_(0), 0.0);
  sgn_c = complex(1.0, 0.0);
  CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
  for (int i = 0; i < 3; i++)
    v(i) += sbp_deriv_.Hinv()*sat_c(i).imag()/ceps;

  int nm1 = num_nodes_-1;
  for (int i = 0; i < 3; i++) {
    q_c(i) = complex(q(nm1,i), ceps*u(3*nm1+i));
    bc_c(i) = complex(bc_right_(i), 0.0);
  }
  area_c = complex(area_(nm1), 0.0);
  sgn_c = complex(-1.0, 0.0);
  CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
  for (int i = 0; i < 3; i++)
    v(3*nm1+i) += sbp_deriv_.Hinv()*sat_c(i).imag()/ceps;
}

// ======================================================================

void Quasi1DEuler::TestJacobianStateProduct() {
  // create a random vector to apply Jacobian to
  InnerProdVector u(3*num_nodes_, 0.0), v(3*num_nodes_, 0.0),
      v_fd(3*num_nodes_, 0.0), q_save(3*num_nodes_, 0.0);
  boost::random::mt19937 gen;
  boost::random::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (int i = 0; i < 3*num_nodes_; i++)
    u(i) = dist(gen);

  // evaluate Jacobian-vector product analytically
  JacobianStateProduct(u, v);

#if 0
  // uncomment to test the JacobianVectorProduct
  MatrixVectorProduct<InnerProdVector>*
      mat_vec = new JacobianVectorProduct(this);
  InnerProdVector u_tmp(u), v_tmp(v);
  (*mat_vec)(u_tmp, v_tmp);
  delete mat_vec;
  v = v_tmp;
#endif

  // evaluate the Jacobian-vector product using backward difference
  q_save = q_;  // save flow state for later

  // evaluate residual and save
  CalcResidual();
  v_fd = res_;

  // perturb flow and re-evaluate residual
  double fd_eps = 1.E-7;
  q_ -= fd_eps*u;
  CalcResidual();
  v_fd -= res_;
  v_fd /= fd_eps;

  // take difference between two products and store in q_ for output
  q_.EqualsAXPlusBY(1.0, v, -1.0, v_fd);
  double L2_error = sbp_deriv_.NormSBP(3, q_);
  cout << "Quasi1DEuler::TestJacobianStateProduct(): "
       << "L2 error between analytical and FD Jacobian-vector product: "
       << L2_error << endl;
}

// ======================================================================

void Quasi1DEuler::JacobianTransposedStateProduct(
    const InnerProdVector & u, InnerProdVector & v) {
  // check for consistent sizes
  if ( (u.size() != 3*num_nodes_) || (v.size() != 3*num_nodes_) ) {
    cerr << "Quasi1DEuler::TransposedJacobianStateProduct(): "
         << "inconsistent sizes.";
    throw(-1);
  }
  v = 0.0;

  // apply the first derivative to u
  InnerProdVector dudx(3*num_nodes_, 0.0);
  sbp_deriv_.ApplyTranspose(3, u, dudx);

  // left multiply dudx by the block diagonal matrix diag(A^{T})
  ublas::matrix<double> flux_jac(3, 3, 0.0);
  for (int i = 0; i < num_nodes_; i++) {
    CalcFluxJacobian(area_(i), q(i), flux_jac);
    ublas::range irange(3*i, 3*(i+1));
    v(irange) += ublas::prod(dudx(irange), flux_jac);
  }

  // add terms corresponding to source term
  InnerProdVector work(num_nodes_, 0.0);
  sbp_deriv_.Apply(1, area_, work);
  work *= (kGamma-1.0);
  for (int i = 0; i < num_nodes_; i++) {
    int indx = 3*i;
    double vel = q(i,1)/q(i,0);
    v(indx)   -= work(i)*0.5*vel*vel*u(indx+1);
    v(indx+1) += work(i)*vel*u(indx+1);
    v(indx+2) -= work(i)*u(indx+1);
  }

  // add terms corresponding to numerical dissipation
  for (int i = 0; i < num_nodes_; i++)
    work(i) = diss_coeff_;
  sbp_diss_.ApplyTranspose(3, u, work, dudx);
  v += dudx;

  // fill in elements corresponding to SAT penalties
  // Here we use the complex-step method
  double ceps = 1.E-30;
  ublas::bounded_vector<complex, 3> sat_c, bc_c, q_c;
  complex area_c, sgn_c;

  // initialize complex values for left end of domain
  for (int i = 0; i < 3; i++) {
    q_c(i) = complex(q(0,i), 0.0);
    bc_c(i) = complex(bc_left_(i), 0.0);
  }
  area_c = complex(area_(0), 0.0);
  sgn_c = complex(1.0, 0.0);

  // loop over variables that we differentiate w.r.t
  for (int i = 0; i < 3; i++) {
    q_c(i) += complex(0.0, ceps); // perturb ith variable
    CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
    for (int j = 0; j < 3; j++)
      v(i) += sbp_deriv_.Hinv()*sat_c(j).imag()*u(j)/ceps;
    q_c(i) -= complex(0.0, ceps); // unperturb ith variable
  }

  // initialize complex values for right end of domain
  int nm1 = num_nodes_-1;
  for (int i = 0; i < 3; i++) {
    q_c(i) = complex(q(nm1,i), 0.0);
    bc_c(i) = complex(bc_right_(i), 0.0);
  }
  area_c = complex(area_(nm1), 0.0);
  sgn_c = complex(-1.0, 0.0);

  // loop over variables that we differentiate w.r.t
  for (int i = 0; i < 3; i++) {
    q_c(i) += complex(0.0, ceps); // perturb ith variable
    CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
    for (int j = 0; j < 3; j++)
      v(3*nm1+i) += sbp_deriv_.Hinv()*sat_c(j).imag()*u(3*nm1+j)/ceps;
    q_c(i) -= complex(0.0, ceps); // unperturb ith variable
  }
}

// ======================================================================

void Quasi1DEuler::TestJacobianTransposedStateProduct() {
  // create a random vector to apply transposed Jacobian to
  InnerProdVector u(3*num_nodes_, 0.0), v(3*num_nodes_, 0.0),
      w(3*num_nodes_, 0.0);
  boost::random::mt19937 gen;
  boost::random::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (int i = 0; i < 3*num_nodes_; i++) {
    u(i) = dist(gen);
    v(i) = dist(gen);
  }
#if 0
  int i = 0;
  ublas::range irange(3*i, 3*(i+1));
  u(irange) = ublas::zero_vector<double>(3);
  v(irange) = ublas::zero_vector<double>(3);
  i = num_nodes_-1;
  irange = ublas::range(3*i, 3*(i+1));
  u(irange) = ublas::zero_vector<double>(3);
  v(irange) = ublas::zero_vector<double>(3);
#endif
  // evaluate Jacobian-vector product and contract with v
  JacobianStateProduct(u, w);
  double forward = InnerProd(v, w); //sbp_deriv_.InnerProductSBP(3, v, w);
  // evaluate the transposed-Jacobian-vector product and contract with u
  JacobianTransposedStateProduct(v, w);
  double backward = InnerProd(u, w); //sbp_deriv_.InnerProductSBP(3, u, w);
  cout << "Difference between forward and backward = "
       << backward - forward << endl;
}

// ======================================================================

void Quasi1DEuler::ResidualHessianProduct(
    const InnerProdVector & psi, InnerProdVector & w,
    InnerProdVector & v) {
  // check for consistent sizes
  if ( (psi.size() != 3*num_nodes_) || (w.size() != 3*num_nodes_) ||
       (v.size() != 3*num_nodes_) ) {
    cerr << "Quasi1DEuler::ResidualHessianProduct(): "
         << "inconsistent sizes.";
    throw(-1);
  }
  v = 0.0;

  // apply the first derivative to psi
  InnerProdVector dpsidx(3*num_nodes_, 0.0);
  sbp_deriv_.ApplyTranspose(3, psi, dpsidx);

  // left multiply dpsidx by the block diagonal matrix diag(w H^{T})
  // where H is the flux Hessian
  ublas::matrix<double> flux_hess(3, 3, 0.0);
  for (int i = 0; i < num_nodes_; i++) {
    ublas::range irange(3*i, 3*(i+1));
    CalcFluxHessianProduct(area_(i), q(i), w(irange), flux_hess);
    //CalcFluxHessianProductHD(area_(i), q(i), w(irange), flux_hess);
    v(irange) += ublas::prod(dpsidx(irange), flux_hess);
  }

  // add terms corresponding to source term
  InnerProdVector work(num_nodes_, 0.0);
  sbp_deriv_.Apply(1, area_, work);
  work *= (kGamma-1.0);
  for (int i = 0; i < num_nodes_; i++) {
    int indx = 3*i;
    double rho = q(i,0);
    double vel = q(i,1)/rho;
    v(indx)   -= work(i)*psi(indx+1)
        *(vel*w(indx+1) - vel*vel*w(indx))/rho;
    v(indx+1) -= work(i)*psi(indx+1)
        *(vel*w(indx) - w(indx+1))/rho;
  }

  // numerical dissipation is a linear operator, so nothing to add for it

  // compute the hessian of the SAT terms using hyperdual numbers

  ublas::bounded_vector<HyperDual, 3> sat_hd, bc_hd, q_hd;
  HyperDual area_hd, sgn_hd;

  // initialize hyper-dual values for left end of domain
  for (int i = 0; i < 3; i++) {
    q_hd(i).setvalues(q(0,i), 0.0, 0.0, 0.0);
    bc_hd(i).setvalues(bc_left_(i), 0.0, 0.0, 0.0);
  }
  area_hd.setvalues(area_(0), 0.0, 0.0, 0.0);
  sgn_hd.setvalues(1.0, 0.0, 0.0, 0.0);

  // loop over the first variable that we differentiate w.r.t.
  for (int i = 0; i < 3; i++) {
    q_hd(i).ipart() = 1.0;
    // loop over the second variable that we differentiate w.r.t
    for (int j = 0; j < 3; j++) {
      q_hd(j).jpart() = 1.0;
      CalcSAT<HyperDual>(bc_hd, area_hd, sgn_hd, q_hd, sat_hd);
      for (int k = 0; k < 3; k++) {
        v(i) += sbp_deriv_.Hinv()*psi(k)*sat_hd(k).ijpart()*w(j);
      }
      q_hd(j).jpart() = 0.0;
    }
    q_hd(i).ipart() = 0.0;
  }

  // initialize hyper-dual values for right end of domain
  int nm1 = num_nodes_-1;
  for (int i = 0; i < 3; i++) {
    q_hd(i).setvalues(q(nm1,i), 0.0, 0.0, 0.0);
    bc_hd(i).setvalues(bc_right_(i), 0.0, 0.0, 0.0);
  }
  area_hd.setvalues(area_(nm1), 0.0, 0.0, 0.0);
  sgn_hd.setvalues(-1.0, 0.0, 0.0, 0.0);

  // loop over the first variable that we differentiate w.r.t.
  for (int i = 0; i < 3; i++) {
    q_hd(i).ipart() = 1.0;
    // loop over the second variable that we differentiate w.r.t
    for (int j = 0; j < 3; j++) {
      q_hd(j).jpart() = 1.0;
      CalcSAT<HyperDual>(bc_hd, area_hd, sgn_hd, q_hd, sat_hd);
      for (int k = 0; k < 3; k++) {
        v(3*nm1+i) += sbp_deriv_.Hinv()*psi(3*nm1+k)
            *sat_hd(k).ijpart()*w(3*nm1+j);
      }
      q_hd(j).jpart() = 0.0;
    }
    q_hd(i).ipart() = 0.0;
  }
}

// ======================================================================

void Quasi1DEuler::TestResidualHessianProduct() {
  // create a random vectors to apply Hessian to
  InnerProdVector psi(3*num_nodes_, 0.0), w(3*num_nodes_, 0.0),
      Hprod(3*num_nodes_, 0.0), Hprod_fd(3*num_nodes_, 0.0),
      psidRdQ(3*num_nodes_, 0.0), q_save(3*num_nodes_, 0.0);
  boost::random::mt19937 gen;
  boost::random::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (int i = 0; i < 3*num_nodes_; i++) {
    psi(i) = dist(gen);
    w(i) = dist(gen);
  }

  // evaluate Residual Hessian product analytically
  ResidualHessianProduct(psi, w, Hprod);

  // evaluate the Jacobian-vector product using backward difference
  q_save = q_;  // save flow state for later

  // evaluate Transposed-Jacobian-vector product and save
  JacobianTransposedStateProduct(psi, Hprod_fd);

  // perturb flow and re-evaluate Transposed-Jacobian-vector product
  double fd_eps = 1.E-7;
  q_ -= fd_eps*w;
  JacobianTransposedStateProduct(psi, psidRdQ);
  Hprod_fd -= psidRdQ;
  Hprod_fd /= fd_eps;

  // take difference between two products and store in q_ for output
  q_.EqualsAXPlusBY(1.0, Hprod, -1.0, Hprod_fd);
  double L2_error = sbp_deriv_.NormSBP(3, q_);
  cout << "Quasi1DEuler::TestResidualHessianProduct(): "
       << "L2 error between analytical and FD product: "
       << L2_error << endl;
  //q_ = q_save;
}

// ======================================================================

void Quasi1DEuler::UnsteadyJacobianStateProduct(const InnerProdVector & u,
                                                InnerProdVector & v,
                                                const bool & plus_time) {
  // check for consistent sizes
  if ( (u.size() != 3*num_nodes_) || (v.size() != 3*num_nodes_) ) {
    cerr << "Quasi1DEuler::UnsteadyJacobianStateProduct(): "
         << "inconsistent sizes.";
    throw(-1);
  }
  v = 0.0;

  InnerProdVector q_mid(3*num_nodes_, 0.0);
  q_mid.EqualsAXPlusBY(0.5, q_, 0.5, q_old_);

  // compute the product of the flux Jacobian matrices with u
  InnerProdVector Au(3*num_nodes_, 0.0);
  ublas::matrix<double> flux_jac(3, 3, 0.0);
  for (int i = 0; i < num_nodes_; i++) {
    ublas::vector_range<ublas::vector<double> >
        q_mid_at_node(q_mid, ublas::range(3*i, 3*i + 3));
    CalcFluxJacobian(area_(i), q_mid_at_node, flux_jac);
    ublas::range irange(3*i, 3*(i+1));
    Au(irange) = ublas::prod(flux_jac, u(irange));
  }
  // apply the SBP first derivative to the vector diag(A)*u
  sbp_deriv_.Apply(3, Au, v);

  // add terms corresponding to source term
  // do nothing since dA/dx = 0.0

  // add terms corresponding to numerical dissipation
  InnerProdVector work(num_nodes_, 0.0);
  for (int i = 0; i < num_nodes_; i++)
    work(i) = diss_coeff_;
  sbp_diss_.Apply(3, u, work, Au);
  v += Au;

  // add terms corresponding to the SAT boundary penalties
  // Here we use the complex-step method
  double ceps = 1.E-30;
  ublas::bounded_vector<complex, 3> sat_c, bc_c, q_c;
  complex area_c, sgn_c;
  for (int i = 0; i < 3; i++) {
    int ptr = i;
    q_c(i) = complex(q_mid(ptr), ceps*u(i));
    bc_c(i) = complex(bc_left_(i), 0.0);
  }
  area_c = complex(area_(0), 0.0);
  sgn_c = complex(1.0, 0.0);
  CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
  for (int i = 0; i < 3; i++)
    v(i) += sbp_deriv_.Hinv()*sat_c(i).imag()/ceps;

  int nm1 = num_nodes_-1;
  for (int i = 0; i < 3; i++) {
    int ptr = 3*nm1 + i;
    q_c(i) = complex(q_mid(ptr), ceps*u(3*nm1+i));
    bc_c(i) = complex(bc_right_(i), 0.0);
  }
  area_c = complex(area_(nm1), 0.0);
  sgn_c = complex(-1.0, 0.0);
  CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
  for (int i = 0; i < 3; i++)
    v(3*nm1+i) += sbp_deriv_.Hinv()*sat_c(i).imag()/ceps;

  // finally, multiply by 0.5*dt and add time term
  v *= (0.5*dt()/dxi());
  if (plus_time)
    v += u;
  else
    v -= u;
}

// ======================================================================

void Quasi1DEuler::UnsteadyApproxJacStateProduct(const InnerProdVector & u,
                                                 InnerProdVector & v,
                                                 const bool & plus_time) {
  // check for consistent sizes
  if ( (u.size() != 3*num_nodes_) || (v.size() != 3*num_nodes_) ) {
    cerr << "Quasi1DEuler::UnsteadyApproxJacStateProduct(): "
         << "inconsistent sizes.";
    throw(-1);
  }
  v = 0.0;

  // build SBP operators, which must be consistent with the operators used
  // in the preconditioner
  SBP1stDerivative deriv(num_nodes_, 2);
  SBPDissipation diss(num_nodes_, 2, 2);

  InnerProdVector q_mid(3*num_nodes_, 0.0);
  q_mid.EqualsAXPlusBY(0.5, q_, 0.5, q_old_);

  // compute the product of the flux Jacobian matrices with u
  InnerProdVector Au(3*num_nodes_, 0.0);
  ublas::matrix<double> flux_jac(3, 3, 0.0);
  for (int i = 0; i < num_nodes_; i++) {
    ublas::vector_range<ublas::vector<double> >
        q_mid_at_node(q_mid, ublas::range(3*i, 3*i + 3));
    CalcFluxJacobian(area_(i), q_mid_at_node, flux_jac);
    ublas::range irange(3*i, 3*(i+1));
    Au(irange) = ublas::prod(flux_jac, u(irange));
  }
  // apply the SBP first derivative to the vector diag(A)*u
  deriv.Apply(3, Au, v);

  // add terms corresponding to source term
  // do nothing since dA/dx = 0.0

  // add terms corresponding to numerical dissipation
  InnerProdVector work(num_nodes_, 0.0);
  for (int i = 0; i < num_nodes_; i++)
    work(i) = diss_coeff_;
  diss.Apply(3, u, work, Au);
  v += Au;

  // add terms corresponding to the SAT boundary penalties
  // Here we use the complex-step method
  double ceps = 1.E-30;
  ublas::bounded_vector<complex, 3> sat_c, bc_c, q_c;
  complex area_c, sgn_c;
  for (int i = 0; i < 3; i++) {
    int ptr = i;
    q_c(i) = complex(q_mid(ptr), ceps*u(i));
    bc_c(i) = complex(bc_left_(i), 0.0);
  }
  area_c = complex(area_(0), 0.0);
  sgn_c = complex(1.0, 0.0);
  CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
  for (int i = 0; i < 3; i++)
    v(i) += deriv.Hinv()*sat_c(i).imag()/ceps;

  int nm1 = num_nodes_-1;
  for (int i = 0; i < 3; i++) {
    int ptr = 3*nm1 + i;
    q_c(i) = complex(q_mid(ptr), ceps*u(3*nm1+i));
    bc_c(i) = complex(bc_right_(i), 0.0);
  }
  area_c = complex(area_(nm1), 0.0);
  sgn_c = complex(-1.0, 0.0);
  CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
  for (int i = 0; i < 3; i++)
    v(3*nm1+i) += deriv.Hinv()*sat_c(i).imag()/ceps;

  // finally, multiply by 0.5*dt and add time term
  v *= (0.5*dt()/dxi());
  if (plus_time)
    v += u;
  else
    v -= u;
}

// ======================================================================

void Quasi1DEuler::UnsteadyJacTransStateProduct(const InnerProdVector & u,
                                                InnerProdVector & v,
                                                const bool & plus_time) {
  // check for consistent sizes
  if ( (u.size() != 3*num_nodes_) || (v.size() != 3*num_nodes_) ) {
    cerr << "Quasi1DEuler::UnsteadyJacTransStateProduct(): "
         << "inconsistent sizes.";
    throw(-1);
  }
  v = 0.0;

  InnerProdVector q_mid(3*num_nodes_, 0.0);
  q_mid.EqualsAXPlusBY(0.5, q_, 0.5, q_old_);

  // apply the first derivative to u
  InnerProdVector dudx(3*num_nodes_, 0.0);
  sbp_deriv_.ApplyTranspose(3, u, dudx);

  // left multiply dudx by the block diagonal matrix diag(A^{T})
  ublas::matrix<double> flux_jac(3, 3, 0.0);
  for (int i = 0; i < num_nodes_; i++) {
    ublas::vector_range<ublas::vector<double> >
        q_mid_at_node(q_mid, ublas::range(3*i, 3*i + 3));
    CalcFluxJacobian(area_(i), q_mid_at_node, flux_jac);
    ublas::range irange(3*i, 3*(i+1));
    v(irange) += ublas::prod(dudx(irange), flux_jac);
  }

  // add terms corresponding to source term
  // do nothing since dA/dx = 0.0

  // add terms corresponding to numerical dissipation
  InnerProdVector work(num_nodes_, 0.0);
  for (int i = 0; i < num_nodes_; i++)
    work(i) = diss_coeff_;
  sbp_diss_.ApplyTranspose(3, u, work, dudx);
  v += dudx;

  // fill in elements corresponding to SAT penalties
  // Here we use the complex-step method
  double ceps = 1.E-30;
  ublas::bounded_vector<complex, 3> sat_c, bc_c, q_c;
  complex area_c, sgn_c;

  // initialize complex values for left end of domain
  for (int i = 0; i < 3; i++) {
    int ptr = i;
    q_c(i) = complex(q_mid(ptr), 0.0);
    bc_c(i) = complex(bc_left_(i), 0.0);
  }
  area_c = complex(area_(0), 0.0);
  sgn_c = complex(1.0, 0.0);

  // loop over variables that we differentiate w.r.t
  for (int i = 0; i < 3; i++) {
    q_c(i) += complex(0.0, ceps); // perturb ith variable
    CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
    for (int j = 0; j < 3; j++)
      v(i) += sbp_deriv_.Hinv()*sat_c(j).imag()*u(j)/ceps;
    q_c(i) -= complex(0.0, ceps); // unperturb ith variable
  }

  // initialize complex values for right end of domain
  int nm1 = num_nodes_-1;
  for (int i = 0; i < 3; i++) {
    int ptr = 3*nm1 + i;
    q_c(i) = complex(q_mid(ptr), 0.0);
    bc_c(i) = complex(bc_right_(i), 0.0);
  }
  area_c = complex(area_(nm1), 0.0);
  sgn_c = complex(-1.0, 0.0);

  // loop over variables that we differentiate w.r.t
  for (int i = 0; i < 3; i++) {
    q_c(i) += complex(0.0, ceps); // perturb ith variable
    CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
    for (int j = 0; j < 3; j++)
      v(3*nm1+i) += sbp_deriv_.Hinv()*sat_c(j).imag()*u(3*nm1+j)/ceps;
    q_c(i) -= complex(0.0, ceps); // unperturb ith variable
  }

  // finally, multiply by 0.5*dt and add/subtract time term
  v *= (0.5*dt()/dxi());
  if (plus_time)
    v += u;
  else
    v -= u;
}

// ======================================================================

void Quasi1DEuler::UnsteadyApproxJacTransStateProduct(const InnerProdVector & u,
                                                      InnerProdVector & v,
                                                      const bool & plus_time) {
    // check for consistent sizes
  if ( (u.size() != 3*num_nodes_) || (v.size() != 3*num_nodes_) ) {
    cerr << "Quasi1DEuler::UnsteadyApproxJacTransStateProduct(): "
         << "inconsistent sizes.";
    throw(-1);
  }
  v = 0.0;

  // build SBP operators, which must be consistent with the operators used
  // in the preconditioner
  SBP1stDerivative deriv(num_nodes_, 2);
  SBPDissipation diss(num_nodes_, 2, 2);

  InnerProdVector q_mid(3*num_nodes_, 0.0);
  q_mid.EqualsAXPlusBY(0.5, q_, 0.5, q_old_);

  // apply the first derivative to u
  InnerProdVector dudx(3*num_nodes_, 0.0);
  deriv.ApplyTranspose(3, u, dudx);

  // left multiply dudx by the block diagonal matrix diag(A^{T})
  ublas::matrix<double> flux_jac(3, 3, 0.0);
  for (int i = 0; i < num_nodes_; i++) {
    ublas::vector_range<ublas::vector<double> >
        q_mid_at_node(q_mid, ublas::range(3*i, 3*i + 3));
    CalcFluxJacobian(area_(i), q_mid_at_node, flux_jac);
    ublas::range irange(3*i, 3*(i+1));
    v(irange) += ublas::prod(dudx(irange), flux_jac);
  }

  // add terms corresponding to source term
  // do nothing since dA/dx = 0.0

  // add terms corresponding to numerical dissipation
  InnerProdVector work(num_nodes_, 0.0);
  for (int i = 0; i < num_nodes_; i++)
    work(i) = diss_coeff_;
  diss.ApplyTranspose(3, u, work, dudx);
  v += dudx;

  // fill in elements corresponding to SAT penalties
  // Here we use the complex-step method
  double ceps = 1.E-30;
  ublas::bounded_vector<complex, 3> sat_c, bc_c, q_c;
  complex area_c, sgn_c;

  // initialize complex values for left end of domain
  for (int i = 0; i < 3; i++) {
    int ptr = i;
    q_c(i) = complex(q_mid(ptr), 0.0);
    bc_c(i) = complex(bc_left_(i), 0.0);
  }
  area_c = complex(area_(0), 0.0);
  sgn_c = complex(1.0, 0.0);

  // loop over variables that we differentiate w.r.t
  for (int i = 0; i < 3; i++) {
    q_c(i) += complex(0.0, ceps); // perturb ith variable
    CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
    for (int j = 0; j < 3; j++)
      v(i) += deriv.Hinv()*sat_c(j).imag()*u(j)/ceps;
    q_c(i) -= complex(0.0, ceps); // unperturb ith variable
  }

  // initialize complex values for right end of domain
  int nm1 = num_nodes_-1;
  for (int i = 0; i < 3; i++) {
    int ptr = 3*nm1 + i;
    q_c(i) = complex(q_mid(ptr), 0.0);
    bc_c(i) = complex(bc_right_(i), 0.0);
  }
  area_c = complex(area_(nm1), 0.0);
  sgn_c = complex(-1.0, 0.0);

  // loop over variables that we differentiate w.r.t
  for (int i = 0; i < 3; i++) {
    q_c(i) += complex(0.0, ceps); // perturb ith variable
    CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
    for (int j = 0; j < 3; j++)
      v(3*nm1+i) += deriv.Hinv()*sat_c(j).imag()*u(3*nm1+j)/ceps;
    q_c(i) -= complex(0.0, ceps); // unperturb ith variable
  }

  // finally, multiply by 0.5*dt and add/subtract time term
  v *= (0.5*dt()/dxi());
  if (plus_time)
    v += u;
  else
    v -= u;
}

// ======================================================================

void Quasi1DEuler::TestUnsteadyJacTransStateProduct() {
  // create a random vector to apply transposed Jacobian to
  InnerProdVector u(3*num_nodes_, 0.0), v(3*num_nodes_, 0.0),
      w(3*num_nodes_, 0.0);
  boost::random::mt19937 gen;
  boost::random::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (int i = 0; i < 3*num_nodes_; i++) {
    u(i) = dist(gen);
    v(i) = dist(gen);
  }
#if 0
  int i = 0;
  ublas::range irange(3*i, 3*(i+1));
  u(irange) = ublas::zero_vector<double>(3);
  v(irange) = ublas::zero_vector<double>(3);
  i = num_nodes_-1;
  irange = ublas::range(3*i, 3*(i+1));
  u(irange) = ublas::zero_vector<double>(3);
  v(irange) = ublas::zero_vector<double>(3);
#endif
  // evaluate Jacobian-vector product and contract with v
  UnsteadyJacobianStateProduct(u, w);
  double forward = InnerProd(v, w); //sbp_deriv_.InnerProductSBP(3, v, w);
  // evaluate the transposed-Jacobian-vector product and contract with u
  UnsteadyJacTransStateProduct(v, w);
  double backward = InnerProd(u, w); //sbp_deriv_.InnerProductSBP(3, u, w);
  cout << "Difference between forward and backward = "
       << backward - forward << endl;
}

// ======================================================================

void Quasi1DEuler::BuildAndFactorPreconditioner() {
  // set all elements to zero (could also use ublas::zero_matrix)
  ublas::banded_matrix<double>::iterator1 itrow;
  ublas::banded_matrix<double>::iterator2 itcol;
  for (itrow = prec_.begin1(); itrow != prec_.end1(); ++itrow) {
    for (itcol = itrow.begin(); itcol != itrow.end(); ++itcol) {
      *itcol = 0.0;
    }
  }

  // fill in elements corresponding to fluxes
  int i;
  ublas::range i_range, nbr_range;
  ublas::matrix<double> flux_jac(3, 3, 0.0);
  for (i = 1; i < num_nodes_-1; i++) {
    CalcFluxJacobian(area_(i), q(i), flux_jac);
    flux_jac *= 0.5;
    i_range = ublas::range(3*i, 3*(i+1));
    // add to nbr on left
    nbr_range = ublas::range(3*(i-1), 3*i);
    prec_(nbr_range, i_range) += flux_jac;
    // subtract from nbr on right
    nbr_range = ublas::range(3*(i+1), 3*(i+2));
    prec_(nbr_range, i_range) -= flux_jac;
  }
  // corrections for left side
  i = 0;
  CalcFluxJacobian(area_(i), q(i), flux_jac);
  i_range = ublas::range(3*i, 3*(i+1));
  prec_(i_range, i_range) -= flux_jac;
  nbr_range = ublas::range(3*(i+1), 3*(i+2));
  prec_(nbr_range, i_range) -= 0.5*flux_jac;
  CalcFluxJacobian(area_(i+1), q(i+1), flux_jac);
  prec_(i_range, nbr_range) += 0.5*flux_jac;
  // corrections for right side
  i = num_nodes_-1;
  CalcFluxJacobian(area_(i), q(i), flux_jac);
  i_range = ublas::range(3*i, 3*(i+1));
  prec_(i_range, i_range) += flux_jac;
  nbr_range = ublas::range(3*(i-1), 3*i);
  prec_(nbr_range, i_range) += 0.5*flux_jac;
  CalcFluxJacobian(area_(i-1), q(i-1), flux_jac);
  prec_(i_range, nbr_range) -= 0.5*flux_jac;

  // fill in elements corresponding to source
  InnerProdVector work(num_nodes_, 0.0);
  sbp_deriv_.Apply(1, area_, work);
  work *= (kGamma-1.0);
  for (i = 0; i < num_nodes_; i++) {
    int indx = 3*i+1;
    double vel = q(i,1)/q(i,0);
    prec_(indx, indx-1) -= work(i)*0.5*vel*vel;
    prec_(indx, indx)   -= work(i)*vel;
    prec_(indx, indx+1) -= work(i);
  }

  // fill in elements corresponding to numerical dissipation
  ublas::matrix<double> diag3x3(3, 3, 0.0);
  for (int j = 0; j < 3; j++)
    diag3x3(j,j) = 0.5; //*diss_coeff_;  // TEMP
  for (i = 1; i < num_nodes_-1; i++) {
    i_range = ublas::range(3*i, 3*(i+1));
    prec_(i_range, i_range) += 2.0*diag3x3;
    nbr_range = ublas::range(3*(i-1), 3*i);
    prec_(i_range, nbr_range) -= diag3x3;
    nbr_range = ublas::range(3*(i+1), 3*(i+2));
    prec_(i_range, nbr_range) -= diag3x3;
  }
  i_range = ublas::range(0, 3);
  prec_(i_range, i_range) += 2.0*diag3x3;
  nbr_range = ublas::range(3, 6);
  prec_(i_range, nbr_range) -= 2.0*diag3x3;
  i = num_nodes_-1;
  i_range = ublas::range(3*i, 3*(i+1));
  prec_(i_range, i_range) += 2.0*diag3x3;
  nbr_range = ublas::range(3*(i-1), 3*i);
  prec_(i_range, nbr_range) -= 2.0*diag3x3;

  // fill in elements corresponding to SAT penalties
  // Here we use the complex-step method
  double ceps = 1.E-30;
  ublas::bounded_vector<complex, 3> sat_c, bc_c, q_c;
  complex area_c, sgn_c;

  // initialize complex values for left end of domain
  for (i = 0; i < 3; i++) {
    q_c(i) = complex(q(0,i), 0.0);
    bc_c(i) = complex(bc_left_(i), 0.0);
  }
  area_c = complex(area_(0), 0.0);
  sgn_c = complex(1.0, 0.0);

  // loop over variables that we differentiate w.r.t
  for (i = 0; i < 3; i++) {
    q_c(i) += complex(0.0, ceps); // perturb ith variable
    CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
    for (int j = 0; j < 3; j++)
      prec_(j, i) += 2.0*sat_c(j).imag()/ceps;
    q_c(i) -= complex(0.0, ceps); // unperturb ith variable
  }

  // initialize complex values for right end of domain
  int indx = num_nodes_-1;
  for (i = 0; i < 3; i++) {
    q_c(i) = complex(q(indx,i), 0.0);
    bc_c(i) = complex(bc_right_(i), 0.0);
  }
  area_c = complex(area_(indx), 0.0);
  sgn_c = complex(-1.0, 0.0);

  // loop over variables that we differentiate w.r.t
  for (i = 0; i < 3; i++) {
    q_c(i) += complex(0.0, ceps); // perturb ith variable
    CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
    for (int j = 0; j < 3; j++)
      prec_(3*indx + j, 3*indx + i) += 2.0*sat_c(j).imag()/ceps;
    q_c(i) -= complex(0.0, ceps); // unperturb ith variable
  }

#if 0
  // uncomment to display matrix entries
  for (itrow = prec_.begin1(); itrow != prec_.end1(); ++itrow) {
    for (itcol = itrow.begin(); itcol != itrow.end(); ++itcol) {
      cout << "itcol.index1() = " << itcol.index1() << " :";
      cout << "itcol.index2() = " << itcol.index2() << " :";
      cout << "*itcol = " << *itcol << endl;
    }
  }
#endif

  // perform LU-factorization (no pivoting, but dissipation should prevent
  // catastrophe...hopefully)
  int err = ublas::lu_factorize(prec_); //, permute);
  if( err != 0 ) {
    cout << "Quasi1DEuler::BuildAndFactorPreconditioner(): "
         << "lu_factorize() returned with err = " << err << endl;
  }
}

// ======================================================================

void Quasi1DEuler::BuildAndFactorUnsteadyPreconditioner(
    const bool & factor) {
  // set all elements to zero (could also use ublas::zero_matrix)
  ublas::banded_matrix<double>::iterator1 itrow;
  ublas::banded_matrix<double>::iterator2 itcol;
  for (itrow = prec_.begin1(); itrow != prec_.end1(); ++itrow) {
    for (itcol = itrow.begin(); itcol != itrow.end(); ++itcol) {
      *itcol = 0.0;
    }
  }

  InnerProdVector q_mid(3*num_nodes_, 0.0);
  q_mid.EqualsAXPlusBY(0.5, q_, 0.5, q_old_);

  // fill in elements corresponding to fluxes
  int i;
  ublas::range i_range, nbr_range;
  ublas::matrix<double> flux_jac(3, 3, 0.0), flux_jac_old(3, 3, 0.0);
  for (i = 1; i < num_nodes_-1; i++) {
    ublas::vector_range<ublas::vector<double> >
        q_mid_at_node(q_mid, ublas::range(3*i, 3*i + 3));
    CalcFluxJacobian(area_(i), q_mid_at_node, flux_jac);
    flux_jac *= 0.5;

    i_range = ublas::range(3*i, 3*(i+1));
    // add to nbr on left
    nbr_range = ublas::range(3*(i-1), 3*i);
    prec_(nbr_range, i_range) += flux_jac;
    // subtract from nbr on right
    nbr_range = ublas::range(3*(i+1), 3*(i+2));
    prec_(nbr_range, i_range) -= flux_jac;
  }
  // corrections for left side
  i = 0;
  ublas::vector_range<ublas::vector<double> >
      q_mid_at_node(q_mid, ublas::range(0, 3));
  CalcFluxJacobian(area_(i), q_mid_at_node, flux_jac);
  i_range = ublas::range(3*i, 3*(i+1));
  prec_(i_range, i_range) -= flux_jac;
  nbr_range = ublas::range(3*(i+1), 3*(i+2));
  prec_(nbr_range, i_range) -= 0.5*flux_jac;
  q_mid_at_node = ublas::vector_range<ublas::vector<double> >
      (q_mid, ublas::range(3, 6));
  CalcFluxJacobian(area_(i+1), q_mid_at_node, flux_jac);
  prec_(i_range, nbr_range) += 0.5*flux_jac;
  // corrections for right side
  i = num_nodes_-1;
  q_mid_at_node = ublas::vector_range<ublas::vector<double> >
      (q_mid, ublas::range(3*i, 3*i + 3));
  CalcFluxJacobian(area_(i), q_mid_at_node, flux_jac);
  i_range = ublas::range(3*i, 3*(i+1));
  prec_(i_range, i_range) += flux_jac;
  nbr_range = ublas::range(3*(i-1), 3*i);
  prec_(nbr_range, i_range) += 0.5*flux_jac;
  i -= 1;
  q_mid_at_node = ublas::vector_range<ublas::vector<double> >
      (q_mid, ublas::range(3*i, 3*i + 3));
  CalcFluxJacobian(area_(i-1), q_mid_at_node, flux_jac);
  prec_(i_range, nbr_range) -= 0.5*flux_jac;

  // fill in elements corresponding to source
  // do nothing because dA/dx = 0

  // fill in elements corresponding to numerical dissipation
  ublas::matrix<double> diag3x3(3, 3, 0.0);
  for (int j = 0; j < 3; j++)
    diag3x3(j,j) = 0.5;
  for (i = 1; i < num_nodes_-1; i++) {
    i_range = ublas::range(3*i, 3*(i+1));
    prec_(i_range, i_range) += 2.0*diag3x3;
    nbr_range = ublas::range(3*(i-1), 3*i);
    prec_(i_range, nbr_range) -= diag3x3;
    nbr_range = ublas::range(3*(i+1), 3*(i+2));
    prec_(i_range, nbr_range) -= diag3x3;
  }
  i_range = ublas::range(0, 3);
  prec_(i_range, i_range) += 2.0*diag3x3;
  nbr_range = ublas::range(3, 6);
  prec_(i_range, nbr_range) -= 2.0*diag3x3;
  i = num_nodes_-1;
  i_range = ublas::range(3*i, 3*(i+1));
  prec_(i_range, i_range) += 2.0*diag3x3;
  nbr_range = ublas::range(3*(i-1), 3*i);
  prec_(i_range, nbr_range) -= 2.0*diag3x3;

  // fill in elements corresponding to SAT penalties
  // Here we use the complex-step method
  double ceps = 1.E-30;
  ublas::bounded_vector<complex, 3> sat_c, bc_c, q_c;
  complex area_c, sgn_c;

  // initialize complex values for left end of domain
  for (i = 0; i < 3; i++) {
    int ptr = i;
    q_c(i) = complex(q_mid(ptr), 0.0);
    bc_c(i) = complex(bc_left_(i), 0.0);
  }
  area_c = complex(area_(0), 0.0);
  sgn_c = complex(1.0, 0.0);

  // loop over variables that we differentiate w.r.t
  for (i = 0; i < 3; i++) {
    q_c(i) += complex(0.0, ceps); // perturb ith variable
    CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
    for (int j = 0; j < 3; j++)
      prec_(j, i) += 2.0*sat_c(j).imag()/ceps;
    q_c(i) -= complex(0.0, ceps); // unperturb ith variable
  }

  // initialize complex values for right end of domain
  int indx = num_nodes_-1;
  for (i = 0; i < 3; i++) {
    int ptr = 3*indx + i;
    q_c(i) = complex(q_mid(ptr), 0.0);
    bc_c(i) = complex(bc_right_(i), 0.0);
  }
  area_c = complex(area_(indx), 0.0);
  sgn_c = complex(-1.0, 0.0);

  // loop over variables that we differentiate w.r.t
  for (i = 0; i < 3; i++) {
    q_c(i) += complex(0.0, ceps); // perturb ith variable
    CalcSAT<complex>(bc_c, area_c, sgn_c, q_c, sat_c);
    for (int j = 0; j < 3; j++)
      prec_(3*indx + j, 3*indx + i) += 2.0*sat_c(j).imag()/ceps;
    q_c(i) -= complex(0.0, ceps); // unperturb ith variable
  }

#if 0
  // uncomment to display matrix entries
  for (itrow = prec_.begin1(); itrow != prec_.end1(); ++itrow) {
    for (itcol = itrow.begin(); itcol != itrow.end(); ++itcol) {
      cout << "itcol.index1() = " << itcol.index1() << " :";
      cout << "itcol.index2() = " << itcol.index2() << " :";
      cout << "*itcol = " << *itcol << endl;
    }
  }
#endif

  // multiply by 0.5*dt/dxi and add the diagonal time term
  prec_ *= (0.5*dt()/dxi());
  for (int j = 0; j < 3; j++)
    diag3x3(j,j) = 1.0;
  for (i = 1; i < num_nodes_-1; i++) {
    i_range = ublas::range(3*i, 3*(i+1));
    prec_(i_range, i_range) += diag3x3;
  }

  if (factor) {
    // perform LU-factorization (no pivoting, but dissipation should prevent
    // catastrophe...hopefully)
    int err = ublas::lu_factorize(prec_); //, permute);
    if( err != 0 ) {
      cout << "Quasi1DEuler::BuildAndFactorPreconditioner(): "
           << "lu_factorize() returned with err = " << err << endl;
    }
  }
}

// ======================================================================

void Quasi1DEuler::Precondition(const InnerProdVector & u,
                                InnerProdVector & v) {
  v = u;
  ublas::lu_substitute(
      static_cast<ublas::banded_matrix<double> >(prec_),v);
  num_primal_precond_++;
}

// ======================================================================

void Quasi1DEuler::PreconditionTransposed(const InnerProdVector & u,
                                          InnerProdVector & v) {
  v = u;
  ublas::lu_substitute(v,
      static_cast<ublas::banded_matrix<double> >(prec_));
  num_adjoint_precond_++;
}

// ======================================================================

void Quasi1DEuler::PreconditionerMultiply(const InnerProdVector & u,
                                          InnerProdVector & v) {
  //v = ublas::prod(static_cast<ublas::banded_matrix<double> >(prec_), u);
  v = static_cast<InnerProdVector>(
      ublas::prod(prec_, static_cast<ublas::vector<double> >(u)));
#if 0
  ublas::banded_matrix<double>::iterator1 itrow;
  ublas::banded_matrix<double>::iterator2 itcol;
  double max_elem = 0.0;
  for (itrow = prec_.begin1(); itrow != prec_.end1(); ++itrow) {
    for (itcol = itrow.begin(); itcol != itrow.end(); ++itcol) {
      max_elem = std::max(fabs(*itcol), max_elem);
    }
  }
  cout << "Max prec_ element = " << max_elem << endl;
#endif
}

// ======================================================================

void Quasi1DEuler::JacobianAreaProduct(const InnerProdVector & u,
                                       InnerProdVector & v) {
  // check for consistent sizes
  if ( (u.size() != num_nodes_) || (v.size() != 3*num_nodes_) ) {
    cerr << "Quasi1DEuler::JacobianAreaProduct(): "
         << "inconsistent sizes.";
    throw(-1);
  }
  v = 0.0;

  // evaluate the pressure, soundspeed, and spectral radius
  CalcAuxiliaryVariables(q_);

  // evaluate terms that depend on Euler flux
  InnerProdVector flux(3*num_nodes_, 0.0);
  for (int i = 0; i < num_nodes_; i++) {
    double u_i = u(i);
    double press = press_(i);
    double rho = q(i,0);
    double vel = q(i,1) / rho;
    double e = q(i,2);
    int ptr = (i*3);
    flux(ptr+0) = rho * vel * u_i;
    flux(ptr+1) = (rho * vel * vel + press) * u_i;
    flux(ptr+2) = vel * ( e + press) * u_i;
  }
  sbp_deriv_.Apply(3, flux, v);

  // add terms that depend on the source
  InnerProdVector work(num_nodes_, 0.0);
  sbp_deriv_.Apply(1, u, work);
  for (int i = 0; i < num_nodes_; i++) {
    int ptr = (i*3) + 1;
    v(ptr) -= work(i)*press_(i);
  }

  // numerical dissipation does not depend on area (no spect_ term)

  // add terms that depend on SAT boundary penalties
  double dx = 1.0/static_cast<double>(num_nodes_-1);
  ublas::bounded_vector<double, 3> sat;
  CalcSAT(bc_left_, area_(0), 1.0, q(0), sat);
  for (int i = 0; i < 3; i++)
    v(i) += sbp_deriv_.Hinv()*sat(i)*u(0)/area_(0);
  int nm1 = num_nodes_-1;
  CalcSAT(bc_right_, area_(nm1), -1.0, q(nm1), sat);
  for (int i = 0; i < 3; i++)
    v(3*nm1 + i) += sbp_deriv_.Hinv()*sat(i)*u(nm1)/area_(nm1);
}

// ======================================================================

void Quasi1DEuler::JacobianTransposedAreaProduct(
    const InnerProdVector & u, InnerProdVector & v) {
  // check for consistent sizes
  if ( (u.size() != 3*num_nodes_) || (v.size() != num_nodes_) ) {
    cerr << "Quasi1DEuler::JacobianTransposedAreaProduct(): "
         << "inconsistent sizes.";
    throw(-1);
  }
  v = 0.0;

  // evaluate the pressure, soundspeed, and spectral radius
  CalcAuxiliaryVariables(q_);

  // add terms that depend on the source
  // Note: do this first so we do not need an extra work vector
  InnerProdVector work(num_nodes_, 0.0);
  for (int i = 0; i < num_nodes_; i++) {
    int ptr = (i*3) + 1;
    work(i) = -press_(i)*u(ptr);
  }
  sbp_deriv_.ApplyTranspose(1, work, v);

  // evaluate terms that depend on Euler flux
  InnerProdVector transDu(3*num_nodes_, 0.0);
  sbp_deriv_.ApplyTranspose(3, u, transDu);
  for (int i = 0; i < num_nodes_; i++) {
    double u_i = u(i);
    double press = press_(i);
    double rho = q(i,0);
    double vel = q(i,1) / rho;
    double e = q(i,2);
    int ptr = (i*3);
    v(i) += rho * vel * transDu(ptr+0);
    v(i) += (rho * vel * vel + press) * transDu(ptr+1);
    v(i) += vel * ( e + press) * transDu(ptr+2);
  }

  // numerical dissipation does not depend on area (no spect_ term)

  // add terms that depend on SAT boundary penalties
  double dx = 1.0/static_cast<double>(num_nodes_-1);
  ublas::bounded_vector<double, 3> sat;
  CalcSAT(bc_left_, area_(0), 1.0, q(0), sat);
  for (int i = 0; i < 3; i++)
    v(0) += sbp_deriv_.Hinv()*sat(i)*u(i)/area_(0);
  int nm1 = num_nodes_-1;
  CalcSAT(bc_right_, area_(nm1), -1.0, q(nm1), sat);
  for (int i = 0; i < 3; i++)
    v(nm1) += sbp_deriv_.Hinv()*sat(i)*u(3*nm1+i)/area_(nm1);
}

// ======================================================================

void Quasi1DEuler::TestJacobianAreaProducts() {
  // create a random vector to apply transposed Jacobian to
  InnerProdVector u(num_nodes_, 0.0), v(3*num_nodes_, 0.0),
      w(3*num_nodes_, 0.0), z(num_nodes_, 0.0);
  boost::random::mt19937 gen;
  boost::random::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (int i = 0; i < num_nodes_; i++)
    u(i) = dist(gen);
  for (int i = 0; i < 3*num_nodes_; i++)
    v(i) = dist(gen);

#if 0
  int i = 0;
  ublas::range irange(3*i, 3*(i+1));
  u(irange) = ublas::zero_vector<double>(3);
  v(irange) = ublas::zero_vector<double>(3);
  i = num_nodes_-1;
  irange = ublas::range(3*i, 3*(i+1));
  u(irange) = ublas::zero_vector<double>(3);
  v(irange) = ublas::zero_vector<double>(3);
#endif
  // evaluate Jacobian-vector product and contract with v
  JacobianAreaProduct(u, w);
  double forward = InnerProd(v, w);
  // evaluate the transposed-Jacobian-vector product and contract with u
  JacobianTransposedAreaProduct(v, z);
  double backward = InnerProd(z, u);
  cout << "Quasi1DEuler::TestJacobianAreaProducts():" << endl;
  cout << "Difference between forward and backward = "
       << backward - forward << endl;
}

// ======================================================================

void Quasi1DEuler::ExplicitEuler(const int & max_iter,
                                 const double & target_cfl,
                                 const double & tol) {
  int iter = 0;
  while (iter < max_iter) {
    CalcResidual();
    double norm = ResidualNorm();
    if ( iter % 100 == 0)
      cout << "iter = " << iter
           << ": L2 norm of residual = " << norm << endl;
    if (norm < tol) {
      cout << "iter = " << iter
           << ": L2 norm of residual = " << norm << endl;
      return;
    }
    double dt = CalcDtFromCFL(target_cfl);
    double dx; // 1.0/static_cast<double>(num_nodes_-1);
    for (int i = 0; i < num_nodes_; i++) {
      dx = met_jac_(i);
      q(i) -= (dt/area_(i))*res(i);
    }
    iter++;
  }
  // if we get here, we failed to converge
  cout << "Quasi1DEuler::ExplicitEuler(): "
       << "failed to converge in " << max_iter << " iterations." << endl;
  //throw(-1);
}

// ======================================================================

int Quasi1DEuler::NewtonKrylov(const int & max_iter,
                                const double & tol) {
  kona::MatrixVectorProduct<InnerProdVector>*
      mat_vec = new JacobianVectorProduct(this);
  kona::Preconditioner<InnerProdVector>*
      precond = new ApproxJacobian(this);

  string filename = "quasi1d_primal_krylov.dat";
  ofstream fout(filename.c_str());

  int iter = 0;
  int precond_calls = 0;
  while (iter < max_iter) {
    // evaluate the residual and its norm
    CalcResidual();
    double norm = ResidualNorm();
    InnerProdVector b(-res_);
    cout << "iter = " << iter
         << ": L2 norm of residual = " << norm << endl;
    if ( (norm < tol) || (norm < 1.e-14) ) {
      cout << "Quasi1DEuler: NewtonKrylov converged" << endl;
      return precond_calls;
    }
    // create the approximate-Jacobian LU preconditioner
    BuildAndFactorPreconditioner();
    // solve for the Newton update dq and add to q
    int m = 10;
    double tol = 1.0e-2;
    InnerProdVector dq(3*num_nodes_, 0.0);
    int krylov_precond_calls;
    try {
      kona::FGMRES(m, tol, b, dq, *mat_vec, *precond,
                   krylov_precond_calls, fout);
    } catch (...) {
      cout << "Quasi1DEuler: FGMRES failed in NewtonKrylov" << endl;
      return -precond_calls;
    }

    q_ += dq;
    precond_calls += krylov_precond_calls;
    iter++;
  }
  // if we get here, we failed to converge
  cout << "Quasi1DEuler::NewtonKrylov(): "
       << "failed to converge in " << max_iter << " iterations." << endl;
  //throw(-1);
  return -precond_calls;
}

// ======================================================================

int Quasi1DEuler::SolveAdjoint(const int & max_iter, const double & tol,
                               const InnerProdVector & dJdQ,
                               InnerProdVector & psi) {
  kona::MatrixVectorProduct<InnerProdVector>*
      mat_vec = new JacobianTransposedVectorProduct(this);
  kona::Preconditioner<InnerProdVector>*
      precond = new ApproxJacobianTransposed(this);

  string filename = "quasi1d_adjoint_krylov.dat";
  ofstream fout(filename.c_str());

  BuildAndFactorPreconditioner();
  psi = 0.0;
  int precond_calls = 0;
  kona::FGMRES(max_iter, tol, dJdQ, psi, *mat_vec, *precond,
               precond_calls, fout);
  return precond_calls;
}

// ======================================================================

int Quasi1DEuler::SolveLinearized(const int & max_iter,
                                  const double & tol,
                                  const InnerProdVector & rhs,
                                  InnerProdVector & dq) {
  kona::MatrixVectorProduct<InnerProdVector>*
      mat_vec = new JacobianVectorProduct(this);
  kona::Preconditioner<InnerProdVector>*
      precond = new ApproxJacobian(this);

  string filename = "quasi1d_linearized_krylov.dat";
  ofstream fout(filename.c_str());

  BuildAndFactorPreconditioner();
  dq = 0.0;
  int precond_calls = 0;
  kona::FGMRES(max_iter, tol, rhs, dq, *mat_vec, *precond,
               precond_calls, fout);
  return precond_calls;
}

// ======================================================================

int Quasi1DEuler::SolveUnsteady(const int & iter, const double & Time,
                                const double & tol, const bool & store,
                                const string & flow_file, const bool & write) {
  kona::MatrixVectorProduct<InnerProdVector>*
      mat_vec = new UnsteadyJacobianVectorProduct(this);
  kona::Preconditioner<InnerProdVector>*
      precond = new ApproxJacobian(this);
  const int max_iter = 15;

  if ( (write) && (!store)) {
    cerr << "Error in Quasi1DEuler::SolveUnsteady(): "
         << "must have store = true if write = true." << endl;
    throw(-1);
  }
  if (iter != src_.size()) {
    cerr << "Error in Quasi1DEuler::SolveUnsteady(): "
         << "number of iterations inconsistent with source size." << endl;
    throw(-1);
  }

  string filename = "quasi1d_primal_krylov.dat";
  ofstream krylov_out(filename.c_str());
  ofstream save_out;
  if (store) {
    // open file to save flow for future adjoint solve
    save_out.open(flow_file.c_str(), ios::out | ios::binary);
    // save header information
    save_out.write(reinterpret_cast<const char*>(&num_nodes_), sizeof(int));
    save_out.write(reinterpret_cast<const char*>(&iter), sizeof(int));
    save_out.write(reinterpret_cast<const char*>(&Time), sizeof(double));
    save_out.flush();
  }

  // loop over the time steps
  int precond_calls = 0;
  dt_ = Time/iter;
  q_ = q_old_;
  if (store) q_.BinaryWrite(save_out);
  for (int n = 0; n < iter; n++) {
    if (write) cout << "iteration = " << n
                    << ": time = " << static_cast< double >(n )*dt() << endl;
    // loop over the Newton iterations
    int newt_iter = 0;
    while (newt_iter < max_iter) {
      // evaluate the unsteady residual and its norm;
      // NOTE: source term must be added separately
      CalcUnsteadyResidual();
      AddUnsteadySource(n);
      //res_(3*src_indx_+1) -= 0.5*(src_(n) + src_(n+1));
      double norm = ResidualNorm();
      InnerProdVector b(-res_);
      if (write) cout << "\titer = " << newt_iter
                      << ": L2 norm of residual = " << norm << endl;
      if ( (norm < tol) || (norm < 1.e-14) ) {
        if (write) cout << "\tQuasi1DEuler: SolveUnsteady Newton converged"
                        << endl;
        break;
        //return precond_calls;
      }
      // create the approximate-Jacobian LU preconditioner
      BuildAndFactorUnsteadyPreconditioner();
      // solve for the Newton update dq and add to q
      int m = 15;
      double tol = 1.0e-2;
      InnerProdVector dq(3*num_nodes_, 0.0);
      int krylov_precond_calls;
      kona::FGMRES(m, tol, b, dq, *mat_vec, *precond,
                   krylov_precond_calls, krylov_out);
      q_ += dq;
      precond_calls += krylov_precond_calls;
      newt_iter++;
    }
    if (newt_iter == max_iter) {
      // if we get here, we failed to converge
      cout << "Quasi1DEuler::SolveUnsteady(): "
           << "failed to converge in " << max_iter << " iterations." << endl;
      //throw(-1);
      //fout.close();
      if (store) save_out.close();
      return -precond_calls;
    }
    q_old_ = q_;
    if (store) q_.BinaryWrite(save_out);
  }
  cout << "Quasi1DEuler::SolveUnsteady(): converged" << endl;
  if (store) save_out.close();
  if (write) {
    filename = "unsteady_flow.dat";
    WriteUnsteadyTecplot(string("save_flow.bin"), filename);
  }
  return precond_calls;
}

// ======================================================================

int Quasi1DEuler::SolveUnsteadyAdjoint(const string & flow_file,
                                       const int & max_krylov,
                                       const double & tol,
                                       const string & dJdQ_file,
                                       const string & psi_file,
                                       const bool & init_adjoint) {
  kona::MatrixVectorProduct<InnerProdVector>*
      mat_vec = new UnsteadyJacTransVectorProduct(this);
  kona::Preconditioner<InnerProdVector>*
      precond = new ApproxJacobianTransposed(this);

  string filename = "quasi1d_adjoint_krylov.dat";
  ofstream krylov_out(filename.c_str());

  // open the solution file, the dJdQ file, and adjoint file
  ifstream fin_q(flow_file.c_str(), ios::in | ios::binary);
  ifstream fin_dJdQ(dJdQ_file.c_str(), ios::in | ios::binary);
  ofstream fout_psi(psi_file.c_str(), ios::out | ios::binary);
  if ( (!fin_q.good()) || (!fin_dJdQ.good()) ) {
    cerr << "Error in Quasi1DEuler::SolveUnsteadyAdjoint(): error opening "
         << flow_file << " or " << dJdQ_file << endl;
    throw(-1);
  }
  // check that num_nodes_ is consistent with value stored
  int nodes;
  fin_q.read(reinterpret_cast<char*>(&nodes), sizeof(int));
  if (nodes != num_nodes_) {
    cerr << "Error in Quasi1DEuler::SolveUnsteadyAdjoint(): "
         << "number of nodes does not match value in saved flow file." << endl;
    cout << "nodes = " << nodes << ": num_nodes_ = " << num_nodes_ << endl;
    throw(-1);
  }
  fin_dJdQ.read(reinterpret_cast<char*>(&nodes), sizeof(int));
  if (nodes != num_nodes_) {
    cerr << "Error in Quasi1DEuler::SolveUnsteadyAdjoint(): "
         << "number of nodes does not match value in saved dJdQ file." << endl;
    cout << "nodes = " << nodes << ": num_nodes_ = " << num_nodes_ << endl;
    throw(-1);
  }
  // get the number of iterations and total time; check for consistency here too
  int iter;
  fin_q.read(reinterpret_cast<char*>(&iter), sizeof(int));
  double total_time;
  fin_q.read(reinterpret_cast<char*>(&total_time), sizeof(double));
  dt_ = total_time/static_cast<double>(iter);
  int check_iter;
  double check_time;
  fin_dJdQ.read(reinterpret_cast<char*>(&check_iter), sizeof(int));
  fin_dJdQ.read(reinterpret_cast<char*>(&check_time), sizeof(double));
  if ( (check_iter != iter) || (check_time != total_time) ) {
    cerr << "Error in Quasi1DEuler::SolveUnsteadyAdjoint(): "
         << "inconsistency between flow and dJdQ files (iter or time)." << endl;
    throw(-1);
  }

  // write header for psi_file
  fout_psi.write(reinterpret_cast<const char*>(&num_nodes_), sizeof(int));
  fout_psi.write(reinterpret_cast<const char*>(&iter), sizeof(int));
  fout_psi.write(reinterpret_cast<const char*>(&total_time), sizeof(double));

  // solve for the nth (i.e. first) adjoint
  unsigned long fptr = 2*sizeof(int) + sizeof(double)
      + iter*(3*num_nodes_*sizeof(double));
  q_.BinaryRead(fin_q, fptr);

  InnerProdVector dJdQ(3*num_nodes_, 0.0);
  dJdQ.BinaryRead(fin_dJdQ, fptr);
  InnerProdVector b(-dJdQ);
  if (init_adjoint)
    b += psi_;

  fptr -= 3*num_nodes_*sizeof(double);
  q_old_.BinaryRead(fin_q, fptr);

  BuildAndFactorUnsteadyPreconditioner();
  psi_ = 0.0;
  int precond_calls = 0;
  int krylov_precond_calls;
  kona::FGMRES(max_krylov, tol, b, psi_, *mat_vec, *precond,
               krylov_precond_calls, krylov_out);
  precond_calls += krylov_precond_calls;
  psi_.BinaryWrite(fout_psi, fptr + 3*num_nodes_*sizeof(double));

  // reverse loop over time
  for (int n = iter-1; n > 0; n--) {
    // reuse q_ before copying in q_old
    UnsteadyJacTransStateProduct(psi_, b, false);
    dJdQ.BinaryRead(fin_dJdQ, fptr);
    b.EqualsAXPlusBY(-1.0, b, -1.0, dJdQ);

    q_ = this->q_old_;
    fptr -= 3*num_nodes_*sizeof(double);
    q_old_.BinaryRead(fin_q, fptr);

    BuildAndFactorUnsteadyPreconditioner();
    psi_ = 0.0;
    kona::FGMRES(max_krylov, tol, b, psi_, *mat_vec, *precond,
                 krylov_precond_calls, krylov_out);
    precond_calls += krylov_precond_calls;
    psi_.BinaryWrite(fout_psi, fptr + 3*num_nodes_*sizeof(double));
  }

  // final adjoint step (the system matrix for psi_ is the identity)
  UnsteadyJacTransStateProduct(psi_, b, false);
  dJdQ.BinaryRead(fin_dJdQ, fptr);
  b.EqualsAXPlusBY(-1.0, b, -1.0, dJdQ);
  psi_ = b;
  psi_.BinaryWrite(fout_psi, fptr);

  fout_psi.close();
  fin_q.close();
  fin_dJdQ.close();
  return precond_calls;
}

// ======================================================================

int Quasi1DEuler::SolveUnsteadyIterLinearized(const int & max_iter,
                                              const double & tol,
                                              const InnerProdVector & rhs,
                                              InnerProdVector & dq) {
  kona::MatrixVectorProduct<InnerProdVector>*
      mat_vec = new UnsteadyJacobianVectorProduct(this);
  kona::Preconditioner<InnerProdVector>*
      precond = new ApproxJacobian(this);

  string filename = "quasi1d_linearized_krylov.dat";
  ofstream fout(filename.c_str());

  BuildAndFactorUnsteadyPreconditioner();
  dq = 0.0;
  int precond_calls = 0;
  kona::FGMRES(max_iter, tol, rhs, dq, *mat_vec, *precond,
               precond_calls, fout);
  return precond_calls;
}

// ======================================================================

int Quasi1DEuler::SolveUnsteadyAdjointIter(const int & max_iter,
                                           const double & tol,
                                           const InnerProdVector & rhs,
                                           InnerProdVector & psi) {
  kona::MatrixVectorProduct<InnerProdVector>*
      mat_vec = new UnsteadyJacTransVectorProduct(this);
  kona::Preconditioner<InnerProdVector>*
      precond = new ApproxJacobianTransposed(this);

  string filename = "quasi1d_adjoint_krylov.dat";
  ofstream fout(filename.c_str());

  BuildAndFactorUnsteadyPreconditioner();
  psi = 0.0;
  int precond_calls = 0;
  kona::FGMRES(max_iter, tol, rhs, psi, *mat_vec, *precond,
               precond_calls, fout);
  return precond_calls;
}

// ======================================================================

void Quasi1DEuler::WriteTecplot(const double & rhoL, const double & aL,
                                const string & filename) {
  CalcAuxiliaryVariables(q_);
  ofstream fout(filename.c_str());
  fout.precision(12);
  fout << "TITLE = \"Quasi-1D-Euler Steady Nozzle Solution\"" << endl;
  fout << "VARIABLES=\"x\",\"area\",\"rho\",\"rho-u\",\"e\",\"press\""
       << ",\"press target\",\"u\",\"Mach\",\"Mach exact\"" << endl;
  fout << "ZONE I=" << num_nodes_ << ", DATAPACKING=POINT" << endl;
  for (int i = 0; i < num_nodes_; i++) {
    fout << x_coord_(i) << " ";
    fout << area_(i) << " ";
    fout << rhoL*q(i,0) << " ";
    fout << rhoL*aL*q(i,1) << " ";
    fout << aL*aL*q(i,2) << " ";
    fout << rhoL*aL*aL*press_(i) << " ";
    fout << rhoL*aL*aL*press_targ_(i) << " ";
    fout << q(i,1)*aL/q(i,0) << " ";
    fout << q(i,1)/(q(i,0)*sndsp_(i)) << " ";
    double mach_exact =
        CalcMachExact<double>(kGamma, 0.8, area_(i), true);
    fout << mach_exact << " ";
#if 0
    // uncomment to calculate pointwise Mach error

    fout << fabs(q(i,1)/(q(i,0)*sndsp_(i)) - mach_exact) << " ";
#endif
    fout << endl;
  }
  fout.close();
}

// ======================================================================

void Quasi1DEuler::WriteUnsteadyTecplot(const int & iter, const double & dt,
                                        ostream & fout) {
  CalcAuxiliaryVariables(q_);
  fout.precision(12);
  if (iter == 0) {
    fout << "TITLE = \"1D-Euler Unsteady Solution\"" << endl;
    fout << "VARIABLES=\"x\",\"area\",\"rho\",\"rho-u\",\"e\",\"press\""
         << ",\"u\",\"Mach\"" << endl;
  }
  double time = static_cast<double>(iter)*dt;
  fout << "ZONE SOLUTIONTIME=" << time << ", I=" << num_nodes_
       << ", DATAPACKING=BLOCK" << endl;
  for (int i = 0; i < num_nodes_; i++)
    fout << x_coord_(i) << " ";
  for (int i = 0; i < num_nodes_; i++)
    fout << area_(i) << " ";
  for (int i = 0; i < num_nodes_; i++)
    fout << q(i,0) << " ";
  for (int i = 0; i < num_nodes_; i++)
    fout << q(i,1) << " ";
  for (int i = 0; i < num_nodes_; i++)
    fout << q(i,2) << " ";
  for (int i = 0; i < num_nodes_; i++)
    fout << press_(i) << " ";
  for (int i = 0; i < num_nodes_; i++)
    fout << q(i,1)/q(i,0) << " ";
  for (int i = 0; i < num_nodes_; i++)
    fout << q(i,1)/(q(i,0)*sndsp_(i)) << " ";
  fout << endl;
}

// ======================================================================

void Quasi1DEuler::WriteUnsteadyTecplot(const string & flow_file,
                                        const string & tec_file) {
  ifstream fin(flow_file.c_str(), ios::in | ios::binary);
  ofstream fout(tec_file.c_str(), ios::out);
  fout.precision(12);
  fout << "TITLE = \"1D-Euler Unsteady Solution\"" << endl;
  fout << "VARIABLES=\"x\",\"area\",\"rho\",\"rho-u\",\"e\",\"press\""
       << ",\"u\",\"Mach\"" << endl;
  //fout << "VARIABLES=\"x\",\"area\",\"rho\",\"rho-u\",\"e\"" << endl;

  // check that num_nodes_ is consistent with value stored
  int nodes;
  fin.read(reinterpret_cast<char*>(&nodes), sizeof(int));
  if (nodes != num_nodes_) {
    cerr << "Error in Quasi1DEuler::WriteUnsteadyTecplot(): "
         << "number of nodes does not match value in saved file." << endl;
    cerr << "nodes = " << nodes << ": num_nodes_ = " << num_nodes_ << endl;
    throw(-1);
  }
  // get the number of iterations and total time
  int iter;
  fin.read(reinterpret_cast<char*>(&iter), sizeof(int));
  double total_time;
  fin.read(reinterpret_cast<char*>(&total_time), sizeof(double));
  double delta_t = total_time/static_cast<double>(iter);

  cout << "iter = " << iter << ": total_time = " << total_time << endl;

  // loop over all times, and write tecplot zones for each
  for (int n = 0; n <= iter; n++) {
    q_.BinaryRead(fin);
    CalcAuxiliaryVariables(q_);
    double time = static_cast<double>(n)*delta_t;
    fout << "ZONE SOLUTIONTIME=" << time << ", I=" << num_nodes_
         << ", DATAPACKING=BLOCK" << endl;
    for (int i = 0; i < num_nodes_; i++)
      fout << x_coord_(i) << " ";
    for (int i = 0; i < num_nodes_; i++)
      fout << area_(i) << " ";
    for (int i = 0; i < num_nodes_; i++)
      fout << q(i,0) << " ";
    for (int i = 0; i < num_nodes_; i++)
      fout << q(i,1) << " ";
    for (int i = 0; i < num_nodes_; i++)
      fout << q(i,2) << " ";
    for (int i = 0; i < num_nodes_; i++)
      fout << press_(i) << " ";
    for (int i = 0; i < num_nodes_; i++)
      fout << q(i,1)/q(i,0) << " ";
    for (int i = 0; i < num_nodes_; i++)
      fout << q(i,1)/(q(i,0)*sndsp_(i)) << " ";
    fout << endl;
  }
  fin.close();
  fout.close();
}

// ======================================================================

void Quasi1DEuler::WriteUnsteadyAdjointTecplot(const string & psi_file,
                                               const string & tec_file) {
  ifstream fin(psi_file.c_str(), ios::in | ios::binary);
  ofstream fout(tec_file.c_str(), ios::out);
  fout.precision(12);
  fout << "TITLE = \"1D-Euler Unsteady Adjoint Solution\"" << endl;
  fout << "VARIABLES=\"x\",\"area\",\"psi-rho\",\"psi-rho-u\",\"psi-e\"" << endl;
  //fout << "VARIABLES=\"x\",\"area\",\"rho\",\"rho-u\",\"e\"" << endl;

  // check that num_nodes_ is consistent with value stored
  int nodes;
  fin.read(reinterpret_cast<char*>(&nodes), sizeof(int));
  if (nodes != num_nodes_) {
    cerr << "Error in Quasi1DEuler::WriteUnsteadyAdjointTecplot(): "
         << "number of nodes does not match value in saved file." << endl;
    cerr << "nodes = " << nodes << ": num_nodes_ = " << num_nodes_ << endl;
    throw(-1);
  }
  // get the number of iterations and total time
  int iter;
  fin.read(reinterpret_cast<char*>(&iter), sizeof(int));
  double total_time;
  fin.read(reinterpret_cast<char*>(&total_time), sizeof(double));
  double delta_t = total_time/static_cast<double>(iter);

  cout << "iter = " << iter << ": total_time = " << total_time << endl;

  // loop over all times in reverse, and write tecplot zones for each
  unsigned long fptr = 2*sizeof(int) + sizeof(double)
      + (iter+1)*(3*num_nodes_*sizeof(double));
  for (int n = iter; n >= 0; n--) {
    fptr -= 3*num_nodes_*sizeof(double);
    psi_.BinaryRead(fin, fptr);
    // multiply psi by inverse SBP norm
    sbp_deriv_.HinvTimesVector(3, psi_, psi_);
    psi_ /= dxi();
    double time = static_cast<double>(iter-n)*delta_t;
    fout << "ZONE SOLUTIONTIME=" << time << ", I=" << num_nodes_
         << ", DATAPACKING=BLOCK" << endl;
    for (int i = 0; i < num_nodes_; i++)
      fout << x_coord_(i) << " ";
    for (int i = 0; i < num_nodes_; i++)
      fout << area_(i) << " ";
    for (int i = 0; i < num_nodes_; i++)
      fout << psi_(3*i) << " ";
    for (int i = 0; i < num_nodes_; i++)
      fout << psi_(3*i+1) << " ";
    for (int i = 0; i < num_nodes_; i++)
      fout << psi_(3*i+2) << " ";
    fout << endl;
  }
  fin.close();
  fout.close();
}

// ======================================================================

void Quasi1DEuler::CalcMachError(const double & area_star,
                                 const bool & subsonic,
                                 double & L2_error, double & max_error) {
  if (!subsonic) {
    cerr << "Quasi1DEuler::CalcMachError(): "
         << "presently not set up for supersonic or transonic cases"
         << endl;
    throw(-1);
  }
  max_error = 0.0;
  InnerProdVector error(num_nodes_, 0.0);
  for (int i = 0; i < num_nodes_; i++) {
    // compute the exact Mach number and find the error at this node
    double mach_exact =
        CalcMachExact<double>(kGamma, area_star, area_(i), subsonic);
    error(i) = fabs(mach_exact - q(i,1)/(q(i,0)*sndsp_(i)));
    if (error(i) > max_error) max_error = error(i);
    error(i) *= error(i);
    //cout << "node " << i << ": Mach error = " << error(i) << endl;
  }
  L2_error = sqrt(sbp_deriv_.InnerProductSBP(1, met_jac_, error));
}

// ======================================================================

double Quasi1DEuler::CalcTotalEnergy(const bool & sbp_quad) {
  InnerProdVector energy(num_nodes_, 0.0);
  for (int i = 0; i < num_nodes_; i++)
    energy(i) = 0.5*q(i,1)*q(i,1)/q(i,0); //q(i,2);
  if (sbp_quad) {
    return sbp_deriv_.InnerProductSBP(1, met_jac_, energy);
  } else {
    if (num_nodes_ % 2 != 1) {
      cerr << "Quasi1DEuler::calcTotalEnergy(): "
           << "sbp_quad == false and num_ndoes_ not odd.";
      throw(-1);
    }
    double total = 0.0;
    for (int i = 0; i < num_nodes_; i++) {
      if ( (i == 0) || (i == num_nodes_-1) ) {
        total += energy(i)*met_jac_(i);
      } else {
        if (i % 2 == 0) {
          total += 2.0*energy(i)*met_jac_(i);
        } else {
          total += 4.0*energy(i)*met_jac_(i);
        }
      }
    }
    total /= 3.0; //(3.0*static_cast<double>(num_nodes_-1));
    return total;
  }
}

// ======================================================================

void Quasi1DEuler::CalcTotalEnergydJdQ(InnerProdVector & dJdQ) {
  if (dJdQ.size() != 3*num_nodes_) {
    cerr << "Quasi1DEuler(CalcTotalEnergydJdQ): "
         << "size of dJdQ is inconsistent with number of nodes"
         << endl;
    throw(-1);
  }
  for (int i = 0; i < num_nodes_; i++) {
    double vel = q(i,1)/q(i,0);
    int ptr = 3*i;
    dJdQ(ptr  ) = -0.5*vel*vel*met_jac_(i);
    dJdQ(ptr+1) =  vel*met_jac_(i);
    dJdQ(ptr+2) =  0.0;
  }
  // scale dJdQ by H norm
  sbp_deriv_.HTimesVector(3, dJdQ, dJdQ);
}

// ======================================================================

void Quasi1DEuler::CalcTotalEnergyd2JdQ2(const InnerProdVector & w,
                                         InnerProdVector & d2JdQ2) {
  if (d2JdQ2.size() != 3*num_nodes_) {
    cerr << "Quasi1DEuler(CalcTotalEnergyd2JdQ2): "
         << "size of d2JdQ2 is inconsistent with number of nodes"
         << endl;
    throw(-1);
  }
  if (w.size() != 3*num_nodes_) {
    cerr << "Quasi1DEuler(CalcTotalEnergyd2JdQ2): "
         << "size of w is inconsistent with number of nodes"
         << endl;
    throw(-1);
  }
  for (int i = 0; i < num_nodes_; i++) {
    double rho = q(i,0);
    double vel = q(i,1)/rho;
    int ptr = 3*i;
    double w1 = w(ptr);
    double w2 = w(ptr+1);
    d2JdQ2(ptr  ) = met_jac_(i)*(vel*vel*w1 - vel*w2)/rho;
    d2JdQ2(ptr+1) = met_jac_(i)*(w2 - vel*w1)/rho;
    d2JdQ2(ptr+2) = 0.0;
  }
  // scale dJdQ by H norm
  sbp_deriv_.HTimesVector(3, d2JdQ2, d2JdQ2);
}

// ======================================================================

double Quasi1DEuler::CalcInverseDesign() {
  CalcAuxiliaryVariables(q_);
  InnerProdVector dpress(num_nodes_, 0.0);
  for (int i = 0; i < num_nodes_; i++)
    dpress(i) = 0.5*(press_(i) - press_targ_(i))
        *(press_(i) - press_targ_(i));
  return sbp_deriv_.InnerProductSBP(1, met_jac_, dpress);
}

// ======================================================================

void Quasi1DEuler::CalcInverseDesigndJdQ(InnerProdVector & dJdQ) {
  if (dJdQ.size() != 3*num_nodes_) {
    cerr << "Quasi1DEuler(CalcInverseDesigndJdQ): "
         << "size of dJdQ is inconsistent with number of nodes"
         << endl;
    throw(-1);
  }
  CalcAuxiliaryVariables(q_);
  for (int i = 0; i < num_nodes_; i++) {
    double dpress = (press_(i) - press_targ_(i))
        *met_jac_(i)*(kGamma-1.0);
    double vel = q(i,1)/q(i,0);
    int ptr = 3*i;
    dJdQ(ptr  ) =  dpress*0.5*vel*vel;
    dJdQ(ptr+1) = -dpress*vel;
    dJdQ(ptr+2) =  dpress;
  }
  // scale dJdQ by H norm
  sbp_deriv_.HTimesVector(3, dJdQ, dJdQ);
}

// ======================================================================

void Quasi1DEuler::CalcInverseDesignd2JdQ2(const InnerProdVector & w,
                                           InnerProdVector & d2JdQ2) {
  if (d2JdQ2.size() != 3*num_nodes_) {
    cerr << "Quasi1DEuler(CalcInverseDesignd2JdQ2): "
         << "size of d2JdQ2 is inconsistent with number of nodes"
         << endl;
    throw(-1);
  }
  if (w.size() != 3*num_nodes_) {
    cerr << "Quasi1DEuler(CalcInverseDesignd2JdQ2): "
         << "size of w is inconsistent with number of nodes"
         << endl;
    throw(-1);
  }
  CalcAuxiliaryVariables(q_);
  const double kGami = kGamma - 1.0;
  for (int i = 0; i < num_nodes_; i++) {
    double met = met_jac_(i);
    double dpress = (press_(i) - press_targ_(i));
    double rho = q(i,0);
    double vel = q(i,1)/rho;
    int ptr = 3*i;
    double w1 = w(ptr);
    double w2 = w(ptr+1);
    double w3 = w(ptr+2);
    double wdpdq = kGami*(0.5*w1*vel*vel - w2*vel + w3);
    d2JdQ2(ptr  ) = met*kGami*(wdpdq*0.5*vel*vel
                         + dpress*(w2*vel/rho - w1*vel*vel/rho));
    d2JdQ2(ptr+1) = met*kGami*(-wdpdq*vel
                         + dpress*(w1*vel/rho - w2/rho));
    d2JdQ2(ptr+2) = met*wdpdq*kGami;
  }
  // scale dJdQ by H norm
  sbp_deriv_.HTimesVector(3, d2JdQ2, d2JdQ2);
}

// ======================================================================

void Quasi1DEuler::Testd2JdQ2(const objective & obj) {
  // create a random vector to apply dJ2dQ2
  InnerProdVector w(3*num_nodes_, 0.0), d2JdQ2(3*num_nodes_, 0.0),
      dJdQ(3*num_nodes_, 0.0), d2JdQ2_fd(3*num_nodes_, 0.0),
      q_save(3*num_nodes_, 0.0);
  boost::random::mt19937 gen;
  boost::random::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (int i = 0; i < 3*num_nodes_; i++)
    w(i) = dist(gen);

  // evaluate d2JdQ2*w product analytically
  switch (obj) {
    case inverse: {// inverse design
      CalcInverseDesignd2JdQ2(w, d2JdQ2);
      break;
    }
    case total_kinetic: {// total kinetic energy
      CalcTotalEnergyd2JdQ2(w, d2JdQ2);
      break;
    }
    default: {
      cerr << "Quasi1DEuler(Testd2JdQ2): "
           << "invalid objective function"
           << endl;
      throw(-1);
      break;
    }
  }

  // evaluate the d2JdQ2*w product using backward difference

  // evaluate and save dJdQ*w
  q_save = q_;  // save flow state for later
  switch (obj) {
    case inverse: {// inverse design
      CalcInverseDesigndJdQ(d2JdQ2_fd);
      break;
    }
    case total_kinetic: {// total kinetic energy
      CalcTotalEnergydJdQ(d2JdQ2_fd);
      break;
    }
    default: {
      cerr << "Quasi1DEuler(Testd2JdQ2): "
           << "invalid objective function"
           << endl;
      throw(-1);
      break;
    }
  }

  // perturb flow and re-evaluate dJdQ
  double fd_eps = 1.E-7;
  q_ -= fd_eps*w;
  switch (obj) {
    case inverse: {// inverse design
      CalcInverseDesigndJdQ(dJdQ);
      break;
    }
    case total_kinetic: {// total kinetic energy
      CalcTotalEnergydJdQ(dJdQ);
      break;
    }
    default: {
      cerr << "Quasi1DEuler(Testd2JdQ2): "
           << "invalid objective function"
           << endl;
      throw(-1);
      break;
    }
  }
  d2JdQ2_fd -= dJdQ;
  d2JdQ2_fd /= fd_eps;

  // take difference between two products and store in q_ for output
  q_.EqualsAXPlusBY(1.0, d2JdQ2, -1.0, d2JdQ2_fd);
  double L2_error = sbp_deriv_.NormSBP(3, q_);
  cout << "Quasi1DEuler::Testd2JdQ2(): "
       << "L2 error between analytical and FD products: "
       << L2_error << endl;
  q_ = q_save;
}

// ======================================================================
#if 0
double Quasi1DEuler::CalcSensor(const double & x_beg, const double & x_end,
                                const double & press_targ,
                                const string & flow_file) {
  ifstream fin(flow_file.c_str(), ios::in | ios::binary);
  if (!fin.good()) {
    cerr << "Error in Quasi1DEuler::CalcSensor(): error opening "
         << flow_file << endl;
    throw(-1);
  }
  // check that num_nodes_ is consistent with value stored
  int nodes;
  fin.read(reinterpret_cast<char*>(&nodes), sizeof(int));
  if (nodes != num_nodes_) {
    cerr << "Error in Quasi1DEuler::CalcSensor(): "
         << "number of nodes does not match value in saved file." << endl;
    cout << "nodes = " << nodes << ": num_nodes_ = " << num_nodes_ << endl;
    throw(-1);
  }
  // get the number of iterations and total time
  int iter;
  fin.read(reinterpret_cast<char*>(&iter), sizeof(int));
  double total_time;
  fin.read(reinterpret_cast<char*>(&total_time), sizeof(double));
  double delta_t = total_time/static_cast<double>(iter);

  // find nodes cooresponding to x_beg and x_end;
  int i_beg = -1;
  int i_end = -1;
  for (int i = 0; i < num_nodes_-1; i++) {
    if ( (x_coord_(i) <= x_beg) && (x_beg < x_coord_(i+1)) ) i_beg = i;
    if ( (x_coord_(i) < x_end) && (x_end <= x_coord_(i+1)) ) i_end = i+2;
  }
  if ( (i_beg == -1) || (i_end == -1) ) {
    cerr << "Error in Quasi1DEuler::CalcSensor(): "
         << "failed to find valid nodes for x_beg and x_end integration region"
         << endl;
    throw(-1);
  }

  q_.BinaryRead(fin);
  InnerProdVector q_mid(3*num_nodes_,0.0);
  int num_pts = i_end - i_beg;
  InnerProdVector dpress(num_pts, 0.0);
  SBP1stDerivative sbp_sensor(num_pts, sbp_deriv_.order());

  // loop over all times, and calculate the sensor objective
  double sensor = 0.0;
  for (int n = 0; n < iter; n++) {
    q_old_ = q_;
    q_.BinaryRead(fin);
    q_mid.EqualsAXPlusBY(0.5, q_, 0.5, q_old_);
    CalcAuxiliaryVariables(q_mid);
    for (int i = i_beg; i < i_end; i++)
      dpress(i-i_beg) = press_(i) - press_targ;
    sensor += delta_t*0.5*dxi()*sbp_sensor.InnerProductSBP(1, dpress, dpress);
  }
  fin.close();
  return sensor;
}
#endif
// ======================================================================

double Quasi1DEuler::CalcSensor(const double & x_center,
                                const double & sigma,
                                const double & press_targ,
                                const double & reg_param,
                                const string & flow_file) {
  ifstream fin(flow_file.c_str(), ios::in | ios::binary);
  if (!fin.good()) {
    cerr << "Error in Quasi1DEuler::CalcSensor(): error opening "
         << flow_file << endl;
    throw(-1);
  }
  // check that num_nodes_ is consistent with value stored
  int nodes;
  fin.read(reinterpret_cast<char*>(&nodes), sizeof(int));
  if (nodes != num_nodes_) {
    cerr << "Error in Quasi1DEuler::CalcSensor(): "
         << "number of nodes does not match value in saved file." << endl;
    cout << "nodes = " << nodes << ": num_nodes_ = " << num_nodes_ << endl;
    throw(-1);
  }
  // get the number of iterations and total time
  int iter;
  fin.read(reinterpret_cast<char*>(&iter), sizeof(int));
  double total_time;
  fin.read(reinterpret_cast<char*>(&total_time), sizeof(double));
  double delta_t = total_time/static_cast<double>(iter);

  q_.BinaryRead(fin);
  InnerProdVector q_mid(3*num_nodes_,0.0);
  InnerProdVector dpress(num_nodes_, 0.0), kernel(num_nodes_, 0.0);

  // loop over all times, and calculate the sensor objective
  double sensor = 0.0;
  double sig2 = sigma*sigma;
  for (int n = 0; n < iter; n++) {
    q_old_ = q_;
    q_.BinaryRead(fin);
    q_mid.EqualsAXPlusBY(0.5, q_, 0.5, q_old_);
    CalcAuxiliaryVariables(q_mid);
    for (int i = 0; i < num_nodes_; i++) {
      double dx = x_coord_(i) - x_center;
      kernel(i) = exp(-dx*dx/sig2);
      dpress(i) = press_(i) - press_targ;
      dpress(i) = dpress(i)*dpress(i);
    }
    sensor += delta_t*0.5*dxi()*sbp_deriv_.InnerProductSBP(1, dpress, kernel);

    // include contributions due to regularization
    // 0.125 is from 0.5*(0.5*(src_(n) + src(n+1)))^2
    //sensor += reg_param*delta_t*0.125*(src_(n) + src_(n+1))*(src_(n) + src_(n+1));
  }
  //cout << "sensor (no reg.)  = " << sensor << endl;

  for (int n = 0; n < iter; n++) {
    sensor += reg_param*delta_t*0.5*src_(n)*src_(n);
  }
  //cout << "sensor (with reg.) = " << sensor << endl;

  fin.close();
  return sensor;
}

// ======================================================================

void Quasi1DEuler::CalcSensordJdQ(const double & x_center, const double & sigma,
                                  const double & press_targ,
                                  const double & reg_param,
                                  const string & flow_file,
                                  const string & dJdQ_file) {
  ifstream fin(flow_file.c_str(), ios::in | ios::binary);
  ofstream fout(dJdQ_file.c_str(), ios::out | ios::binary);
  if (!fin.good()) {
    cerr << "Error in Quasi1DEuler::CalcSensordJdQ(): error opening "
         << flow_file << endl;
    throw(-1);
  }
  // check that num_nodes_ is consistent with value stored
  int nodes;
  fin.read(reinterpret_cast<char*>(&nodes), sizeof(int));
  if (nodes != num_nodes_) {
    cerr << "Error in Quasi1DEuler::CalcSensordJdQ(): "
         << "number of nodes does not match value in saved file." << endl;
    cout << "nodes = " << nodes << ": num_nodes_ = " << num_nodes_ << endl;
    throw(-1);
  }
  // get the number of iterations and total time
  int iter;
  fin.read(reinterpret_cast<char*>(&iter), sizeof(int));
  double total_time;
  fin.read(reinterpret_cast<char*>(&total_time), sizeof(double));
  double delta_t = total_time/static_cast<double>(iter);

  // write header for dJdQ_file
  fout.write(reinterpret_cast<const char*>(&num_nodes_), sizeof(int));
  fout.write(reinterpret_cast<const char*>(&iter), sizeof(int));
  fout.write(reinterpret_cast<const char*>(&total_time), sizeof(double));

  q_.BinaryRead(fin);
  InnerProdVector q_mid(3*num_nodes_,0.0);
  InnerProdVector dJdQ(3*num_nodes_, 0.0), dJdQ_old(3*num_nodes_, 0.0);

  // loop over all times, and calculate the sensor objective derivative
  double sig2 = sigma*sigma;
  for (int n = 0; n < iter; n++) {
    q_old_ = q_;
    q_.BinaryRead(fin);
    q_mid.EqualsAXPlusBY(0.5, q_, 0.5, q_old_);
    CalcAuxiliaryVariables(q_mid);
    for (int i = 0; i < num_nodes_; i++) {
      int ptr = 3*i;
      double dx = x_coord_(i) - x_center;
      double dpress = delta_t*dxi()*exp(-dx*dx/sig2)*(kGamma-1.0)
          *(press_(i) - press_targ);
      double vel = q_mid(ptr+1)/q_mid(ptr);
      dJdQ(ptr)   =  dpress*0.5*vel*vel;
      dJdQ(ptr+1) = -dpress*vel;
      dJdQ(ptr+2) =  dpress;
    }
    sbp_deriv_.HTimesVector(3, dJdQ, dJdQ);
    dJdQ_old.EqualsAXPlusBY(0.5, dJdQ_old, 0.5, dJdQ);
    dJdQ_old.BinaryWrite(fout);
    dJdQ_old = dJdQ;
    dJdQ = 0.0;
  }
  dJdQ_old *= 0.5;
  dJdQ_old.BinaryWrite(fout);
  fout.close();
  fin.close();
}

// ======================================================================

void Quasi1DEuler::TestSensordJdQ(const double & x_center, const double & sigma,
                                  const double & press_targ,
                                  const double & reg_param,
                                  const string & flow_file) {
  // w will be a random vector at each time step
  InnerProdVector w(3*num_nodes_, 0.0);
  boost::random::mt19937 gen;
  boost::random::uniform_real_distribution<double> dist(-1.0, 1.0);
  double fd_eps = 1.E-7;

  // first, compute the analytical derivative
  CalcSensordJdQ(x_center, sigma, press_targ, reg_param, flow_file,
                 string("save_dJdQ.bin"));

  // evaluate the sensor at the unperturbed state
  double sensor = CalcSensor(x_center, sigma, press_targ, reg_param, flow_file);

  ifstream fin_q(flow_file.c_str(), ios::in | ios::binary);
  ifstream fin_dJdQ("save_dJdQ.bin", ios::in | ios::binary);
  ofstream fout("flow_pert.bin", ios::out | ios::binary);
  if ( (!fin_q.good()) || (!fin_dJdQ.good()) ) {
    cerr << "Error in Quasi1DEuler::TestSensordJdQ(): error opening "
         << flow_file << " or " << "save_dJdQ.bin" << endl;
    throw(-1);
  }
  // check that num_nodes_ is consistent with value stored
  int nodes;
  fin_q.read(reinterpret_cast<char*>(&nodes), sizeof(int));
  if (nodes != num_nodes_) {
    cerr << "Error in Quasi1DEuler::TestSensordJdQ(): "
         << "number of nodes does not match value in saved flow file." << endl;
    cout << "nodes = " << nodes << ": num_nodes_ = " << num_nodes_ << endl;
    throw(-1);
  }
  fin_dJdQ.read(reinterpret_cast<char*>(&nodes), sizeof(int));
  if (nodes != num_nodes_) {
    cerr << "Error in Quasi1DEuler::TestSensordJdQ(): "
         << "number of nodes does not match value in saved dJdQ file." << endl;
    cout << "nodes = " << nodes << ": num_nodes_ = " << num_nodes_ << endl;
    throw(-1);
  }

  // get the number of iterations and total time; check for consistency here too
  int iter;
  fin_q.read(reinterpret_cast<char*>(&iter), sizeof(int));
  double total_time;
  fin_q.read(reinterpret_cast<char*>(&total_time), sizeof(double));
  double delta_t = total_time/static_cast<double>(iter);
  int check_iter;
  double check_time;
  fin_dJdQ.read(reinterpret_cast<char*>(&check_iter), sizeof(int));
  fin_dJdQ.read(reinterpret_cast<char*>(&check_time), sizeof(double));
  if ( (check_iter != iter) || (check_time != total_time) ) {
    cerr << "Error in Quasi1DEuler::TestSensordJdQ(): "
         << "inconsistency between flow and dJdQ files (iter or time)." << endl;
    throw(-1);
  }

  // write header for dJdQ_file
  fout.write(reinterpret_cast<const char*>(&num_nodes_), sizeof(int));
  fout.write(reinterpret_cast<const char*>(&iter), sizeof(int));
  fout.write(reinterpret_cast<const char*>(&total_time), sizeof(double));

#if 0
  // loop over all times, and calculate the sensor objective derivative
  InnerProdVector dJdQ(3*num_nodes_, 0.0);
  double prod = 0.0;
  for (int n = 0; n <= iter; n++) {
    q_.BinaryRead(fin_q);
    dJdQ.BinaryRead(fin_dJdQ);
    // make a random vector
    for (int i = 0; i < 3*num_nodes_; i++) {
      w(i) = dist(gen);
    }
    // take product with dJdQ
    prod += InnerProd(dJdQ, w);
    // perturb flow
    q_ += fd_eps*w;
    q_.BinaryWrite(fout);
  }
#endif

  // loop over all times in reverse, and calculate the sensor derivative
  InnerProdVector dJdQ(3*num_nodes_, 0.0);
  double prod = 0.0;
  for (int n = iter; n >= 0; n--) {
    unsigned long fptr = 2*sizeof(int) + sizeof(double)
        + n*(3*num_nodes_*sizeof(double));
    q_.BinaryRead(fin_q, fptr);
    dJdQ.BinaryRead(fin_dJdQ, fptr);
    // make a random vector
    for (int i = 0; i < 3*num_nodes_; i++) {
      w(i) = dist(gen);
    }
    // take product with dJdQ
    prod += InnerProd(dJdQ, w);
    // perturb flow
    q_ += fd_eps*w;
    q_.BinaryWrite(fout, fptr);
  }

  fout.close();
  fin_q.close();
  fin_dJdQ.close();

  double sensor_pert = CalcSensor(x_center, sigma, press_targ, reg_param,
                                  string("flow_pert.bin"));
  double prod_fd = (sensor_pert - sensor)/fd_eps;
  cout << "TestSensordJdQ: prod = " << prod << ": prod_fd = " << prod_fd
       << endl;
  cout << "|prod - prod_fd| = " << fabs(prod - prod_fd) << endl;
}

// ======================================================================

void Quasi1DEuler::CalcGradient(const objective & obj,
                                Nozzle & nozzle_shape,
                                InnerProdVector & dJdX) {
  InnerProdVector dJdQ(3*num_nodes_, 0.0), dJdA(num_nodes_, 0.0);
  switch (obj) {
    case inverse: {// inverse design
      CalcInverseDesigndJdQ(dJdQ);
      break;
    }
    case total_kinetic: {// total kinetic energy
      CalcTotalEnergydJdQ(dJdQ);
      break;
    }
    default: {
      cerr << "Quasi1DEuler(CalcGradient): "
           << "invalid objective function"
           << endl;
      throw(-1);
      break;
    }
  }
  // "move" dJdQ to rhs, and solve adjoint system
  dJdQ *= -1.0;
  SolveAdjoint(100, 1.e-12, dJdQ, psi_);

#if 0
  // uncomment to output adjoint
  sbp_deriv_.HinvTimesVector(3, psi_, psi_);
  set_q(psi_);
  WriteTecplot(1.0, 1.0);
  throw(-1);
#endif

  // differentiate psi^R(q, Area) w.r.t. area, and contract with dArea/dX
  JacobianTransposedAreaProduct(psi_, dJdA);
  dJdX = 0.0;
  dJdX = nozzle_shape.AreaReverseDerivative(get_x_coord(), dJdA);
}

// ======================================================================

void Quasi1DEuler::CalcSensorGradient(const double & x_center,
                                      const double & sigma,
                                      const double & press_targ,
                                      const double & reg_param,
                                      const int & max_krylov,
                                      const string & flow_file,
                                      InnerProdVector & dJdX) {
  // first calculate the rhs for the adjoint equation
  CalcSensordJdQ(x_center, sigma, press_targ, reg_param, flow_file,
                 string("save_dJdQ.bin"));
  // next, solve for the adjoint variables
  SolveUnsteadyAdjoint(flow_file, max_krylov, 1.e-14, string("save_dJdQ.bin"),
                       string("save_adjoint.bin"));
  // contract adjoint with derivative of residual w.r.t. source
  dJdX = 0.0;
  UnsteadyJacTransSourceProduct(string("save_adjoint.bin"), dJdX);

  // add contributions due to regularization term
  AddSensordJdSource(reg_param, dJdX);
}

// ======================================================================

void Quasi1DEuler::AddSensordJdSource(const double & reg_param,
                                      InnerProdVector & dJdX) {
  int iter = src_.size();
  for (int n = 0; n < iter; n++) {
    double dsrc = reg_param*dt()*src_(n);
    dJdX(n) += dsrc;
  }
}

// ======================================================================

void Quasi1DEuler::UnsteadyJacobianSourceProduct(const InnerProdVector & u,
                                                 const string & prod_file) {
  ofstream fout(prod_file.c_str(), ios::out | ios::binary);
  // write header for product file
  int iter = src_.size();
  double total_time = iter*dt();
  fout.write(reinterpret_cast<const char*>(&num_nodes_), sizeof(int));
  fout.write(reinterpret_cast<const char*>(&iter), sizeof(int));
  fout.write(reinterpret_cast<const char*>(&total_time), sizeof(double));
  double fac = dt();
  if (u.size() != iter) {
    cerr << "Error in Quasi1DEuler::UnsteadyJacobianSourceProduct(): "
         << "number of iterations is inconsistent with size of u." << endl;
    throw(-1);
  }
  // the source kernel is constant for all time, so precompute
  InnerProdVector kernel(3*num_nodes_, 0.0), tmp(3*num_nodes_, 0.0);
  for (int i = 0; i < num_nodes_; i++) {
    double dx = x_coord_(i) - src_x_;
    kernel(3*i+1) = dt()*exp(-dx*dx/src_sig2_);
  }

  // initial residual is just R = u - u(x,0), so no source dependence
  tmp.BinaryWrite(fout);

  // loop over time iterations
  for (int n = 0; n < iter; n++) {
    tmp = kernel;
    tmp *= -u(n); // neg sign moves source to left side
    tmp.BinaryWrite(fout);
  }
  fout.close();
}

// ======================================================================

void Quasi1DEuler::UnsteadyJacTransSourceProduct(const string & psi_file,
                                                 InnerProdVector & dJdX) {
  // open the solution file, the dJdQ file, and adjoint file
  ifstream fin(psi_file.c_str(), ios::in | ios::binary);
  if (!fin.good()) {
    cerr << "Error in Quasi1DEuler::UnsteadyJacTransSourceProduct(): "
         << "error opening " << psi_file << endl;
    throw(-1);
  }
  // check that num_nodes_ is consistent with value stored
  int nodes;
  fin.read(reinterpret_cast<char*>(&nodes), sizeof(int));
  if (nodes != num_nodes_) {
    cerr << "Error in Quasi1DEuler::UnsteadyJacTransSourceProduct(): "
         << "number of nodes does not match value in saved psi file." << endl;
    cout << "nodes = " << nodes << ": num_nodes_ = " << num_nodes_ << endl;
    throw(-1);
  }
  // get the number of iterations and total time
  int iter;
  fin.read(reinterpret_cast<char*>(&iter), sizeof(int));
  double total_time;
  fin.read(reinterpret_cast<char*>(&total_time), sizeof(double));
  dt_ = total_time/static_cast<double>(iter);
  double fac = dt();
  if (dJdX.size() != iter) {
    cerr << "Error in Quasi1DEuler::UnsteadyJacTransSourceProduct(): "
         << "number of iterations is inconsistent with size of dJdX." << endl;
    throw(-1);
  }

  // first adjoint corresponds to R = u - u(0,x), which has no source dependence
  // so discard this adjoint
  psi_.BinaryRead(fin);

  // can go in forward direction here
  for (int n = 0; n < iter; n++) {
    psi_.BinaryRead(fin);
    for (int i = 0; i < num_nodes_; i++) {
      double dx = x_coord_(i) - src_x_;
      dJdX(n) -= psi_(3*i+1)*fac*exp(-dx*dx/src_sig2_);
    }
  }
  fin.close();
}

// ======================================================================

void Quasi1DEuler::TestSensorGradient() {

  const double length = 1.0;
  const double x_min = 0.0;
  const double x_max = x_min + length;
  const double T = 2.0;
  const int iter = 200; //200;

  // set the initial conditions
  InnerProdVector q_init(3*num_nodes_, 1.0);
  const double x_c = 0.3*length;
  for (int i = 0; i < num_nodes_; i++) {
    q_init(3*i) = 1.0 + 0.01*exp(-500.0*(x_coord_(i)-x_c)*(x_coord_(i)-x_c));
  }
  InitialCondition(q_init);

  SolveUnsteady(iter, T, 1.e-14, true);
  const string flow_file("save_flow.bin");

  // get a random vector to contract with
  InnerProdVector w(iter, 0.0);
  boost::random::mt19937 gen;
  boost::random::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (int n = 0; n < iter; n++)
    w(n) = dist(gen);

  // compute the analytical gradient
  double x_center = x_min + 0.8*length;
  double sigma = 0.05*length;
  double press_targ = 0.2;
  double reg_param = 1.0e-6;
  int max_krylov = 30;
  InnerProdVector dJdX(iter, 0.0);
  CalcSensorGradient(x_center, sigma, press_targ, reg_param, max_krylov,
                     flow_file, dJdX);

  // compute "exact" inner product with w
  double prod = InnerProd(dJdX, w);

  // compute the gradient using a finite-difference approximation
  double sensor = CalcSensor(x_center, sigma, press_targ, reg_param, flow_file);
  double fd_eps = 1.E-6;
  src_ += fd_eps*w;

  InitialCondition(q_init);
  SolveUnsteady(iter, T, 1.e-14, true);
  double sensor_pert = CalcSensor(x_center, sigma, press_targ, reg_param,
                                  flow_file);
  double prod_fd = (sensor_pert - sensor)/fd_eps;

  cout << "Quasi1DEuler::TestSensorGradient(): " << endl;
  cout << "prod (analytical) = " << prod << endl;
  cout << "prod (FD approx.) = " << prod_fd << endl;
  cout << "|prod difference| = " << fabs(prod - prod_fd) << endl;
  cout << "relative error    = " << fabs(prod - prod_fd)/fabs(prod) << endl;
}

// ======================================================================

double Quasi1DEuler::EstimateGradientError(const norm_type & norm,
    const objective & obj, Nozzle & nozzle_shape,
    const InnerProdVector & dJdX) {
  unsigned int num_design = dJdX.size();
  InnerProdVector rhs(3*num_nodes_, 0.0), w(3*num_nodes_, 0.0),
      lambda(3*num_nodes_, 0.0), dNormdGrad(num_design, 0.0);

  double grad_norm, grad_fac;
  int index_grad_inf = -1;
  switch (norm) {
    case L2: {
      grad_norm = dJdX.Norm2();
      dNormdGrad = dJdX;
      grad_fac = 1.0/grad_norm;
      break;
    }
    case inf: {
      grad_norm = ublas::norm_inf(dJdX);
      index_grad_inf = ublas::index_norm_inf(dJdX);
      dNormdGrad(index_grad_inf) = 1.0;
      grad_fac = kona::sign(1.0, dJdX(index_grad_inf));
      break;
    }
    default: {
      cerr << "Quasi1DEuler(EstimateGradientError): invalid norm."
           << endl;
      throw(-1);
      break;
    }
  }

  // Step 1: solve the first of two adjoint problems
  InnerProdVector dArea(num_nodes_, 0.0);
  dArea = nozzle_shape.AreaForwardDerivative(get_x_coord(), dNormdGrad);
  JacobianAreaProduct(dArea, rhs);
  rhs *= -grad_fac;
  SolveLinearized(100, 1.e-12, rhs, w);

#if 0
  // uncomment to output the adjoint w
  set_q(w);
  WriteTecplot(1.0, 1.0);
  throw(-1);
#endif

  // Step 2: solve the second adjoint problem

  // Note: for the term \partial^{2} L /\partial u \partial x,
  // we contract with dArea, and make use of the fact that the residual
  // depends linearly on the area (so we set area_ = dArea)
  InnerProdVector save_area(area_);
  set_area(dArea);
  double save_diss_coeff = diss_coeff_;
  diss_coeff_ = 0.0;
  JacobianTransposedStateProduct(psi_, rhs);
  set_area(save_area);
  diss_coeff_ = save_diss_coeff;
  rhs *= -grad_fac;

  InnerProdVector work_vec(3*num_nodes_, 0.0);
  switch (obj) {
    case inverse: {// inverse design
      CalcInverseDesignd2JdQ2(w, work_vec);
      break;
    }
    case total_kinetic: {// total kinetic energy
      CalcTotalEnergyd2JdQ2(w, work_vec);
      break;
    }
    default: {
      cerr << "Quasi1DEuler(EstimateGradientError): "
           << "invalid objective function"
           << endl;
      throw(-1);
      break;
    }
  }
  rhs -= work_vec;

  ResidualHessianProduct(psi_, w, work_vec);
  rhs -= work_vec;
  SolveAdjoint(100, 1.e-12, rhs, lambda);
#if 0
  // uncomment to output the adjoint lambda
  sbp_deriv_.HinvTimesVector(3, lambda, lambda);
  set_q(lambda);
  WriteTecplot(1.0, 1.0);
  throw(-1);
#endif

  // evaluate the gradient using higher-order SBP operator

  // remove norm from adjoints that left-multiplied residuals
  sbp_deriv_.HinvTimesVector(3, psi_, psi_);
  sbp_deriv_.HinvTimesVector(3, lambda, lambda);

  InnerProdVector dJdA(num_nodes_, 0.0), dJdX_ho(num_design, 0.0);
  int save_order = sbp_deriv_.order();
  int save_diss_order = sbp_diss_.order();
  sbp_deriv_.Define(num_nodes_, save_order + 1);
  sbp_diss_.Define(num_nodes_, save_order + 1, save_diss_order + 1);

  // apply new norm to adjoints that left-multiplied residuals
  sbp_deriv_.HTimesVector(3, psi_, psi_);
  sbp_deriv_.HTimesVector(3, lambda, lambda);

  sbp_deriv_.Apply(1, x_coord_, met_jac_); // compute HO metric Jac
  JacobianTransposedAreaProduct(psi_, dJdA);
  dJdX_ho = 0.0;
  dJdX_ho = nozzle_shape.AreaReverseDerivative(get_x_coord(), dJdA);

  double grad_err;
  switch (norm) {
    case L2: {
      grad_err = grad_norm - dJdX_ho.Norm2();
      break;
    }
    case inf: {
      grad_err = grad_norm - norm_inf(dJdX_ho);
      break;
    }
    default: {
      cerr << "Quasi1DEuler(EstimateGradientError): invalid norm."
           << endl;
      throw(-1);
      break;
    }
  }
  cout << "grad_norm - grad_norm_ho = " << grad_err << endl;

  // evaluate the residual-error terms

  CalcResidual();
  grad_err -= InnerProd(res_, lambda);
  cout << "residual correction: lambda*Res(q) = "
       << InnerProd(res_, lambda) << endl;

  EvaluateAdjointResidual(obj, res_);
  grad_err -= InnerProd(res_, w);
  cout << "residual correction: w*(A^T psi_ + dJdQ) = "
       << InnerProd(res_, w) << endl;

  // return SBP operators and metrics to lower order
  sbp_deriv_.Define(num_nodes_, save_order);
  sbp_diss_.Define(num_nodes_, save_order, save_diss_order);
  sbp_deriv_.Apply(1, x_coord_, met_jac_);

  return grad_err;
}

// ======================================================================

double Quasi1DEuler::EstimateGradientErrorFD(const norm_type & norm,
    const objective & obj, Nozzle & nozzle_shape,
    const InnerProdVector & dJdX) {
  unsigned int num_design = dJdX.size();
  InnerProdVector rhs(3*num_nodes_, 0.0), w(3*num_nodes_, 0.0),
      lambda(3*num_nodes_, 0.0), dNormdGrad(num_design, 0.0);

  double grad_norm, grad_fac;
  int index_grad_inf = -1;
  switch (norm) {
    case L2: {
      grad_norm = dJdX.Norm2();
      dNormdGrad = dJdX;
      grad_fac = 1.0/grad_norm;
      break;
    }
    case inf: {
      grad_norm = ublas::norm_inf(dJdX);
      index_grad_inf = ublas::index_norm_inf(dJdX);
      dNormdGrad(index_grad_inf) = 1.0;
      grad_fac = kona::sign(1.0, dJdX(index_grad_inf));
      break;
    }
    default: {
      cerr << "Quasi1DEuler(EstimateGradientErrorFD): invalid norm."
           << endl;
      throw(-1);
      break;
    }
  }

  // Step 1: solve the first of two adjoint problems
  InnerProdVector dArea(num_nodes_, 0.0);
  dArea = nozzle_shape.AreaForwardDerivative(get_x_coord(), dNormdGrad);
  JacobianAreaProduct(dArea, rhs);
  rhs *= -grad_fac;
  SolveLinearized(100, 1.e-12, rhs, w);

#if 0
  // uncomment to output the adjoint w
  set_q(w);
  WriteTecplot(1.0, 1.0);
  throw(-1);
#endif

  // Step 2: solve the second adjoint problem

  // Compute the first part of the rhs using finite-differences
  EvaluateAdjointResidual(obj, res_);
  InnerProdVector init_res(res_);

  // perturb the design variables and re-evaluate residual
  InnerProdVector save_design(num_design, 0.0),
      pert_design(num_design, 0.0);
  nozzle_shape.GetCoeff(save_design);
  double eps = kona::CalcEpsilon(save_design.Norm2(),
                                 grad_fac*dNormdGrad.Norm2());
  cout << "eps for FD is = " << eps << endl;

  pert_design = dNormdGrad;
  pert_design *= (eps*grad_fac);
  pert_design += save_design;
  nozzle_shape.SetCoeff(pert_design);
  set_area(nozzle_shape.Area(get_x_coord()));
  EvaluateAdjointResidual(obj, rhs);
  rhs -= res_;
  rhs /= -eps;

  // unperturb design
  nozzle_shape.SetCoeff(save_design);
  set_area(nozzle_shape.Area(get_x_coord()));

  // Compute the second part of the rhs using finite-differences
  InnerProdVector save_q(q_);
  eps = kona::CalcEpsilon(save_q.Norm2(), w.Norm2());
  q_.EqualsAXPlusBY(1.0, save_q, eps, w);
  EvaluateAdjointResidual(obj, res_);
  res_ -= init_res;
  res_ /= -eps;
  rhs += res_;

  SolveAdjoint(100, 1.e-12, rhs, lambda);
#if 0
  // uncomment to output the adjoint lambda
  sbp_deriv_.HinvTimesVector(3, lambda, lambda);
  set_q(lambda);
  WriteTecplot(1.0, 1.0);
  throw(-1);
#endif

  // evaluate the gradient using higher-order SBP operator

  // remove norm from adjoints that left-multiplied residuals
  sbp_deriv_.HinvTimesVector(3, psi_, psi_);
  sbp_deriv_.HinvTimesVector(3, lambda, lambda);

  InnerProdVector dJdA(num_nodes_, 0.0), dJdX_ho(num_design, 0.0);
  int save_order = sbp_deriv_.order();
  int save_diss_order = sbp_diss_.order();
  sbp_deriv_.Define(num_nodes_, save_order + 1);
  sbp_diss_.Define(num_nodes_, save_order + 1, save_diss_order + 1);

  // apply new norm to adjoints that left-multiplied residuals
  sbp_deriv_.HTimesVector(3, psi_, psi_);
  sbp_deriv_.HTimesVector(3, lambda, lambda);

  sbp_deriv_.Apply(1, x_coord_, met_jac_); // compute HO metric Jac
  JacobianTransposedAreaProduct(psi_, dJdA);
  dJdX_ho = 0.0;
  dJdX_ho = nozzle_shape.AreaReverseDerivative(get_x_coord(), dJdA);

  double grad_err;
  switch (norm) {
    case L2: {
      grad_err = grad_norm - dJdX_ho.Norm2();
      break;
    }
    case inf: {
      grad_err = grad_norm - norm_inf(dJdX_ho);
      break;
    }
    default: {
      cerr << "Quasi1DEuler(EstimateGradientErrorFD): invalid norm."
           << endl;
      throw(-1);
      break;
    }
  }
  cout << "grad_norm - grad_norm_ho = " << grad_err << endl;

  // evaluate the residual-error terms

  CalcResidual();
  grad_err -= InnerProd(res_, lambda);
  cout << "residual correction: lambda*Res(q) = "
       << InnerProd(res_, lambda) << endl;

  InnerProdVector dJdQ(3*num_nodes_, 0.0);
  switch (obj) {
    case inverse: {// inverse design
      CalcInverseDesigndJdQ(dJdQ);
      break;
    }
    case total_kinetic: {// total kinetic energy
      CalcTotalEnergydJdQ(dJdQ);
      break;
    }
    default: {
      cerr << "Quasi1DEuler(EstimateGradientErrorFD): "
           << "invalid objective function"
           << endl;
      throw(-1);
      break;
    }
  }
  // "move" dJdQ to rhs, and solve adjoint system
  //dJdQ *= -1.0;

  JacobianTransposedStateProduct(psi_, res_);
  res_ += dJdQ;
  grad_err -= InnerProd(res_, w);
  cout << "residual correction: w*(A^T psi_ + dJdQ) = "
       << InnerProd(res_, w) << endl;

  // return SBP operators and metrics to lower order
  sbp_deriv_.Define(num_nodes_, save_order);
  sbp_diss_.Define(num_nodes_, save_order, save_diss_order);
  sbp_deriv_.Apply(1, x_coord_, met_jac_);

  return grad_err;
}

// ======================================================================

void Quasi1DEuler::EvaluateAdjointResidual(const objective & obj,
                                           InnerProdVector & adj_res) {
  InnerProdVector dJdQ(3*num_nodes_, 0.0);
  CalcAuxiliaryVariables(q_);
  switch (obj) {
    case inverse: {// inverse design
      CalcInverseDesigndJdQ(dJdQ);
      break;
    }
    case total_kinetic: {// total kinetic energy
      CalcTotalEnergydJdQ(dJdQ);
      break;
    }
    default: {
      cerr << "Quasi1DEuler(EvaluateAdjointResidual): "
           << "invalid objective function"
           << endl;
      throw(-1);
      break;
    }
  }
  JacobianTransposedStateProduct(psi_, adj_res);
  adj_res += dJdQ;
}

// ======================================================================

void JacobianVectorProduct::operator()(const InnerProdVector & u,
                                       InnerProdVector & v) {
  solver->JacobianStateProduct(u, v);
}

// ======================================================================

void ApproxJacobian::operator()(InnerProdVector & u,
                                InnerProdVector & v) {
  solver->Precondition(u, v);
}

// ======================================================================

void JacobianTransposedVectorProduct::operator()(
    const InnerProdVector & u, InnerProdVector & v) {
  solver->JacobianTransposedStateProduct(u, v);
}

// ======================================================================

void ApproxJacobianTransposed::operator()(InnerProdVector & u,
                                          InnerProdVector & v) {
  solver->PreconditionTransposed(u, v);
}

// ======================================================================

void UnsteadyJacobianVectorProduct::operator()(const InnerProdVector & u,
                                               InnerProdVector & v) {
  solver->UnsteadyJacobianStateProduct(u, v);
}

// ======================================================================

void UnsteadyJacTransVectorProduct::operator()(const InnerProdVector & u,
                                               InnerProdVector & v) {
  solver->UnsteadyJacTransStateProduct(u, v, true);
}

// ======================================================================

void Quasi1DEuler::CalcEulerFlux(const InnerProdVector & q_var,
                                 InnerProdVector & flux) {
  for (int i = 0; i < num_nodes_; i++) {
    double area = area_(i);
    double press = press_(i);
    int ptr = 3*i;
    double rho = q_var(ptr);
    double u = q_var(ptr+1) / rho;
    double e = q_var(ptr+2);

    flux(ptr+0) = rho * u * area;
    flux(ptr+1) = (rho * u * u + press) * area;
    flux(ptr+2) = u * ( e + press) * area;
  }
}

// ======================================================================

void Quasi1DEuler::CalcAuxiliaryVariables(const InnerProdVector & q_var) {
  for (int i = 0; i < num_nodes_; i++) {
    int ptr = 3*i;
    double rho = q_var(ptr);
    double u = q_var(ptr+1)/rho;
    double e = q_var(ptr+2);
    press_(i) = (kGamma-1.0)*(e - 0.5*rho*u*u);
    sndsp_(i) = sqrt(kGamma*press_(i)/rho);
    spect_(i) = area_(i)*(fabs(u) + sndsp_(i));
  }
}

// ======================================================================

void Quasi1DEuler::CalcFluxJacobian(
    const double & area,
    const ublas::vector_range<ublas::vector<double> > q,
    ublas::matrix<double> & flux_jac) const {
  // calculate primative variables
  double rho = q(0);
  double u = q(1) / rho;
  double e = q(2);
  // first row
  flux_jac(0,0) = 0.0;
  flux_jac(0,1) = area;
  flux_jac(0,2) = 0.0;
  // second row
  flux_jac(1,0) = 0.5*(kGamma-3.0)*u*u*area;
  flux_jac(1,1) = (3.0-kGamma)*u*area;
  flux_jac(1,2) = (kGamma-1.0)*area;
  // third row
  flux_jac(2,0) = u*((kGamma-1.0)*u*u - kGamma*e/rho)*area;
  flux_jac(2,1) = (kGamma*e/rho - 1.5*(kGamma-1.0)*u*u)*area;
  flux_jac(2,2) = kGamma*u*area;
}

// ======================================================================

void Quasi1DEuler::CalcFluxHessianProduct(
    const double & area,
    const ublas::vector_range<ublas::vector<double> > q,
    const ublas::vector_range<ublas::vector<double> > w,
    ublas::matrix<double> & flux_hess) const {
  // calculate primative variables
  double rho = q(0);
  double u = q(1) / rho;
  double e = q(2);
  const double kGami = (kGamma - 1.0);
  double w1 = w(0);
  double w2 = w(1);
  double w3 = w(2);
  // first row
  flux_hess(0,0) = 0.0;
  flux_hess(0,1) = 0.0;
  flux_hess(0,2) = 0.0;
  // second row
  flux_hess(1,0) = -(kGamma - 3.0)*u*u*w1/rho - (3.0 - kGamma)*u*w2/rho;
  flux_hess(1,1) = (kGamma - 3.0)*u*w1/rho + (3.0 - kGamma)*w2/rho;
  flux_hess(1,2) = 0.0;
  // third row
  flux_hess(2,0) = -2.0*kGami*u*u*(u*w1 - 1.5*w2)/rho
      - kGami*u*u*u*w1/rho - kGamma*e*(w2 - u*w1)/(rho*rho)
      + kGamma*e*u*w1/(rho*rho) - kGamma*u*w3/rho;
  flux_hess(2,1) = 2.0*kGami*u*(u*w1 - 1.5*w2)/rho + kGami*u*u*w1/rho
      - kGamma*e*w1/(rho*rho) + kGamma*w3/rho;
  flux_hess(2,2) = kGamma*(w2 - u*w1)/rho;
  flux_hess *= area;
}

// ======================================================================

void Quasi1DEuler::CalcFluxHessianProductHD(
    const double & area,
    const ublas::vector_range<ublas::vector<double> > q,
    const ublas::vector_range<ublas::vector<double> > w,
    ublas::matrix<double> & flux_hess) const {
  ublas::bounded_vector<HyperDual,3> q_hd, flux_hd;
  HyperDual area_hd;
  for (int i = 0; i < 3; i++) {
    q_hd(i).setvalues(q(i), 0.0, 0.0, 0.0);
  }
  area_hd.setvalues(area, 0.0, 0.0, 0.0);
  flux_hess = ublas::zero_matrix<double>(3,3);
  // loop over the first variable to differentiate w.r.t
  for (int i = 0; i < 3; i++) {
    q_hd(i).ipart() = 1.0;
    // loop over the first variable to differentiate w.r.t
    for (int j = 0; j < 3; j++) {
      q_hd(j).jpart() = 1.0;
      CalcFlux(area_hd, q_hd, flux_hd);
      for (int k = 0; k < 3; k++)
        flux_hess(k,i) += flux_hd(k).ijpart()*w(j);
      q_hd(j).jpart() = 0.0;
    }
    q_hd(i).ipart() = 0.0;
  }
}

// ======================================================================
#if 0
void Quasi1DEuler::CalcSAT(const InnerProdVector & bc,
                           const double & area, const double & sgn,
                           ublas::vector_range<ublas::vector<double> > q,
                           ublas::bounded_vector<double,3> & sat) const {
  // calculate primative variables
  double rho = q(0);
  double u = q(1)/rho;
  double E = q(2);

  double p = (kGamma - 1.0)*(E - 0.5*rho*u*u);
  double a = sqrt(kGamma*p/rho);
  double H = (E + p)/rho;
  double phi = 0.5*u*u;

  // calculate the wave speeds
  double lam1 = 0.5*area*(fabs(u + a) + sgn*(u + a));
  double lam2 = 0.5*area*(fabs(u - a) + sgn*(u - a));
  double lam3 = 0.5*area*(fabs(u) + sgn*(u));

#if 0
  double spec = fabs(u) + a;
  double lam1 = 0.5*area*(max(fabs(u + a),kVn*spec) + sgn*(u + a));
  double lam2 = 0.5*area*(max(fabs(u - a),kVn*spec) + sgn*(u - a));
  double lam3 = 0.5*area*(max(fabs(u),kVl*spec) + sgn*(u));
#endif

  // calculate the differences
  double dq1 = q(0) - bc(0);
  double dq2 = q(1) - bc(1);
  double dq3 = q(2) - bc(2);

  sat(0) = lam3*dq1;
  sat(1) = lam3*dq2;
  sat(2) = lam3*dq3;

  // temporary vectors
  ublas::bounded_vector<double, 3> E1dq, E2dq;

  // get E1 times dq
  E1dq(0) = phi*dq1 - u*dq2 + dq3;
  E1dq(1) = E1dq(0)*u;
  E1dq(2) = E1dq(0)*H;

  // get E2 times dq
  E2dq(0) = 0.0;
  E2dq(1) = -u*dq1 + dq2;
  E2dq(2) = E2dq(1)*u;

  // add to pen
  double tmp1 = 0.5*(lam1 + lam2) - lam3;
  double tmp2 = (kGamma-1.0)/(a*a);
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
  sat += tmp1*(E1dq + (kGamma-1.0)*E2dq);
}
#endif
// ======================================================================
#if 0
void Quasi1DEuler::CalcSAT_Complex(
    const ublas::bounded_vector<complex,3> & bc,
    const complex & area, const complex & sgn,
    const ublas::bounded_vector<complex,3> & q,
    ublas::bounded_vector<complex,3> & sat) const {
  // calculate primative variables
  complex rho = q(0);
  complex u = q(1)/rho;
  complex E = q(2);

  complex p = (kGamma - 1.0)*(E - 0.5*rho*u*u);
  complex a = sqrt(kGamma*p/rho);
  complex H = (E + p)/rho;
  complex phi = 0.5*u*u;

  // calculate the wave speeds
  complex lam1 = 0.5*area*(fabs(u + a) + sgn*(u + a));
  complex lam2 = 0.5*area*(fabs(u - a) + sgn*(u - a));
  complex lam3 = 0.5*area*(fabs(u) + sgn*(u));

#if 0
  complex spec = fabs(u) + a;
  complex lam1 = 0.5*area*(max(fabs(u + a),kVn*spec) + sgn*(u + a));
  complex lam2 = 0.5*area*(max(fabs(u - a),kVn*spec) + sgn*(u - a));
  complex lam3 = 0.5*area*(max(fabs(u),kVl*spec) + sgn*(u));
#endif

  // calculate the differences
  complex dq1 = q(0) - bc(0);
  complex dq2 = q(1) - bc(1);
  complex dq3 = q(2) - bc(2);

  sat(0) = lam3*dq1;
  sat(1) = lam3*dq2;
  sat(2) = lam3*dq3;

  // temporary vectors
  ublas::bounded_vector<complex, 3> E1dq, E2dq;

  // get E1 times dq
  E1dq(0) = phi*dq1 - u*dq2 + dq3;
  E1dq(1) = E1dq(0)*u;
  E1dq(2) = E1dq(0)*H;

  // get E2 times dq
  E2dq(0) = complex(0.0, 0.0);
  E2dq(1) = -u*dq1 + dq2;
  E2dq(2) = E2dq(1)*u;

  // add to pen
  complex tmp1 = 0.5*(lam1 + lam2) - lam3;
  complex tmp2 = (kGamma-1.0)/(a*a);
  sat += tmp1*(tmp2*E1dq + E2dq);

  // get E3 times dq
  E1dq(0) = -u*dq1 + dq2;
  E1dq(1) = E1dq(0)*u;
  E1dq(2) = E1dq(0)*H;

  // get E4 times dq
  E2dq(0) = complex(0.0, 0.0);
  E2dq(1) = phi*dq1 - u*dq2 + dq3;
  E2dq(2) = E2dq(1)*u;

  // add to sat
  tmp1 = 0.5*(lam1 - lam2)/a;
  sat += tmp1*(E1dq + (kGamma-1.0)*E2dq);
}
#endif
