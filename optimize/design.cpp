/**
 * \file design.cpp
 * \brief driver for the unit test for the Quasi1DEuler optimization
 * \version 1.0
 */

#include <ostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#include <user_memory.hpp>
#include <kona.hpp>
#include <boost/math/constants/constants.hpp>

#include "../inner_prod_vector.hpp"
#include "../exact_solution.hpp"
#include "../nozzle.hpp"
#include "../quasi_1d_euler.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::ofstream;

const double pi = boost::math::constants::pi<double>();

// global variables...but only to this file
static int num_design_vec = -1; // number of design vectors
static int num_state_vec = -1;  // number of state vectors
static vector<InnerProdVector> design;    // design vectors
static vector<InnerProdVector> state;     // state vectors
static InnerProdVector press_targ;        // target pressure
 
// domain parameters
static const double length = 1.0;
static const double x_min = 0.0;
static const double x_max = x_min + length;

// design parameters
static int num_design;
static const double area_left = 1.5; //2.0;
static const double area_right = 1.25; //1.5;

// discretization parameters
static const int nodes = 201;
static const int num_var = 3*nodes;
static const int order = 3;
static const bool sbp_quad = true; // use SBP norm for quadrature

// solution parameters
static const double kAreaStar = 0.8; // for the exact solution
static const bool subsonic = true; // for the exact solution
static const double kTempStag = 300.0;
static const double kPressStag = 100000;
static const double kRGas = 287.0;
static double rho_R, rho_u_R, e_R;

static Quasi1DEuler solver(nodes, order);
static BsplineNozzle nozzle_shape;

double MeshCoord(const double & length, const int & num_nodes,
                 const int & i);

double TargetNozzleArea(const double & x);

double InitNozzleArea(const double & x);

int userFunc(int request, int leniwrk, int *iwrk, int lendwrk,
	     double *dwrk);

static vector<double> path;
void WritePath();

// ======================================================================

int main(int argc, char *argv[]) {

  if (argc != 2) {
    cerr << "Error in design: expect exactly one command line variable "
         << "(number of B-spline control points defining nozzle area)"
         << endl;
  } else {
    num_design = atoi(argv[1]);
    cout << "Running design with " << num_design << " design vars." << endl;
  }

  // define the x-coordinates
  InnerProdVector x_coord(nodes, 0.0);
  for (int i = 0; i < nodes; i++)
    x_coord(i) = MeshCoord(length, nodes, i);
  
  // define the target nozzle shape and set area 
  // NOTE: the left and right nozzle areas are fixed and are not
  // included in the number of design variables
  InnerProdVector area_targ(nodes, 0.0);
  for (int i = 0; i < nodes; i++)
    area_targ(i) = TargetNozzleArea(x_coord(i));
  
  nozzle_shape.SetAreaAtEnds(area_left, area_right);
  
  // define reference and boundary conditions
  double rho, rho_u, e;
  CalcFlowExact(kGamma, kRGas, kAreaStar, area_left, true, 
                kTempStag, kPressStag, rho, rho_u, e);
  double rho_ref = rho;
  double press = (kGamma - 1.0)*(e - 0.5*rho_u*rho_u/rho);
  double a_ref = sqrt(kGamma*press/rho_ref);
  double rho_L = 1.0;
  double rho_u_L = rho_u/(a_ref*rho_ref);
  double e_L = e/(rho_ref*a_ref*a_ref);
  CalcFlowExact(kGamma, kRGas, kAreaStar, area_right, true,
                kTempStag, kPressStag, rho, rho_u, e);
  rho_R = rho/rho_ref;
  rho_u_R = rho_u/(a_ref*rho_ref);
  e_R =  e/(rho_ref*a_ref*a_ref);  

  // define the target pressure
  press_targ.resize(nodes);
  for (int i = 0; i < nodes; i++) {
    CalcFlowExact(kGamma, kRGas, kAreaStar, area_targ(i), true,
                  kTempStag, kPressStag, rho, rho_u, e);
    rho /= rho_ref;
    rho_u /= (a_ref*rho_ref);
    e /= (rho_ref*a_ref*a_ref);
    press_targ(i) = (kGamma - 1.0)*(e - 0.5*rho_u*rho_u/rho);
  }

  // create the solver
  solver.set_x_coord(x_coord);
  solver.set_area(area_targ);

  // set boundary and initial conditions
  solver.set_bc_left(rho_L, rho_u_L, e_L);
  solver.InitialCondition(rho_R, rho_u_R, e_R);
  solver.set_bc_right(rho_R, rho_u_R, e_R);
  
  // define any discretization and solver paramters
  solver.set_diss_coeff(0.04);

  // find the discrete target pressure
  solver.NewtonKrylov(100, 1.e-10);
  solver.set_press_targ(solver.get_press());
  int solver_precond_calls = solver.TotalPreconditionerCalls();
  solver.ClearPreconditionerCallCounts();

  map<string,string> optns;
  //optns["inner.lambda_init"] = "0.8";
  KonaOptimize(userFunc, optns);
  int kona_precond_calls = solver.TotalPreconditionerCalls();
    
  cout << "Total solver precond calls: " << solver_precond_calls << endl;
  cout << "Total Kona precond calls:   " << kona_precond_calls << endl;
  cout << "Ratio of precond calls:     " 
       << ( static_cast<double>(kona_precond_calls)
            /static_cast<double>(solver_precond_calls) ) << endl;


#if 0
  // solve for the flow with original coeff;
  //coeff(0) = 0.45;
  //nozzle_shape.SetCoeff(coeff);
  //solver.set_area(nozzle_shape.Area(solver.get_x_coord()));
  solver.NewtonKrylov(100, 1.e-14);
  int cold_precond_calls = solver.TotalPreconditionerCalls();

  // solve for the flow with perturbed coeff
  for (int j = 0; j < coeff.size(); j++) 
    coeff(j) += 1.e-7;
  nozzle_shape.SetCoeff(coeff);
  solver.set_area(nozzle_shape.Area(solver.get_x_coord()));
  solver.ClearPreconditionerCallCounts();
  solver.NewtonKrylov(100, 1.e-14);
  int warm_precond_calls = solver.TotalPreconditionerCalls();
  cout << "Cold solver precond calls: " << cold_precond_calls << endl;
  cout << "Warm solver precond calls: " << warm_precond_calls << endl;
  cout << "Ratio of precond calls:     " 
       << ( static_cast<double>(warm_precond_calls)
            /static_cast<double>(cold_precond_calls) ) << endl;
#endif

#if 0
  // solve for the adjoint variables
  InnerProdVector dJdQ(num_var, 0.0), psi(num_var, 0.0);
  solver.CalcInverseDesigndJdQ(dJdQ);
  dJdQ *= -1.0;
  solver.SolveAdjoint(100, 1.e-2, dJdQ, psi);
#endif

  // output the optimized nozzle and flow
  //solver.WriteTecplot(rho_ref, a_ref);
  solver.WriteTecplot(1.0, 1.0, "flow_opt.dat");
  double error_L2, error_inf;
  solver.CalcMachError(kAreaStar, subsonic, error_L2,
                       error_inf);
  cout << "error L2 = " << error_L2 << endl;
  cout << "error_inf = " << error_inf << endl;

  WritePath();
}

// ======================================================================

double MeshCoord(const double & length, const int & num_nodes,
                 const int & i) {
  double xi = static_cast<double>(i)/static_cast<double>(num_nodes-1);
  // uniform spacing
  return length*xi;
  // simple exponential mesh spacing
  //return (exp(4.0*xi)-1.0)/(exp(4.0)-1.0);
  // mesh ponts clustered near center
  //return (xi + (cos(pi*xi) - 1.0)/pi)/(1.0 -2.0/pi);
}

// ======================================================================

double TargetNozzleArea(const double & x) {
  // cubic polynomial nozzle
  const double area_mid = 1.0;
  double a = area_left;
  double b = 4.0*area_mid - 5.0*area_left + area_right;
  double c = -4.0*(area_right -2.0*area_left + area_mid);
  double d = 4.0*(area_right - area_left);
  return a + x*(b + x*(c + x*d));
}

// ======================================================================

double InitNozzleArea(const double & x) {
  // linear nozzle
  return area_left + (area_right - area_left)*x;
}

// ======================================================================

int userFunc(int request, int leniwrk, int *iwrk, int lendwrk,
	     double *dwrk) {
  static int opt_iter = 0;
  switch (request) {
    case kona::allocmem: {// allocate iwrk[0] design vectors and
      // iwrk[1] state vectors
      if (num_design_vec >= 0) { // free design memory first
        if (design.size() == 0) {
          cerr << "userFunc: "
               << "design array is empty but num_design_vec > 0" << endl;
          throw(-1);
        }
        design.clear();
      }        
      if (num_state_vec >= 0) { // free state memory first
        if (state.size() == 0) {
          cerr << "userFunc: "
               << "state array is empty but num_state_vec > 0" << endl;
          throw(-1);
        }
        state.clear();
      }
      num_design_vec = iwrk[0];
      num_state_vec = iwrk[1];
      assert(num_design_vec >= 0);
      assert(num_state_vec >= 0);
      design.resize(num_design_vec);
      for (int i = 0; i < num_design_vec; i++)
        design[i].resize(num_design);
      state.resize(num_state_vec);
      for (int i = 0; i < num_state_vec; i++)
        state[i].resize(num_var);
      break;
    }
    case kona::axpby_d: {// using design array set
      // iwrk[0] = dwrk[0]*iwrk[1] + dwrk[1]*iwrk[2]
      int i = iwrk[0];
      int j = iwrk[1];
      int k = iwrk[2];
      assert((i >= 0) && (i < num_design_vec));
      assert(j < num_design_vec);
      assert(k < num_design_vec);
      
      double scalj = dwrk[0];
      double scalk = dwrk[1];
      if (j == -1) {
        if (k == -1) { // if both indices = -1, then all elements = scalj
          design[i] = scalj;
          
        } else { // if just j = -1 ...
          if (scalk == 1.0) {
            // direct copy of vector k with no scaling
            design[i] = design[k];
            
          } else {
            // scaled copy of vector k
            design[i] = design[k];
            design[i] *= scalk;
          }
        }
      } else if (k == -1) { // if just k = -1 ...
        if (scalj == 1.0) {
          // direct copy of vector j with no scaling
          design[i] = design[j];
          
        } else {
          // scaled copy of vector j
          design[i] = design[j];
          design[i] *= scalj;
        }
      } else { // otherwise, full axpby
        design[i].EqualsAXPlusBY(scalj, design[j], scalk, design[k]);
      }
      break;
    }
    case kona::axpby_s: {// using state array set
      // iwrk[0] = dwrk[0]*iwrk[1] + dwrk[1]*iwrk[2]
      int i = iwrk[0];
      int j = iwrk[1];
      int k = iwrk[2];
      assert((i >= 0) && (i < num_state_vec));
      assert(j < num_state_vec);
      assert(k < num_state_vec);

      double scalj = dwrk[0];
      double scalk = dwrk[1];
      if (j == -1) {
        if (k == -1) { // if both indices = -1, then all elements = scalj
          state[i] = scalj;

        } else { // if just j = -1 ...
          if (scalk == 1.0) {
            // direct copy of vector k with no scaling
            state[i] = state[k];

          } else {
            // scaled copy of vector k
            state[i] = state[k];
            state[i] *= scalk;
          }
        }
      } else if (k == -1) { // if just k = -1 ...
        if (scalj == 1.0) {
          // direct copy of vector j with no scaling
          state[i] = state[j];

        } else {
          // scaled copy of vector j
          state[i] = state[j];
          state[i] *= scalj;
        }
      } else { // otherwise, full axpby
        state[i].EqualsAXPlusBY(scalj, state[j], scalk, state[k]);
      }
      break;
    }
    case kona::innerprod_d: {// using design array set
      // dwrk[0] = (iwrk[0])^{T} * iwrk[1]
      int i = iwrk[0];
      int j = iwrk[1];
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= 0) && (j < num_design_vec));
      dwrk[0] = InnerProd(design[i], design[j]);
      break;
    } 
    case kona::innerprod_s: {// using state array set
      // dwrk[0] = (iwrk[0])^{T} * iwrk[1]
      int i = iwrk[0];
      int j = iwrk[1];
      assert((i >= 0) && (i < num_state_vec));
      assert((j >= 0) && (j < num_state_vec));
      dwrk[0] = InnerProd(state[i], state[j]);
      break;
    }
    case kona::eval_obj: {// evaluate the objective
      int i = iwrk[0];
      int j = iwrk[1];      
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= -1) && (j < num_state_vec));
      nozzle_shape.SetCoeff(design[i]);
      solver.set_area(nozzle_shape.Area(solver.get_x_coord()));
      if (j == -1) {
        // need to solve for the state first
	solver.InitialCondition(rho_R, rho_u_R, e_R);
        iwrk[0] = solver.NewtonKrylov(100, 1.e-6);
      } else {
        iwrk[0] = 0; // no precondition calls
        solver.set_q(state[j]);
      }
      dwrk[0] = solver.CalcInverseDesign();
      break;
    }
    case kona::eval_pde: {// evaluate PDE at (design,state) =
                         // (iwrk[0],iwrk[1])
      int i = iwrk[0];
      int j = iwrk[1];
      int k = iwrk[2];
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= 0) && (j < num_state_vec));
      assert((k >= 0) && (k < num_state_vec));
      nozzle_shape.SetCoeff(design[i]);
      solver.set_area(nozzle_shape.Area(solver.get_x_coord()));
      solver.set_q(state[j]);
      solver.CalcResidual();
      state[k] = solver.get_res();
      break;
    }
    case kona::jacvec_d: {// apply design component of the Jacobian-vec
      int i = iwrk[0];
      int j = iwrk[1];
      int k = iwrk[2];
      int m = iwrk[3];
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= 0) && (j < num_state_vec));
      assert((k >= 0) && (k < num_design_vec));
      assert((m >= 0) && (m < num_state_vec));
      nozzle_shape.SetCoeff(design[i]);
      solver.set_area(nozzle_shape.Area(solver.get_x_coord()));
      solver.set_q(state[j]);

      InnerProdVector dArea(num_design, 0.0);
      dArea = nozzle_shape.AreaForwardDerivative(solver.get_x_coord(),
                                                 design[k]);
      solver.JacobianAreaProduct(dArea, state[m]);
      break;
    }
    case kona::jacvec_s: {// apply state component of the Jacobian-vector
      // product to vector iwrk[2]
      // recall; iwrk[0], iwrk[1] denote where Jacobian is evaluated
      int i = iwrk[0];
      int j = iwrk[1];
      int k = iwrk[2];
      int m = iwrk[3];
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= 0) && (j < num_state_vec));
      assert((k >= 0) && (k < num_state_vec));
      assert((m >= 0) && (m < num_state_vec));
      nozzle_shape.SetCoeff(design[i]);
      solver.set_area(nozzle_shape.Area(solver.get_x_coord()));
      solver.set_q(state[j]);
      solver.JacobianStateProduct(state[k], state[m]);
      break;
    }
    case kona::tjacvec_d: {// apply design component of Jacobian to adj
      int i = iwrk[0];
      int j = iwrk[1];
      int k = iwrk[2];
      int m = iwrk[3];
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= 0) && (j < num_state_vec));
      assert((k >= 0) && (k < num_state_vec));
      assert((m >= 0) && (m < num_design_vec));
      nozzle_shape.SetCoeff(design[i]);
      solver.set_area(nozzle_shape.Area(solver.get_x_coord()));
      solver.set_q(state[j]);

      InnerProdVector dArea(nodes, 0.0);
      solver.JacobianTransposedAreaProduct(state[k], dArea);
      design[m] = nozzle_shape.AreaReverseDerivative(
          solver.get_x_coord(), dArea);
      break;
    }
    case kona::tjacvec_s: {// apply state component of Jacobian to adj
      int i = iwrk[0];
      int j = iwrk[1];
      int k = iwrk[2];
      int m = iwrk[3];
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= 0) && (j < num_state_vec));
      assert((k >= 0) && (k < num_state_vec));
      assert((m >= 0) && (m < num_state_vec));
      nozzle_shape.SetCoeff(design[i]);
      solver.set_area(nozzle_shape.Area(solver.get_x_coord()));
      solver.set_q(state[j]);
      solver.JacobianTransposedStateProduct(state[k], state[m]);
      break;
    }
    case kona::eval_precond: {// build the preconditioner if necessary
      int i = iwrk[0];
      int j = iwrk[1];
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= 0) && (j < num_state_vec));
      nozzle_shape.SetCoeff(design[i]);
      solver.set_area(nozzle_shape.Area(solver.get_x_coord()));
      solver.set_q(state[j]);
      solver.BuildAndFactorPreconditioner();
      break;
    } 
    case kona::precond_s: {// apply primal preconditioner to iwrk[2]
      // recall; iwrk[0], iwrk[1] denote where preconditioner is
      // evaluated and in this case, they are not needed
      int i = iwrk[0];
      int j = iwrk[1];
      int k = iwrk[2];
      int m = iwrk[3];
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= 0) && (j < num_state_vec));
      assert((k >= 0) && (k < num_state_vec));
      assert((m >= 0) && (m < num_state_vec));
      solver.Precondition(state[k], state[m]);
      iwrk[0] = 1; // one preconditioner application
      break;
    }
    case kona::tprecond_s: {// apply adjoint preconditioner to iwrk[2]
      int i = iwrk[0];
      int j = iwrk[1];
      int k = iwrk[2];
      int m = iwrk[3];
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= 0) && (j < num_state_vec));
      assert((k >= 0) && (k < num_state_vec));
      assert((m >= 0) && (m < num_state_vec));
      solver.PreconditionTransposed(state[k], state[m]);
      iwrk[0] = 1; // one preconditioner application
      break;
    }
    case kona::grad_d: {// design component of objective gradient
      int i = iwrk[0];
      int j = iwrk[1];
      int k = iwrk[2];
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= 0) && (j < num_state_vec));
      assert((k >= 0) && (k < num_design_vec));
      design[k] = 0.0;
      break;
    }
    case kona::grad_s: {// state component of objective gradient
      int i = iwrk[0];
      int j = iwrk[1];
      int k = iwrk[2];
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= 0) && (j < num_state_vec));
      assert((k >= 0) && (k < num_state_vec));
      nozzle_shape.SetCoeff(design[i]);
      solver.set_area(nozzle_shape.Area(solver.get_x_coord()));
      solver.set_q(state[j]);
      solver.CalcInverseDesigndJdQ(state[k]);
      //cout << "dJdQ.Norm2() = " << state[k].Norm2() << endl;
      break;
    }
    case kona::initdesign: {// initialize the design variables
      int i = iwrk[0];
      assert((i >= 0) && (i < num_design_vec));
      //design[i] = 0.0; // all coefficients set to zero

      // in case the nozzle has not been initiated
      design[i] = 1.0;
      nozzle_shape.SetCoeff(design[i]);
      // fit a b-spline nozzle to a given shape
      InnerProdVector x_coord(nodes, 0.0), area(nodes, 0.0);
      for (int j = 0; j < nodes; j++) {
        x_coord(j) = MeshCoord(length, nodes, j);
        area(j) = InitNozzleArea(x_coord(j)/length);
      }
      nozzle_shape.FitNozzle(x_coord, area);
      nozzle_shape.GetCoeff(design[i]);

      cout << "kona::initdesign design coeff = ";
      for (int n = 0; n < num_design; n++)
        cout << design[i](n) << " ";
      cout << endl;
      
      break;
    }
    case kona::solve: { // solve the primal equations
      int i = iwrk[0];
      int j = iwrk[1];
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= 0) && (j < num_state_vec));
      nozzle_shape.SetCoeff(design[i]);
      solver.set_area(nozzle_shape.Area(solver.get_x_coord()));
      solver.WriteTecplot(1.0, 1.0);
      solver.InitialCondition(rho_R, rho_u_R, e_R);
      iwrk[0] = solver.NewtonKrylov(20, 1.e-10);
      state[j] = solver.get_q();
      break;
    }
    case kona::adjsolve: {// solve the adjoint equations
      int i = iwrk[0];
      int j = iwrk[1];
      int k = iwrk[2];
      assert((i >= 0) && (i < num_design_vec));
      assert((j >= 0) && (j < num_state_vec));
      assert((k >= 0) && (k < num_state_vec));
      nozzle_shape.SetCoeff(design[i]);
      solver.set_area(nozzle_shape.Area(solver.get_x_coord()));
      solver.set_q(state[j]);
      InnerProdVector dJdQ(num_var, 0.0);
      solver.CalcInverseDesigndJdQ(dJdQ);
      dJdQ *= -1.0;
      iwrk[0] = solver.SolveAdjoint(100, 1.e-8, dJdQ, state[k]);
      break;
    }
    case kona::info: {// supplies information to user
      // current design is in iwrk[0]
      // current pde solution is in iwrk[1]
      // current adjoint solution is in iwrk[2]
      int i = iwrk[0];
      int j = iwrk[1];
      int k = iwrk[2];
      int iter = iwrk[3];
      nozzle_shape.SetCoeff(design[i]);

#if 0
      // uncomment to list B-spline coefficients
      cout << "before kona::info set_area: design coeff = ";
      for (int n = 0; n < num_design; n++)
        cout << design[i](n) << " ";
      cout << endl;
#endif
      
#if 0
      cout << "total preconditioner calls (solver says) = " 
           << solver.TotalPreconditionerCalls() << endl;
#endif
      
#if 0
      solver.set_area(nozzle_shape.Area(solver.get_x_coord()));
      cout << "before kona::info NewtonKrylov..." << endl;
      solver.NewtonKrylov(100, 1.e-8);
      //if (iter > 0) solver.set_q(state[j]);

      solver.WriteTecplot(1.0, 1.0);
#endif

      path.push_back(design[i](0));
#if 0
      cout << "path.back() = " << path.back()
           << ": design[i](0) = " << design[i](0) << endl;
#endif

#if 0
      // need rho_ref and a_ref for WriteTecplot
      double rho, rho_u, e;
      CalcFlowExact(kGamma, kRGas, kAreaStar, area_left, true, 
                    kTempStag, kPressStag, rho, rho_u, e);
      double rho_ref = rho;
      double press = (kGamma - 1.0)*(e - 0.5*rho_u*rho_u/rho);
      double a_ref = sqrt(kGamma*press/rho_ref);
#endif
#if 0
      string filename("flow_at_opt_iter");
      std::stringstream ss;
      ss << opt_iter;
      filename.append(ss.str());
      filename.append(".dat");
      //solver.WriteTecplot(rho_ref, a_ref, filename);
      solver.WriteTecplot(1.0, 1.0, filename);
      opt_iter++;
#endif

#if 0
      cout << "(design), (state), (adjoint) = (" 
           << design[i][0] << "), (" 
           << state[j][0] << "," << state[j][1] << "), ("
           << state[k][0] << "," << state[k][1] << ")" << endl;
#endif
      break;
    }
    default: {
      cerr << "userFunc: "
           << "unrecognized request value: request = "
           << request << endl;
      throw(-1);
      break;
    }
  }
  return 0;
}

void WritePath() {
  string filename = "path.dat";
  ofstream fout(filename.c_str());
  fout.precision(12);
  fout << "TITLE = \"Quasi-1D-Euler Solution Path\"" << endl;
  fout << "VARIABLES=\"lambda\",\"coeff\"" << endl;
  fout << "ZONE I=" << path.size() << ", DATAPACKING=POINT" << endl;
  for (int i = 0; i < path.size(); i++) {
    fout << static_cast<double>(i) 
        / static_cast<double>(path.size()-1) << " ";
    fout << path[i] << " ";
    fout << endl;
  }
  fout.close();
}
