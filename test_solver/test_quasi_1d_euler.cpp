/**
 * \file test_quasi_1d_euler.cpp
 * \brief unit test for the Quasi1DEuler class
 * \version 1.0
 */

#include "../quasi_1d_euler.hpp"

using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char *argv[]) {

  const int nodes = 81;
  const int order = 4;
  const double length = 10.0;

  Quasi1DEuler solver(nodes, order);

  // define and set the x-coordinates and nozzle area
  InnerProdVector x_coord(nodes, 0.0);
  InnerProdVector area(nodes, 0.0);
  for (int i = 0; i < nodes; i++) {
    x_coord(i) = static_cast<double>(i*length) 
        / static_cast<double>(nodes-1);
    if (x_coord(i) < 0.5*length) {
      area(i) = 1.0 + 1.5*(1.0 - x_coord(i)/5.0)*(1.0 - x_coord(i)/5.0);
    } else {
      area(i) = 1.0 + 0.5*(1.0 - x_coord(i)/5.0)*(1.0 - x_coord(i)/5.0);
    }
  }
  solver.set_x_coord(x_coord);
  solver.set_area(area);

  // define and set the left-end boundary conditions
  double rho_ref = 1.14091202011454;
  double p = 9.753431315656936E4;
  double rho = rho_ref;
  double rho_u = rho*65.45103620864865;
  double a_ref = sqrt(kGamma*p/rho);
  double e = p/(kGamma-1.0) + 0.5*rho_u*rho_u/rho;
  // nondimensionalize
  rho_u /= (a_ref*rho_ref);
  e /= (rho_ref*a_ref*a_ref);
  rho = 1.0; // i.e. non-dimensionalize density by itself
  solver.set_bc_left(rho, rho_u, e);

  // set the initial flow to the left boundary
  solver.InitialCondition(rho, rho_u, e);

  // define and set the right-end boundary conditions
  p = 9.277211161772544E4;
  rho = 1.10083845647356;
  rho_u = rho * 113.0560581652928;
  e = p/(kGamma-1.0) + 0.5*rho_u*rho_u/rho;
  // nondimensionalize
  rho /= rho_ref;
  rho_u /= (rho_ref*a_ref);
  e /= (rho_ref*a_ref*a_ref);
  solver.set_bc_right(rho, rho_u, e);

  // define any discretization and solver paramters
  solver.set_diss_coeff(0.04);

  //solver.ExplicitEuler(10000, 1.0, 1e-10);
  solver.NewtonKrylov(100, 1.e-10);
  //solver.CalcResidual();

  solver.WriteTecplot(rho_ref, a_ref);

  const double area_star = 0.8;
  const bool subsonic = true;
  double L2_error, max_error;
  solver.CalcMachError(area_star, subsonic, L2_error, max_error);
  cout << "L2 error in Mach number   = " << L2_error << endl;
  cout << "Lmax error in Mach number = " << max_error << endl;

#if 0
  // test the building of the preconditioner
  solver.BuildAndFactorPreconditioner();
#endif

#if 0
  // test the Jacobian-vector product
  solver.TestJacobianStateProduct();
  solver.WriteTecplot(rho_ref, a_ref);
#endif

#if 0
  // test the transponsed Jacobian-vector product
  solver.TestJacobianTransposedStateProduct();
#endif

#if 0
  // test the Jacobian (w.r.t. area) vector products
  solver.TestJacobianAreaProducts();
#endif

#if 0
  // test the dPress/dQ vector products
  solver.TestDPressDQProducts();
#endif
}
