#include <capd/ddeshelper/DDEHelperRigorous.h>
#include <capd/ddeshelper/DDEHelperNonrigorous.h>
#include <iomanip>
#include <algorithm>

#include "equation.h"

using Eq = ElNinoRig<double>;
using Setup = capd::ddeshelper::NonrigorousHelper<Eq, 1>;
using Solution = Setup::Solution;
using Solver = Setup::Solver;
using ParamType = Setup::ParamType;
using ParamsVector = Setup::ParamsVector;

using REq = ElNinoRig<capd::DInterval>;
using RSetup = capd::ddeshelper::RigorousHelper<REq, 1>;
using RSolution = RSetup::Solution;
using RSolver = RSetup::Solver;
using RParamType = RSetup::ParamType;
using RParamsVector = RSetup::ParamsVector;

int main() {

  // NONRIGOROUS

  const ParamType alpha = {1.0};
  const ParamType beta = {3.0};
  const ParamType kappa = {11.0};
  const ParamType tau = {1.2};
  const Setup::ParamsVector params = {alpha, beta, kappa, tau};
  int pna1 = 800;
  const int p = tau * pna1; // number of points in the interval [-tau, 0]
  const int n = 6;   // degree of taylor polynomials
  Setup setup(params, p, n);
  const Setup::Grid& grid = setup.grid();
  const auto dtau = grid.point(-p);
  const auto zero = grid.point(0);
  std::cout << "dtau = " << dtau << std::endl;
  
  Solution initial(grid, dtau, zero, n, {0.5});
  auto X = initial;
  auto Y = initial;
  auto solution = initial;
  solution = setup.integrate(100*p, initial, X);
  setup.drawSolution(".", "d-initial", initial);
  setup.drawSolution(".", "d-solution", solution);

  solution = setup.integrate(5 * pna1, X, Y);
  setup.drawSolution(".", "d-Ysolution", solution);
  setup.drawSolution(".", "d-X", X);
  setup.drawSolution(".", "d-Y", Y);

  // RIGOROUS

  const RParamType ralpha = {1.0};
  const RParamType rbeta = {3.0};
  const RParamType rkappa = {11.0};
  const RParamType rtau = {1.2};
  const RSetup::ParamsVector rparams = {ralpha, rbeta, rkappa, rtau};
  RSetup rsetup(p, n, rparams, 0, 0, 200);
  const RSetup::Grid& rgrid = rsetup.grid();
  const auto rdtau = rgrid.point(-p);
  const auto rzero = rgrid.point(0);
  RSolution rinitial(rgrid, rdtau, rzero, n, {0.5}, 0);
  capd::IVector r0(rsetup.M());
  for(auto& x: r0) x = {-1.0e-4, 1.0e-4};
  std::cout << "x" << std::endl;
  rinitial.set_x(capd::IVector(X.get_x()));
  std::cout << "Cr0" << std::endl;
  rinitial.set_Cr0(capd::IMatrix::Identity(rsetup.M()), r0);

  capd::IVector Xi(rsetup.p() * rsetup.d());
  for(auto& xi: Xi) xi = {-1.0e-1, 1.0e-1};
  rinitial.set_Xi(Xi);
  /*std::cout << "B" << std::endl;*/
  /*rinitial.set_B(capd::IMatrix::Identity(2305));*/
  /*std::cout << "Binv" << std::endl;*/
  /*rinitial.set_Binv(capd::IMatrix::Identity(2305));*/
  /*std::cout << "r" << std::endl;*/
  /*rinitial.set_r(capd::IVector(1+(n+1)*p));*/

  // RSolution rX = rsetup.timemap(rinitial, p);
  RSolution rX = rinitial;
  for(auto& x: r0) x = {0., 0.};
  rinitial.set_r0(r0);
  RSolution rY = rsetup.timemap(rinitial, 5 * pna1);
  auto hull = rY.hull();
  capd::vectalg::split(hull, r0);
  rinitial = rX;
  rinitial.set_r0(r0 * 0.1);
  rX.set_r0(r0 * 2);
  rY = rsetup.timemap(rinitial, 5 * pna1);

  rsetup.drawSolution(".", "rsolution", rinitial);
  rsetup.drawSolution(".", "rX", rX);
  rsetup.drawSolution(".", "rY", rY);
  int subset = 0;
  auto rXhull = rX.hull();
  auto rYhull = rY.hull();
  for(int i = 0; i < rXhull.dimension(); ++i) {
    if(!rXhull[i].contains(rYhull[i])) {
      std::cout << std::setprecision(16);
      std::cout << i << "\n\t" << rXhull[i] << "\n\t" << rYhull[i] << "\n";
    }else{
        subset += 1;
    }
  }
  std::cout << subset << "/" << rsetup.M() << "\n";
}
