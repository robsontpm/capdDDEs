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
  const int p = 256; // number of points in the interval [-tau, 0]
  const int n = 4; // degree of taylor polynomials
  Setup setup(params, p, n);
  const Setup::Grid& grid = setup.grid();
  const auto dtau = grid.point(-p);
  const auto zero = grid.point(0);
  
  Solution initial(grid, dtau, zero, n, {0.5});
  auto X = initial;
  auto Y = initial;
  initial = setup.integrate(75*p, initial, initial);
  initial = setup.integrate(25*p, initial, X);
  initial = setup.integrate(25*p, initial, Y);

  setup.drawSolution(".", "initial", initial);
  setup.drawSolution(".", "X", X);
  setup.drawSolution(".", "Y", Y);

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
  RSolution rinitial(rgrid, rdtau, rzero, n, {{0.5-1.0e-12, 0.5+1.0e-12}}, 0);
  capd::IVector r0(1+(n+1)*p);
  for(auto& x: r0) {
    x = {-1.0e-8, 1.0e-8};
  }
  std::cout << "x" << std::endl;
  rinitial.set_x(capd::IVector(X.get_x()));
  std::cout << "Cr0" << std::endl;
  rinitial.set_Cr0(capd::IMatrix::Identity(1+(n+1)*p), r0);
  /*std::cout << "B" << std::endl;*/
  /*rinitial.set_B(capd::IMatrix::Identity(2305));*/
  /*std::cout << "Binv" << std::endl;*/
  /*rinitial.set_Binv(capd::IMatrix::Identity(2305));*/
  /*std::cout << "r" << std::endl;*/
  /*rinitial.set_r(capd::IVector(1+(n+1)*p));*/

  RSolution rX = rsetup.timemap(rinitial, p);
  RSolution rY = rsetup.timemap(rinitial, p);

  rsetup.drawSolution(".", "rX", rX);
  rsetup.drawSolution(".", "rY", rY);
  bool subset = true;
  auto rXhull = rX.hull();
  auto rYhull = rY.hull();
  for(int i = 0; i < rXhull.dimension(); ++i) {
    if(!rXhull[i].contains(rYhull[i])) {
      std::cout << std::setprecision(16);
      std::cout << i << "\n" << rXhull[i] << "\n" << rYhull[i] << "\n";
      subset = false;
      break;
    }
  }
  std::cout << subset << "\n";
}
