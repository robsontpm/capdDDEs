#include "equation.h"
#include <cmath>

// Aliases for the function in the equation and the helper
using Eq = DelayedSpring<double>;
using Setup = capd::ddeshelper::NonrigorousHelper<Eq>;

int main(int argc, char** argv) {
  // Initialize the vector of parameters
  // params[0] - \alpha in the DelayedSpring
  // params[1] - the value of the delay
  Setup::Vector params(2);
  params[0] = 1.0;
  params[1] = 2*M_PI;

  // degree of the taylor polynomials
  int n = 4;

  // number of points in the inverval [t, t + \tau]
  int p = 128;

  // number of iterations
  // if iters < 0, then it will integrate the solution up to
  // t = |iters| * \tau
  int iters = -10;

  // initialization of the helper
  Setup setup(params, p, n);

  // initial function on the interval [-\tau, 0]
  capd::DMap map("var: x; fun: sin(x), cos(x);", 10);
  Setup::Solution initial(setup.grid(), setup.grid().point(-p), setup.grid().point(0), n, map);
  
  // solving the equation
  Setup::Solution solution = setup.integrate(iters, initial);

  // setting up paths to store the results
  // the results will be stored in the current directory
  // with prefix main-*
  capd::ddeshelper::PathConfig paths{ ".", "main" };

  // plots the solution (x_1(t)) to the equation
  setup.drawSolution(paths, solution);

  // plots the (x_1(t), x_1(t - \tau)) phasespace
  setup.drawDelayMap(paths, solution);
}
