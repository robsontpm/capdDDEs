#include <iostream>
#include "equation.h"
#include <cmath>

using Eq = DelayedStringEq<double>;
using Setup = capd::ddeshelper::NonrigorousHelper<Eq>;

int main(int argc, char** argv) {
  std::cout.precision(15);
  Setup::Vector params(2);
  params[0] = 1.0;
  params[1] = 2*M_PI;
  int n = 4;
  int p = 128;
  int iters = -10;

  capd::DMap map("var: x; fun: sin(x), cos(x);", 10);

  Setup setup(params, p, n);
  Setup::Solution initial(setup.grid(), setup.grid().point(-p), setup.grid().point(0), n, map);
  
  Setup::Solution x(setup.grid(), p, n, capd::DVector{1.0, 0.0});
  Setup::Solution solution = setup.integrate(iters, initial, x);
  capd::ddeshelper::PathConfig paths{ ".", "main" };
  setup.drawSolution(paths, solution);
}
