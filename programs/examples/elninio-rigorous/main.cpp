#include <capd/ddeshelper/DDEHelperRigorous.h>

#include "equation.h"

/*using Eq = ElNinoRig<capd::intervals::Interval<double>>;*/
using Eq = ElNinoRig<capd::DInterval>;
using Setup = capd::ddeshelper::RigorousHelper<Eq, 1>;

int main() {
  Setup::ParamsVector params(4);
  params[0] = {1.0}; // alpha
  params[1] = {3.0}; // beta
  params[2] = {11.0}; // kappa
  params[3] = {0.62}; // delay
  
  int p = 128; // number of points in the interval [-tau, 0]
  int n = 4; // degree of taylor polynomials
  //Setup setup(p, n, params /*, reqSteps = 0, maxSteps = 0, maxOrder = 10 */);
  Setup setup(p, n, params, 0, 0, 100);
  capd::IMap initial_map("var:t;fun:t;", 10);
  /*capd::map::Map<capd::vectalg::Matrix<capd::intervals::Interval<double>, 0, 0>> initial_map("var:t;fun:1", 10);*/
  Setup::Solution initial(setup.grid(), setup.grid().point(-p), setup.grid().point(0), n, initial_map /*, N0=0*/);

  Setup::Solution X = setup.timemap(initial, 20*128 /*, epsilon= 0. */);
  setup.drawSolution(".", "solution", initial);
  
  Setup::Section section(1, 0);
  capd::poincare::CrossingDirection dir = capd::poincare::CrossingDirection::Both;
  capd::DInterval z = {0};
  Setup::Solution Y = X;
  setup.poincare(section, dir, X, z, Y);
  setup.drawSolution(".", "prepoincare", X);
  setup.drawSolution(".", "postpoincare", Y);
}
