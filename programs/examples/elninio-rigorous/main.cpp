#include <capd/ddeshelper/DDEHelperRigorous.h>

#include "equation.h"

/*using Eq = ElNinoRig<capd::intervals::Interval<double>>;*/
using Eq = ElNinoRig<capd::DInterval>;
using Setup = capd::ddeshelper::RigorousHelper<Eq, 1>;

int main() {
  Setup::ParamsVector params(4);
  params[0] = {1.0}; // alpha
  params[1] = {1.0}; // gamma
  params[2] = {1.0}; // kappa
  params[3] = {1.0}; // delay
  
  int p = 128; // number of points in the interval [-tau, 0]
  int n = 4; // degree of taylor polynomials
  Setup setup(p, n, params /*, reqSteps = 0, maxSteps = 0, maxOrder = 10 */);
  capd::IMap initial_map("var:t;fun:1;", 10);
  /*capd::map::Map<capd::vectalg::Matrix<capd::intervals::Interval<double>, 0, 0>> initial_map("var:t;fun:1", 10);*/
  Setup::Solution initial(setup.grid(), setup.grid().point(-p), setup.grid().point(0), n, initial_map /*, N0=0*/);

  setup.timemap(initial, 10*128 /*, epsilon= 0. */);
  setup.drawSolution(".", "solution", initial);
}
