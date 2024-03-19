/**
 * original code of the Example for R\"ossler ODE taken from:
 *
 *    http://capd.sourceforge.net/capdDynSys/docs/html/examples_RosslerChaoticDynamics.html
 *
 * See documentation there.
 *
 * I have removed ConeCondition check, as the capdDDEs cannot do that for now.
 */

#include <iostream>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;

// y - Coordinates of the trapping region of the attractor.
double g_bottom = 0.028;
double g_top = 0.034;

// z - Coordinates of the trapping region of the attractor.
double g_left = -10.7;
double g_right = -2.3;

// z - Coordinates of sets on which symbolic dynamics is observed.
// y - Coordinates are the same.
// M = [-8.4,-7.6]x[0.028,0.034]
// N = [-5.7,-4.6]x[0.028,0.034]
double g_leftM  = -8.4;
double g_rightM = -7.6;
double g_leftN  = -5.7;
double g_rightN = -4.6;

/*
 * A generic routine that checks
 * if the image under some iteration of Poincare map of the set
 *   [y1,y2]x[g_bottom,g_top]
 * satisfies condition 'c'.
 */
template<class Condition>
bool checkCondition(IPoincareMap& pm, double y1, double y2, int N, Condition c, int iteration = 2) {
  bool result = true;
  interval p = (interval(y2) - interval(y1)) / interval(N);
  for (int i = 0; i < N; ++i) {
    IVector x = {interval(0.), interval(y1) + interval(i,i+1) * p, interval(g_bottom, g_top)};
    C0HOTripletonSet s(x);
    IVector y = pm(s, iteration);
    result = result and c(y);
  }
  return result;
}

int main() {
  cout << boolalpha;
  try {
    int order = 20;
    interval a = interval(57) / interval(10);
    interval b = interval(2) / interval(10);
    IMap vf("par:a,b;var:x,y,z;fun:-(y+z),x+b*y,b+z*(x-a);");
    vf.setParameter("a", a);
    vf.setParameter("b", b);

    IOdeSolver solver(vf, order);
    ICoordinateSection section(3, 0);
    IPoincareMap pm(solver, section, poincare::MinusPlus);

	// Lambda functions that check some inequalities.
	auto mappedLeft   = [] (IVector u) { return u[1] < g_leftM; };
	auto mappedRight  = [] (IVector u) { return u[1] > g_rightN; };
	auto mappedIn     = [] (IVector u) { return u[2] > g_bottom and u[2] < g_top and u[1] > g_left and u[1]<g_right; };

    // Here we check if [g_left,g_right]x[g_bottom,g_top] is mapped into itself by Poincare map.
    // From these computations we also obtain that the sets N and M are mapped across the horizontal strip.
    // This is one of the required conditions for the covering relations.
    cout << "Existence of attractor: " << checkCondition(pm, g_left, g_right, 200, mappedIn, 1) << endl;

    // Remaining inequalities for the covering relations N=>N, N=>M, M=>M, M=>N.
    cout << "P^2( Left (M) ) < Left (M): " << checkCondition( pm, g_leftM,  g_leftM,  1, mappedLeft  ) << endl;
    cout << "P^2( Right(M) ) > Right(N): " << checkCondition( pm, g_rightM, g_rightM, 1, mappedRight ) << endl;
    cout << "P^2( Left (N) ) > Right(N): " << checkCondition( pm, g_leftN,  g_leftN,  1, mappedRight ) << endl;
    cout << "P^2( Right(N) ) < Left (M): " << checkCondition( pm, g_rightN, g_rightN, 1, mappedLeft  ) << endl;

  } catch (exception& e) {
    cout << "\n\nException caught: " << e.what() << endl;
  }
  return 0;
}

// END
