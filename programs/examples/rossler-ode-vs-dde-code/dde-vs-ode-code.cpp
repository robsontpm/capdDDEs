/*
 * Here we reproduce the results from the following source:
 *
 *     [1] http://capd.sourceforge.net/capdDynSys/docs/html/examples_RosslerChaoticDynamics.html
 *
 * See there for the definition of two-parameter R\"ossler ODE system.
 *
 * For parameters a=5.7, b=0.2, we do the computer assisted proof of two assertions:
 * (1) There exists trapping region A for the system where the attractor exists
 * (2) There exists a subset of A on which the dynamics is chaotic
 *     (a conjugacy to the shift on 2 symbols for a second iteration of a Poincare map to the section x=0)
 *
 * We also show that DDE code can be used to integrate ODEs almost in the same way as in CAPD,
 * when the number of delays is 0 (no delays).
 *
 * The program is a little bit expanded, to provide more output and to produce some pictures for gnuplot.
 * The original code for the topological conditions is in the file 'capd-original-code.cpp'
 *
 * NOTE: The code is very similar to ODE as we are dealing with ODEs and we does not need to compare
 *       Solution Segments. The code is more involved when dealing with solution segments, but this
 *       cannot be avoided (similarly in PDEs you cannot avoid treatment of the infinite dimensional
 *       Tail part, see works of D. Wilczak and P. Zgliczyński about CAPD use in PDEs for reference).
 */

// comment or uncomment the following for extra options
// DDES_ALLOW_SYSTEM 		when uncommented, the program can produce some pictures from the proof,
//							you need gnuplot installed fpr this to work, and system() command available in C++
// CHECK_TRAPPING_REGION 	checking that the trapping region exists is a long run (many subdivisions),
//							if you only interested in the proof of (2) then you can uncomment it
// EXTRA_OUTPUT				uncommenting this will output more information to the standard output

#define DDES_ALLOW_SYSTEM
#define CHECK_TRAPPING_REGION
//#define EXTRA_OUTPUT

// ====================== DO NOT MODIFY ANYTHING BELOW THIS LINE ===========================================
#include "setup.h"
using namespace capd;

// z - coordinates of the trapping region of the attractor, taken from [1]
const double g_bottom = 0.021;
const double g_top = 0.041;

// y - coordinates of the trapping region, taken from [1]
const double g_left = -10.7;
const double g_right = -2.3;

// y - coordinates of sets on which symbolic dynamics is observed.
// z - coordinates are the same as trapping region
// M = [-8.4,-7.6]x[0.028,0.034]
// N = [-5.7,-4.6]x[0.028,0.034]
const double g_leftM  = -8.4;
const double g_rightM = -7.6;
const double g_leftN  = -5.7;
const double g_rightN = -4.6;

// how much subdivisions to make when proving trapping region
const int SLICE_N = 200;

// this needs to be here, as I use the global variable g_* there. TODO: Not the best programming practice...
#include "gnuplot.h"

/*
 * A generic routine that checks if the image under some iteration of Poincare map of the set
 *
 *   [y1,y2]x[g_bottom,g_top]
 *
 * satisfies condition (function \R^3 \to {0, 1}) given in parameter c.
 *
 * This version uses original ODE codes (http://capd.sourceforge.net/capdDynSys/docs/html/examples_RosslerChaoticDynamics.html).
 */
template<class Condition>
std::pair<IVector, IVector>  checkConditionODE(std::ostream& out, IPoincareMap& pm, double y1, double y2, int i, int N, Condition c, bool& result, int iteration = 2) {
	interval p = (interval(y2) - interval(y1)) / interval(N);
	IVector x = {interval(0.), interval(y1) + interval(i,i+1) * p, interval(g_bottom, g_top)};
	C1Rect2Set s(x);
	interval t = 0;
	IVector y = pm(s, t, iteration);
	out << "    ODE: " << y << endl;
	result = c(y) and result;
	return std::make_pair(x, y);
}

/*
 * A generic routine that checks if the image under some iteration of Poincare map of the set
 *
 *   [y1,y2]x[g_bottom,g_top]
 *
 * satisfies a given condition (some function \R^3 \to {0, 1}) , that is given as the parameter c.
 *
 * This version uses DDE codes.
 */
template<class Condition>
std::pair<IVector, IVector> checkConditionDDE(std::ostream &out, Grid& grid, DDEPoincare& pm, double y1, double y2, int i, int N, Condition c, bool& result, int iteration = 2) {
	interval p = (interval(y2) - interval(y1)) / interval(N);							// same as ODE
	IVector x = {interval(0.), y1 + interval(i,i+1) * p, interval(g_bottom, g_top)};	// same as ODE
	DDESolution X(grid, grid(0), x); 		// compared to ODE, we need to make a set in form of a Solution over a time line divided into a Grid.
	interval t = 0;							// same as ODE
	DDESolution Y = pm(X, t, iteration);	// almost the same as ODE, we only store result as a Solution
	IVector	y = Y;							// and we need to cast it to Vector
	out << "    DDE: " << y << endl;		// same as ODE
	result = c(y) and result;				// same as ODE
	return std::make_pair(x, y);			// same as ODE
}

int main(int, char**){
	// set common parameters for both codes
	int p = 32;			// number of grid points in the basic interval. This will translate to the step size of the method := tau/p.
	int max_order = 4;  // we are doing low order method (original ODE code was order = 20, but this is the test for possible chaos in DDE-perturbed ODE, so we should keep this low (see papers and 'Long Enough Integration time')
	interval a = interval(57) / interval(10);	// parameter a * 10
	interval b = interval(2) / interval(10);	// parameter b * 10
	interval h = interval(1.0) / interval(p);	// the step size for both methods (CAPD ODES fixed step and capdDDEs grid size)

	// setup ODE CAPD equation
	IMap vf("par:a,b;var:x,y,z;fun:-(y+z),x+b*y,b+z*(x-a);");
	vf.setParameter("a", a);
	vf.setParameter("b", b);
	// setup ODE CAPD solver and Poincare Map
	IOdeSolver ode_solver(vf, max_order);
	ICoordinateSection ode_section(3, 0);
	IPoincareMap ode_pm(ode_solver, ode_section, poincare::MinusPlus);
	// extra setup, but similar to what we do for DDEs
	ode_pm.setMaxReturnTime(20);

	// setup the delayed equation.
	// Eq does not have fancy setParameter() functions, but id does not need ones, as it is tailored
	// for this specific Equation, not a general parser of equations as in CAPD.
	// Therefore, it suffices to pass them in the constructor.
	// The DDE-- classes declarations are in setup.h file, see how its done. You can use other names for those classes.
	Grid grid(h);
	Eq eq(a, b);
	DDEq dde(eq);
	// setup DDE solver and Poincare map
	DDESolver dde_solver(dde, max_order);
	DDESection dde_section(3, 0);
	DDEPoincare dde_pm(dde_solver, dde_section, poincare::MinusPlus);
	// extra setup, only for DDEs
	dde_pm.setRequiredSteps(0);	// ODE should not care about this, but we must disable this manually. This is important for true DDEs, but it is controlled by the DDE code itself to be safe.
	dde_pm.setMaxSteps(1000);	// this is equivalent to ode_pm.setMaxReturnTime(), but less convenient probably

	// setting up extra output - if needed
	#ifdef EXTRA_OUTPUT
	std::ostream& out = std::cout;
	#else
	std::ostringstream devnull;
	std::ostream& out = devnull;
	#endif

	// Lambda functions that check some inequalities - they are the same as in original code
	auto mappedLeft   = [] (IVector u) { return u[1] < g_leftM; };
	auto mappedRight  = [] (IVector u) { return u[1] > g_rightN; };
	auto mappedIn     = [] (IVector u) { return u[2] > g_bottom and u[2] < g_top and u[1] > g_left and u[1]<g_right; };

	// here we will store results of the tests
	bool resultDDE = true, resultODE = true;

	// Here we check if [g_left,g_right]x[g_bottom,g_top] is mapped into itself by Poincare map.
	// From these computations we also obtain that the sets N and M are mapped across the horizontal strip.
	// This is one of the required conditions for the covering relations.
	#ifdef CHECK_TRAPPING_REGION

	int N = SLICE_N;
	std::vector<IVector> x0_slices;
	// here we will store data for nice plots
	std::vector<IVector> Px_ODE, Px_DDE;
	for (int i = 0; i < N; ++i){
		ode_pm.setStep(h.rightBound()); // we are doing low order method, without this ODE code will struggle with the step control
		ode_pm.turnOffStepControl();    // we are doing low order method, without this ODE code will struggle with the step control
		auto ddePx = checkConditionDDE(out, grid, dde_pm, g_left, g_right, i, N, mappedIn, resultDDE, 1); // the difference here is that DDEs use Grid to define solutions.
		auto odePx = checkConditionODE(out,       ode_pm, g_left, g_right, i, N, mappedIn, resultODE, 1); // ODEs code does not need grid, it can even have various time steps, see comment about step control above
		// push data for draw later
		x0_slices.push_back(odePx.first);
		Px_ODE.push_back(odePx.second);
		Px_DDE.push_back(ddePx.second);
	}
	cout << "DDES Existence of attractor: " << resultDDE << endl;
	cout << "CAPD Existence of attractor: " << resultODE << endl;

	// output data needed for nice pictures
	double box[] = {g_left, g_right, g_bottom, g_top};
	plotTrappingRegion(box, "ddeplot", x0_slices, Px_DDE);
	plotTrappingRegion(box, "odeplot", x0_slices, Px_ODE);

	#endif // CHECK_TRAPPING_REGION

	resultDDE = true; resultODE = true;
	// Remaining inequalities for the covering relations N=>N, N=>M, M=>M, M=>N.
	checkConditionDDE(out, grid,  dde_pm, g_leftM,  g_leftM,  0, 1, mappedLeft,  resultDDE);
	checkConditionODE(out,        ode_pm, g_leftM,  g_leftM,  0, 1, mappedLeft,  resultODE);
	cout << "DDES P^2( Left (M) ) < Left (M): " << resultDDE << endl;
	cout << "CAPD P^2( Left (M) ) < Left (M): " << resultODE << endl;

	resultDDE = true; resultODE = true;
	checkConditionDDE(out, grid,  dde_pm, g_rightM, g_rightM, 0, 1, mappedRight, resultDDE);
	checkConditionODE(out,        ode_pm, g_rightM, g_rightM, 0, 1, mappedRight, resultODE);
	cout << "DDES P^2( Right(M) ) > Right(N): " << resultDDE << endl;
	cout << "CAPD P^2( Right(M) ) > Right(N): " << resultODE << endl;

	resultDDE = true; resultODE = true;
	checkConditionDDE(out, grid,  dde_pm, g_leftN,  g_leftN,  0, 1, mappedRight, resultDDE);
	checkConditionODE(out,        ode_pm, g_leftN,  g_leftN,  0, 1, mappedRight, resultODE);
	cout << "DDES P^2( Left (N) ) > Right(N): " << resultDDE << endl;
	cout << "CAPD P^2( Left (N) ) > Right(N): " << resultODE << endl;

	resultDDE = true; resultODE = true;
	checkConditionDDE(out, grid,  dde_pm, g_rightN, g_rightN, 0, 1, mappedLeft,  resultDDE);
	checkConditionODE(out,        ode_pm, g_rightN, g_rightN, 0, 1, mappedLeft,  resultODE);
	cout << "DDES P^2( Right(N) ) < Left (M): " << resultDDE << endl;
	cout << "CAPD P^2( Right(N) ) < Left (M): " << resultODE << endl;

	return 0;
}

