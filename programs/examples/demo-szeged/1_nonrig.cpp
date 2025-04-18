// see this file for the initialization of the library for this particular example
#include "setup.h"

using namespace std;
// you can uncomment this and remove capd:: everywhere in this file
// using namespace capd;

///////////////////////////////////////////
///////////////////////////////////////////
int main(){
	// you can switch this to false to get
	// figures saved to .png instead of showing on screen
	bool show_on_screen = true;
	// the parameters for our DDE
	// (*) x'(t) = -a*x(t) + b*(x(t-tau))^2  =: f(x(t), x(t-tau))
	double delay = 1.0;
	double a = 3.1;
	double b = 7.95;
	// initial value to be used later, close to the fixed point a/b
	double x0 = a/b + 0.1;
	// the technical parameters of the representation
	// of solution segments to be used in computations.
	int p = 64, n = 4;	// n will be the order of representation (approximate order of error of the method)
	// grid step size is basic delay divided by the steps per full delay
	double h = delay / p;
	// make a grid and compute important time points on the grid
	// Note: we are using DS:: to access the data types defined in a helper class NonrigorousSetup defined in "setup.h"
	DS::Grid grid(h);
	DS::TimePoint tau = grid(p); 	// tau must be TimePoint from the grid! Not double!
	DS::TimePoint zero = grid(0);	// tau must be TimePoint from the grid! Not double!
	// define f in r.h.s. of (*) with specified parameters (no delay here!)
	DS::Eq eq(a, b);			// define f(x, y) with a, b set
	DS::DDEq dde(eq, tau);		// say that x = x(t), y = x(t-tau) in the above
	// basic class to extend the solution over intervals of length h
	DS::Solver dde_solver(dde);

	// this does method of steps, but on intervals of length h
	// (better approximation). We will do 4 * p steps, so it gives time = 4 * tau
	cout << "A) SIMPLE INTEGRATION" << endl;
	// define constant solution on interval [-tau, 0], with order of representation n
	// the value must be vector, even in scalar case (so we use {x0} to make a DVector).
	DS::Solution solution_A(-tau, zero, n, {x0});
	for (int i = 0; i < 4*p; i++){
		dde_solver(solution_A);
	}
	// the followinf function is from ddeshelper - a collection of useful small routines
	// in this case we are using gnuplot to generate values
	// of solution_A(t_i) for t_i = i*0.1*h, such that t_i \in Domain(solution_A) = [-tau, 4*tau]
	// note, the plotter is getting better graphical representation by evaluating
	// the solution_A in points between the grid points!
	capd::ddeshelper::plot_value("nonrig-A-", 0.1*h, solution_A, show_on_screen);
	// NOTE: by default, if you have gnuplot installed, the above command
	// should produce a gnuplot window with a plot. In case it does not,
	// it might be that you a) do not have gnuplot, or b) have compiled the library without a DDES_ALLOW_SYSTEM flag.
	// in the latter case, you ight add the following line:
	// int sys = system("gnuplot -p nonrig-A-ddes-plot.gp");
	// if that works, you might need to add this kind of line after each plot_value().

	cout << "B) INTEGRATION WITH BUILD-IN time map" << endl;
	// we are re-using the existing variable to store initial data
	solution_A = DS::Solution (-tau, zero, n, {x0});
	dde_solver(solution_A, 4*tau); // we need second argument to be TimePoint on the grid!
	capd::ddeshelper::plot_value("nonrig-B-", 0.1*h, solution_A, show_on_screen);

	cout << "C) FUNCTION AS AN INITIAL SEGMENT" << endl;
	// this defines a function y_0(t) = ((t + tau)/tau)*c
	// i.e. linear from y_0(-tau) = 0 to y_0(0) = c
	// the last argument n+1 tells CAPD to generate code to differentiate the function w.r.t. t
	// up to n+1 derivatives. This value should be bigger than n used to define segments, so we supply n+1.
	capd::DMap y_0("var:t;par:c,tau;fun:((t + tau)/tau)*c;", n+1);
	y_0.setParameter("c", 1.7*a/b);
	y_0.setParameter("tau", tau);
	// define segment with the function y_0(t) given above
	// here we do not use { } as the y_0 is a function that returns vectors.
	DS::Solution solution_B(-tau, zero, n, y_0);
	dde_solver(solution_B, 4*tau);
	capd::ddeshelper::plot_value("nonrig-C-", 0.1*h, solution_B, show_on_screen);

	// we want to find the periodic solution, which has at least one unstable direction
	// we will use Newton operator for that

	// First, define Poincare map to the Section
	DS::Section S(1, 0, a/b); // a/b is a fixed point, we are looking for a periodic orbit around this point.
	// this says, we are looking for P(x_0) = x_t, such that x_t \in S, i.e. x_t(0) = a/b
	// capd::poincare::MinusPlus adds a requirement that x'_t(0) > 0 (so x(t+s) goes from Minus to Plus sign
	// for s \in small neigbourhood [-delta, delta] around t)
	DS::PoincareMap P(dde_solver, S, capd::poincare::MinusPlus);
	P.setRequiredSteps(0);	// this is technical, it says that P is not required to find just the first intersection with the section.
	P.setMaxSteps(4*p);		// this is technical, it prevents infinite loop in case of the section cnnot be reached (e.g. solution diverges to infinity)

	cout << "D) INITIAL POINCARE MAP" << endl;
	// define initial segment, as before
	DS::Solution X(-tau, zero, n, y_0);
	// we draw initial segment x_0
	capd::ddeshelper::plot_value("nonrig-D-X-", 0.1*h, X, show_on_screen);
	// compute the image under the Poincare map P.
	auto PX = P(X);
	// we draw the P(x_0)
	capd::ddeshelper::plot_value("nonrig-D-PX-", 0.1*h, PX, show_on_screen);

	// we did one Poincare to get the initial guess
	// now we will use newton to refine this.
	// we are using ddeshelper::refinePeriodic for it
	// if you are interested how it works, check its source code.
	cout << "E) REFINE" << endl;
	DHelper helper({a, b, tau}, p, n, 0, 4*p);
	helper.setCrossingDirection(capd::poincare::MinusPlus);
	// the matrix V will be the monodromy matrix, the DP will be the Jacobian matrix DP(x_0)
	capd::DMatrix V(X.dimension(), X.dimension()), DP = V;
	// it does one step of Newton method to find the solution of P(x) = x
	// it returns the numerical approximation (finite dimensional) of \frac{\partial}{\partial x0} \varphi(t_P, x_0) = V
	// and the Jacobian of Poincare matrix: \frac{d}{d x0} \varphi(t_P(x_0), x_0) = DP (it takes into account the t_P(x_0)!)
	auto diff = helper.refinePeriodic(cerr, S, X, V, DP);
	// the smaller the diff, the better approximation, but more steps might be needed
	double epsilon = 1e-10;
	while (std::abs(diff) > epsilon){
		cout << "diff " << diff << endl;
		diff = helper.refinePeriodic(cerr, S, PX, V, DP);
	}
	// copy PX to X and then draw the X and P(X)
	// we are re-using variables here.
	X = PX; PX = P(X);
	capd::ddeshelper::plot_value("nonrig-E-X-ref-", 0.1*h, X, show_on_screen);
	capd::ddeshelper::plot_value("nonrig-E-PX-ref-", 0.1*h, PX, show_on_screen);

	// we draw the periodic orbit candidate over a very long interval
	// it should diverge after maybe 20-30 delays, depending on epsilon
	cout << "F) TEST THE PERIODIC CANDIDATE" << endl;
	auto solution_X = X;
	for (int i = 0; i < 50*p; i++){
		dde_solver(solution_X);
	}
	capd::ddeshelper::plot_value("nonrig-F-TEST-", 0.1*h, solution_X, show_on_screen);

	// now we can see how the solution segments are defined:
	cout << "= X ========================================================" << endl;
	cout << X << endl << endl;
	// this produces a lot of output, so uncomment only if you want to see it.
	//	cout << "= solution_X ===============================================" << endl;
	//	cout << solution_X << endl << endl;
	cout << "= X as a vector ============================================" << endl;
	cout << (capd::DVector)X << endl << endl;

	cout << "X.dimension() = " << X.dimension() << endl;
	cout << "(this is ,,ambient space'' dimension, it is one, because the equation is scalar)" << endl;

	cout << "dimension of representation of X = " << capd::DVector(X).dimension() << endl;
	cout << "(this is actual number of parameters that defines the shape of the segment X)" << endl;

	// we save the candidate to be used later
	// first, we convert to interval vector
	capd::IVector icandidate { capd::DVector(X) };
	// first, we save as a regular text file.
	// this file is portable between operating systems, and is somewhat human-readable
	// it takes more space. If we read it, we might get a slightly different values than those presented.
	ofstream text_out("x0.ivector.txt");
	text_out.precision(16); // more than 16 is too much, less than 15 is to small for accuracy
	text_out << icandidate << endl;
	text_out.close();
	// then we save as a binary file
	// it is not human readable, and not portable between systems
	// but it saves exact values, not some approximations.
	// the name of the file is arbitrary. I like to use extension .ivector.bin
	// to show that it is ivector and in binary format
	capd::ddeshelper::saveBinary("x0.ivector.bin", icandidate);

	cout << "\n\n" << "EVERYTHING DONE. Check 'nonrig-...' files generated in the directory for extra information!" << endl;

	return 0;

}
///////////////////////////////////////////
///////////////////////////////////////////

