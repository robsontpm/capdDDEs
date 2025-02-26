// see this file for the initialization of the library for this particular example
#include "setup.h"

#include <capd/ddeshelper/ddeshelperlib.h>

#include <capd/ddeshelper/DDEHelperDrawing.hpp>
#include <capd/ddeshelper/DDEHelperNonrigorous.hpp>

using namespace std;
// you can uncomment this and remove capd:: everywhere in this file
// using namespace capd;

///////////////////////////////////////////
///////////////////////////////////////////
int main(int argc, char** argv){
	// the parameters for our DDE
	// (*) x'(t) = -a*x(t) + b*(x(t-tau))^2  =: f(x(t), x(t-tau))
	double delay = 1.0;
	double a = 3.1;
	double b = 7.95;
	double x0 = a/b + 0.1;
	int p = 16, n = 4;
	double h = delay / p;

	DS::Grid grid(h);
	DS::TimePoint tau = grid(p); 	// tau must be TimePoint from the grid! Not double!
	DS::TimePoint zero = grid(0);	// tau must be TimePoint from the grid! Not double!
	DS::Eq eq(a, b);			// define f(x, y) with a, b set
	DS::DDEq dde(eq, tau);		// say that x = x(t), y = x(t-tau) in the above
	DS::Solver dde_solver(dde);

	capd::DMap y_0("var:t;fun:sin(t)+cos(2*t);", n+1);
	DS::Solution X(grid, -tau, zero, n, y_0);

	// those operations are for the nonrigorous version,
	// but **should** be available also in rigorous.
	// if there are some operations missing in rigorous version please let me know!

	cout << "EVALUATION AT A GIVEN TIME" << endl;
	// TODO: do it

	cout << "TAKING DERIVATIVES" << endl;
	try {
		for (int i = 0; i < 2*n; ++i){
			// you can compute the representations of the i-th derivative
			// but to do so, we need to provide somehow the value at the end point
			// one way is to supply the DDE equations to automatically evaluate x^{(i)}
			// at t = current time.
			auto Y = X.dt(dde, i);
			std::ostringstream name; name << "nonrig-computed-Xdt" << i << "-";
			capd::ddeshelper::plot_value(name.str(), 0.1*h, Y, false);
		}
	} catch (std::exception& e){
		// we do the try {} catch {} block, as when system cannot
		// evalueate proper representation of the derivative it throws an exception.
		cout << "EXCEPTION: " << e.what() << endl;
	}
	try {
		for (int i = 0; i < 2*n; ++i){
			// the other option is just to supply a fixed value.
			auto Y = X.dt({0.}, i);
			std::ostringstream name; name << "nonrig-fixed-Xdt" << i << "-";
			capd::ddeshelper::plot_value(name.str(), 0.1*h, Y, false);
		}
	} catch (std::exception& e){
		cout << "EXCEPTION:" << e.what() << endl;
	}
	cout << "derivatives has the same representation, so you can evaluate them: " << endl;
	auto Y = X.dt(dde, 2);
	cout << "Y(-0.33) = " << Y.eval(-0.33) << endl;

	return 0;

}
///////////////////////////////////////////
///////////////////////////////////////////

