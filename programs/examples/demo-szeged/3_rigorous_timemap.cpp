// see this file for the initialization of the library for this particular example
#include "setup.h"

#include <capd/ddeshelper/ddeshelperlib.h>

#include <capd/ddeshelper/DDEHelperDrawing.hpp>
#include <capd/ddeshelper/DDEHelperNonrigorous.hpp>

using namespace std;
using namespace capd;

///////////////////////////////////////////
///////////////////////////////////////////
int main(int argc, char** argv){
	bool show_on_screen = true;
	// instead of double we use interval
	interval delay = 1.0;
	interval a = 3.1;
	interval b = 7.95;
	int p = 64, n = 4;
	interval h = delay / p;

	IS::Grid grid(h);
	IS::TimePoint tau = grid(p);
	IS::TimePoint zero = grid(0);
	IS::Eq eq(a, b);
	IS::DDEq dde(eq, tau);
	IS::Solver dde_solver(dde, 3*n);

	IHelper helper(p, n, {a, b, tau}, n, 4*p);
	// we use helper to get the size of the representation vectors
	// we must create large enough vector to store data,
	// as the size is not saved with the vector.
	IVector candidate(helper.M());
	// we read binary, you can also read from txt (below)
	//capd::ddeshelper::readBinary("x0.ivector.bin", candidate);
	ifstream input("x0.ivector.txt"); input >> candidate; input.close();

	IMatrix C(helper.M(), helper.M()); C.setToIdentity();
	IVector r0(helper.M()), Xi(p);
	for (int i = 0; i < helper.M(); ++i) r0[i] = interval(-1, 1) * 1e-3;
	for (int i = 0; i < p; ++i) Xi[i] = interval(-1, 1) * 0.5;

	// make an affine set in the space (a big cube around the candidate).
	IS::Solution X(-tau, zero, n, { interval(0.) });
	X.set_x(candidate);
	X.set_Cr0(C, r0);
	X.set_Xi(Xi);

	// draw it
	capd::ddeshelper::plot_value("rigorous-timemap-X0-", h, X, show_on_screen);

	// integration, all methods below are equivalent.

	// integrate forward in time - 1st version
	// (explicitly calling the Lohner algorithm - procedure move on the set)
	auto X1 = X; // copy the data
	for (int i = 0; i < 4*p; ++i)
		X1.move(dde_solver);
	capd::ddeshelper::plot_value("rigorous-timemap-X1-", h, X1, show_on_screen);

	// integrate forward in time - 2nd version
	// (similar to nonrigorous version, we move by the given number of steps, each step is h from the grid)
	auto X2 = X; // copy the data
	dde_solver(X2, 4*p);
	capd::ddeshelper::plot_value("rigorous-timemap-X2-", h, X2, show_on_screen);

	// integrate forward in time - 3rd version
	// (similar to nonrigorous version, we can move by a given time that is multiple of basic time step)
	// this looks the closest to what mathematician would like to write :)
	auto X3 = X; // copy the data
	dde_solver(X3, 4*tau);
	capd::ddeshelper::plot_value("rigorous-timemap-X3-", h, X3, show_on_screen);

	// draw a subsurve of length tau (last segment of the solution)
	// TODO: add a method so that X.segment(t) is the equivalent to $X_t$...
	auto Y = X.subcurve(X.currentTime() - tau);
	capd::ddeshelper::plot_value("rigorous-timemap-X-subcurve-", 0.1*h, Y, show_on_screen);

	// starting from an explicitely given initial function
	capd::IMap sinus("var:t;par:a,b;fun:a*sin(b*t);");
	sinus.setDegree(n+2);					 // this must be set bigger than n+1 for n used below to construct a segment.
	sinus.setParameter("a", 0.25);			 // set the parameter
	sinus.setParameter("b", interval::pi()); // rigorous representation of pi
	IS::Solution Z(-tau, zero, n, sinus);
	dde_solver(Z, 20*tau);
	capd::ddeshelper::plot_value("rigorous-timemap-Z-", h, Z, show_on_screen);

	// test how far I can push the single trajectory close to
	// the periodic solution before it explodes
	r0 *= 0.;	// make it a vector of width 0
	Xi *= 0.;	// make it a vector of width 0
	X.set_r0(r0);
	X.set_Xi(Xi);
	dde_solver(X, 40*tau);
	capd::ddeshelper::plot_value("rigorous-timemap-X-thin-", h, X, show_on_screen);

	return 0;

}
///////////////////////////////////////////
///////////////////////////////////////////

