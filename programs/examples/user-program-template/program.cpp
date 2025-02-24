// see this file for the initialization of the library for this particular example
#include "setup.h"

#include <capd/ddeshelper/DDEHelperDrawing.hpp>

using namespace std;
using namespace capd;

///////////////////////////////////////////
///////////////////////////////////////////
int main(int argc, char** argv){
	bool show_on_screen = true;
	// instead of double we use interval
	interval delay = 1.0;
	int p = 64, n = 4;
	interval h = delay / p;

	IS::Grid grid(h);
	IS::TimePoint tau = grid(p);
	IS::TimePoint zero = grid(0);
	IS::Eq eq;
	IS::DDEq dde(eq, tau);
	IS::Solver dde_solver(dde, 4*n);

	// starting from an explicitely given initial function
	capd::IMap ivp("var:t;par:a;fun:1.1;");
	ivp.setDegree(n+2);	// this must be set bigger than n+1 for n used below to construct a segment.
	ivp.setParameter("a", interval::pi() * 2);
	IS::Solution X(-tau, zero, n, ivp);
	capd::ddeshelper::plot_value("X0-", h, X, show_on_screen);

	auto T = 2*tau;
	dde_solver(X, T);
	capd::ddeshelper::plot_value("XT-", h, X, show_on_screen);
	if (T >= tau)
		capd::ddeshelper::plot_value("XT-last-segment-", h, X.subcurve(X.t0()-tau), show_on_screen);

	return 0;

}
///////////////////////////////////////////
///////////////////////////////////////////

