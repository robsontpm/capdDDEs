#include "setup.h"
using namespace std;

int main(){
	int p = 128, n = 4;
	double delay = 1.0, a = 1.64, b = 3.4116, h = delay/p;

	DS::Grid grid(h);
	DS::TimePoint tau = grid(p), t_0 = grid(0);
	DS::Eq f {a};
	DS::DDEq dde(f, tau);
	DS::Solver ndsolve(dde);

	capd::DMap y_0("var:t;par:a,b;fun:b*exp(-a*(t+1))-a;", n+1);
	y_0.setParameter("a", a);
	y_0.setParameter("b", b);
	DS::Solution x(t_0-tau, t_0, n, y_0);

	cout << "Eval: " << x(t_0) << " " << x(-tau) << " " << x(-0.1) << "\n";
	try { x(0.1); } catch (std::exception& e) { cout << e.what() << "\n"; }

	ndsolve(x, 20*tau);
	cout << x.eval(x.rightDomain()) << "\n";
	capd::ddeshelper::plot_value("numeric-solution", h, x, false);

	auto dx = x.dt();
	capd::ddeshelper::plot_value("numeric-derivative", h, dx, false);

	cout << "Sample jet: " << x.jet(t_0) << "\nzeros: \n";
	for (auto& a: x) {
		auto v0 = a->evalAtDelta(0)[0], vh = a->evalAtDelta(h)[0];
		if (v0 * vh < 0.){
			auto dv = a->evalCoeffAtDelta(1, 0)[0];
			auto tz = double(a->t0()) - v0/dv;
			cout << tz << " " << a->eval(tz) << endl;
		}
	}
}

