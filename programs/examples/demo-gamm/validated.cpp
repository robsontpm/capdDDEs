#include "setup.h"
using namespace std;
using namespace capd;

int main(){
	int p = 128, n = 4; double a = 1.64;
	interval delay = 1.0, b {3.4115, 3.4116}, h = delay/interval(p);

	IS::Grid grid(h);
	IS::TimePoint tau = grid(p), t_0 = grid(0);
	IS::Eq f {a};
	IS::DDEq dde(f, tau);
	IS::Solver ndsolve(dde, 2*n);

	IMap y_0("var:t;par:a,b;fun:b*exp(-a*(t+1))-a;", n+2);
	y_0.setParameter("a", a);
	y_0.setParameter("b", b);
	IS::Solution x(t_0-tau, t_0, n, y_0);
	cout << "Eval: " << x.eval(t_0) << " " << x.eval(-tau) << "\n";

	ndsolve(x, 20*tau);
	cout << x.eval(x.rightDomain()) << "\n";
	ddeshelper::plot_value("validated-solution", h, x, false);

	cout << "Sample jet: " << x.jet(t_0) << "\nzeros: \n";
	for (auto& a: x) {
		auto v0 = a->evalAtDelta(0)[0], vh = a->evalAtDelta(h)[0];
		if (v0 * vh < 0.){
			auto dv = a->evalCoeffAtDelta(1, interval(0, 1)*h)[0];
			auto tz = interval(a->t0()) - v0/dv;
			cout << tz << " " << a->eval(tz) << endl;
		}
	}
}

