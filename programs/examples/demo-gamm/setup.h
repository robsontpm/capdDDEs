#include <iostream>
#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>

// represents a r.h.s. f : \R^m \to \R^d,
// d=imageDimension(), m=dimension()
class DemoGamm{
public:
	double a;
	static int imageDimension() { return 1; }
	static int dimension() { return 2; }

	template<typename Time, typename Vars, typename Out>
	void operator()(const Time& t, const Vars& x, Out& fx) const {
		auto& x_t = x[0];		// x(t)
		auto& x_tau = x[1];		// x(t-tau)
		//f(x(t), x(t-tau)) = -a x(t) + b x(t-tau)
		fx[0] = -a*x_t + 2*a*x_tau + x_tau*x_tau;
	}
};

typedef capd::ddes::NonrigorousSetup<DemoGamm> DS;
typedef capd::ddes::RigorousSetup<DemoGamm> IS;
