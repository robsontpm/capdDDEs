#ifndef DEMO_SZEGED_SETUP_H_
#define DEMO_SZEGED_SETUP_H_

#include <iostream>
#include <iomanip>
#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>

/**
 * A model map of what I expect from a Map R^m -> R^n
 *
 * x'(t) = -a*x(t) + b * (x(t-tau))^2
 */
template<typename ScalarSpec, typename ParamSpec = ScalarSpec>
class DemoSzeged{
public:
	typedef ScalarSpec ScalarType;
	typedef unsigned int size_type;
	typedef ScalarType RealType;
	typedef ParamSpec ParamType;
	typedef capd::vectalg::Vector<ParamSpec, 0> ParamsVectorType;
	typedef capd::vectalg::Vector<ScalarType, 0> VectorType;

	DemoSzeged(ParamType a=-1, ParamType b=2): a(a), b(b) {}
	DemoSzeged(DemoSzeged const & orig): a(orig.a), b(orig.b) {}
	DemoSzeged(capd::vectalg::Vector<ParamSpec, 0> const & params): a(params[0]), b(params[1]) {}

	DemoSzeged& operator=(DemoSzeged const & orig){ a = orig.a; b = orig.b; return *this; }

	/** output dimension of the internal map */
	static size_type imageDimension() { return 1; }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	static size_type dimension() { return 2; }
	/** number of parameters to fully configure equation */
	static size_type getParamsCount() { return 2; }

	/**
	 * We require that solution class has a templated operator of this signature
	 * the program will pass in x the m values where m = dimension()
	 * in fx there will be reference to already initialized vector of
	 * dimension d = imageDimension().
	 */
	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		typedef typename OutVectorSpec::ScalarType OutScalarType;
		auto x_t = x[0];
		auto x_tau = x[1];
		fx[0] = -a * x_t + b * x_tau * x_tau;
	}

	static std::string show(){
		return "Some exemplary system";
	}

protected:
	ParamType a;
	ParamType b;
};

typedef DemoSzeged<double, double> DEq;
typedef DemoSzeged<capd::interval, capd::interval> IEq;

typedef capd::ddes::NonrigorousSetup<DEq> DS;
typedef capd::ddes::RigorousSetup<IEq> IS;

#endif /* DEMO_ELNINIO_NONRIG_SETUP_H_ */
