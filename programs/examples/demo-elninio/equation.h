#ifndef DEMOS_ELNINIO_H_
#define DEMOS_ELNINIO_H_

#include <iostream>
#include <iomanip>
#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>
#include <capd/ddeshelper/DDEHelperNonrigorous.hpp>

// Initial setting of parameters. You can change some of them
// by the program params. Then, other programs use values stored in run files!
const int NUM_PARAMS = 6;
const double DEFAULT_PAR_ALPHA = 1.0;
const double DEFAULT_PAR_BETA = 1.0;
const double DEFAULT_PAR_KAPPA = 1.0;
const double DEFAULT_PAR_GAMMA = 1.0;
const double DEFAULT_PAR_TAU1 = 2.0;
const double DEFAULT_PAR_TAU2 = 1.0;
const int DEFAULT_P = 128;
const int DEFAULT_N = 4;
const int DEFAULT_STEPS = -10;
const int DEFAULT_REQUIRED_STEPS = -1;
const int DEFAULT_MAX_STEPS = -10;

/*** definition must be in .h, as this is template. */
template<typename T> inline T get_pi(){ return T::pi(); }
/*** definition must be in .h, as this is inline function. */
template<> inline double get_pi(){ return 3.14159265358979323846; }

/*** used in some DDEs from the article */
template<typename T> inline T tanh(T const& x){
	return (exp(2 * x) - 1) / (exp(2 * x) + 1);
}

/**
 * A model map of what I expect from a Map R^m -> R^n
 * to be good for use with DDE codes.
 *
 */
template<typename ScalarSpec, typename ParamSpec = ScalarSpec>
class ElNino{
public:
	typedef ScalarSpec ScalarType;
	typedef unsigned int size_type;
	typedef ScalarType RealType;
	typedef ParamSpec ParamType;
	typedef capd::vectalg::Vector<ParamSpec, 0> ParamsVectorType;
	typedef capd::vectalg::Vector<ScalarType, 0> VectorType;

	ElNino(): alpha(1), beta(1), gamma(1), kappa(1) {}
	ElNino(ParamType kappa): alpha(1), beta(1), gamma(1), kappa(kappa) {}
	ElNino(ParamType alpha, ParamType beta, ParamType gamma, ParamType kappa): alpha(alpha), beta(beta), gamma(gamma), kappa(kappa) {}
	ElNino(ElNino const & orig): alpha(orig.alpha), beta(orig.beta), gamma(orig.gamma), kappa(orig.kappa) {}
	ElNino(capd::vectalg::Vector<ParamSpec, 0> const & params): alpha(params[0]), beta(params[1]), gamma(params[2]), kappa(params[3]) {}

	ElNino& operator=(ElNino const & orig){ alpha = orig.alpha; beta = orig.beta; gamma = orig.gamma; kappa = orig.kappa; return *this; }

	/** output dimension of the internal map */
	static size_type imageDimension() { return 1; }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	static size_type dimension() { return 2; }
	/** number of parameters to fully configure equation */
	static size_type getParamsCount() { return 4; }

	/**
	 * We require that solution class has a templated operator of this signature
	 * the program will pass in x the m values where m = dimension()
	 * in fx there will be reference to already initialized vector of
	 * dimension d = imageDimension().
	 * The dimensions are as follows: if d = imageDimension() then usually
	 * m = d * (number of delayed terms + 1).
	 * +1 is for the current term, which is always present, and always
	 * assumed to be stored in x[0]. (so if your equation is not dependent on
	 * this, then you should not use it in your formulas).
	 * In this example (Mackey-Glass Eq.), we have: d = 1 (scalar), m = 2
	 * (two terms: current value at 0 and delayed term at t = t - \tau)
	 *
	 * Note that tau is not explicitely present in the equation. The concrete
	 * value of the delay is defined when constructing DiscreteDelaysFunctionalMap
	 * from this template function. See docs there.
	 */
	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		typedef typename OutVectorSpec::ScalarType OutScalarType;
		auto Tt = x[0];
		auto Ttau = x[1];
		fx[0] = -alpha * tanh(kappa * Ttau) + gamma * cos(2 * get_pi<ScalarType>() * t);
	}

	static std::string show(){
		return "El Nino equation as autonomous system: $T'(t) = \\gamma \\cos(2 \\pi s) - \\alpha \\tanh( \\kappa T(t-\\tau))$.";
	}

protected:
	ParamType alpha;
	ParamType beta;
	ParamType gamma;
	ParamType kappa;
};


#endif /* DEMOS_ELNINIO_H_ */

