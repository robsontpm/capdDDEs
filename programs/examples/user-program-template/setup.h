#ifndef PROGRAM_SETUP_H_
#define PROGRAM_SETUP_H_

#include <iostream>
#include <iomanip>
#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>

template<typename ScalarSpec, typename ParamSpec = ScalarSpec>
class ExemplaryDDE {
public:
	typedef ScalarSpec ScalarType;
	typedef ScalarType RealType;
	typedef ParamSpec ParamType;
	typedef capd::vectalg::Vector<ScalarType, 0> VectorType;
	typedef capd::vectalg::Vector<ParamType, 0> ParamsVectorType;
	typedef unsigned int size_type;

	ExemplaryDDE(ParamType mu=0.01, ParamType ro=1.1, ParamType gamma=1.0): mu(mu), ro(ro), gamma(gamma) {}

	static size_type imageDimension() { return 1; }
	static size_type dimension() { return 2; }
	static size_type getParamsCount() { return 3; }

	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec& x, OutVectorSpec& fx) const {
		auto& x_t = x[0];
		auto& x_tau = x[1];	
		//fx[0] = -mu * x_t + gamma * ro * x_tau; //ro*exp(gamma*x_tau);
		fx[0] = -mu * x_t + ro*exp(gamma*x_tau);

	}
	static std::string show(){ return "Lasota-Wazewski DDE (similar to Mackey-Glass)"; }

protected:
	ParamType mu, ro, gamma;
};

// this version will be used in nonrigorous code
typedef ExemplaryDDE<double, double> DEq;
// this version will be used in rigorous code
typedef ExemplaryDDE<capd::interval, capd::interval> IEq;

// this is a basic class that creates a collection of classes
// (type names) to be used in your code.
// you might define classes by hand in a similar fashion it is done
// in *Setup classes, but it is a little bit messy.
typedef capd::ddes::NonrigorousSetup<DEq> DS;
// this is the setup for rigorous computations.
typedef capd::ddes::RigorousSetup<IEq> IS;


#endif /* PROGRAM_SETUP_H_ */
