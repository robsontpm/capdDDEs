#ifndef DEMO_SZEGED_SETUP_H_
#define DEMO_SZEGED_SETUP_H_

#include <iostream>
#include <iomanip>
#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>

/**
 * A model map of what I expect from a r.h.s. of a DDE equation.
 * In our case we are dealing with:
 *
 * $x'(t) = -a*x(t) + b * (x(t-tau))^2 = f(x(t), x(t-tau))$
 *
 * This class defines the r.h.s. of this equation, i.e. the function f
 * evaluated at some x, y, f(x, y), such that x, y \in \R^d.
 * The information on delays (i.e. that y=x(t-tau) will be supplied elsewhere later.
 * NOTE: we always expect the first argument of f is x = x(t) (value at the current time)
 * even in cases when the r.h.s. does not depend on x(t) (e.g. x'(t) = f(x(t-tau)))
 *
 * The implemented function is a function $\R^{dimension()} \to \R^{imageDimension()}$,
 * see below.
 */
template<typename ScalarSpec, typename ParamSpec = ScalarSpec>
class DemoSzeged{
public:
	// those are the standard things. They should be more or less the same for any
	// function that fits the capdDDEs code. In this case we use templates
	// to define the ScalarType and ParamType, so we do not need to
	// write two versions later: one for double and one for capd::interval
	// in case you only use one version, you can skip templates and specify the
	// explicit type here.
	typedef ScalarSpec ScalarType;
	typedef ScalarType RealType;
	typedef ParamSpec ParamType;
	typedef capd::vectalg::Vector<ParamSpec, 0> ParamsVectorType;
	typedef capd::vectalg::Vector<ScalarType, 0> VectorType;
	typedef unsigned int size_type;

	/** default constructor, also alows to set the parameters. */
	DemoSzeged(ParamType a=-1, ParamType b=2): a(a), b(b) {}
	/** copy constructor. standard thing in C++ */
	DemoSzeged(DemoSzeged const & orig): a(orig.a), b(orig.b) {}
	/** a constructor that is required with RigorousHelper and NonrigorousHelper classes from capd::ddeshelper namespace */
	DemoSzeged(capd::vectalg::Vector<ParamSpec, 0> const & params): a(params[0]), b(params[1]) {}

	/** output dimension of the internal map f (in other words the dimension of x(t)) */
	static size_type imageDimension() { return 1; }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another. */
	static size_type dimension() { return 2; }
	/** number of parameters to fully configure equation (needed by *Helper classes) */
	static size_type getParamsCount() { return 2; }

	/**
	 * We require that solution class has a templated operator of this signature
	 * the program will pass in the x variable the m values where m = dimension()
	 * in fx there will be a reference to a vector of dimension d = imageDimension().
	 * You do not need to worry about those, they are used internally, you
	 * just need to follow the convention.
	 *
	 * NOTE: x[0],...,x[d] is the value of the solution at the current time x(t) \in \R^d
	 * In our case d = 1 (scalar equation).
	 *
	 * NOTE: pay attention to & and const&. You might skip the reference for
	 * two first arguments, but in the current form it is a little bit faster.
	 */
	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec& x, OutVectorSpec& fx) const {
		// you can rename your variables if you want:
		auto& x_t = x[0];	// x(t)
		auto& x_tau = x[1];	// x(t-tau)
		// then compute the value of f
		fx[0] = -a * x_t + b * x_tau * x_tau;  //f(x(t)) = -a x(t) + b x(t-tau)
		// NOTE: there is no 'return fx;' statement!
	}

	/** this should print the info on your system. It is used by *Helper classes to generate human-readable infos. */
	static std::string show(){
		return "Some exemplary system";
	}

protected:
	ParamType a;
	ParamType b;
};

// this version will be used in nonrigorous code
typedef DemoSzeged<double, double> DEq;
// this version will be used in rigorous code
typedef DemoSzeged<capd::interval, capd::interval> IEq;

// this is a basic class that creates a collection of classes
// (type names) to be used in your code.
// you might define classes by hand in a similar fashion it is done
// in *Setup classes, but it is a little bit messy.
typedef capd::ddes::NonrigorousSetup<DEq> DS;
// this is the setup for rigorous computations.
typedef capd::ddes::RigorousSetup<IEq> IS;

// those are optional libraries (capd::ddeshelper namespace)
// we use ploting subrotines that generates gnuplot compatible figures
// and the *Helper classes with some subroutines usually needed
// to do computer assisted proofs (CAP) with Covering Relations.
// the helper part is not neccessary to do computations and proofs, but it helps :)
#include <capd/ddeshelper/ddeshelperlib.h>
#include <capd/ddeshelper/DDEHelperDrawing.hpp>
#include <capd/ddeshelper/DDEHelperNonrigorous.hpp>
#include <capd/ddeshelper/DDEHelperRigorous.hpp>

// this is a helper to handle approximate solutions and maybe generate data for proofs
typedef capd::ddeshelper::NonrigorousHelper<DEq, 1, capd::DMatrix, capd::DVector> DHelper;
// this is a helper to assist in rigorous computations.
typedef capd::ddeshelper::RigorousHelper<IEq, 1, capd::IMatrix, capd::IVector> IHelper;





//// TODO: add it to helper!
//#include <capd/mpcapdlib.h>
//
//using namespace capd;
//using namespace multiPrec;
//
//MpInterval interval_to_mpi(const interval& a){
//	return MpInterval(a.leftBound(), a.rightBound());
//}
//
//MpInterval interval_to_mpi(const Interval& a, int precission){
//	int old_precision = MpReal::getDefaultPrecision();
//	MpReal::setDefaultPrecision(precission);
//	return MpInterval(a.leftBound(), a.rightBound());
//	MpReal::setDefaultPrecision(old_precision);
//}
//
//template<typename MatrixSpec>
//void to_mpi_matrix(const MatrixSpec& inA, MpIMatrix& outA, int precission = 512){
//	int old_precision = MpReal::getDefaultPrecision();
//	MpReal::setDefaultPrecision(precission);
//	int size = inA.numberOfColumns(); // assume square matrix
//	outA = MpIMatrix(size, size);
//	for (int i = 0; i < size; ++i)
//		for (int j = 0; j < size; ++j)
//			outA[i][j] = MpInterval(leftBound(inA[i][j]), rightBound(inA[i][j]));
//	MpReal::setDefaultPrecision(old_precision);
//}
//
//template<typename MatrixSpec>
//void from_mpi_matrix(const MpIMatrix& inA, MatrixSpec& outA){
//	int size = inA.numberOfColumns(); // assume square matrix
//	outA = MatrixSpec(size, size);
//	for (int i = 0; i < size; ++i)
//		for (int j = 0; j < size; ++j)
//			outA[i][j] = mpi_to_interval(inA[i][j]);
//}

#endif /* DEMO_ELNINIO_NONRIG_SETUP_H_ */
