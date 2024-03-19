/*
 * This file constitutes part of DDEs rigorous integration framework developed
 * in PhD Thesis under supervision of prof. Piotr Zgliczynski:
 *
 * 		"Rigorous Integration of Delay Differential Equations", Jagiellonian University, 2015
 *
 * When using in scientific work please consider citing my articles (preffered),
 * PhD thesis and/or my webpage. For the publications describing the source code
 * please refer to http://scirsc.org/p/papers. Most notable paper up to date is:
 *
 *     Szczelina, R.; Zgliczyński, P.; "Algorithm for rigorous integration of 
 *     Delay Differential Equations and the computer-assisted proof of periodic 
 *     orbits in the Mackey-Glass equation", http://dx.doi.org/10.1007/s10208-017-9369-5
 *     Foundations of Computational Mathematics (2018), Vol. 18, Iss 6, Pages 1299--1332
 *
 * This work would not be possible without aid and expertise of people involved in
 * CAPD developing library (Computer Assisted Proofs in Dynamics).
 * Please refer to http://capd.ii.uj.edu.pl and consider citing also this library 
 * when using those codes in any scientific work.
 *
 * Author: Robert Szczelina, PhD
 * Faculty of Mathematics and Computer Science, Jagiellonian University AND
 * (former) Małopolska Center of Biotechnology, Jagiellonian University
 * email: 	robert.szczelina@uj.edu.pl
 * www: 	scirsc.org
 *
 * This source code is provided under GNU GPL license 
 * (v.2 or whatever compatible with CAPD license)
 */

#ifndef _CAPD_DDES_SAMPLEEQNS_H_
#define _CAPD_DDES_SAMPLEEQNS_H_

#include <capd/ddes/DDECommon.h>
#include <stdexcept>
#include <string>

namespace capd{
namespace ddes{

/**
 * For now, you need to use class similar to those below to define f(x, y, ..)
 * in your equation. Then you define delays in 'DiscreteDelaysFunctionalMap'
 * or some other capd::ddes::FunctionalMap derived classes. Please
 * remember that usually 'DiscreteDelaysFunctionalMap' must specify
 * the Solution template parameter to be able to determine how to handle solution curves.
 *
 * The setup is much complicated than in the case of ODEs, as the equations must take into account
 * the problem with the grid over which the solution is defined (as of now). See the problem
 * with the continuity of solution at grid points and long enough integration time described
 * in the papers. Because of those problems, the grid, the solution and the equation
 * (DiscreteDelaysFunctionalMap) must be aware to some extent of each other.
 * TODO: (FAR FUTURE) try to reduce dependency only to Grid (i.e. DiscreteDelaysFunctionalMap and Solution depend on the Grid,
 * TODO: (FAR FUTURE) but not on each other). It could be done, if DiscreteDelaysFunctionalMap accepts all solutions that can
 * TODO: (FAR FUTURE) produce concrete Interface of a Jet at Grid::TimePoints)
 *
 * Example:
 * 		typedef BasicDoubleton<IMatrix> SetType;
 * 		typedef DDESolutionCurve<SetType> SolutionType;
 * 		typedef MackeyGlass<typename SolutionType::ScalarType> F;
 * 		typedef DiscreteDelaysFunctionalMap<F> RHS;
 * 		typedef typename SolutionType::GridType GridType;
 * 		int p = 128;
 * 		GridType grid(2.0 / p);
 * 		F f(2.0, 1.0, 8.0);
 * 		RHS rhs(f, grid(2 * p));
 *
 * TODO: (FAR FUTURE) make a class that can use integral differential equation to computeDDECoefficients()
 *
 * TODO: (FAR FUTURE) Unfortunatelly, as of today, capd::map::Map is not good for my use.
 * TODO: (FAR FUTURE) Need to consult Daniel Wilczak on the matter.
 * TODO: (FAR FUTURE) Need to understand the Node and the Graph of execution used by Daniel Wilczak
 * TODO: (FAR FUTURE) to implement my version of computeDDECoefficients with capd::map::Map (seems doable)
 */

/**
 * A model map of what I expect from a Map R^m -> R^n
 * to be good for use with DDE codes.
 *
 */
template<typename ScalarSpec, typename ParamSpec = ScalarSpec>
class MackeyGlass{
public:
	typedef ScalarSpec ScalarType;
	typedef unsigned int size_type;
	typedef ScalarType RealType;
	typedef ParamSpec ParamType;
	typedef capd::vectalg::Vector<ParamSpec, 0> ParamsVectorType;
	typedef capd::vectalg::Vector<ScalarType, 0> VectorType;

	MackeyGlass(): beta(2), gamma(1), n(6) {}
	MackeyGlass(ParamType n): beta(2), gamma(1), n(n) {}
	MackeyGlass(ParamType beta, ParamType gamma, ParamType n): beta(beta), gamma(gamma), n(n) {}
	MackeyGlass(MackeyGlass const & orig): beta(orig.beta), gamma(orig.gamma), n(orig.n) {}
	MackeyGlass(capd::vectalg::Vector<ParamSpec, 0> const & params): beta(params[1]), gamma(params[0]), n(params[2]) {}

	MackeyGlass& operator=(MackeyGlass const & orig){ beta = orig.beta; gamma = orig.gamma; n = orig.n; return *this; }

	/** output dimension of the internal map */
	static size_type imageDimension() { return 1; }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	static size_type dimension() { return 2; }
	/** number of parameters to fully configure equation */
	static size_type getParamsCount() { return 3; }

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
		// x[0] is value at a current time, x[1] is at t-\tau
		// we use identity:
		// (x(t-tau))^n = (e^log(x(t-tau))) ^ n = e^(log(x(t)-tau) * n)
		OutScalarType x1_to_n = exp(log(x[1]) * OutScalarType(n));
		fx[0] = -gamma * x[0] + ( (beta * x[1]) / (1.0 + x1_to_n) );
	}

	static std::string show(){
		return "Mackey-Glass equation: $x'(t) = -gamma x(t) + \\beta {{ x(t-\\tau) }\\over{ 1 + (x(t-\\tau))^n }}$.";
	}

protected:
	ParamType beta;		/** classical parameter, default: 2 */
	ParamType gamma;	/** classical parameter, default: 1 */
	ParamType n;		/** classical parameter, used for investigating biffurcations */
};


/**
 * Scalar eq. x(t) * x(t-tau).
 *
 * Use with 'DiscreteDelaysFunctionalMap' to define tau
 */
class ToyModel{
public:
	typedef long size_type;

	ToyModel(){}
	ToyModel(ToyModel const & orig){}
	ToyModel& operator=(ToyModel const & orig){ return *this; }

	/** output dimension of the internal map */
	size_type imageDimension() const { return 1; }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	size_type dimension() const { return 2; }

	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		fx[0] = x[0]*x[1];
	}

	static std::string show(){
		return "Test: TODO: description";
	}
};

/**
 * Scalar eq. Depends only on the past, f(x(t), x(t-1)) = 2*x(t-1) + (x(t-1))^2. Used for tests.
 *
 * Use with 'DiscreteDelaysFunctionalMap' to define tau
 */
class ToyModelSq{
public:
	typedef long size_type;

	ToyModelSq(){}
	ToyModelSq(ToyModelSq const & orig){}
	ToyModelSq& operator=(ToyModelSq const & orig){ return *this; }

	/** output dimension of the internal map */
	size_type imageDimension() const { return 1; }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	size_type dimension() const { return 2; }

	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		fx[0] = 2. * x[1] + x[1] * x[1];
	}

	static std::string show(){
		return "Test: TODO: description";
	}
};

/**
 * Scalar eq. Depends only on the past, f(x(t), x(t-1)) = -a*x(t-1) + c * x(t-1) + (x(t-1))^2.
 * For c = 0  it has Hopf biffurcation at a^* = (5*pi) / (3*sqrt(3))
 * It has 2 fixed points: 0 and a, a is repelling with one unstable to infty and to 0, 0 is stable
 * Param c is used to move the fixed point a to 0, set c = 2a
 * Use with 'DiscreteDelaysFunctionalMap' to define tau. Tau should be 1
 */
template<typename ScalarSpec, typename ParamSpec>
class ToyModelSqA{
public:
	typedef unsigned int size_type;
	typedef ParamSpec ParamType;
	typedef ScalarSpec ScalarType;
	typedef ScalarType RealType;
	typedef capd::vectalg::Vector<ParamSpec, 0> ParamsVectorType;
	typedef capd::vectalg::Vector<ScalarSpec, 0> VectorType;
	ParamSpec a, c;

	ToyModelSqA(): a(2), c(0) {}
	ToyModelSqA(ParamType a): a(a), c(0) {}
	ToyModelSqA(ToyModelSqA const & orig): a(orig.a), c(orig.c) {}
	ToyModelSqA(ParamsVectorType const & params): a(params[0]), c(params[1]) {}

	ToyModelSqA& operator=(ToyModelSqA const & orig){ a = orig.a; c = orig.c; return *this; }

	/** output dimension of the internal map */
	static size_type imageDimension() { return 1; }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	static size_type dimension() { return 2; }
	/** number of parameters to fully configure equation */
	static size_type getParamsCount() { return 2; }

	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		fx[0] = -a * x[0] + c * x[1] + x[1] * x[1];
	}

	static std::string show(){
		return "Test: TODO: description";
	}
};


/**
 * Vector linear equation of the form:
 *
 * 		x'(t) = A * x(t) + B * x(t-tau) + b
 *
 * Use with 'DiscreteDelaysFunctionalMap' to define tau
 */
template<typename MatrixSpec, typename ParamSpec = typename MatrixSpec::ScalarType>
class LinearMap{
public:
	typedef MatrixSpec MatrixType;
	typedef typename MatrixSpec::RowVectorType VectorType;
	typedef typename VectorType::ScalarType ScalarType;
	typedef typename VectorType::size_type size_type;
	typedef ScalarType RealType;
	typedef ParamSpec ParamType;

	void reinitialize(){
		if (A.numberOfRows() != B.numberOfRows()) throw std::logic_error("LinearMap: A and B has different number of rows");
		if (A.numberOfColumns() != B.numberOfColumns()) throw std::logic_error("LinearMap: A and B has different number of cols");
		if (A.numberOfColumns() != A.numberOfRows()) throw std::logic_error("LinearMap: A and B must be square");
		if (A.numberOfColumns() != b.dimension()) throw std::logic_error("LinearMap: b must be compatible with A and B");
		AB = MatrixType(imageDimension(), dimension());
		for (size_type i = 0; i < imageDimension(); ++i){
			for (size_type j = 0; j < imageDimension(); ++j){
				AB[i][j] = A[i][j];
				AB[i][j+imageDimension()] = B[i][j];
			}
		}
	}
	LinearMap(): A(0), B(0), b(0) {}
	LinearMap(size_type d): A(d, d), B(d, d), b(d), AB(0, 0) { reinitialize(); }
	LinearMap(MatrixType const & A, MatrixType const & B, VectorType const & b): A(A), B(B), b(b) { reinitialize(); }
	LinearMap(LinearMap const & orig): A(orig.A), B(orig.B), b(orig.b), AB(orig.AB) {}
	LinearMap& operator=(LinearMap const & orig){ A = orig.A; B = orig.B; b = orig.b; reinitialize(); return *this; }

	/** output dimension of the internal map */
	size_type imageDimension() const { return A.numberOfRows(); }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	size_type dimension() const { return 2 * A.numberOfRows(); }

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
		// simulation of matrix multiplication by vector-vector multiplication
		InVectorSpec ABrow(AB.numberOfColumns());
		for (size_type i = 0; i < AB.numberOfRows(); ++i){
			for (size_type j = 0; j < ABrow.dimension(); ++j) ABrow[j] = AB[i][j];
			fx[i] = ABrow * x + b[i];
		}
	}

protected:
	MatrixType A;
	MatrixType B;
	VectorType b;
	MatrixType AB;

	static std::string show(){
		return "Simple linear: TODO: description";
	}
};

/**
 * Scalar linear equation of the form:
 *
 * 		x'(t) = a * x(t) + b * x(t-tau) + c
 *
 * Use with 'DiscreteDelaysFunctionalMap' to define tau
 */
template<typename ScalarSpec, typename ParamSpec = ScalarSpec>
class ScalarLinearMap{
public:
	typedef ScalarSpec ScalarType;
	typedef ScalarType RealType;
	typedef ParamSpec ParamType;
	typedef int size_type;

	ScalarLinearMap() {}
	ScalarLinearMap(ScalarType const & a, ScalarType const & b, ScalarType const & c): a(a), b(b), c(c) { }
	ScalarLinearMap(ScalarLinearMap const & orig): a(orig.a), b(orig.b), c(orig.c) {}
	ScalarLinearMap& operator=(ScalarLinearMap const & orig){ a = orig.a; b = orig.b; c = orig.c; return *this; }

	/** output dimension of the internal map */
	size_type imageDimension() const { return 1; }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	size_type dimension() const { return 2; }

	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		fx[0] = a * x[0] + b * x[1] + c;
	}

	static std::string show(){
		return "Simple linear: TODO: description";
	}

protected:
	ParamType a;
	ParamType b;
	ParamType c;
};

/** Mainly for testing how DDE code compares to ODE in CAPD */
class ODEPendulum{
public:
	typedef long size_type;

	ODEPendulum(){}
	ODEPendulum(ODEPendulum const & orig){}
	ODEPendulum& operator=(ODEPendulum const & orig){ return *this; }

	/** output dimension of the internal map */
	size_type imageDimension() const { return 2; }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	size_type dimension() const { return 2; }

	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		fx[0] = x[1];
		fx[1] = -x[0];
	}

	static std::string show(){
		return "Simple oscillator: $x'(t) = y, y'(t) = -x(t)$";
	}
};


/**
 * A Toy Example for proof of chaotic behaviour in DDE
 * (chaotic ODE perturbed by a small DDE term - delay might be big, but perturbation should be small)
 * The equation is:
 *
 *     v'(t) = f(v(t)) + \epsilon f(v(t-\tau))
 *
 * where v = (x, y, z) \in \R^3 and f is the r.h.s. of original R\"ossler ODE:
 *
 * f(x, y, z) = ( -(y+z), (x + ay), (b + z * (x - c)) )
 *
 * The parameters are: a, b, c (as in R\"ossler) and \epsilon.
 */
template<typename ScalarSpec, typename ParamSpec = ScalarSpec>
class RosslerDelay{
public:
	typedef ScalarSpec ScalarType;
	typedef unsigned int size_type;
	typedef ScalarType RealType;
	typedef ParamSpec ParamType;

	RosslerDelay(): a(ParamSpec(4) / ParamSpec(10)), b(2.), c(4.), epsi(0.) {}
	RosslerDelay(ParamType a, ParamType b, ParamType c): a(a), b(b), c(c), epsi(0.) {}
	RosslerDelay(RosslerDelay const & orig): a(orig.a), b(orig.b), c(orig.c), epsi(orig.epsi) {}
	template<typename VecSpec>
	RosslerDelay(VecSpec const& params): a(params[0]), b(params[1]), c(params[2]), epsi(params[3]) {}
	RosslerDelay& operator=(RosslerDelay const & orig){ a = orig.a; b = orig.b; c = orig.c; epsi = orig.epsi; return *this; }

	/** output dimension of the internal map */
	static size_type imageDimension() { return 3; }
	static size_type dimension() { return 6; } // (1 + one delay) * imageDimension()!
	/** number of parameters to fully configure equation */
	static size_type getParamsCount() { return 4; }

	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		// "classic" rossler
		fx[0] = -(x[1]+x[2]);
		fx[1] = x[0]+a*x[1];
		fx[2] = b+x[2]*(x[0]-c);

		// "delayed" term
		fx[0] += epsi * (-(x[4]+x[5]));
		fx[1] += epsi * (x[3]+a*x[4]);
		fx[2] += epsi * (b+x[5]*(x[3]-c));
	}

	static std::string show(){
		return "Rossler Toy Delay Model: $x'(t) = f(x(t)) + \\epsilon f(x(t-\\tau))$, where $f$ is r.h.s of classic Rossler ODE in 3D with a,b,c parameters.";
	}

protected:
	/** classic rossler parameters, see scholarpedia. */
	ParamType a;
	/** classic rossler parameters, see scholarpedia. */
	ParamType b;
	/** classic rossler parameters, see scholarpedia. */
	ParamType c;
	/** extra parameter for a size of the delayed term */
	ParamType epsi;
};




/**
 * x'(t) = a*sin(x(t-1))
 */
template<typename ScalarSpec, typename ParamSpec = ScalarSpec>
class Wischert{
public:
	typedef ScalarSpec ScalarType;
	typedef unsigned int size_type;
	typedef ScalarType RealType;
	typedef ParamSpec ParamType;
	typedef capd::vectalg::Vector<ParamSpec, 0> ParamsVectorType;
	typedef capd::vectalg::Vector<ScalarType, 0> VectorType;

	Wischert(ParamType a = 2.): a(a) {}
	Wischert(Wischert const & orig): a(orig.a) {}
	Wischert(capd::vectalg::Vector<ParamSpec, 0> const & params): a(params[0]) {}

	Wischert& operator=(Wischert const & orig){ a = orig.a; return *this; }

	/** output dimension of the internal map */
	static size_type imageDimension() { return 1; }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	static size_type dimension() { return 2; }
	/** number of parameters to fully configure equation */
	static size_type getParamsCount() { return 1; }

	/** evaluation of rhs */
	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		fx[0] = a * sin(x[1]);
	}

	static std::string show(){
		return "Wischert equation: $x'(t) = a sin(x(t-\\tau))$";
	}

protected:
	ParamType a;		/** classical parameter, default: 2 */
};



/**
 * Lasota-Wazewska equation for chaos as mentioned in scholarpedia article on Mackey-Glass equation.
 *
 * x'(t) = -gamma * x(t) + beta * x(t-tau)^n * e^(-x(t-tau))
 */
template<typename ScalarSpec, typename ParamSpec = ScalarSpec>
class LasotaWazewska{
public:
	typedef ScalarSpec ScalarType;
	typedef unsigned int size_type;
	typedef ScalarType RealType;
	typedef ParamSpec ParamType;
	typedef capd::vectalg::Vector<ParamSpec, 0> ParamsVectorType;
	typedef capd::vectalg::Vector<ScalarType, 0> VectorType;

	LasotaWazewska(): beta(2), gamma(1), n(6) {}
	LasotaWazewska(ParamType n): beta(2), gamma(1), n(n) {}
	LasotaWazewska(ParamType beta, ParamType gamma, ParamType n): beta(beta), gamma(gamma), n(n) {}
	LasotaWazewska(LasotaWazewska const & orig): beta(orig.beta), gamma(orig.gamma), n(orig.n) {}
	LasotaWazewska(capd::vectalg::Vector<ParamSpec, 0> const & params): beta(params[1]), gamma(params[0]), n(params[2]) {}

	LasotaWazewska& operator=(LasotaWazewska const & orig){ beta = orig.beta; gamma = orig.gamma; n = orig.n; return *this; }

	/** output dimension of the internal map */
	static size_type imageDimension() { return 1; }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	static size_type dimension() { return 2; }
	/** number of parameters to fully configure equation */
	static size_type getParamsCount() { return 3; }

	/** Evaluate r.h.s of the equation. */
	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		typedef typename OutVectorSpec::ScalarType OutScalarType;
		// x[0] is value at a current time, x[1] is at t-\tau
		// we use identity:
		// (x(t-tau))^n = (e^log(x(t-tau))) ^ n = e^(log(x(t)-tau) * n)
		OutScalarType x1_to_n = exp(log(x[1]) * OutScalarType(n));
		fx[0] = -gamma * x[0] + beta * x1_to_n * exp(-x[1]);
	}

	static std::string show(){
		return "Lasota-Ważewska equation: $x'(t) = -\\gamma x(t) + \\beta x(t-\\tau)^n e^x(t-\\tau)$";
	}

protected:
	ParamType beta;
	ParamType gamma;
	ParamType n;
};


/**
 * Perez-Malta-Coutinho equation for chaos as mentioned in scholarpedia article on Mackey-Glass equation.
 *
 * x'(t) = -gamma * x(t) + (b0 - b1 * x(t-tau)) * x(t-tau)
 */
template<typename ScalarSpec, typename ParamSpec = ScalarSpec>
class PerezMaltaCoutinho{
public:
	typedef ScalarSpec ScalarType;
	typedef unsigned int size_type;
	typedef ScalarType RealType;
	typedef ParamSpec ParamType;
	typedef capd::vectalg::Vector<ParamSpec, 0> ParamsVectorType;
	typedef capd::vectalg::Vector<ScalarType, 0> VectorType;

	PerezMaltaCoutinho(): gamma(1), b0(2), b1(6) {}
	PerezMaltaCoutinho(ParamType b1): gamma(1), b0(2), b1(b1) {}
	PerezMaltaCoutinho(ParamType gamma, ParamType b0, ParamType b1): gamma(gamma), b0(b0), b1(b1) {}
	PerezMaltaCoutinho(PerezMaltaCoutinho const & orig): gamma(orig.gamma), b0(orig.b0), b1(orig.b1) {}
	PerezMaltaCoutinho(capd::vectalg::Vector<ParamSpec, 0> const & params): gamma(params[0]), b0(params[1]), b1(params[2]) {}

	PerezMaltaCoutinho& operator=(PerezMaltaCoutinho const & orig){ b0 = orig.b0; gamma = orig.gamma; b1 = orig.b1; return *this; }

	/** output dimension of the internal map */
	static size_type imageDimension() { return 1; }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	static size_type dimension() { return 2; }
	/** number of parameters to fully configure equation */
	static size_type getParamsCount() { return 3; }

	/** Evaluate r.h.s of the equation. */
	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		fx[0] = -gamma * x[0] + x[1] * (b0 - b1 * x[1]);
	}

	static std::string show(){
		return "Perez Malta Coutinho model: $x'(t) = -\\gamma x(t) + (b_0 - b_1 x(t-\\tau)) x(t-\\tau)$";
	}

protected:
	ParamType gamma;
	ParamType b0, b1;
};



/**
 * A model map of what I expect from a Map R^m -> R^n
 * to be good for use with DDE codes.
 *
 */
template<typename ScalarSpec, typename ParamSpec = ScalarSpec>
class CubicIkeda{
public:
	typedef ScalarSpec ScalarType;
	typedef unsigned int size_type;
	typedef ScalarType RealType;
	typedef ParamSpec ParamType;
	typedef capd::vectalg::Vector<ParamSpec, 0> ParamsVectorType;
	typedef capd::vectalg::Vector<ScalarType, 0> VectorType;

	CubicIkeda(): a(1.) {}
	CubicIkeda(ParamType a): a(a) {}
	CubicIkeda(CubicIkeda const & orig): a(orig.a) {}
	CubicIkeda(capd::vectalg::Vector<ParamSpec, 0> const & params): a(params[0]) {}

	CubicIkeda& operator=(CubicIkeda const & orig){ a = orig.a; return *this; }

	/** output dimension of the internal map */
	static size_type imageDimension() { return 1; }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	static size_type dimension() { return 2; }
	/** number of parameters to fully configure equation */
	static size_type getParamsCount() { return 1; }

	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		fx[0] = a * (x[1] - x[1]*x[1]*x[1]);
	}

	static std::string show(){
		return "Cubic Ikeda equation: $x'(t) = a(x(t - \\tau) + x(t-\\tau)^3)$ (a = 1 = classic).";
	}

protected:
	ParamType a;		/** classical parameter, default: 1 */
};



} // namespace ddes
} // namespace capd


#endif /* _CAPD_DDES_SAMPLEEQNS_H_ */
