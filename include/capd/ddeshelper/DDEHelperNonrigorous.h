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

#ifndef _CAPD_DDEHELPERNONRIGOROUS_H_
#define _CAPD_DDEHELPERNONRIGOROUS_H_

#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/DDEHelperCommon.h>
#include <capd/ddeshelper/DDECoordinateFrameHelper.h>

// TODO: (IMPORTANT) code below assume delaysSpec=1, it will not work for many delays, FIX!

// TODO: (NOT URGENT): we have NonrigorousHelper but it is in DDEHelperNonrigorous.h/.hpp,... refactor...?

namespace capd {
namespace ddeshelper {

/**
 * A helper class that made all the fuzzy boilerplate
 * for you, so you do not need to know various things by yourself.
 *
 * It provides functions like integrate and poincare
 * to compute images of the initial functions under your own
 * equations. You only need to supply the equation and the number
 * of delays in the system.
 *
 * The Equation (template @param EqSpec) must be a class of special form.
 * See programs/examples for usage and examples and/or include/capd/ddes/SampleEqns.h for more information.
 * A general form is as follows:
 *
 *    struct MyEquation{
 *    	   typedef unsigned int size_type;								  // this is required to know what the three static functions return. This is usually unsigned int. At the cost of some warnings, you can use int as well.
 *    	   typedef double ParamType;									  // the type to hold parameters. I usually make it a template parameter to be able to use in both rigorous and non-rigorous setting.
 *    	   typedef capd::vectalg::Vector<ParamType, 0> ParamsVectorType;  // this defines the type to hold a vector of parameters
 *
 *    	   MyEquation(ParamsVectorType const & params); // we require that the equation has a canstructor that gets parameters vector.
 *    	   MyEquation();						// default constructor, it must be defined, as we have another constructor (with params vector). All other constructors are entirely up to you.
 *    	   static size_type imageDimension();   // output dimension of the internal map (usually, the dimension d of the space for x(t) variable)
 *    	   static size_type dimension(); 		// input dimension of the internal map (usually, imageDimension() * (number_of_delays+1))
 *    	   static size_type getParamsCount();   // number of parameters to fully configure the equation, but without delays! Depends on the equation, e.g. Mackey-Glass has 3 parameters, except the delay.
 *    	   static std::string show();			// this is used in debug output, it can return empty string if you wish.
 *
 *    	   // this is the ost important part
 *    	   // it takes t -time, x - the collection of values at 0 and consecutive delays, then
 *    	   // it should compute the value of R.h.s f(x(0), x(t_1), ... ) in fx.
 *    	   // the dimension of x is given by dimension(), the
 *    	   // dimension of fx is given by imageDimension();
 *    	   template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
 *    	   void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const;
 *    };
 *
 * You might also supply Matrix and Vector type, but the dfault versions works ok.
 *
 * As of now, it can only handle discrete delay differential equations in the form of
 *
 * x'(t) = f(x(t), x(t-\tau_1), ... , x(t-tau_m))
 *
 * where tau_i's matches time points of the integration grid.
 *
 * C1-versions allows to do nonrigorous computations of the variational equation associated to the
 * original equation, allowing e.g. to apply Newton method to find better candidates.
 *
 * TODO: (IMPORTANT) Currently it only support one delay!!!!!!
 * TODO: (URGENT?) DEV: make sure all const qualifiers are put in the right places!
 * TODO: (FUTURE) move implementation to hpp file to shorten this for end user
 * TODO: (FUTURE): the normal and C^1 versions are totaly separated. Rethink if they can reuse some code, because there is a lot of repetition...
 */
template<
	typename EqSpec,
	int delaysSpec=1,
	typename MatrixSpec = capd::vectalg::Matrix<typename EqSpec::ParamType, 0, 0>,
	typename VectorSpec = capd::vectalg::Vector<typename EqSpec::ParamType, 0>
>
class NonrigorousHelper{
public:
	typedef EqSpec Eq;
	typedef typename Eq::ParamType ParamType;
	typedef VectorSpec Vector;
	typedef MatrixSpec Matrix;
	typedef typename Matrix::ScalarType Scalar;
	typedef typename Matrix::ScalarType Real; // TODO: (FUTURE) Rethink? What if scalar is Complex?
	typedef typename Eq::ParamsVectorType ParamsVector;
	typedef capd::ddes::DiscreteTimeGrid<Real> Grid;
	typedef typename Grid::TimePointType TimePoint;
	typedef capd::ddes::GenericJet<TimePoint, Vector, Vector, Matrix> Jet;
	typedef capd::ddes::GenericJet<TimePoint, capd::ddes::VectorWithJacData<Vector, Matrix>, Vector, Matrix> C1Jet;
	typedef capd::ddes::DDEPiecewisePolynomialCurve<Grid, Jet> Solution;
	typedef capd::ddes::DDEPiecewisePolynomialCurve<Grid, C1Jet> C1Solution;
	typedef Jet CurvePiece;
	typedef capd::ddes::BasicDiscreteDelaysFunctionalMap<Eq, C1Solution> C1DDEq;
	typedef capd::ddes::BasicDiscreteDelaysFunctionalMap<Eq, Solution> DDEq;
	typedef capd::ddes::DDENonrigorousTaylorSolver<C1DDEq> C1Solver;
	typedef capd::ddes::DDENonrigorousTaylorSolver<DDEq> Solver;
	typedef typename C1Solver::VariableStorageType Variables;
	typedef typename C1Solver::JacobianStorageType Jacobians;
	typedef typename C1Solver::ValueStorageType Values;
	typedef typename C1Solver::size_type size_type;
	typedef int step_type; // we need to have signed values here
	typedef capd::ddes::DDEJetSection<C1Solution> C1JetSection;
	typedef typename C1JetSection::JetType C1SecJet;
	typedef capd::ddes::DDEBasicPoincareMap<C1Solver, C1JetSection> C1PoincareMap;
	typedef capd::ddes::DDEJetSection<Solution> JetSection;
	typedef capd::ddes::DDEBasicPoincareMap<Solver, JetSection> PoincareMap;
	static const size_type PARAMS_COUNT;
	static const size_type DIMENSION;

	typedef capd::ddeshelper::CoordinateFrame<Matrix, Vector, Scalar> CoordinateFrame;

	/**
	 * you need to call at least one of the constructors and then use the setup variable
	 *
	 * This constructor takes parameters vector:
	 *
	 * (q_1, ..., q_k, tau_1, ... tau_m)
	 *
	 * where  q_i's are the parameters of the equation (r.h.s f : \R^{d*(m+1)} \to \R^d of the DDE)
	 * and tau_i's are the delays, so that the equation is
	 *
	 * x'(t) = f(x(t), x(t-tau_1), ...)
	 *
	 * p, n - are parameters of the (p,n)-representation of the functions.
	 * p is the number of grid points in the basic interval [-tau, 0], tau = max(tau_i)
	 * n is the order of the jet at each point.
	 * See papers for more details.
	 */
	NonrigorousHelper(ParamsVector params, size_type p, size_type n, step_type reqSteps=0, step_type maxSteps=0, size_type maxOrder=10):
			m_params(params),
			m_p(p), m_n(n),
			m_reqSteps(reqSteps), m_maxSteps(maxSteps), m_maxOrder(maxOrder),
			crossingDirection(capd::poincare::MinusPlus),
			m_experimentalRenormalizeVariational(false)
	{ updateGrid(); updateSteps(); }

	/** similar to the first one, but reads the parameters and p and n from a given file. */
	NonrigorousHelper(std::string filepath, step_type reqSteps=0, step_type maxSteps=0, size_type maxOrder=10):
			m_params(PARAMS_COUNT),
			m_p(1), m_n(1),
			m_reqSteps(reqSteps), m_maxSteps(maxSteps), m_maxOrder(maxOrder),
			crossingDirection(capd::poincare::MinusPlus),
			m_experimentalRenormalizeVariational(false)
	{ loadSetup(filepath); }

	/** similar to the second one, but user gives istream instead of filepath */
	NonrigorousHelper(std::istream& input, step_type reqSteps=0, step_type maxSteps=0, size_type maxOrder=10):
			m_params(PARAMS_COUNT),
			m_p(1), m_n(1),
			m_reqSteps(reqSteps), m_maxSteps(maxSteps), m_maxOrder(maxOrder),
			crossingDirection(capd::poincare::MinusPlus),
			m_experimentalRenormalizeVariational(false)
	{
		rawLoadSetup(input, m_p, m_n, m_params);
		updateGrid();
		updateSteps();
	}

	/**
	 * this controls how many steps are needed before the section crossing
	 * would be detected when computing poincare maps. This is important
	 * more in the rigorous setting, but to have consistent computations you might
	 * want to keep it the same. The helper will set the right (long enough as in papers)
	 * number of steps for you. Please make it at least p steps if you want to change it.
	 */
	void setRequiredSteps(step_type reqSteps, bool control_steps = true){ m_reqSteps = reqSteps; updateSteps(control_steps); }
	/** sets the maximum allowed steps. This prevents the infinite loops when the solution goes astray */
	void setMaximumSteps(step_type maxSteps, bool control_steps = true){ m_maxSteps = maxSteps; updateSteps(control_steps); }
	/**
	 * sets the maximum allowed order of the expanded representation as described in the 2023 FOCM paper.
	 * it can be as large as you want. But please remember that the expansions occur by one every p steps (full delay).
	 * Also, the size makes computations expensive very fast.
	 */
	void setMaximumOrder(step_type maxOrder, bool control_steps = true){ m_maxOrder = maxOrder; }
	/** @see docs for setter */
	step_type getRequiredSteps() const { return m_reqSteps; }
	/** @see docs for setter */
	step_type getMaximumSteps() const { return m_maxSteps; }
	/** returns  */
	step_type getMaximumOrder() const { return m_maxOrder; }

	/** returns the parameters of the r.h.s of the DDE set for this helper */
	ParamsVector params() const { return m_params; }
	/** the intrinsic dimension of the (d,p,n)-representation. This is the constant M = M(d, p, n) from papers. It is usually a big number = O(d*p*n) */
	size_type M() const { return DIMENSION * (1 + m_p * (m_n +1)); }
	/** the grid size, i.e. into how many subintervals we divide basic delay */
	size_type p() const { return m_p; }
	/** the order of the representation of the solutions used in this helper */
	size_type n() const { return m_n; }
	/** dimension of the 'ambient' space, i.e. the dimension of x(t) for a fixed t */
	size_type d() const { return DIMENSION; }
	/** returns basic grid size step h := tau/p */
	Real h() const { return getBasicIntervalLength() / p(); }
	/** returns the basic interval (the longest one) */
	TimePoint tau() const { return m_grid(p()); }
	/** returns the size of the longest interval but as a Real value, not TimePoint */
	Real getBasicIntervalLength() const { return m_params[m_params.dimension()-1]; /* TODO: (IMPORTANT): make it more general, many delays */ }
	/** returns the i-th grid point, where e.g. t(-p) = -tau, t(0) = 0 */
	TimePoint t(int i) const { return m_grid(i); }
	/** returns the internal grid used by this helper */
	const Grid& grid() const { return m_grid; }

	/** sets the crossing directions for PoincareMaps used in this helper subroutines */
	void setCrossingDirection(capd::poincare::CrossingDirection const& d) { this->crossingDirection = d; }

	/** This creates an object representing the equation for the current value of parameters set in the helper. */
	DDEq makeEquation(){ return makeEquationTemplate<DDEq>(); }
	/** T@see the other. This is to be used with C1*** types.  */
	C1DDEq makeC1Equation(){ return makeEquationTemplate<C1DDEq>(); }

	/** sets the value of the current params (but not the delays!) */
	void setParams(ParamsVector const& new_params) {
		int param_count = new_params.dimension();
		if (m_params.dimension() - delaysSpec < param_count)
			param_count = m_params.dimension() - delaysSpec < param_count;
		for (int i = 0; i < param_count; i++){
			m_params[i] = new_params[i];
		}
	}

	/** sets the value of some parameter. Can't change value of the delay! Returns the previous value of the param. */
	ParamType setParam(size_type i, ParamType const& value) {
		if (i < m_params.dimension() - delaysSpec){
			ParamType prev_val = m_params[i];
			m_params[i] = value;
			return prev_val;
		} else {
			throw std::logic_error("NonrigorousHelper::setParam(): cannot change value of the delays after creation!");
		}
	}

	/**
	 * This creates a raw solver, if you want to do nonstandard tasks and have an Equation made on the side.
	 * If you want just to integrate initial values, consider using iterate() or poincare() instead.
	 */
	Solver makeSolver(DDEq const& eq){ return Solver(eq, m_maxOrder); }
	/** It uses the other version with the current equation (@see makeEquation()). */
	Solver makeSolver(){ return makeSolver(makeEquation()); }
	/**
	 * T@see makeSolver().
	 *
	 * This version can be used to compute monodromy matrix along the solution.
	 */
	C1Solver makeC1Solver(C1DDEq const& eq){ return C1Solver(eq, m_maxOrder); }
	/** It uses the other version with the current equation (@see makeEquation()). */
	C1Solver makeC1Solver(){ return makeC1Solver(makeC1Equation()); }

	/**
	 * This creates a PoincareMap.
	 *
	 * Warning: the solver and section are passed as references, so they need
	 * to be non-temporary variables. You cannot do makePoincare(makeSolver(), makeSection());
	 * You need to have:
	 * auto solver = makeSolver();
	 * auto section = makeSection();
	 * auto pm = makePoincareMap(solver, section);
	 */
	PoincareMap makePoincareMap(Solver& solver, JetSection& section){
		PoincareMap pm(solver, section);
		pm.setDirection(this->crossingDirection);
		pm.setRequiredSteps(m_reqSteps);
		pm.setMaxSteps(m_maxSteps);
		return pm;
	}

	/**
	 * This creates a PoincareMap. See the other version for more docs.
	 *
	 * This version creates a Poincare map that can also compute Jacobian and Monodromy matrix (approximations).
	 * Those are computationally expensive, so if you do not need them, then the other version is much faster.
	 */
	C1PoincareMap makeC1PoincareMap(C1Solver& solver, C1JetSection& section){
		C1PoincareMap pm(solver, section);
		pm.setDirection(this->crossingDirection);
		pm.setRequiredSteps(m_reqSteps);
		pm.setMaxSteps(m_maxSteps);
		return pm;
	}

	/**
	 * makes a segment of solution over [-tau, 0] of order n,
	 * with a given value. If the dim(v) == d (dimension of the 'ambient' space),
	 * then the constant solution will be made. If dim(v) = M(), then
	 * a solution with the given data is produced.
	 *
	 * Without the v parameter, it creates a zero solution (v = {0,..,0}).
	 * Without n, it creates solution with order m_n as set in this Helper.
	 */
	Solution makeSegment(int n, Vector v=Vector({})) const { return makeSegmentTemplate<Solution>(n, v); }
	/** just a rename of makeSegment, for backward compatibility. */
	Solution makeSolution(int n, Vector v=Vector({})) const { return makeSegment(n, v); }
	/** Creates a constant value segment over [-tau, 0] with representation order n as set in this helper. See other @method makeSegment() */
	Solution makeSegment(Vector v=Vector({})) const { return makeSegment(m_n, v); }
	/** just a rename of makeSegment, for backward compatibility. */
	Solution makeSolution(Vector v=Vector({})) const { return makeSegment(m_n, v); }
	/** @deprecated use makeSegment() instead! */
	Solution vectorToSolution(Vector const& x) const { return makeSegment(x); }
	/**
	 * This version is suitable for approximate C^1 computations - computing Jacobian of a Poincare/Time map.
	 *
	 * makes a segment of solution over [-tau, 0] of order n,
	 * with a given value. If the dim(v) == d (dimension of the 'ambient' space),
	 * then the constant solution will be made. If dim(v) = M(), then
	 * a solution with the given data is produced.
	 *
	 * Without the v parameter, it creates a zero solution (v = {0,..,0}).
	 * Without n, it creates solution with order m_n as set in this Helper.
	 */
	C1Solution makeC1Segment(int n, Vector v=Vector({})) const { return makeSegmentTemplate<C1Solution>(n, v); }
	/** just a rename of makeSegment, for backward compatibility. */
	C1Solution makeC1Solution(int n, Vector v=Vector({})) const { return makeC1Segment(n, v); }
	/** Creates a segment over [-tau, 0] with order n as set in this helper. See other @method makeSegment() */
	C1Solution makeC1Segment(Vector v=Vector({})) const { return makeC1Segment(m_n, v); }
	/** just a rename of makeSegment, for backward compatibility. */
	C1Solution makeC1Solution(Vector v=Vector({})) const { return makeC1Segment(m_n, v); }
	/** @deprecated use makeJacSegment() instead! */
	C1Solution vectorToC1Solution(Vector const& x) const{ return makeC1Segment(x); }
	/**
	 * makes a section and stores it in out_section
	 *
	 * if dim(s) == d, then we setup section in the coordinate x(0) . s = c
	 * otherwise we set full-space section s . x_0 = c
	 */
	JetSection makeSection(Vector const& s, Scalar const& c){ return makeSectionTemplate<JetSection>(s, c); }
	/**
	 * same as makeSection(), but returns a datatype needed for computation of poincare map
	 * with monodromy and jacobian of the poincare map, @see the relevant makePoincareMap().
	 */
	C1JetSection makeC1Section(Vector const& s, Scalar const& c){ return makeSectionTemplate<C1JetSection>(s, c); }

	/**
	 * Just integrate the solution.
	 *
	 * Integrates for time T = iters * tau_max / p = iters * h.
	 * The value h is the step size of the method / grid
	 *
     * If @param iter < 0 then it will be converted to full delays! The DDEs cannot be integrated backward in time!
	 *
	 * @return the solution on the whole integration time.
	 *
	 * The initial solution segment must have the same grid as the Helper, otherwise the exception
	 * will be thrown. It is best to make a solution with makeSolution() methods of the helper.
	 *
	 * If @param use_extension is set to true, then the more accurate results will be produced due to use
	 * of extension algorithm described in FoCM 2023 paper. But then, the returned full Solution
	 * will have varying degree of the representation at various jets.
	 *
	 * In @param return it returns the last segment of the solution in the @param result. The structure
	 * of result will be preserved, i.e. if extension is used, and the structure of result is the same as
	 * initial, then the result still have the same structure as initial. This might be useful for
	 * iterating maps.
	 *
	 * TODO: (FUTURE, IMPORTANT?) As in CAPD, create Timemap class to handle this in general...
	 */
	Solution integrate(int iters, const Solution& initial, Solution& result, bool use_extension=false){
		try {
			checkGrid(initial); checkGrid(result);
		} catch (std::logic_error& e){
			throw capd::ddes::rethrow("NonrigorousHelper::integrate(int, Solution, Solution, bool): initial and result must have the same grid!", e);
		}
		if (iters < 0) iters = -m_p * iters;
		auto tau = m_grid.point(m_p);
		DDEq dde(Eq(m_params), tau);
		Solver solver(dde, use_extension ? m_maxOrder : m_n);
		Solution XXX = initial;
		for (int i = 0; i < iters; i++)
			solver(XXX);
		 Solution X = XXX.subcurve(XXX.currentTime() - tau);
		 capd::ddeshelper::copyReduce(X, result);
		 return XXX;
	}
	/** @see integrate(int, const Solution&, Solution&, bool) */
	Solution integrate(int iters, const Solution& initial, bool use_extension=false){
		try{
			Solution dump = initial;
			return integrate(iters, initial, dump, use_extension);
		} catch (std::logic_error &e) {
			throw capd::ddes::rethrow("NonrigorousHelper::integrate(int, Solution, bool):", e);
		}
	}

	/**
	 * just integrate the solution and return the solution
	 *
	 * initial can be either of (1) d-dimension (\R^d) or (2) the M(d,p,n)-dimension (see papers)
	 * In the (1) case it integrates initial function x_0(s) == initial for all s \in [-tau, 0]
	 * In the (2) case it constructs initial function from data stored in x0 (the order of coefficients
	 * in this vector is described elsewhere (see papers for example). But suppling it by hand
	 * is very cumbersome. It is best to use output of other functions o get initials in this
	 * form.
	 *
	 * Integrates for time T = iters * tau_max / p = iters * h. The value h is the step size of the method / grid
	 *
	 * In the result the last segment x_T from the solution in the vector format.
	 * Can be used as an input to the next integrate() procedure, to create an initial solution
	 * curve, or in other functions from the helper (e.g. poincare).
	 *
	 * @deprecated It is better to use function based on Solution parameters
	 */
	Solution integrate(int iters, const Vector& initial, Vector& result, bool use_extension=false){
		try{
			Solution X = makeSegment(initial);
			Solution TX = X;
			Solution XXX = integrate(iters, X, TX, use_extension);
			result = TX;
			return XXX;
		} catch (std::logic_error &e) {
			throw capd::ddes::rethrow("NonrigorousHelper::integrate(int, Vector, Vector, bool):", e);
		}
	}
	/**
	 * just integrate the solution and return the solution
	 *
	 * initial can be either of (1) d-dimension (\R^d) or (2) the M(d,p,n)-dimension (see papers)
	 * In the (1) case it integrates initial function x_0(s) == initial for all s \in [-tau, 0]
	 * In the (2) case it constructs initial function from data stored in x0 (the order of coefficients
	 * in this vector is described elsewhere (see papers for example). But suppling it by hand
	 * is extremely cumbersome. It is best to use output of other functions o get initials in this
	 * form.
	 *
	 * @deprecated It is better to use function based on Solution parameters
	 */
	Solution integrate(int iters, const Vector& initial, bool use_extension=false){
		try{
			Vector dump;
			return integrate(iters, initial, dump, use_extension);
		} catch (std::logic_error &e) {
			throw capd::ddes::rethrow("NonrigorousHelper::integrate(int, Vector, bool):", e);
		}
	}
	/**
	 * integrate, save output to files and draw the solution. See other versions
	 * for more information on input and output.
	 *
	 * The variable outconfig holds a pair of strings that define paths to where store
	 * the results of the computations. See documentation there.
	 *
	 * @deprecated For backward compatibility...
	 */
	Solution integrate(int iters, const Solution& initial, Solution& result, PathConfig const& outconfig){
		Solution solution = integrate(iters, initial, result);
		drawMap(solution, outconfig);
		return solution;
	}
	/** @deprecated For backward compatibility... */
	Solution integrate(int iters, const Vector& initial, Vector& result, PathConfig const& outconfig){
		Solution solution = integrate(iters, initial, result);
		drawMap(solution, outconfig);
		return solution;
	}
//	/**
//	 * Detects the crossing direction of the section just by integration.
//	 * Can be helpful, when section is weird and you cannot decide.
//
//	TODO: (IMPORTANT): finish it...
//	 */
//	capd::poincare::CrossingDirection detectCorssingDirection(JacJetSection section, Vector const& x){
//		return detectCorssingDirectionTemplate<JacJetSection>(section, x)
//		auto tau = m_grid.point(m_p);
//		auto t_0 = m_grid.point(0);
//		JacSolution X(m_grid, -tau, t_0, m_n, Vector(DIMENSION));
//		X.set_x(x);
//		JacDDEq dde(Eq(m_params), tau);
//		JacSolver solver(dde, m_maxOrder);
//		JacPoincareMap pm(solver, section);
//		return pm.detectCrossingDirection(X);
//	}

	/**
	 * Computes a poincare map, starting from point x, until section s is crossed
	 * in the direction set for this Helper. The result will be in Px.
	 * The (p, \eta)-structure of Px is preserved (see papers, references).
	 * the time to reach section is in reachTime.
	 *
	 * The reachTime will be at least base grid step h * requiredSteps, which are set for this Helper.
     */
	Solution poincare(
				JetSection section, Solution const& x,
				double& reachTime, Solution& Px){
		try {
			checkGrid(x); checkGrid(Px);
		} catch (std::logic_error& e){
			throw capd::ddes::rethrow("NonrigorousHelper::poincare(Section, Solution, double, Solution): initial and result must have the same grid!", e);
		}

		Solution X = x;

		Solver solver = makeSolver();
		PoincareMap pm = makePoincareMap(solver, section);

		pm(X, Px, reachTime);
		return X; // X contains full trajectory
	}
	/** @see poincare(JetSection, Solution const&, double, Solution&) */
	Solution poincare(JetSection section, Solution const& x, Solution& Px){
		double dumpTime;
		return poincare(section, x, dumpTime, Px);
	}
	/**
	 * The simples computation of poincare map, without extra data.
	 * This version operates on raw vectors (they need to comply
	 * with the dimension M() of the Helper.
	 */
	Solution poincare(
				JetSection section, Vector const& x,
				double& reachTime, Vector& Px){
		Solution X = makeSegment(x);
		Solution PX = X;

		auto trajectory = poincare(section, X, reachTime, PX);
		Px = PX.get_x();

		return trajectory;
	}
	/** the simplest computation of poincare map, without extra data (even without return time) */
	Solution poincare(JetSection section, Vector const& x, Vector& Px){
		double dumpReachTime;
		return poincare(section, x, dumpReachTime, Px); // X contains full trajectory
	}

	/**
	 * if you do not know what initV is, then
	 * you should probably stick to using
	 * the other C1Solution poincare() method (without initV).
	 */
	C1Solution poincare(
				C1JetSection section,
				Vector const& x, Matrix& initV,
				double& reachTime, int& steps, Vector& Px, Vector& fPx,
				Matrix& V, Matrix& DP){

//		// TODO: (NOT URGENT) refactor it to be more DRY.
//		auto tau = m_grid.point(m_p);
//		auto t_0 = m_grid.point(0);
//		C1Solution X(m_grid, -tau, t_0, m_n, Vector(DIMENSION));
//		X.set_x(x);
//
//		C1DDEq dde(Eq(m_params), tau);
//		C1Solver solver(dde, m_maxOrder);
//
//		C1PoincareMap pm(solver, section);
//		pm.setDirection(this->crossingDirection);
//		pm.setRequiredSteps(m_reqSteps);
//		pm.setMaxSteps(m_maxSteps);
//
//		capd::vectalg::EuclLNorm<Vector, Matrix> euclNorm;
//		Px = Vector(M()); fPx = Vector(M());
//		V = Matrix(M(), M()); DP = Matrix(M(), M());
//		Matrix Id(X.storageDimension(), X.storageDimension()); Id.setToIdentity();
//		pm.setNormalizeVariational(m_experimentalRenormalizeVariational);
//		pm.setInitialV(X, initV);
//		C1Solution PX(X); PX *= 0.;
//		Vector dump = x;
//		pm(X, PX, reachTime, dump, Px, fPx, V, DP);
//		steps = pm.getLastStepsAfterSection();
//		return X; // X contains full trajectory


		C1Solution X0 = makeC1Segment(x);
		C1Solution PX = X0; PX *= 0.;
		C1Solution trajectory = poincare(section, X0, initV, reachTime, steps, PX, fPx, V, DP);
		Px = Vector(PX);
		return trajectory; // returns full trajectory from the beginning
	}
	/**
	 * computes poincare map and the (approximate)
	 * solution to the variational equation on the coefficients
	 * in V it returns the monodromy matrix, that is D_x \varph(t_P(x0), x)
	 * in DP you get D \varphi(t_p(x), x)
	 * (sic! the difference in t_p argument)
	 * in fPx you get the value of the "vector field" at P(x0)
	 * this is not straightforward as in ODE
	 * that is for ODE you have fPx = f(P(x0)), where f is r.h.s. of the ODE.
	 * But for DDE you do not have r.h.s for all of the points in [-tau, 0]
	 * But it can be shown, that if $t_p \ge tau$ then $fPx = (Px)'$
	 * (note Px : [-tau,0] \to \R^d is a function of time, so it can be differentiated)
	 */
	C1Solution poincare(
				C1JetSection section, Vector const& x0,
				double& reachTime, Vector& Px, Vector& fPx,
				Matrix& V, Matrix& DP){
		Matrix Id(x0.dimension(), x0.dimension()); Id.setToIdentity();
		int steps;
		return 	poincare(section, x0, Id, reachTime, steps, Px, fPx, V, DP);
	}

	/**
	 * if you do not know what initV is, then
	 * you should probably stick to using
	 * the other C1Solution poincare() method (without initV).
	 */
	C1Solution poincare(
				C1JetSection section,
				C1Solution const& X0, Matrix& initV,
				double& reachTime, int& steps, C1Solution& PX, Vector& fPX,
				Matrix& V, Matrix& DP){

		C1Solution X(X0); // make a copy for modification
		size_type M = X.storageDimension();

		C1Solver solver = makeC1Solver();
		C1PoincareMap pm = makeC1PoincareMap(solver, section);
		pm.setDirection(this->crossingDirection);
		pm.setRequiredSteps(m_reqSteps);
		pm.setMaxSteps(m_maxSteps);

		capd::vectalg::EuclLNorm<Vector, Matrix> euclNorm;
		Vector x(M), Px(M); fPX = Vector(M);
		V = Matrix(M, M); DP = Matrix(M, M);
		// pm.setNormalizeVariational(m_experimentalRenormalizeVariational);
		pm.setInitialV(X, initV);
		PX = X; PX *= 0.;
		pm(X, PX, reachTime, x, Px, fPX, V, DP);
		steps = pm.getLastStepsAfterSection();
		return X; // X contains full trajectory
	}

	/**
	 * TODO: (URGENT!!!) docs
	 */
	C1Solution poincare(
				C1JetSection section,
				C1Solution const& X0,
				double& reachTime,
				C1Solution& PX,
				Matrix& V, Matrix& DP){
		size_type M = X0.storageDimension();
		Matrix Id(M, M); Id.setToIdentity();
		int dump_steps; Vector dump_fPx(M);
		return poincare(section, X0, Id, reachTime, dump_steps, PX, dump_fPx, V, DP);
	}

	/**
	 * TODO: (URGENT!!!) docs
	 */
	C1Solution poincare(
				C1JetSection section,
				C1Solution const& X0,
				C1Solution& PX,
				Matrix& V, Matrix& DP){
		double dump_time;
		return poincare(section, X0, dump_time, PX, V, DP);
	}

	/**
	 * TODO: (URGENT!!!) docs
	 */
	C1Solution poincare(
				C1JetSection section, C1Solution const& X0,
				C1Solution& PX, Matrix& DP){
		size_type M = X0.storageDimension();
		Real dump_reachTime; Matrix dump_V(M, M);
		return poincare(section, X0, dump_reachTime, PX, dump_V, DP);
	}

//	/**
//	 * TODO: (URGENT!!!) docs!
//	 */
//	JacSolution poincare(
//				JetSection section, Solution const& x0,
//				double& reachTime, Solution& Px, Vector& fPx,
//				Matrix& V, Matrix& DP){
//		Matrix Id(x0.dimension(), x0.dimension()); Id.setToIdentity();
//		int steps;
//		JacJetSection jac_section = makeJacSection(section.get_s(), section.get_c());
//		return 	poincare(section, x0, Id, reachTime, steps, Px, fPx, V, DP);
//	}

	/**
	 * Helper function to refine a candidate periodic orbit with a Newton method.
	 *
	 * First parameter is to get some text info back during the process, you can pass std::cout there or some std::ostringstream to ignore it.
	 *
	 * For technical reasons section must be of type C1JetSection, not JetSection.
	 *
	 * The V and DP matrices are the monodromy and Jacobian of Poincare map, respectively.
	 * They might be helpful in finding good candidate for section coordinates.
	 * @see periodicCoordinates(), where an example procedure described in 2023 FOCM paper is implemented.
	 *
	 * It returns the difference between candidate x and P(x) (in basic euclidean norm on the vector representations!).
	 *
	 * You should iterate the method for desired number of iterates or until desired accuracy is obtained.
	 *
	 */
	Real refinePeriodic(std::ostream& info, C1JetSection& section, Vector& x, Matrix& V, Matrix&DP){
		info << "# refinePeriodic START" << std::endl;

		auto X = makeC1Segment(x);
		auto solver = makeC1Solver();
		auto pm = makeC1PoincareMap(solver, section);

		auto vsec = section(X);
		info << "# section(X) = " << vsec << " (should be small)" << std::endl;

		capd::vectalg::EuclLNorm<Vector, Matrix> euclNorm;
		Real reachTime;
		Vector Px, fPx, dx(X.storageDimension());
		Matrix Id(X.storageDimension(), X.storageDimension()); Id.setToIdentity();
		pm.setInitialV(X);
		C1Solution PX(X); PX *= 0.;
		pm(X, PX, reachTime, x, Px, fPx, V, DP);

		double diff = euclNorm(Px-x);
		info << "# reachTime = " << reachTime << std::endl;
		info << "# ||PX - X|| = " << diff << std::endl;
		info << "# Solving for Newton..." << std::flush;
		capd::matrixAlgorithms::gauss(DP-Id, x-Px, dx);
		info << "# DONE!" << std::endl;
		x += dx;
		info << "# refinePeriodic END" << std::endl;
		return diff;
	}
	/**
	 * refines, by one Newton iteration, the X if it is a candidate for periodic orbit on the given section.
	 *
	 * If X is not on the section or far from the true solution, then the method might fail.
	 *
	 * You should iterate the method for desired number of iterates or until desired accuracy is obtained.
	 *
	 * The V is the Monodromy matrix and DP is the Jacobian of the Poincare Map. They can be used to obtain good coordinates on the section.
	 *
	 * TODO: (FUTURE) this will only work if the segment X is of uniform order n(). This should be corrected to work for any...
	 */
	Real refinePeriodic(std::ostream& info, JetSection& section, Solution& X, Matrix& V, Matrix&DP){
		auto x = X.get_x();
		auto jacSection = makeC1Section(Vector(section), section.get_c());
		auto diff = refinePeriodic(info, jacSection, x, V, DP);
		X.set_x(x);
		return diff;
	}

	// TODO: (FUTURE) copy / implement
	// void periodicCoordinates(std::ostream& info, Matrix const& V, Matrix const& DP, Matrix& coords);

	/** TODO: (URGENT): docs! */
	template<typename VectorIteratorSpec>
	void saveData(std::string filepath, VectorIteratorSpec start, VectorIteratorSpec end, std::string extraComment=""){
		std::ofstream outf(filepath); outf.precision(15);
		outf << PARAMS_COUNT << " " << m_params << std::endl;
		outf << DIMENSION << " " << m_p << " " << m_n << std::endl;
		outf << int(end - start) << std::endl;
		for (; start != end; ++start)
			outf << *start << std::endl;
		outf << "# INFO on format: " << std::endl;
		outf << "# first line:      number_of_params (" << PARAMS_COUNT << "), then params." << std::endl; // TODO: doac zeby rownania potrafily mowic jakie maja prarametry
		outf << "# second line:     d p n, as in (p,n)-representations." << std::endl;
		outf << "# third line:      count, the number of data for solutions. How the programs will use solutions, it depends on the program." << std::endl;
		outf << "# following lines: the solutions that fit into (p,n)-representation. " << std::endl;
		outf << std::endl << extraComment << std::endl;
		outf.close();
	}

	/** TODO: (URGENT): docs! */
	void loadSetup(std::string filepath);
	/** load data only if compatible with setup */
	std::vector<Vector> loadData(std::string filepath){
		size_type dp, dn; Vector dumppar(PARAMS_COUNT);
		std::ifstream in(filepath);
		rawLoadSetup(in, dp, dn, dumppar);
		if (dp != m_p || dn != m_n)
			throw std::logic_error("NonrigorousHelper::loadData(): incompatible dimensions.");
		int count = 0;
		in >> count;
		std::vector<Vector> result(count, Vector(M()));
		if (in.good() && count > 0){
			for (int i = 0; i < count; ++i)
				in >> result[i];
		}
		in.close();
		return result;
	}

	void drawMap(Solution const& solution, PathConfig const& outconfig){
		drawSolution(outconfig.dirPath, outconfig.prefix, solution);
		Vector start = solution.subcurve(solution.pastTime(), m_grid(0)).get_x();
		Vector iterated = solution.subcurve(solution.currentTime() - m_grid(m_p)).get_x();
		saveData(outconfig.filepath("start"), &start, &start+1);
		saveData(outconfig.filepath("iterated"), &iterated, &iterated+1);
	}
	// TODO: (FUTURE): Why those are templated?
	// TODO: (FUTURE): some code repeats what is in plot_value), consider refactor DRY
	// TODO: (FUTURE): the code doeas not test for system() or respect DDES_ALLOW_SYSTEM flag.
	template<typename SolutionT>
	std::string drawSolution(std::string dirpath, std::string filename, SolutionT const& X, double tshift = 0.) const {
		{
			std::ostringstream cmd;
			cmd << "mkdir -p '" << dirpath << "' ";
			capd::ddeshelper::runSystemCommand(cmd.str());
		}
		std::ofstream outfX(dirpath + "/" + filename + ".dat"); outfX.precision(15);
		capd::ddeshelper::value_to_gnuplot(outfX, X.pastTime(), X.currentTime(), X.getStep(), X);
		outfX.close();
		std::ofstream outg(dirpath + "/" + filename + ".gp");
		outg << "set terminal png size 800,600" << std::endl;
		outg << "set output '" << filename << ".png" << "'" << std::endl;
		outg << "plot ";
		std::ostringstream plotcmd;
		if (tshift == 0.)
			plotcmd << "'" << filename << ".dat" << "' using 1:3 with lines"; // 1, 3 because we have thickness 0 as the 2 and 4 column, see value_to_gnuplot()
		else
			plotcmd << "'" << filename << ".dat" << "' using ($1+" << tshift <<"):3 with lines"; // same here as above
		outg << plotcmd.str();
		outg.close();
		{
			std::ostringstream cmd;
			cmd << "cd '" << dirpath << "' && gnuplot '" << filename << ".gp'";
			capd::ddeshelper::runSystemCommand(cmd.str());
		}
		return plotcmd.str();
	}

	std::string drawSolution(std::string dirpath, std::string filename, Vector const& x, double tshift = 0.) const {
		if (x.dimension() != M())
			throw std::logic_error("NonrigorousHelper::drawSolution(Vector): bad x dimension");
		auto tau = m_grid.point(m_p);
		auto t_0 = m_grid.point(0);
		Vector v(DIMENSION);
		Solution X(m_grid, -tau, t_0, m_n, v);
		X.set_x(x);
		return drawSolution(dirpath, filename, X, tshift);
	}

	// TODO: (IMPORTANT): Why those are templated?
	template<typename SolutionT>
	void drawSolution(std::string filename, SolutionT const& X, double tshift = 0.) const {
		drawSolution(".", filename, X, tshift);
	}

	// TODO: (IMPORTANT): Why those are templated?
	template<typename SolutionT>
	void drawSolution(PathConfig const& paths, SolutionT const& X, double tshift = 0.) const {
		drawSolution(paths.dirpath(), paths.filename("solution"), X, tshift);
	}

	// TODO: (IMPORTANT): Why those are templated?
	template<typename SolutionT>
	void drawDelayMap(std::string dirpath, std::string filename, SolutionT const& X) const {
		{
			std::ostringstream cmd;
			cmd << "mkdir -p '" << dirpath << "' ";
			capd::ddeshelper::runSystemCommand(cmd.str());
		}
		std::ofstream outfX(dirpath + "/" + filename + ".dat"); outfX.precision(15);
		auto past = X.pastTime();
		auto present = m_grid(0);
		auto step = X.getStep();
		while (present <= X.currentTime()){
			outfX << X.eval(present)[0] << " " << X.eval(past)[0] << "\n";
			past += step; present += step;
		}
		outfX.close();
		std::ofstream outg(dirpath + "/" + filename + ".gp");
		outg << "set terminal png size 800,600" << std::endl;
		outg << "set output '" << filename << ".png" << "'" << std::endl;
		outg << "plot '" <<  filename << ".dat" << "' using 1:2 with lines";
		outg.close();
		{
			std::ostringstream cmd;
			cmd << "cd '" << dirpath << "' && gnuplot '" << filename << ".gp'";
			capd::ddeshelper::runSystemCommand(cmd.str().c_str());
		}
	}

	// TODO: (IMPORTANT): Why those are templated?
	template<typename SolutionT>
	void drawDelayMap(std::string filename, SolutionT const& X) const {
		drawDelayMap(".", filename, X);
	}

	// TODO: (IMPORTANT): Why those are templated?
	template<typename SolutionT>
	void drawDelayMap(PathConfig const& paths, SolutionT const& X) const {
		drawDelayMap(paths.dirpath(), paths.filename("phasespace"), X);
	}

	/** This is DEVELOPMENTAL, pleas do not use */
	NonrigorousHelper& setExperimentalRenormalizeVariational(bool v) { m_experimentalRenormalizeVariational = v; return *this; }

private:
	ParamsVector m_params;
	size_type m_p;
	size_type m_n;
	step_type m_reqSteps;
	step_type m_maxSteps;
	step_type m_maxOrder;
	Grid m_grid;
	capd::poincare::CrossingDirection crossingDirection; // TODO: rename m_... (cosistency)

	// experimental setups...
	bool m_experimentalRenormalizeVariational;

	void updateGrid() { m_grid = Grid(getBasicIntervalLength() / m_p); };
	void rawLoadSetup(std::istream& in, size_type &p, size_type &n, Vector &params);
	void updateSteps(bool control_steps = true) {
		m_reqSteps = (m_reqSteps < 0 ? -m_reqSteps * m_p : m_reqSteps);
		if (control_steps)
			m_reqSteps = (m_reqSteps < m_p ? m_p : m_reqSteps);
		m_maxSteps = (m_maxSteps < 0 ? -m_maxSteps * m_p : m_maxSteps);
		if (control_steps)
			m_maxSteps = (m_maxSteps < m_reqSteps ? 10 * m_reqSteps : m_maxSteps);
	}

	template<typename SegSpec>
	SegSpec makeSegmentTemplate(size_type n, Vector v) const {
		if (!v.dimension()) v = Vector(d());
		auto zero = m_grid.point(0);
		auto tau = m_grid.point(m_p);
		Vector v0(d());
		if (v.dimension() == d()) v0 = v;
		SegSpec X(m_grid, -tau, zero, n, v0);
		if (v.dimension() != d() && v.dimension() != Vector(X).dimension()){
			std::ostringstream info;
			info << "NonrigorousHelper::makeSegment/Solution(): dim(v) = " << v.dimension() << ", ";
			info << "should be either " << d() << " or " << Vector(X).dimension();
			throw std::logic_error(info.str());
		}
		if (v.dimension() != d()) X.set_x(v);
		return X;
	}

	template<typename SecionSpec>
	SecionSpec makeSectionTemplate(Vector const& s, Scalar const& c){
		Vector v;
		if (s.dimension() == d()){
			v = Vector(M());
			for (size_t i = 0; i < d(); ++i)
				v[i] = s[i];
		}else if (s.dimension() == M()){
			v = s;
		}else {
			throw std::logic_error("NonrigorousHelper::makeSection(): vector must be either DIMENSION (for const) or M dimensional (for any initial)");
		}
		return SecionSpec(d(), p(), n(), v, c);
	}

	template<typename DDEEqSpec>
	DDEEqSpec makeEquationTemplate(){
		auto tau = m_grid.point(m_p);
		DDEEqSpec dde(Eq(m_params), tau);
		return dde;
	}

	void checkGrid(Solution const& segment){
		if (segment.grid() != m_grid){
			throw std::logic_error("NonrigorousHelper::checkGrid(Solution): Solution grid not compatible with this helper!");
		}
	}
};

template<typename EqSpec, int delaysSpec, typename MatrixSpec, typename VectorSpec>
const typename NonrigorousHelper<EqSpec, delaysSpec, MatrixSpec, VectorSpec>::size_type
NonrigorousHelper<EqSpec, delaysSpec, MatrixSpec, VectorSpec>::PARAMS_COUNT = EqSpec::getParamsCount() + delaysSpec;

template<typename EqSpec, int delaysSpec, typename MatrixSpec, typename VectorSpec>
const typename NonrigorousHelper<EqSpec, delaysSpec, MatrixSpec, VectorSpec>::size_type
NonrigorousHelper<EqSpec, delaysSpec, MatrixSpec, VectorSpec>::DIMENSION = EqSpec::imageDimension();

} // namespace ddeshelper
} // namespace capd

#endif /* _CAPD_DDEHELPERNONRIGOROUS_H_ */
