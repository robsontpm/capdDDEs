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

// TODO: (NOT URGENT): we have NonrigorousHelper but it is in DDEHelperNonrigorous.h/.hpp,... refactor...

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
 * You might also supply Matrix and Vector type.
 *
 * It can only handle discrete delay differential equations in the form of
 *
 * x'(t) = f(x(t), x(t-\tau_1), ... , x(t-tau_m))
 *
 * TODO: (IMPORTANT) Currently it only support one delay!!!!!!
 * TODO: (URGENT?) DEV: make sure all const qualifiers are put in the right places!
 */
template<
	typename EqSpec,
	int delaysSpec=1,
	typename MatrixSpec = capd::vectalg::Matrix<typename EqSpec::ParamType, 0, 0>,
	typename VectorSpec = capd::vectalg::Vector<typename EqSpec::ParamType, 0>,
	typename PoliciesSpec=capd::dynset::C11Rect2Policies
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
	typedef capd::ddes::GenericJet<TimePoint, capd::ddes::VectorWithJacData<Vector, Matrix>, Vector, Matrix> JetWithJacobian;
	typedef capd::ddes::DDEPiecewisePolynomialCurve<Grid, Jet> Solution;
	typedef capd::ddes::DDEPiecewisePolynomialCurve<Grid, JetWithJacobian> JacSolution;
	typedef Jet CurvePiece;
	typedef capd::ddes::BasicDiscreteDelaysFunctionalMap<Eq, JacSolution> JacDDEq;
	typedef capd::ddes::BasicDiscreteDelaysFunctionalMap<Eq, Solution> DDEq;
	typedef capd::ddes::DDENonrigorousTaylorSolver<JacDDEq> JacSolver;
	typedef capd::ddes::DDENonrigorousTaylorSolver<DDEq> Solver;
	typedef typename JacSolver::VariableStorageType Variables;
	typedef typename JacSolver::JacobianStorageType Jacobians;
	typedef typename JacSolver::ValueStorageType Values;
	typedef typename JacSolver::size_type size_type;
	typedef int step_type; // we need to have signed values here
	typedef capd::ddes::DDEJetSection<JacSolution> JacJetSection;
	typedef typename JacJetSection::JetType JacSecJet;
	typedef capd::ddes::DDEBasicPoincareMap<JacSolver, JacJetSection> JacPoincareMap;
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
	 * This creates a raw solver, if you want to do nonstandard tasks
	 * If you want just to integrate initial values, consider using iterate() or poincare() instead.
	 * TODO: add functions to create more elements
	 */
	Solver makeSolver(){
		auto tau = m_grid.point(m_p);
		DDEq dde(Eq(m_params), tau);
		Solver solver(dde, m_maxOrder);
		return solver;
	}
	/**
	 * just integrate the solution and return the solution
	 *
	 * initial can be either of (1) d-dimension (\R^d) or (2) the M(d,p,n)-dimension (see papers)
	 * In the (1) case it integrates initial function x_0(s) == initial for all s \in [-tau, 0]
	 * In the (2) case it constructs initial function from data stored in x0 (the order of coefficients
	 * in this vector is described elsewhere (see papers for example). But suppling it by hand
	 * is extremally cumbersome. It is best to use output of other functions o get initials in this
	 * form.
	 *
	 * Integrates for time T = iters * tau_max / p = iters * h. The value h is the step size of the method / grid
	 *
	 * In the result the last segment x_T from the solution in the vector format.
	 * Can be used as an input to the next integrate() procedure, to create an initial solution
	 * curve, or in other functions from the helper (e.g. poincare).
	 *
	 * TODO: use_extension was used for backward compatibility of some programs
	 *       it should work with use_extension=true in the default version
	 *       and user should control it with setting m_maxOrder with getter/setter
	 */
	Solution integrate(int iters, const Vector& initial, Vector& result, bool use_extension=false){
		if (initial.dimension() != DIMENSION && initial.dimension() != M()){
			throw std::logic_error("NonrigorousHelper::initialIteration(): vector must be either DIMENSION (for const) or M dimensional (for any initial)");
		}
		if (iters < 0) iters = -m_p * iters;
		auto tau = m_grid.point(m_p);
		auto t_0 = m_grid.point(0);
		Vector v(DIMENSION); if (initial.dimension() == DIMENSION) v = initial;
		Solution X(m_grid, -tau, t_0, m_n, v);
		Solution XXX(m_grid, -tau, t_0, m_n, v);
		if (initial.dimension() == X.storageDimension()){ X.set_x(initial); XXX.set_x(initial); }
		DDEq dde(Eq(m_params), tau);
		Solver solver(dde, use_extension ? m_maxOrder : m_n);
		for (int i = 0; i < iters; i++)
			solver(XXX);
		 result = XXX.subcurve(XXX.currentTime() - tau).get_x();
		 return XXX;
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
	 */
	Solution integrate(int iters, const Vector& initial){
		Vector dump;
		return integrate(iters, initial, dump);
	}
	/**
	 * integrate, save output to files and draw the solution. See other versions
	 * for more information on input and output.
	 *
	 * The variable outconfig holds a pair of strings that define paths to where store
	 * the results of the computations. See documentation there.
	 */
	Solution integrate(int iters, const Vector& initial, Vector& result, PathConfig const& outconfig){
		Solution solution = integrate(iters, initial, result);
		drawSolution(outconfig.dirPath, outconfig.prefix, solution);
		Vector start = solution.subcurve(solution.pastTime(), m_grid(0)).get_x();
		Vector iterated = solution.subcurve(solution.currentTime() - m_grid(m_p)).get_x();
		saveData(outconfig.filepath("start"), &start, &start+1);
		saveData(outconfig.filepath("iterated"), &iterated, &iterated+1);
		return solution;
	}

	Real iteratePoincare(std::ostream& info, int iters, JacJetSection& section, Vector &x, Matrix& V, Matrix& DP){
		info << "# iteratePoincare START" << std::endl;
		auto tau = m_grid.point(m_p);
		auto t_0 = m_grid.point(0);
		JacSolution X(m_grid, -tau, t_0, m_n, Vector(DIMENSION));
		X.set_x(x);

		JacDDEq dde(Eq(m_params), tau);
		JacSolver solver(dde, m_maxOrder);
		JacPoincareMap pm(solver, section);
		pm.setDirection(this->crossingDirection);
		pm.setRequiredSteps(m_reqSteps);
		pm.setMaxSteps(m_maxSteps);

		Real reachTime; Jacobians Jac;
		capd::vectalg::EuclLNorm<Vector, Matrix> euclNorm;
		JacSolution on_section(X); on_section *= 0.;
		pm(X, on_section, reachTime);
		X = on_section; on_section *= 0.;
		info << "# time of X on section: " << X.getT0() << " should be 0" << std::endl;
		info << "# reqSteps: " << pm.getRequiredSteps() << std::endl;
		info << "# maxSteps: " << pm.getMaximumSteps() << std::endl;
		double diff;
		for (int i = 0; i < iters; i++){
			Vector xx = X.get_x();
			pm(X, on_section, reachTime);
			Vector yy = on_section.get_x();
			diff = euclNorm(xx-yy);
			X = on_section; on_section *= 0.;
			info << "# reqSteps: " << pm.getRequiredSteps() << std::endl;
			info << "# maxSteps: " << pm.getMaximumSteps() << std::endl;
			info << "# time of X on section: " << X.getT0() << " should be 0" << std::endl;
			info << "# reachTime = " << reachTime << std::endl;
			info << "# ||PX - X|| = " << diff << std::endl;
		}
		info << "# Final poincare: " << std::flush;
		pm.setRequiredSteps(m_reqSteps);			info << "#" << std::flush;
		pm.setMaxSteps(m_maxSteps);					info << "#" << std::flush;
		x = X.get_x();								info << "#" << std::flush;
		Vector Px, fPx;								info << "#" << std::flush;
		pm.setInitialV(X);							info << "#" << std::flush;
		JacSolution PX(X); PX *= 0.; 				info << "#" << std::flush;
		pm(X, PX, reachTime, x, Px, fPx, V, DP); 	info << "#" << std::flush;
		diff = euclNorm(x - PX.get_x()); 			info << " DONE" << std::endl;
		info << "# iteratePoincare END" << std::endl;
		return diff;
	}

	/**
	 * makes a section and stores it in out_section
	 *
	 * if dim(s) == d, then we setup section in the coordinate x(0) . s = c
	 * otherwise we set full-space section s . x_0 = c
	 */
	template<typename SecionType>
	void makeSection(Vector const& s, Scalar const& c, SecionType& out_section){
		Vector v;
		if (s.dimension() == d()){
			v = Vector(M());
			for (size_t i = 0; i < d(); ++i)
				v[i] = s[i];
		}else if (s.dimension() == M()){
			v = s;
		}else {
			throw std::logic_error("NonrigorousHelper::initialIteration(): vector must be either DIMENSION (for const) or M dimensional (for any initial)");
		}
		out_section = SecionType(d(), p(), n(), v, c);
	}

	/**
	 * helper function to refine a candidate periodic orbit with a Newton method.
	 *
	 * First parameter is to get some text info back during the process, you can pass std::cout there.
	 *
	 * Section you need to set-up, it will be s * x_0 = x, s \in \R^{M(d,p,n), c \in \R}
	 */
	Real refinePeriodic(std::ostream& info, JacJetSection& section, Vector& x, Matrix& V, Matrix&DP){
		info << "refinePeriodic START" << std::endl;
		auto tau = m_grid.point(m_p);
		auto t_0 = m_grid.point(0);
		JacSolution X(m_grid, -tau, t_0, m_n, Vector(DIMENSION));
		X.set_x(x);

		JacDDEq dde(Eq(m_params), tau);
		JacSolver solver(dde, m_maxOrder);

		JacPoincareMap pm(solver, section);
		pm.setDirection(this->crossingDirection);
		pm.setRequiredSteps(m_reqSteps);
		pm.setMaxSteps(m_maxSteps);

		capd::vectalg::EuclLNorm<Vector, Matrix> euclNorm;
		Real reachTime;
		Vector Px, fPx, dx(X.storageDimension());
		Matrix Id(X.storageDimension(), X.storageDimension()); Id.setToIdentity();
		pm.setInitialV(X);
		JacSolution PX(X); PX *= 0.;
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
	capd::poincare::CrossingDirection detectCorssingDirection(JacJetSection section, Vector const& x){
		auto tau = m_grid.point(m_p);
		auto t_0 = m_grid.point(0);
		JacSolution X(m_grid, -tau, t_0, m_n, Vector(DIMENSION));
		X.set_x(x);
		JacDDEq dde(Eq(m_params), tau);
		JacSolver solver(dde, m_maxOrder);
		JacPoincareMap pm(solver, section);
		return pm.detectCrossingDirection(X);
	}
	/**
	 * if you do not know what initV is, then
	 * you should probably stick to using
	 * the other JacSolution poincare() method (without initV).
	 */
	JacSolution poincare(
				JacJetSection section,
				Vector const& x, Matrix& initV,
				double& reachTime, int& steps, Vector& Px, Vector& fPx,
				Matrix& V, Matrix& DP){
		auto tau = m_grid.point(m_p);
		auto t_0 = m_grid.point(0);
		JacSolution X(m_grid, -tau, t_0, m_n, Vector(DIMENSION));
		X.set_x(x);

		JacDDEq dde(Eq(m_params), tau);
		JacSolver solver(dde, m_maxOrder);

		JacPoincareMap pm(solver, section);
		pm.setDirection(this->crossingDirection);
		pm.setRequiredSteps(m_reqSteps);
		pm.setMaxSteps(m_maxSteps);

		capd::vectalg::EuclLNorm<Vector, Matrix> euclNorm;
		Px = Vector(M()); fPx = Vector(M());
		V = Matrix(M(), M()); DP = Matrix(M(), M());
		Matrix Id(X.storageDimension(), X.storageDimension()); Id.setToIdentity();
		pm.setNormalizeVariational(m_experimentalRenormalizeVariational);
		pm.setInitialV(X, initV);
		JacSolution PX(X); PX *= 0.;
		Vector dump = x;
		pm(X, PX, reachTime, dump, Px, fPx, V, DP);
		steps = pm.getLastStepsAfterSection();
		return X; // X contains full trajectory
	}
	/**
	 * compuytes poincare map and the (approximate)
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
	JacSolution poincare(
				JacJetSection section, Vector const& x0,
				double& reachTime, Vector& Px, Vector& fPx,
				Matrix& V, Matrix& DP){
		Matrix Id(x0.dimension(), x0.dimension()); Id.setToIdentity();
		int steps;
		return 	poincare(section, x0, Id, reachTime, steps, Px, fPx, V, DP);
	}
	/**
	 * the simples computation of poincare map, without extra data
	 */
	Solution poincare(
				JetSection section, Vector const& x,
				double& reachTime, Vector& Px){
		auto tau = m_grid.point(m_p);
		auto t_0 = m_grid.point(0);
		Solution X(m_grid, -tau, t_0, m_n, Vector(DIMENSION));
		X.set_x(x);

		DDEq dde(Eq(m_params), tau);
		Solver solver(dde, m_maxOrder);

		PoincareMap pm(solver, section);
		pm.setDirection(this->crossingDirection);
		pm.setRequiredSteps(m_reqSteps);
		pm.setMaxSteps(m_maxSteps);

		capd::vectalg::EuclLNorm<Vector, Matrix> euclNorm;
		Solution PX(X); PX *= 0.;
		pm(X, PX, reachTime);
		Px = PX.get_x();
		return X; // X contains full trajectory
	}
	// void computeCoordinates(std::ostream& info, Matrix const& V, Matrix const& DP, Matrix& coords);

	Solution vectorToSolution(Vector const& x) const{
		auto tau = m_grid.point(m_p);
		auto t_0 = m_grid.point(0);
		Solution X(m_grid, -tau, t_0, m_n, Vector(DIMENSION));
		X.set_x(x);
		return X;
	}
	JacSolution vectorToJacSolution(Vector const& x) const{
		auto tau = m_grid.point(m_p);
		auto t_0 = m_grid.point(0);
		JacSolution X(m_grid, -tau, t_0, m_n, Vector(DIMENSION));
		X.set_x(x);
		return X;
	}

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
			plotcmd << "'" << filename << ".dat" << "' using 1:2 with lines";
		else
			plotcmd << "'" << filename << ".dat" << "' using ($1+" << tshift <<"):2 with lines";
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
	template<typename SolutionT>
	void drawSolution(PathConfig const& paths, SolutionT const& X, double tshift = 0.) const {
		drawSolution(paths.dirpath(), paths.filename("solution"), X, tshift);
	}
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
	template<typename SolutionT>
	void drawDelayMap(PathConfig const& paths, SolutionT const& X) const {
		drawDelayMap(paths.dirpath(), paths.filename("phasespace"), X);
	}
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

	void setRequiredSteps(step_type reqSteps, bool control_steps = true){ m_reqSteps = reqSteps; updateSteps(control_steps); }
	void setMaximumSteps(step_type maxSteps, bool control_steps = true){ m_maxSteps = maxSteps; updateSteps(control_steps); }
	step_type getRequiredSteps() const { return m_reqSteps; }
	step_type getMaximumSteps() const { return m_maxSteps; }
	step_type getMaximumOrder() const { return m_maxOrder; }

	ParamsVector params() const { return m_params; }
	size_type M() const { return DIMENSION * (1 + m_p * (m_n +1)); }
	size_type p() const { return m_p; }
	size_type n() const { return m_n; }
	size_type d() const { return DIMENSION; }
	Real h() const { return getBasicIntervalLength() / p(); }
	Real getBasicIntervalLength() const { return m_params[m_params.dimension()-1]; /* TODO: (IMPORTANT): make it more general, many delays */ }

	void setCrossingDirection(capd::poincare::CrossingDirection const& d) { this->crossingDirection = d; }

	const Grid& grid() const { return m_grid; }

	NonrigorousHelper& setExperimentalRenormalizeVariational(bool v) { m_experimentalRenormalizeVariational = v; return *this; }

	/** returns old params */
	ParamsVector setParams(ParamsVector const& new_params) {
		auto old_params = m_params;
		m_params = new_params;
		return old_params;
	};

	/** returns old param value; TODO: add index checking */
	ParamType setParam(size_type index, ParamType const& new_param) {
		auto old_param = m_params[index];
		m_params[index] = new_param;
		return old_param;
	};

	TimePoint t(int i) const { return m_grid(i); }

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
};

template<typename EqSpec, int delaysSpec, typename MatrixSpec, typename VectorSpec, typename PoliciesSpec>
const typename NonrigorousHelper<EqSpec, delaysSpec, MatrixSpec, VectorSpec, PoliciesSpec>::size_type
NonrigorousHelper<EqSpec, delaysSpec, MatrixSpec, VectorSpec, PoliciesSpec>::PARAMS_COUNT = EqSpec::getParamsCount() + delaysSpec;

template<typename EqSpec, int delaysSpec, typename MatrixSpec, typename VectorSpec, typename PoliciesSpec>
const typename NonrigorousHelper<EqSpec, delaysSpec, MatrixSpec, VectorSpec, PoliciesSpec>::size_type
NonrigorousHelper<EqSpec, delaysSpec, MatrixSpec, VectorSpec, PoliciesSpec>::DIMENSION = EqSpec::imageDimension();

} // namespace ddeshelper
} // namespace capd

#endif /* _CAPD_DDEHELPERNONRIGOROUS_H_ */
