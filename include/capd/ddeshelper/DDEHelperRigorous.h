/*
 * DDEHelperRigorous.h
 *
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

#ifndef _CAPD_DDEHELPERRIGOROUS_H_
#define _CAPD_DDEHELPERRIGOROUS_H_

#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/DDEHelperCommon.h>
#include <capd/ddeshelper/DDECompareHelper.h>
#include <capd/ddeshelper/DDECubeSet.h>

namespace capd {
namespace ddeshelper {

/**
 * A helper class to facilitate programmer using capdDDEs.
 * It takes the burden of __defining__ the final library classes
 * from templates. Basically, the CAPD code is written in such a way
 * that templates do not need to be tailored for specyfic
 * equations, as RHS of ODE is given in explicit form (a string).
 * Therefore CAPD can have ready types such as DSolver, DPoincareMap, etc.
 * for the user to use right away.
 *
 * In capdDDEs we supply RHS as a templated class for Automatic Differentiation
 * to work. Therefore we need to generate a code specyfic to that DDE.
 * In order to do that, we need to generate code for Solver, Poincare, Section, etc.
 * There are like 10 classes to configure in this way, but the way is usually the same
 * each time. Therefore I have gathered this process in this helper class to
 * simplify things for the user. Simply define a class MyEquation (see examples for docs)
 * then use (numDelays is an integer):
 *
 * 		class RigorousHelper<MyEquation, numDelays>;
 *
 * To generate your code. You can add:
 *
 * 		typedef RigorousHelper<MyEquation, numDelays> Setup;
 *
 * And then use Setup::Solver, Setup::Poincare, Setup::Solution, etc.
 * classes as you would do in CAPD (e.g. instead of IOdeSolver, IPoincare, etc.).
 *
 * The default values of the template arguments MatrixSpec and VectorSpec
 * are tested to work with many equations, and with 'complicated' functions
 * in RHS, like log, exp, sin, tanh, etc. Hovewer, you can use
 * IVector and IMatrix to get maybe a little faster computations by
 * using:
 *
 *		typedef RigorousHelper<MyEquation, numDelays, capd::IMatrix, capd::IVector> Setup;
 *		class RigorousHelper<MyEquation, numDelays, capd::IMatrix, capd::IVector>;
 *
 * but it might not work with some functions in the RHS. This might be
 * corrected in the future, when CAPD resolves a problem with autodiff
 * for capd::Interval type. The setup will definitely work for polynomial RHS.
 *
 * TODO: (NOT URGENT): work more on the docs below
 *
 * Besides, this class also allows to supply good coordinates and compute
 * in those coordinates (also, to give data in those coordinates).
 * So for example if you have coordinates C such that
 * the are given by eigenvectors of Jacobian of poincare map
 * then you can slice your data by 1 along any vector j by specifying
 * dx = {0,...,0,1,0,...}, where 1 is on j-th place.
 *
 * The coordinate change P to good coordinates is:
 * B(dx) = reference + C * dx
 * (so that B(0) == reference)
 * So the backward change is
 * B^{-1} (y) = C^{-1} * (y - reference)
 * Backward change is used in poincare to produce result in 'good' coordinates
 * That is we compute B^{-1}(P(B(.))) applied on a set (dx, r0, Xi)
 * And we get (Pdx, Pr0, PXi) for an easy comparison (coord by coord)
 *
 *
 * 	TODO: (IMPORTANT!!! IN THIS CASE!) docs - GOOD COMMENTS TO ALL ROUTINES!
 */
template<
	typename EqSpec,
	int delaysSpec=1,
	typename MatrixSpec = capd::vectalg::Matrix<typename EqSpec::ParamType, 0, 0>,
	typename VectorSpec = capd::vectalg::Vector<typename EqSpec::ParamType, 0>,
	typename PoliciesSpec=capd::dynset::C11Rect2Policies
>
class RigorousHelper{
public:
	typedef EqSpec Eq;
	typedef typename Eq::ParamType ParamType;
	typedef VectorSpec Vector;
	typedef MatrixSpec Matrix;
	typedef typename Matrix::ScalarType Scalar;
	typedef typename Matrix::ScalarType Real; // TODO: (FUTURE) Rethink? What if scalar is Complex?
	typedef typename Eq::ParamsVectorType ParamsVector;
	typedef PoliciesSpec Policies;
	typedef capd::ddes::SharedDoubleton<Matrix, Policies> SetType;
	typedef capd::ddes::DDESolutionCurve<SetType> Solution;
	typedef typename Solution::GridType Grid;
	typedef typename Solution::TimePointType TimePoint;
	typedef typename Solution::CurvePieceType CurvePiece;
	typedef capd::ddes::DiscreteDelaysFunctionalMap<Eq, Solution> DDEq;
	typedef capd::ddes::DDETaylorSolver<DDEq> Solver;
	typedef typename Solver::VariableStorageType Variables;
	typedef typename Solver::JacobianStorageType Jacobians;
	typedef typename Solver::ValueStorageType Values;
	typedef typename Solver::size_type size_type;
	typedef capd::ddes::DDEJetSection<Solution> Section;
	typedef typename Section::JetType SectionJet;
	typedef capd::ddes::DDEPoincareMap<Solver, Section> PoincareMap;
	typedef DDECompareHelper<Vector> Comparator;
	static const int PARAMS_COUNT;
	static const int DIMENSION;

	/**
	 * Internal class to hold a affine coordinate change
	 * on a (p,n)-fset.
	 */
	class CoordinateSystem {
	public:
		int M;
		Vector reference;
		Vector vsection;
		Scalar secvalue;
		capd::poincare::CrossingDirection crossingDirection;
		Matrix coords;
		Matrix inverseCoords;
		std::string filepath;

		CoordinateSystem(Vector x, Vector sec, Scalar secval, capd::poincare::CrossingDirection d, Matrix C, Matrix invC):
			M(x.dimension()), reference(x), vsection(sec), secvalue(secval), crossingDirection(d), coords(C), inverseCoords(invC)
		{}

		CoordinateSystem(Vector x, Matrix C, Matrix invC, capd::poincare::CrossingDirection d = capd::poincare::MinusPlus):
			M(x.dimension()), reference(x), vsection(C.column(0)), secvalue(C.column(0) * x), crossingDirection(d), coords(C), inverseCoords(invC)
		{}

		explicit CoordinateSystem(int M):
			M(M), reference(M), vsection(M), crossingDirection(capd::poincare::MinusPlus), coords(M, M), inverseCoords(M, M)
		{
			coords.setToIdentity();
			inverseCoords.setToIdentity();
			vsection[0] = 1.0;
		}

		void setCrossingDirection(capd::poincare::CrossingDirection d){ crossingDirection = d; }

		friend std::istream& operator>>(std::istream& in, CoordinateSystem& coords){
			int dir;
			in >> dir;
			if (dir < 0) coords.crossingDirection = capd::poincare::MinusPlus;
			else if (dir > 0) coords.crossingDirection = capd::poincare::PlusMinus;
			else if (dir == 0) coords.crossingDirection = capd::poincare::Both;
			in >> coords.secvalue;
			relativeReadData(in, coords.vsection);
			relativeReadData(in, coords.reference);
			relativeReadData(in, coords.coords);
			relativeReadData(in, coords.inverseCoords);
			return in;
		}

		friend bool operator==(CoordinateSystem& first, CoordinateSystem& second){
			return (first.reference == second.reference) &&
					(first.coords == second.coords);
		}
		friend bool operator!=(CoordinateSystem& first, CoordinateSystem& second){
			return !(first == second);
		}

		Vector x0(){ return reference; }
		Matrix C(){ return coords; }
		Matrix invC(){ return inverseCoords; }
		Vector sectionVector(){ return vsection; }
		Scalar sectionValue(){ return secvalue; }
	};

	RigorousHelper(std::string filepath, int reqSteps=0, int maxSteps=0, int maxOrder=10):
			m_params(PARAMS_COUNT),
			m_p(1), m_n(1),
			m_reqSteps(reqSteps), m_maxSteps(maxSteps), m_maxOrder(maxOrder),
			m_coords(M())
	{ loadSetup(filepath); m_coords.filepath = filepath; }

	RigorousHelper(int p, int n, ParamsVector const& params, int reqSteps=0, int maxSteps=0, int maxOrder=10):
			m_params(PARAMS_COUNT),
			m_p(p), m_n(n),
			m_reqSteps(reqSteps), m_maxSteps(maxSteps), m_maxOrder(maxOrder),
			m_coords(M())
	{
		if (params.dimension() != PARAMS_COUNT)
			throw std::logic_error("RigorousHelper::__construct__(): Bad number of parameters.");

		m_params = params;
		m_coords.coords.setToIdentity();
		m_coords.inverseCoords.setToIdentity();

		updateGrid();
		updateSteps();
	}

	/**
	 * Makes a set of the form (reference + C * dx) + C * r0
	 * It is best to supply r0 as zero centered vector, but
	 * the program makes an adjustment to assure that.
	 */
	Solution dataToSolution(CoordinateSystem const& coords, Vector dx, Vector r0, Vector Xi) const {
		Vector rr = r0;
		capd::vectalg::split(rr, r0);
		auto tau = m_grid.point(m_p);
		auto t_0 = m_grid.point(0);
		Solution X(m_grid, -tau, t_0, m_n, Vector(DIMENSION));
		X.set_x(coords.reference + coords.coords * dx + coords.coords * rr);
		X.set_Cr0(coords.coords, r0);
		X.set_Xi(Xi);
		return X;
	}
	/**
	 * Makes a set of the form (reference + C * dx) + C * r0
	 * It is best to supply r0 as zero centered vector, but
	 * the program makes an adjustment to assure that.
	 */
	Solution dataToSolution(Vector dx, Vector r0, Vector Xi) const {
		return dataToSolution(this->m_coords, dx, r0, Xi);
	}

	/** makes a constant solution over basic delay interval [-tau, 0] */
	Solution constantInitialSolution(Vector const& vx) const {
		return constantInitialSolution(vx, m_n);
	}

	/** makes a constant solution over basic delay interval [-tau, 0] */
	Solution constantInitialSolution(Vector const& vx, size_type order) const {
		if (vx.dimension() != DIMENSION)
			throw std::logic_error("DDERigorousHelper::constantInitialSolution(): incompatible dimension in value.");
		auto tau = m_grid.point(m_p);
		auto t_0 = m_grid.point(0);
		Solution X(m_grid, -tau, t_0, order, vx);
		if (order == m_n)
			X.set_Cr0(coords().coords, Vector(coords().coords.numberOfColumns()));
		else {
			auto dim = X.get_x().dimension();
			X.set_Cr0(Matrix::Identity(dim), Vector(dim));
		}
		return X;
	}

	/**
	 * makes a solution based on the function values over basic delay interval [-tau, 0]
	 * The function must be capd::map::Map type, IMAP works great.
	 * Function must be R -> R^dimension, otherwise - exception!
	 */
	template<typename AnyMatrixSpec>
	Solution functionToSolution(capd::map::Map<AnyMatrixSpec> f) const {
		if (f.dimension() != 1 || f.imageDimension() != DIMENSION){
			std::ostringstream info;
			info << "DDERigorousHelper::functionToSolution(): f should be R -> R^" << DIMENSION << ", ";
			info << "is R^" << f.dimension() << " -> R^" << f.imageDimension() << ".";
			throw std::logic_error(info.str());
		}
		auto tau = m_grid.point(m_p);
		auto t_0 = m_grid.point(0);
		Solution X(m_grid, -tau, t_0, m_n, f);
		X.set_Cr0(coords().coords, Vector(coords().coords.numberOfColumns()));
		return X;
	}

	/**
	 * makes a solution based on the function values over basic delay interval [-tau, 0]
	 * The function must be capd::map::Map type, IMAP works great.
	 * Function must be R -> R^dimension, otherwise - exception!
	 */
	template<typename AnyMatrixSpec>
	Solution functionToSolution(capd::map::Map<AnyMatrixSpec> f, int starti) const {
		if (f.dimension() != 1 || f.imageDimension() != DIMENSION){
			std::ostringstream info;
			info << "DDERigorousHelper::functionToSolution(): f should be R -> R^" << DIMENSION << ", ";
			info << "is R^" << f.dimension() << " -> R^" << f.imageDimension() << ".";
			throw std::logic_error(info.str());
		}
		auto t_tau = m_grid.point(-m_p + starti);
		auto t_0 = m_grid.point(starti);
		Solution X(m_grid, t_tau, t_0, m_n, f);
		X.set_Cr0(coords().coords, Vector(coords().coords.numberOfColumns()));
		return X;
	}

	/**
	 * Makes a set of the form reference + C * r0
	 * It is best to supply r0 as zero centered vector, but
	 * the program makes an adjustment to assure that.
	 */
	Solution r0ToSolution(CoordinateSystem const& coords, Vector r0, Vector Xi) const {
		Vector d0(M());
		return dataToSolution(coords, d0, r0, Xi);
	}
	/**
	 * Makes a set of the form reference + C * r0
	 * It is best to supply r0 as zero centered vector, but
	 * the program makes an adjustment to assure that.
	 */
	Solution r0ToSolution(Vector r0, Vector Xi) const {
		Vector d0(M());
		return dataToSolution(m_coords, d0, r0, Xi);
	}

	/**
	 * TODO: Not implemented yet!
	 */
	capd::poincare::CrossingDirection detectCrossingDirection(){
//		auto tau = m_grid.point(m_p);
//		auto t_0 = m_grid.point(0);
//		Solution X(m_grid, -tau, t_0, m_n, Vector(DIMENSION));
//		X.set_x(x);
//		DDEq dde(Eq(m_params), tau);
//		Solver solver(dde, m_maxOrder);
//		PoincareMap pm(solver, section);
//		return pm.detectCrossingDirection(X);
		throw std::logic_error("DDERigorousHelper::detectCrossingDirection(): Not Implemented Yet");
		return capd::poincare::MinusPlus; // TODO (IMPORTANT) NAPISZ!
	}

	Solver makeSolver(){
		auto tau = m_grid.point(m_p);
		DDEq dde(Eq(m_params), tau);
		Solver solver(dde, m_maxOrder);
		return solver;
	}

	/**
	 * computes timemap using the solution as initial data
	 * WARNING: it will extend the solution given in X by a given number of steps!
	 * returns: it returns the last segment of the solution in a shape of the
	 *          same representation as the initial segment of X.
	 *          if the segment cannot be produced (due to continuity class issues)
	 *          the exception will be thrown.
	 *          if epsilon is given as != 0. then the procedure will try to
	 *          do the epsilonStep procedure.
	 * */
	Solution timemap(Solution& X, const size_type& steps, const Real& epsilon = 0.){
		auto tau = m_grid.point(m_p);
		DDEq dde(Eq(m_params), tau);
		Solver solver(dde, m_maxOrder);

		Solution TX(X.subcurve(0, m_p)); TX *= 0.;
		size_type count = epsilon != 0. ? steps + 1 : steps;
		for (size_type step = 0; step < count; ++step) X.move(solver);
		if (epsilon != 0.){
			X.epsilonShift(solver, epsilon, TX);
		} else {
			X.subcurve(X.length() - m_p).reduce(TX);
		}
		return TX;
	}

	Solution timemap(
			CoordinateSystem const& in_coords,
			Vector const& dx, Vector const& r0, Vector const& Xi,
			const size_type& steps, const Real& epsilon,
			CoordinateSystem const& out_coords,
			Vector& Tdx, Vector& Tr0, Vector& TXi){
		Solution X = dataToSolution(in_coords, dx, r0, Xi);

		auto tau = m_grid.point(m_p);
		DDEq dde(Eq(m_params), tau);
		Solver solver(dde, m_maxOrder);

		// TODO: (NOT URGENT) refactor: DRY, use timemap(Solution, steps, epsilon) implementation
		Tdx = Vector(M());
		Solution TX(X); TX *= 0.;
		size_type count = epsilon != 0. ? steps + 1 : steps;
		for (size_type step = 0; step < count; ++step) X.move(solver);
		if (epsilon != 0.){
			X.epsilonShift(solver, epsilon, TX);
		} else {
			X.subcurve(X.length() - m_p).reduce(TX);
		}
		// we compute transformation: ref + C^-1 * (PX - ref)
		Tdx = (out_coords.inverseCoords * (TX.get_x() - out_coords.reference)) + 	// shift middle by reference and then apply coords change on the remainder
			  (out_coords.inverseCoords * TX.get_C()) * TX.get_r0() +				// apply coords change to Cr0 part
			  (out_coords.inverseCoords * TX.get_B()) * TX.get_r();					// apply coords change to Br part
		capd::vectalg::split(Tdx, Tr0);
		TXi = TX.get_Xi();
		return X; // X contains full trajectory
	}

	Solution timemap(
			Vector const& dx, Vector const& r0, Vector const& Xi,
			const size_type& steps, const Real& epsilon,
			Vector& Tdx, Vector& Tr0, Vector& TXi){
		return timemap(coords(), dx, r0, Xi, steps, epsilon, coords(), Tdx, Tr0, TXi);
	}

	/**
	 * TODO: (IMPORTANT!!! IN THIS CASE!) docs,
     */
	Solution poincare(
			Section section,
			capd::poincare::CrossingDirection crossing_direction,
			Solution const& X,
			Real& reachTime, size_type& steps, Real& epsilon,
			Solution& PX){
		if (X.length() != m_p){
			std::ostringstream info;
			info << "RigorousHelper::poincare(section, X, ...): ";
			info << "X must be of length exactly " << m_p << ", is of length " << X.length() << std::endl;
			throw std::logic_error(info.str());
		}

		auto tau = m_grid.point(m_p);
		DDEq dde(Eq(m_params), tau);
		Solver solver(dde, m_maxOrder);
		PoincareMap pm(solver, section);
		pm.setDirection(crossing_direction);
		pm.setRequiredSteps(m_reqSteps);
		pm.setMaxSteps(m_maxSteps);

		Solution TX = X;
		PX *= 0.;
		pm(TX, reachTime, PX);
		epsilon = pm.getLastEpsilonTime();
		steps = TX.length() - m_p - 1; // TODO: this is not accurate... but for now it must suffice
		return TX; // TX contains full trajectory
	}

	/**
	 * TODO: (IMPORTANT!!! IN THIS CASE!) docs,
     */
	Solution poincare(
			Section section,
			capd::poincare::CrossingDirection crossing_direction,
			Solution const& X,
			Real& reachTime,
			Solution& PX){
		Real eps; size_type steps;
		return poincare(section, crossing_direction, X, reachTime, steps, eps, PX);
	}

	/**
	 * TODO: (IMPORTANT!!! IN THIS CASE!) docs,
     */
	Solution poincare(
			CoordinateSystem const& in_coords,
			Vector const& dx, Vector const& r0, Vector const& Xi,
			CoordinateSystem const& out_coords,
			Real& reachTime, size_type& steps, Real& epsilon,
			Vector& Pdx, Vector& Pr0, Vector& PXi){
		Solution X = dataToSolution(in_coords, dx, r0, Xi);

		auto tau = m_grid.point(m_p);
		DDEq dde(Eq(m_params), tau);
		Solver solver(dde, m_maxOrder);

		Section section(d(), p(), n(), out_coords.vsection, out_coords.secvalue);
		PoincareMap pm(solver, section);
		pm.setDirection(out_coords.crossingDirection);
		pm.setRequiredSteps(m_reqSteps);
		pm.setMaxSteps(m_maxSteps);

		Pdx = Vector(M());
		Solution PX(X); PX *= 0.;
		pm(X, reachTime, PX);
		epsilon = pm.getLastEpsilonTime();
		steps = X.length() - m_p - 1; // TODO: this is not accurate... but for now it must suffice
		// we compute transformation: ref + C^-1 * (PX - ref)
		Pdx = (out_coords.inverseCoords * (PX.get_x() - out_coords.reference)) + 	// shift middle by reference and then apply coords change on the remainder
			  (out_coords.inverseCoords * PX.get_C()) * PX.get_r0() +				// apply coords change to Cr0 part
			  (out_coords.inverseCoords * PX.get_B()) * PX.get_r();					// apply coords change to Br part
		capd::vectalg::split(Pdx, Pr0);
		PXi = PX.get_Xi();
		return X; // X contains full trajectory
	}

	void toCoords(Solution const& PX, CoordinateSystem const& coords, Vector& Pdx, Vector& Pr0, Vector& PXi){
		// TODO: check dimensions of PX!
		Pdx = Vector(M());
		Pdx = (coords.inverseCoords * (PX.get_x() - coords.reference)) + 	// shift middle by reference and then apply coords change on the remainder
			  (coords.inverseCoords * PX.get_C()) * PX.get_r0() +			// apply coords change to Cr0 part
			  (coords.inverseCoords * PX.get_B()) * PX.get_r();				// apply coords change to Br part
		capd::vectalg::split(Pdx, Pr0);
		PXi = PX.get_Xi();
	}

	/** out_coords are given by the coords of this Helper */
	Solution poincare(
			CoordinateSystem const& in_coords,
			Vector const& dx, Vector const& r0, Vector const& Xi,
			Real& reachTime, Vector& Pdx, Vector& Pr0, Vector& PXi){
		size_type dump_s;
		Real dump_e;
		return poincare(in_coords, dx, r0, Xi, m_coords, reachTime, dump_s, dump_e, Pdx, Pr0, PXi);
	}

	/** in_ and out_ coords are given by the coords of this Helper */
	Solution poincare(
			Vector const& dx, Vector const& r0, Vector const& Xi,
			Real& reachTime, Vector& Pdx, Vector& Pr0, Vector& PXi){
		return poincare(m_coords, dx, r0, Xi, reachTime, Pdx, Pr0, PXi);
	}

//	template<typename VectorIteratorSpec>
//	void saveData(std::string filepath, VectorIteratorSpec start, VectorIteratorSpec end, std::string extraComment=""){
//		std::ofstream outf(filepath); outf.precision(15);
//		outf << PARAMS_COUNT << " " << m_params << std::endl;
//		outf << DIMENSION << " " << m_p << " " << m_n << std::endl;
//		outf << int(end - start) << std::endl;
//		for (; start != end; ++start)
//			outf << *start << std::endl;
//		outf << "# INFO on format: " << std::endl;
//		outf << "# first line:      number_of_params (" << PARAMS_COUNT << "), then params." << std::endl; // TODO: doac zeby rownania potrafily mowic jakie maja prarametry
//		outf << "# second line:     d p n, as in (p,n)-representations." << std::endl;
//		outf << "# third line:      count, the number of data for solutions. How the programs will use solutions, it depends on the program." << std::endl;
//		outf << "# following lines: the solutions that fit into (p,n)-representation. " << std::endl;
//		outf << std::endl << extraComment << std::endl;
//		outf.close();
//	}

	void loadSetup(std::string filepath);
	void setRequiredSteps(int reqSteps, bool ensure_long_enough = true){ m_reqSteps = reqSteps; updateSteps(ensure_long_enough); }
	int getRequiredSteps() const { return m_reqSteps; }
	void setMaximumSteps(int maxSteps, bool ensure_long_enough = true){ m_maxSteps = maxSteps; updateSteps(ensure_long_enough); }
	int getMaximumSteps() const { return m_maxSteps; }
	void setMaximumOrder(int maxOrder){ m_maxOrder = maxOrder; }
	int getMaximumOrder() const { return m_maxOrder; }

	ParamsVector params() const { return m_params; }
	int M() const { return DIMENSION * (1 + m_p * (m_n +1)); }
	int p() const { return m_p; }
	int n() const { return m_n; }
	int d() const { return DIMENSION; }
	TimePoint h() const { return m_grid.point(1); }
	Real H() const { return Real(m_grid.point(1)) * Real(0, 1); }
	Real getBasicIntervalLength() const { return m_params[m_params.dimension()-1]; /* TODO: (IMPORTANT): make it more general, many delays */ }
	TimePoint tau() const { return m_grid.point(m_p); }
	Grid& grid() { return m_grid; }
	const Grid& grid() const { return m_grid; }

	capd::poincare::CrossingDirection getCrossingDirection() { return m_coords.crossingDirection; }
	void setCrossingDirection(capd::poincare::CrossingDirection d) { m_coords.crossingDirection = d; }

	CoordinateSystem coords() const { return m_coords; }
	void setCoords(CoordinateSystem const& c) { m_coords = c; }
	Vector reference(){ return coords().reference; }
	Vector x0(){ return coords().reference; }
	Matrix C(){ return coords().coords; }
	Matrix invC(){ return coords().inverseCoords; }
	Vector sectionVector(){ return coords().vsection; }
	Scalar sectionValue(){ return coords().sectionValue(); }
	Section section(){ return Section(d(), p(), n(), coords().vsection, coords().secvalue); }

	/**
	 * this is to dump data as text format (vectors as text representation)
	 * Not the best to use in the proofs, but can be of use in preparation of data.
	 */
	void dumpData(
			std::ostream& out, CoordinateSystem const& coords,
			Vector const& dx, Vector const& r0, Vector const& Xi){
		if (coords.filepath == ""){
			// TODO: save the CoordsSystem along with the data
			throw std::logic_error("RigorousHelper::dumpData(): Empty filepath data in coordinate system - cannot handle that kind of situation.");
		}
		out.precision(15);
		out << coords.filepath << std::endl;
		out << "RELATIVE" << std::endl;
		out << dx << std::endl;
		out << r0 << std::endl;
		out << Xi << std::endl;
	}

	void saveData(
			PathConfig const& paths, CoordinateSystem const& coords,
			Vector const& dx, Vector const& r0, Vector const& Xi,
			std::string extraMsg=""){
		if (coords.filepath == ""){
			// TODO: save the CoordsSystem along with the data
			throw std::logic_error("RigorousHelper::saveData(): Empty filepath data in coordinate system. (NOT IMPLEMENTED STATE YET)");
		}
		mkdir_p(paths.dirPath);
		std::string binpath = paths.dirPath + "/" + paths.prefix + "-bin-initial-rig.txt";
		std::string txtpath = paths.dirPath + "/" + paths.prefix + "-txt-initial-rig.txt";
		std::ofstream outbin(binpath);
		std::ofstream outtxt(txtpath);
		outtxt.precision(15);
		outtxt << coords.filepath << std::endl;
		outtxt << "RELATIVE" << std::endl;
		outtxt << dx << std::endl;
		outtxt << r0 << std::endl;
		outtxt << Xi << std::endl;
		outtxt << std::endl;
		outtxt << "# This file is a human-readable version of '" << binpath << "'" << std::endl;
		outtxt << "# You should use binary version for actual rigorous computations! " << std::endl;
		outtxt << extraMsg << std::endl;

		std::string dxbinpath = paths.dirPath + "/" + paths.prefix + "-initial-dx.bin";
		std::string r0binpath = paths.dirPath + "/" + paths.prefix + "-initial-r0.bin";
		std::string Xibinpath = paths.dirPath + "/" + paths.prefix + "-initial-Xi.bin";

		outbin << coords.filepath << std::endl;
		outbin << "RELATIVE" << std::endl;
		outbin << dxbinpath << std::endl;
		outbin << r0binpath << std::endl;
		outbin << Xibinpath << std::endl;
		outbin << std::endl;
		outbin << "# You should use THIS FILE for actual rigorous computations! " << std::endl;
		outbin << "# This file is binary (exact) version of '" << txtpath << "' if you need to look at human readable version" << std::endl;
		outbin << extraMsg << std::endl;

		saveBinary(dxbinpath, dx);
		saveBinary(r0binpath, r0);
		saveBinary(Xibinpath, Xi);
	}

	CoordinateSystem loadData(std::istream& input, Vector& dx, Vector &r0, Vector &Xi){
		std::string type, coordspath;
		input >> coordspath;
		std::string dump;
		input >> dump; // TODO: old pomysl, zeby nie zmieniac na raie plikow, potem wywalic i naprawic pliki, zeby nie bylo tej linii RELATIVE/ABSOLITE
		RigorousHelper tmp(coordspath);
		relativeReadData(input, dx);
		relativeReadData(input, r0);
		relativeReadData(input, Xi);
		return tmp.coords();
	}

	CoordinateSystem loadData(std::string filepath, Vector& dx, Vector &r0, Vector &Xi){
		std::ifstream input(filepath);
		return loadData(input, dx, r0, Xi);
	}

	template<typename SolutionT>
	void drawSolution(std::string dirpath, std::string filename, SolutionT const& X){
		{
			std::ostringstream cmd;
			cmd << "mkdir -p '" << dirpath << "' ";
			capd::ddeshelper::runSystemCommand(cmd.str());
		}
		std::ofstream outfX(dirpath + "/" + filename + ".dat"); outfX.precision(15);
		capd::ddeshelper::value_to_gnuplot(outfX, X.pastTime(), X.currentTime(), X.getStep(), X, H());
		outfX.close();
		std::ofstream outg(dirpath + "/" + filename + ".gp");
		outg << "set terminal png size 800,600" << std::endl;
		outg << "set output '" << filename << ".png" << "'" << std::endl;
		outg << "plot '" <<  filename << ".dat" << "' using 1:3 with lines, \\\n";
		outg << "     '" <<  filename << ".dat" << "' using 1:($3-$4):($3+$4) with filledcu notitle";
		outg.close();
		{
			std::ostringstream cmd;
			cmd << "cd '" << dirpath << "' && gnuplot '" << filename << ".gp'";
			capd::ddeshelper::runSystemCommand(cmd.str());
		}
	}

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

private:
	ParamsVector m_params;
	int m_p;
	int m_n;
	int m_reqSteps;
	int m_maxSteps;
	int m_maxOrder;
	Grid m_grid;
	CoordinateSystem m_coords;

	void updateGrid() { m_grid = Grid(getBasicIntervalLength() / m_p); };
	void rawLoadSetup(
			std::istream& in,
			int &p, int &n, Vector &params,
			CoordinateSystem& coords);
	void updateSteps(bool ensure_long_enough = true) {
		m_reqSteps = (m_reqSteps < 0 ? -m_reqSteps * m_p : m_reqSteps);
		if (ensure_long_enough)
			m_reqSteps = (m_reqSteps < (m_n+1) * m_p ? (m_n+1) * m_p : m_reqSteps); // long enough steps
		m_maxSteps = (m_maxSteps < m_reqSteps ? m_reqSteps + 1 : m_maxSteps);
	}
};

template<typename EqSpec, int delaysSpec, typename MatrixSpec, typename VectorSpec, typename PoliciesSpec>
const int RigorousHelper<EqSpec, delaysSpec, MatrixSpec, VectorSpec, PoliciesSpec>::PARAMS_COUNT = EqSpec::getParamsCount() + delaysSpec;
template<typename EqSpec, int delaysSpec, typename MatrixSpec, typename VectorSpec, typename PoliciesSpec>
const int RigorousHelper<EqSpec, delaysSpec,  MatrixSpec, VectorSpec, PoliciesSpec>::DIMENSION = EqSpec::imageDimension();

} // namespace ddeshelper
} // namespace capd



#endif /* _CAPD_DDEHELPERRIGOROUS_H_ */
