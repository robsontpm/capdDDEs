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

#ifndef _CAPD_DDES_DDE_BASICPOINCARE_MAP_H_
#define _CAPD_DDES_DDE_BASICPOINCARE_MAP_H_

#include <capd/ddes/DDECommon.h>
#include <capd/poincare/BasicPoincareMap.h>

namespace capd{
namespace ddes{

using namespace poincare;

#define DEFAULT_BINSEARCH_EPSILON 1e-14
#define DEFAULT_MAX_ITERATION 1000

/**
 * Implementation of Poincare Map for use with nonrigorous DDESolvers.
 *
 * In theory should work with any compatible DDEDynSys and Section.
 *
 * TODO: (NOT URGENT): check the interface is compatible with CAPD Poincare Map
 * TODO: (NOT URGENT): Especially operators and  integrateUntilSectionCrossing, if they need extra n (numnber of iterates) parameter
 * TODO: (FUTURE, NOT URGENT, RETHINK): should I rewrite the library so that normal PoncareMap from CAPD can be used with my system? It may not be feasible, as CAPD assumes not (semi) in dynamical system.
 */
template<typename DynSysSpec, typename SectionSpec>
class DDEBasicPoincareMap {
public:
	typedef DDEBasicPoincareMap<DynSysSpec, SectionSpec> Class;
	typedef DynSysSpec DynSysType;
	typedef SectionSpec SectionType;
	typedef typename DynSysSpec::CurveType CurveType;
	typedef typename CurveType::TimePointType TimePointType;
	typedef typename CurveType::VectorType VectorType;
	typedef typename CurveType::MatrixType MatrixType;
	typedef typename CurveType::ScalarType ScalarType;
	typedef typename CurveType::RealType RealType;
	typedef typename DynSysType::JacobianStorageType JacobianStorageType;
	typedef typename CurveType::size_type size_type;

private:
	typedef void (Class::*SingleStepFn)(CurveType&);
	void singleStep(CurveType& curve){ m_dynsys(curve); }
	void singleStepWithVariational(CurveType& curve){
		m_dynsys(curve, m_variational);
		// if (m_normalizeVariational && m_steps % m_storedInitialShape.size() == 0){
		if (m_normalizeVariational && m_steps > 0){
			std::cout << "Doing renormalization at " << m_steps << std::endl;
			MatrixType fullV, redV;
			extractVariationalMatrix(curve, fullV, redV, m_storedInitialShape, m_storedInitialShape.size());
			std::cout << fullV.numberOfColumns() << " " << redV.numberOfColumns() << std::endl;
			if (fullV.numberOfColumns() != redV.numberOfColumns())
				throw std::logic_error("DDEBasicPoincareMap::singleStepWithVariational(): (EXPERIMENTAL) fullV has diff shape than redV (num cols) ");
			for (size_type j = 0; j < fullV.numberOfColumns(); ++j){
				ScalarType norm = sqrt(redV.column(j) * redV.column(j)); // TODO: use EuclNorm from CAPD
				fullV.column(j) /= norm;
				// std::cout << "   Norm at " << j << " : " << norm << std::endl;
			}
			setCurrentV(curve, fullV);
		}
	}

public:

	DDEBasicPoincareMap(DDEBasicPoincareMap const& other):
		m_dynsys(other.m_dynsys), m_section(other.m_section),
		m_direction(other.m_direction),
		m_requiredSteps(other.m_requiredSteps), m_maxSteps(other.m_maxSteps), m_steps(other.m_steps),
		m_binsearchEpsilon(other.m_binsearchEpsilon),
		m_normalizeVariational(false)
	{}

	DDEBasicPoincareMap(
			DynSysType& dynsys,
			SectionType& section,
			CrossingDirection direction = CrossingDirection::Both,
			int reqSteps = 0,
			int maxSteps = -1,
			double binsearchEpsilon = DEFAULT_BINSEARCH_EPSILON):
		m_dynsys(dynsys), m_section(section),
		m_direction(direction),
		m_requiredSteps(reqSteps), m_maxSteps(maxSteps), m_steps(0),
		m_binsearchEpsilon(binsearchEpsilon),
		m_normalizeVariational(false)
	{}

	/**
	 * This is probably the simplest version to be used by the end user.
	 *
	 * It is simply applied like:
	 *
	 * 		CurveType X(...); // init to what you want
	 * 		auto PX = P(X);   // PX will be the same shape as X
	 * 		// so you can e.g. compute difference between representations
	 * 		capd::vectalg::EuclNorm<VectorType, MatrixType> euclNorm;
	 * 		auto diff = euclNorm((VectorType)PX - (VectorType)X)
	 *
	 * to be compatible with standard notion and CAPD interface,
	 * It makes following assumptions:
	 *  * X and PX will have the same structure
	 *    (i.e. over same grid points and jets of the same order).
	 *
	 * NOTE: you can retrieve information on reach time and epsilon step
	 *       (see getLast***() functions),
	 *       but you cannot acquire the solution on the whole time.
	 *       If you need solution over the full integration time,
	 *       then use other operator().
	 */
	CurveType operator()(CurveType const& X) {
		CurveType solution(X);
		CurveType PX(X); PX *= 0.; // convenient way to have X and PX of the same structure, but PX \equiv 0
		RealType dump;
		this->operator()(solution, PX, dump);
		return PX;
	}

	/**
	 * First parameter will change! There will be a curve that reaches first full step after the section,
	 * define dover the whole integration time. So if you want to save it, use a copy in curve parameter.
	 *
	 * In on_section there will be the image on the section, as close to section as possible.
	 *
	 * The parameter on_section must be initialized with the desired "length" (time interval)
	 * and "order" of the representation. It also must be compatible (length, order) with the section.
	 *
	 * Both parameters should be defined on the same Grid.
	 *
	 * In out_approachTime you will get the return time to section.
	 */
	void operator()(CurveType& curve, CurveType& on_section, RealType& out_approachTime) {
		integrateUntilSectionCrossing(curve, m_lastTimeBeforeSection, &Class::singleStep);
		findCrossingTime(curve, on_section, m_lastEpsilonTime, -1);
		curve.epsilonShift(m_dynsys, m_lastEpsilonTime, on_section);
		out_approachTime = m_lastTimeBeforeSection + m_lastEpsilonTime;
	}

	/**
	 * For other parameters and base description, see the other operator()
	 *
	 * This is more complicated than the case of ODEs, but for technical reasons
	 * it is most desirable to work with vector/matrix representations of the objects
	 * when dealing with variational equation on the coefficients.
	 *
	 * NOTE: CurveType must be defined with DataType that allows to handle extra Matrix data.
	 * 	     I have special structures for this. CHeck examples and helper classes.
	 *       Thanks to C++ lazy evaluation of templates this also works on the basic data structure
	 *       as long as you do not call the function in the program.
	 *
	 * out_x is vector representation of input curve
	 * out_Px is vector representation of the image on the section (Pcurve)
	 * out_fPx is the value of the "Vector Field" of the equation as defined in the Banach space
	 *         in practice, if c = Px is a curve in C([-tau, 0], \R^d), then out_fPx = c'.
	 *         This is theoretically correct when c is on the C^1 solutions manifold, i.e. c(0^-) == f(c)
	 *         where f is the r.h.s. of the DDE. This is true for any initial curve, if the reachTime > delay.
	 *         See works for more explanation.
	 * out_V is the variational equation solution on the coefficients, i.e. d\varphi(t_p, x)/dx
	 *       (in sense of the coefficients defining x) therefore it is of dimension
	 *       out_x.dimension() x out_x.dimension()
	 * out_DP is the dP/dx (x)  (derivative of P w.r.t. initial data, computed at point x)
	 *        out_DP is obtained from out_V by some correction on the reach time t_p(x) in the neighbourhood of x.
	 *        i.e. dP/dx (x) = d(\varphi(t_p(x), x))/dx =
	 *                       = \partial(\varphi)/\partial x (t_p(x), x) +
	 *                       + \partial(\varphi)/\partial t (t_p(x), x) * d(t_p)/dx (x)
	 *                       = V + f(Px) * dt_p/dx (x)
	 *        Please note that f(Px) is "vertical" and dt_p/dx (x) is "horizontal" vector
	 *        so their product is a full matrix of desired dimension! We see that
	 *        f(Px) is in out_fPx, V in out_V, then the last term is obtained
	 *        by differentiation of the section condition equation:
	 *
	 *        	s(\varphi(t_p(x), x)) == 0
	 *
	 *        so we get (* means scalar product in \R^M where needed, M = out_x.dimension())
	 *
	 *        	\gradient s * (V + f(Px) * dt_p/dx (x)) == 0
	 *        so
	 *        	dt_p/dx (x) = - ({\gradient s} * V) / ({\gradient s} * f(Px))
	 *
	 *        Note, that  ({\gradient s} * f(Px)) is simply a scalar, and it needs to be not equal 0
	 *        for the formula to make sense. This is usual notion of the section to be transversal
	 *        to the flow at x.
	 *
	 *  Finally, note that x, Px, fPx can be used in .set_x() to
	 *  a curve of a proper structure (important!) to get interesting data.
	 *  The proper structure is stored in Pcurve, i.e. usually, one can do:
	 *  curve_x = Pcurve; curve_x.set_x(out_x); 		// guaranteed: curve_x == curve.subcurve([-tau, 0])
	 *  curve_Px = Pcurve; curve_Px.set_x(out_Px); 		// guaranteed: curve_Px == Pcurve
	 *  curve_fPx = Pcurve; curve_fPx.set_x(out_fPx); 	// guaranteed: curve_fPx == (Pcurve)'
	 *
	 *  IMPORTANT: curve_fPx will be of one order higher than simply trying to compute Pcurve.dt()!
	 *  IMPORTANT: as .dt() loses one order. So Pcurve.dt() is not a proper way to get "force field"
	 *  IMPORTANT: on curve_Px, as the dimensions will be different (orders different).
	 *  IMPORTANT: Therefore, use curve_fPx to get the value of the "force field" on curve_Px.
	 */
	void operator()(
				CurveType& in_out_curve, CurveType& in_out_Pcurve,
				RealType& out_approachTime,
				VectorType& out_x, VectorType& out_Px, VectorType& out_fPx,
				MatrixType& out_V, MatrixType& out_DP) {
		// TODO: initialize curve derivatives? (JacData); (NOW IT IS AWKWARD FOR THE USER)

		// This is just a naming convention change (refactor)
		auto& curve = in_out_curve;
		auto& Pcurve = in_out_Pcurve;

		// we need to take this before we start to move things! To store its dimension etc.
		out_x = curve.get_x();
		m_storedInitialShape.clear();  // we store shape because I need it in one experimental procedure to propagate vectors through Variational equation.
		for (auto ijet = curve.rbegin(); ijet != curve.rend(); ++ijet)
			m_storedInitialShape.push_back((*ijet)->end() - (*ijet)->begin());

		integrateUntilSectionCrossing(curve, m_lastTimeBeforeSection, &Class::singleStepWithVariational);
		CurveType on_section = Pcurve.increasedOrder(); // we need extra order to compute fPx
		findCrossingTime(curve, on_section, m_lastEpsilonTime, -1);
		curve.epsilonShift(m_dynsys, m_lastEpsilonTime, on_section);
		out_approachTime = m_lastTimeBeforeSection + m_lastEpsilonTime;

		Pcurve = on_section.decreasedOrder();					// convert back to the good format (from increased order)
		CurveType fPcurve = on_section.dt(m_dynsys.getMap());	// get the value of the "force field", lose one order
		// now both Pcurve and fPcurve has the same order
		out_Px = Pcurve.get_x();
		out_fPx = fPcurve.get_x();
		out_V = MatrixType(out_Px.dimension(), out_x.dimension());

		extractVariationalMatrix(Pcurve, out_V);
//		JacobianStorageType V; // it will be easier to iterate over the complete structure to join it into V matrix later.
//		V.push_back(Pcurve.getValueAtCurrent().getMatrixData());
//		for (auto ijet = Pcurve.rbegin(); ijet != Pcurve.rend(); ++ijet)
//			for (auto icoeff = (*ijet)->begin(); icoeff != (*ijet)->end(); ++icoeff)
//				V.push_back(icoeff->getMatrixData());
//
//		size_type I = 0, J = 0;
//		size_type i = 0, j = 0;
//		for (auto a = V.begin(); a != V.end(); ++a){
//			J = 0;
//			for (auto b = a->begin(); b != a->end(); ++b){
//				for (i = 0; i < b->numberOfRows(); ++i)
//					for (j = 0; j < b->numberOfColumns(); ++j)
//						out_V[I+i][J+j] = (*b)[i][j];
//				J += j;
//			}
//			I += i;
//		}

		// makes a Vector compatible with the structure of fPcurve
		// (note that section might not be defined for higher orders in curve!)
		VectorType gradS = m_section.getGradient(Pcurve);
		ScalarType factor = 1.0 / (gradS * out_fPx);
		// we combine
		// 		DP = V + f(Px) * dt_p/dx (x)
		// with
		// 		dt_p/dx (x) = - ({\gradient s} * V) / ({\gradient s} * f(Px))
		// to get
		// 		V - ((({\gradient s} * f(Px))^-1) * (f(Px) * {\gradient s})) * V
		//      V - factor * C * V;
		size_type rows = out_fPx.dimension();
		size_type cols = gradS.dimension();
		MatrixType C(rows, cols);
		for (size_type i = 0; i < rows; ++i)
			for (size_type j = 0; j < cols; ++j)
				C[i][j] = out_fPx[i] * gradS[j];
		out_DP = out_V - factor * (C * out_V);
	}

	/**
	 * helper function, sets to Id the Variational matrix w.r.t. initial coefficients in curve
	 * TODO: DRY
	 */
	void setInitialV(CurveType& in_out_curve){
		size_type xcount = in_out_curve.storageDimension();
		size_type d = in_out_curve.dimension();
		MatrixType initialD(d, xcount);
		size_type I = 0;
		for (size_type i = 0; i < d; ++i) initialD[i][I+i] = 1.0;
		in_out_curve.getValueAtCurrent().setMatrix(initialD);

		for (auto ij = in_out_curve.rbegin(); ij != in_out_curve.rend(); ++ij){
			for (auto ik = (*ij)->begin(); ik != (*ij)->end(); ++ik){
				for (size_type i = 0; i < d; ++i) initialD[i][I+i] = 0.0;
				I += d;
				for (size_type i = 0; i < d; ++i) initialD[i][I+i] = 1.0;
				ik->setMatrix(initialD);
			}
		}
	}
	/**
	 * sets initial V to a given Matrix. Matrix must be of a good shape
	 * TODO: DRY
	 */
	void setInitialV(CurveType& in_out_curve, MatrixType const& V){
		size_type xcount = in_out_curve.storageDimension();
		size_type d = in_out_curve.dimension();
		MatrixType initialD(d, xcount);
		size_type I = 0;
		for (size_type i = 0; i < d; ++i)
			for (size_type j = 0; j < xcount; ++j)
				initialD[i][j] = V[I + i][j];

		in_out_curve.getValueAtCurrent().setMatrix(initialD);
		for (auto ij = in_out_curve.rbegin(); ij != in_out_curve.rend(); ++ij){
			for (auto ik = (*ij)->begin(); ik != (*ij)->end(); ++ik){
				I += d;
				for (size_type i = 0; i < d; ++i)
					for (size_type j = 0; j < xcount; ++j)
						initialD[i][j] = V[I + i][j];
				ik->setMatrix(initialD);
			}
		}
	}

	/**
	 * Experimental...
	 */
	void setCurrentV(CurveType& in_out_curve, MatrixType const& V){
		size_type d = in_out_curve.dimension();
		if (V.numberOfRows() < d) return;
		size_type xcount = V.numberOfColumns();
		MatrixType initialD(d, xcount);
		size_type I = 0;
		for (size_type i = 0; i < d; ++i)
			for (size_type j = 0; j < xcount; ++j)
				initialD[i][j] = V[I + i][j];

		in_out_curve.getValueAtCurrent().setMatrix(initialD);
		for (auto ij = in_out_curve.rbegin(); ij != in_out_curve.rend(); ++ij){
			for (auto ik = (*ij)->begin(); ik != (*ij)->end(); ++ik){
				I += d;
				if (V.numberOfRows() < I + d)
					continue;
				for (size_type i = 0; i < d; ++i)
					for (size_type j = 0; j < xcount; ++j)
						initialD[i][j] = V[I + i][j];
				ik->setMatrix(initialD);
			}
		}
	}

	DDEBasicPoincareMap& setDirection(CrossingDirection direction) { m_direction = direction; return *this; }
	CrossingDirection getDirection() { return m_direction; }

	DDEBasicPoincareMap& setMaxSteps(int maxSteps) { m_maxSteps = maxSteps; return *this; }
	int getMaxSteps() { return getMaximumSteps(); }
	int getMaximumSteps() { return m_maxSteps; }

	DDEBasicPoincareMap& setRequiredSteps(int requiredSteps) { m_requiredSteps = requiredSteps; return *this; }
	int getRequiredSteps() { return m_requiredSteps; }

	/**
	 * it extends the solution until the crossing
	 * with the section. Also, it takes into consideration the
	 * required number of steps (should be set-up before)
	 * It also moves away of section is initially curve crosses the section.
	 *
	 * WARNING: it does not check for sanity of argument as of now, so it might cause runtime-errors if the proof is not well prepared.
	 */
	void integrateUntilSectionCrossing(CurveType& curve, RealType& timeBeforeSection, SingleStepFn stepFn){
		ScalarType sign = m_section(curve);
		// m_dynsys(curve);
		(this->*stepFn)(curve);
		sign = m_section(curve);
		m_steps = 1;
		try { checkSteps(); }
		catch (std::logic_error& e) { throw rethrow("DDEBasicPoincareMap::approachSection: While leaving section: ", e); }

		//first we go to the right place - we want to leave section,
		// or if we are on the opposite side, then go slightly after section, then proceed
		while( (m_requiredSteps > 0 && m_steps < m_requiredSteps) ||  // we did not do enough steps
				(sign == 0.0) ||                				   // We are on the section
				!(													   // or section sign is not correct
						(m_direction == CrossingDirection::Both)
						|| ((m_direction == CrossingDirection::MinusPlus) && (sign < 0.0))
						|| ((m_direction == CrossingDirection::PlusMinus) && (sign > 0.0))
				)
				){
			timeBeforeSection = curve.t0();
			// m_dynsys(curve);
			(this->*stepFn)(curve);
			sign = m_section(curve);
			++m_steps;
			try { checkSteps(); }
			catch (std::logic_error& e) { throw rethrow("DDEBasicPoincareMap::approachSection: While leaving section: ", e); }
		}
		//do the iteration until section is crossed (we have already done enough steps)
		while(m_direction * sign > 0.0){
			timeBeforeSection = curve.t0();
			// m_dynsys(curve);
			(this->*stepFn)(curve);
			sign = m_section(curve);
			++m_steps;
			try { checkSteps(); }
			catch (std::logic_error& e) { throw rethrow("DDEBasicPoincareMap::approachSection: While leaving section: ", e); }
		}
	}

	/**
	 * testDirection = -1 for before section
	 * testDirection = +1 for after section
	 */
	void findCrossingTime(CurveType const& curve, CurveType& requested, RealType &crossingTime, int testDirection){
		/** binsearch for now */
		RealType lo = 0.0;
		RealType up = curve.getStep();
		RealType& epsi = crossingTime;
		while (up - lo > m_binsearchEpsilon){
			requested *= 0;
			epsi = 0.5 * (lo + up);
			curve.epsilonShift(m_dynsys, epsi, requested);
			// TODO: (NOT URGENT) this is so bad looking...
			if (m_section(requested) * m_direction * testDirection < 0.0){
				if (testDirection < 0)
					lo = epsi;
				else
					up = epsi;
			}else{
				if (testDirection < 0)
					up = epsi;
				else
					lo = epsi;
			}
		} // end while epsi
	} // end findCrossingTime

	// TODO: (FUTURE) implement findCrossingTime using Newtons algorithm (and compare with bin-search!)

	/** warning: might be time consuming, it (might) integrate the solution... it assumes curve is on section. */
	CrossingDirection detectCrossingDirection(CurveType const& curve){
		CurveType x(curve);
		m_dynsys(x);
		ScalarType sx = m_section(x);
		while (sx == 0.){
			m_dynsys(x);
			sx = m_section(x);
		}
		return sx < 0. ? CrossingDirection::PlusMinus : CrossingDirection::MinusPlus;
	}

	RealType getLastEpsilonTime() const { return m_lastEpsilonTime; }
	RealType getLastTimeBeforeSection() const { return m_lastTimeBeforeSection; }
	RealType getLastReachTime() const { return m_lastTimeBeforeSection + m_lastEpsilonTime; }
	JacobianStorageType getLastVariational() const { return m_variational; }

	int getLastStepsAfterSection() { return m_steps; }

	bool isNormalizeVariational() const { return m_normalizeVariational; }
	Class& setNormalizeVariational(bool value) { m_normalizeVariational = value; return *this; }

protected:
	DynSysType& m_dynsys;
	SectionType& m_section;
	CrossingDirection m_direction;
	int m_requiredSteps;
	int m_maxSteps;
	int m_steps;
	double m_binsearchEpsilon;
	RealType m_lastTimeBeforeSection;
	RealType m_lastEpsilonTime;
	JacobianStorageType m_variational;
	bool m_normalizeVariational;
	std::vector<size_type> m_storedInitialShape;

	void checkSteps(){
		if (m_maxSteps > 0 && m_steps > m_maxSteps){
			std::ostringstream info;
			info << "DDEBasicPoincareMap::checkSteps(): Made a maximum of " << m_steps << ". Integration stopped.";
			throw std::logic_error(info.str());
		}
	}

	void extractVariationalMatrix(
				CurveType const& curve,
				MatrixType& fullV, MatrixType& reducedV,
				std::vector<size_type> reducedShape, int p_howFar = -1) const {
		JacobianStorageType V, R; // it will be easier to iterate over the complete structure to join it into V matrix later.
		bool doReduced = reducedShape.size() != 0;
		size_type howFar = (p_howFar >= 0 ? p_howFar : curve.rend() - curve.rbegin());

		V.push_back(curve.getValueAtCurrent().getMatrixData());
		if (doReduced) R.push_back(curve.getValueAtCurrent().getMatrixData());

		size_type pMax = reducedShape.size();
		size_type pCount = 0;
		for (auto ijet = curve.rbegin(); (ijet != curve.rend()) && (pCount < howFar); ++ijet, ++pCount){
			size_type nCount = 0;
			for (auto icoeff = (*ijet)->begin(); icoeff != (*ijet)->end(); ++icoeff, ++nCount){
				V.push_back(icoeff->getMatrixData());
				if (doReduced && pCount < pMax && nCount < reducedShape[pCount])
					R.push_back(icoeff->getMatrixData());
			}
		}

		size_type allRowsCount = V.size() * curve.dimension();
		size_type redRowsCount = R.size() * curve.dimension();
		size_type columnsCount = V[0].size() * curve.dimension();

		size_type I = 0, J = 0;
		size_type i = 0, j = 0;
		fullV = MatrixType(allRowsCount, columnsCount);
		for (auto a = V.begin(); a != V.end(); ++a){
			J = 0;
			for (auto b = a->begin(); b != a->end(); ++b){
				for (i = 0; i < b->numberOfRows(); ++i)
					for (j = 0; j < b->numberOfColumns(); ++j)
						fullV[I+i][J+j] = (*b)[i][j];
				J += j;
			}
			I += i;
		}

		if (redRowsCount){
			if (&fullV == &reducedV)
				throw std::logic_error("Same storage variable used for full and reduced V. Might result is a crash.");
			reducedV = MatrixType(redRowsCount, columnsCount);
			size_type I = 0, J = 0;
			size_type i = 0, j = 0;
			for (auto a = R.begin(); a != R.end(); ++a){
				J = 0;
				for (auto b = a->begin(); b != a->end(); ++b){
					for (i = 0; i < b->numberOfRows(); ++i)
						for (j = 0; j < b->numberOfColumns(); ++j)
							reducedV[I+i][J+j] = (*b)[i][j];
					J += j;
				}
				I += i;
			}
		}

	} // END of void extractVariationalMatrix()

	void extractVariationalMatrix(CurveType const& curve, MatrixType& fullV) const {
		std::vector<size_type> noShape;
		extractVariationalMatrix(curve, fullV, fullV, noShape);
	}
};

} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_DDE_BASICPOINCARE_MAP_H_ */
