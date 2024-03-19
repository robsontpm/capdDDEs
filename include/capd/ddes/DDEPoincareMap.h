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

#ifndef _CAPD_DDES_DDE_POINCARE_MAP_H_
#define _CAPD_DDES_DDE_POINCARE_MAP_H_

#include <capd/ddes/DDECommon.h>
#include <capd/poincare/BasicPoincareMap.h>

namespace capd{
namespace ddes{

using namespace poincare;

#define DEFAULT_BINSEARCH_EPSILON 1e-14
#define DEFAULT_MAX_ITERATION 1000

/**
 * Implementation of Rigorous Poincare Map for use with DDESolvers.
 *
 * In theory should work with any compatible DDEDynSys and Section.
 */
template<typename DynSysSpec, typename SectionSpec>
class DDEPoincareMap {
public:
	typedef DDEPoincareMap<DynSysSpec, SectionSpec> Class;
	typedef DynSysSpec DynSysType;
	typedef SectionSpec SectionType;
	typedef typename DynSysSpec::CurveType CurveType;
	typedef typename CurveType::TimePointType TimePointType;
	typedef typename CurveType::VectorType VectorType;
	typedef typename CurveType::MatrixType MatrixType;
	typedef typename CurveType::ScalarType ScalarType;
	typedef typename CurveType::RealType RealType;
	typedef unsigned int size_type;

	/** copy constructor */
	DDEPoincareMap(DDEPoincareMap const& other):
		m_dynsys(other.m_dynsys), m_section(other.m_section),
		m_direction(other.m_direction),
		m_requiredSteps(other.m_requiredSteps), m_maxSteps(other.m_maxSteps), m_steps(other.m_steps),
		m_binsearchEpsilon(other.m_binsearchEpsilon)
	{}

	/**
	 * PoincareMap cannot be constructed without two main ingredients:
	 *  - (semi)dynamical system to push solutions forward (for now, we can use only @ref DDETaylorSolver for this)
	 *    This dynsys contains the information of the DDE we are solving. See documentation there.
	 *  - section to which we want to go (you do not need to be on it with your initial set)
	 *    the section is some codim-1 curve, usually a hyperplane in the space of C^n_p functions.
	 *    See the examples and docs of @ref DDEJetSection.
	 *
	 * The other settings are as in CAPD - direction of the crossing we want to check
	 * (from 'left to right' or 'right to left', or 'both' (any)). The Poincare will detect
	 * set crossing the section if it is going in the right direction. The flags are taken
	 * from capd::poincare::CorssingDirection: one of capd::poincare::MinusPlus, capd::poincare::PlusMinus,
	 * capd::poincare::Both.
	 *
	 * @param reqSteps - see @ref setRequiredSteps - this is an important parameter! But the class
	 * can safely handle it for you!
	 *
	 * @param maxSteps - see @setMaxSteps - this is self explanatory
	 *
	 * @param binsearchEpsilon - this controls a very basic algorithm of finding the crossing.
	 * The smaller the parameter the better results, but longer computations. Do not set this below 1e-14.
	 *
	 * The set is given later as input to the @ref operator().
	 */
	DDEPoincareMap(
			DynSysType& dynsys,
			SectionType& section,
			CrossingDirection direction = CrossingDirection::Both,
			int reqSteps = 0,
			int maxSteps = -1,
			double binsearchEpsilon = DEFAULT_BINSEARCH_EPSILON):
		m_dynsys(dynsys), m_section(section),
		m_direction(direction),
		m_requiredSteps(reqSteps), m_maxSteps(maxSteps), m_steps(0),
		m_binsearchEpsilon(binsearchEpsilon)
	{}

	/**
	 * in curve there is a curve that reaches just after the section
	 * in on_section there is the image on the section, as close to section as possible
	 * on_section must be initialized with the desired "length" (time interval) and "order" of the representation.
	 */
	void operator()(CurveType& curve, RealType& out_approachTime, CurveType& on_section, size_type n=1) {
		this->computePoincareMap(curve, out_approachTime, on_section, n);
	}

	/**
	 * To be compatible with standard notion and CAPD interface.
	 *
	 * It makes following assumptions:
	 *  * X and PX will have the same structure
	 *    (i.e. over same grid points and jets of the same order).
	 *
	 * NOTE: you can retrieve information on reach time and epsilon step,
	 *       but you cannot acquire the solution on the whole time.
	 *       If you need solution over the full integration time,
	 *       then use other operator().
	 */
	CurveType operator()(CurveType const& X, RealType& out_approachTime, size_type n=1) {
		CurveType solution(X);
		CurveType PX(X); PX *= 0.; // convenient way to have X and PX of the same structure, but PX \equiv 0
		RealType dump;
		this->computePoincareMap(solution, out_approachTime, PX, n);
		return PX;
	}

	/**
	 * See @ref operator() with time. This version does not give the return time.
	 */
	CurveType operator()(CurveType const& X, size_type n=1) {
		RealType dump_time;
		return this->operator()(X, dump_time, n);
	}

	/** as in CAPD Poincare */
	DDEPoincareMap& setDirection(CrossingDirection direction) { m_direction = direction; return *this; }
	/** as in CAPD Poincare */
	CrossingDirection getDirection() { return m_direction; }

	/**
	 * sets the maximal number of steps to be made. This is in contrast to CAPD, as there you only sets max time.
	 * Here, since we have fixed size steps, it is better to prescirbe maximal steps.
	 *
	 * TODO: (FUTURE) add setMaxReturnTime(Real maxReturnTime) as in capd::Poincare.
	 */
	DDEPoincareMap& setMaxSteps(int maxSteps) { m_maxSteps = maxSteps; return *this; }
	/** getter for the upper bound of the steps before an exception is thrown */
	int getMaxSteps() { return m_maxSteps; }

	/**
	 * THIS IS IMPORTANT SETTING. It is set usually automatically in the rigorous code
	 * to the 'Long Enough Integration time' (see papers), which guarantees that the enclosures
	 * are correct. Otherwise the image on the section for a given initial set might be ill-defined
	 * for some points and thus useless from CAP point of view.
	 *
	 * Use with caution!
	 */
	DDEPoincareMap& setRequiredSteps(int requiredSteps) { m_requiredSteps = requiredSteps; return *this; }
	/** see @setRequiredSteps */
	int getRequiredSteps() { return m_requiredSteps; }

	/**
	 * returns the last not full step made to reach the section.
	 * Usually it is an interval 0 < [\epsi_1, \epsi_2] < h := tau/p.
	 * It is valid for the last computation of the poincare map.
	 */
	RealType getLastEpsilonTime() const { return m_lastEpsilonTime; }
	/**
	 * returns the time before the crossing of the section begun.
	 * It is a multiple of the grid size h.
	 * It is valid for the last computation of the poincare map.
	 */
	RealType getLastTimeBeforeSection() const { return m_lastTimeBeforeSection; }
	/** returns the total time needed to return to section for the last computation of the Poincare map. */
	RealType getLastReachTime() const { return m_lastTimeBeforeSection + m_lastEpsilonTime; }

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

	/**
	 * in curve there is a curve that reaches just after the section
	 * in on_section there is the image on the section, as close to section as possible
	 * on_section must be initialized with the desired "length" (time interval) and "order" of the representation.
	 */
	void computePoincareMap(CurveType& curve, RealType& out_approachTime, CurveType& on_section, size_type n) {
		// TODO: (!!!URGENT, RETHINK) check if the section is not crossed between steps. For this I need to have a method to obtain enclosures over intervals of the curve and this should take section into consideration.
		// TODO: (!!!URGENT, RETHINK)  My first bet is to put code in Section and require from Solution to be able to produceenclosures (after the step is made - becouse then I can simply produce estimate over current time)
		// TODO: (!!!URGENT, RETHINK) this would go nicely with the imporvement to hold enclosures already computed in the representation to make better estimates.
		bool isExtraStepNeeded = false;
		RealType epsiLo, epsiUp;
		integrateUntilSectionCrossing(curve, m_lastTimeBeforeSection, isExtraStepNeeded, n);
		findCrossingTime(curve, on_section, epsiLo, -1);

		// TODO: (NOT URGENT): for some reason intervalHull from capd does not work on default capd::interval (using filib), therefore I have decided to write it myself; if capd repaired - you can chenge this to capd::intervals::intervalHull()
		auto myIntervalHull = [](RealType& a, RealType b) {
			auto lo = min(a.leftBound(), b.leftBound());
			auto up = max(a.rightBound(), b.rightBound());
			return RealType(lo, up);
		};

		if (isExtraStepNeeded){
			// TODO: (NOT URGENT) STRASZNIE SKOMPLIKOWANE ODTAD - Refactor later
			m_lastEpsilonTime = epsiLo.leftBound();
			RealType step = curve.getStep();
			RealType fullStep = RealType(0.0, 1.0) * step;
			epsiLo = RealType(epsiLo.leftBound(), RealType(curve.getStep()).rightBound());
			curve.epsilonShift(m_dynsys, epsiLo, on_section);
			CurveType extraImage = on_section;
			bool done = false;
			while (!done){
				extraImage *= 0.;
				curve.move(m_dynsys);
				RealType sign = m_section(curve);
				if (m_direction * sign < 0.){
					done = true;
					findCrossingTime(curve, extraImage, epsiUp, +1);
					epsiUp = RealType(0., epsiUp.rightBound());
					curve.epsilonShift(m_dynsys, epsiUp, extraImage);
				}else{
					curve.epsilonShift(m_dynsys, fullStep, extraImage);
					m_lastEpsilonTime += fullStep;
				}
				// no matter if we did intermediate full step, or if we done last (up) step
				// we need to hull intermediate result with what we have from lo.
				capd::ddes::hull(on_section, extraImage, on_section);
			}
			//m_lastEpsilonTime += capd::intervals::intervalHull(epsiLo, epsiUp); // TODO: (check with Daniel) filib interval does not have this???
			m_lastEpsilonTime += myIntervalHull(epsiLo, epsiUp);
			// TODO: (NOT URGENT) DOTAD REFACTOR
		}else{
			findCrossingTime(curve, on_section, epsiUp, +1);
			// m_lastEpsilonTime = capd::intervals::intervalHull(epsiLo, epsiUp);   // TODO: (check with Daniel) as before
			m_lastEpsilonTime = myIntervalHull(epsiLo, epsiUp);
			curve.epsilonShift(m_dynsys, m_lastEpsilonTime, on_section);
		}
		out_approachTime = m_lastTimeBeforeSection + m_lastEpsilonTime;
	}

	/**
	 * helper used to check if maximum number of steps was made.
	 * Throws std::logic_error in such a case
	 */
	void checkSteps() {
		if (m_maxSteps > 0 && m_steps > m_maxSteps){
			std::ostringstream info;
			info << "Made a maximum of " << m_steps << ". Integration stopped.";
			throw std::logic_error(info.str());
		}
	}

	/**
	 * it extends the solution until the crossing
	 * with the section. Also, it takes into consideration the
	 * required number of steps (should be set-up before)
	 * It also moves away of section is initially curve crosses the section.
	 *
	 * WARNING: it does not check for sanity of argument as of now, so it might cause runtime-errors if the proof is not well prepared.
	 */
	void integrateUntilSectionCrossing(CurveType& curve, RealType& timeBeforeSection, bool& isExtraStepNeeded, size_type n=1){
		ScalarType sign = m_section(curve);
		curve.move(m_dynsys);
		sign = m_section(curve);
		//std::cout << "at t = " << ScalarType(curve.t0()) << " step: " << int(curve.t0()) << " sign: " << sign << std::endl;
		m_steps = 1;
		try { checkSteps(); }
		catch (std::logic_error& e) { throw rethrow("DDEPoincareMap::approachSection: While initial leaving section: ", e); }

		for (size_type count = 0; count < n; ++count){
			//first we go to the right place - we want to leave section,
			// or if we are on the opposite side, then go slightly after section, then proceed
			while( (m_requiredSteps > 0 && m_steps < m_requiredSteps) ||  // we did not do enough steps
					(sign.contains(0.0)) ||                				   // We are on the section
					!(													   // or section sign is not correct
							(m_direction == CrossingDirection::Both)
							|| ((m_direction == CrossingDirection::MinusPlus) && (sign < 0.0))
							|| ((m_direction == CrossingDirection::PlusMinus) && (sign > 0.0))
					)
					){
				timeBeforeSection = curve.t0();
				curve.move(m_dynsys);
				sign = m_section(curve);
				//std::cout << "at t = " << ScalarType(curve.t0()) << " step: " << int(curve.t0()) << " sign: " << sign << std::endl;
				++m_steps;
				try { checkSteps(); }
				catch (std::logic_error& e) { throw rethrow("DDEPoincareMap::approachSection: While approaching section: ", e); }
			}
			//std::cout << "Section approached at " << m_steps << " with sign: " << sign << " direction is " << m_direction << std::endl;
			//do the iteration until section is crossed (we have already done enough steps)
			while(m_direction * sign > 0.0){
				timeBeforeSection = curve.t0();
				curve.move(m_dynsys);
				sign = m_section(curve);
				//std::cout << "at t = " << ScalarType(curve.t0()) << " step: " << int(curve.t0()) << " sign: " << sign << std::endl;
				++m_steps;
				try { checkSteps(); }
				catch (std::logic_error& e) { throw rethrow("DDEPoincareMap::approachSection: While crossing section after approach: ", e); }
			}
			isExtraStepNeeded = sign.contains(0.);
		}
	}

	/**
	 * testDirection = -1 for before section
	 * testDirection = +1 for after section
	 */
	void findCrossingTime(CurveType const& curve, CurveType& requested, RealType &crossingTime, int testDirection){
		/** binsearch for now */
		crossingTime = RealType(0.0, RealType(curve.getStep()).rightBound());
		RealType epsi;
		while (crossingTime.rightBound() - crossingTime.leftBound() > m_binsearchEpsilon){
			requested *= 0;
			epsi = 0.5 * (crossingTime.leftBound() + crossingTime.rightBound());
			curve.epsilonShift(m_dynsys, epsi, requested);
			// TODO: (NOT URGENT) this is so bad looking...
			if (m_section(requested) * m_direction * testDirection < 0.0)
				// I am on the good side of the section after step epsi (good - I wanted to be on this side)
				crossingTime = testDirection < 0 ? RealType(epsi.leftBound(), crossingTime.rightBound()) : RealType(crossingTime.leftBound(), epsi.rightBound());
			else
				// I am on section or on the bad side of it
				crossingTime = testDirection < 0 ? RealType(crossingTime.leftBound(), epsi.rightBound()) : RealType(epsi.leftBound(), crossingTime.rightBound());
		}
	}

	/** warning: might be time consuming, it (might) integrate the solution... it assumes curve is on section. */
	CrossingDirection detectCrossingDirection(CurveType const& curve){
		CurveType x(curve);
		x.move(m_dynsys);
		ScalarType sx = m_section(x);
		while (sx.contains(0.)){
			x.move(m_dynsys);
			sx = m_section(x);
		}
		return sx < 0. ? CrossingDirection::PlusMinus : CrossingDirection::MinusPlus;
	}

	/** TODO: (FUTURE): assign operator hidden now for DDEPoincareMap */
	Class& operator=(Class const& other){
		throw std::logic_error("DDEPoincareMap::operator=(): Not Implemented Yet");
	}
};

} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_DDE_POINCARE_MAP_H_ */
