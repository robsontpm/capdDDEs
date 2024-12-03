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

#ifndef _CAPD_DDES_DISCRETEDELAYSFUNCTIONALMAP_H_
#define _CAPD_DDES_DISCRETEDELAYSFUNCTIONALMAP_H_

#include <capd/ddes/DDECommon.h>
#include <capd/ddes/FunctionalMap.h>
#include <capd/basicalg/TypeTraits.h>
#include <stdexcept>
#include <string>
#include <capd/fadbad/fadbad.h>
#include <capd/fadbad/tadiff.h>
#include <capd/fadbad/fadiff.h>

namespace capd{
namespace ddes{

/**
 * A class to represent right hand side of the DDE of the form:
 *
 *     x' = f(t, x(t), x(t-tau_1), ..., x(t-tau_m))
 *
 * NOTE: This is a demo version of what I want to achieve with the DDEs code interface
 *       I assume it should be similar to ODEs in original CAPD.
 *       I assume that any RHS map f of the equation x'(t) = f(x_\tau)
 *       should be representable with similar interface (i.e. the function
 *       accepts some SolutionCurve and can produce Jet at current time in the
 *       solution defined with the equation (function computeDDECoefficients)
 *       I do not assume anything about SolutionCurve except, it can produce
 *       Forward Taylor Jets of itself at a given point (that is used in the equation)
 *       and that it can return its value at a current time.
 *
 *       In the implementation below, I assume that I can ask about SOlutionCUrve and its Jet
 *       at $t-\tau$, $t$ is current time. Then I construct an AD jet of x(t-\tau) and pass
 *       it to standard CAPD map.
 *
 * NOTE: I do not know if integral equations (i.e. x'(t) = F(\int_{t-\tau}^t g(x(s)) ds) could be described in this manner...
 *       But I hope it could be done, at least if g is linear in x (e.g. g = Id).
 * NOTE: in fact, functional map could have operators that are templates, but we assume that
 *       the argument curve might have some requirements when computing, for example
 *       it could determine only the delays at the specific grid points (because of the loss of regularity,
 *       we are not able to allways give jets over some general intervals, see new notes).
 *       Therefore, we supply the Map with appropriate CurveType that it can manipulate.
 */
template<
	typename FinDimMapSpec,
	typename SolutionCurveSpec,
	typename JetSpec=typename SolutionCurveSpec::JetType
>
class DiscreteDelaysFunctionalMap : public DDERigorousFunctionalMap<SolutionCurveSpec, JetSpec> {
public:
	typedef DDERigorousFunctionalMap<SolutionCurveSpec, JetSpec> BaseClass;
	typedef DiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec> Class;
	typedef SolutionCurveSpec CurveType;
	typedef JetSpec JetType;
	typedef FinDimMapSpec MapType;
	typedef typename BaseClass::RealType RealType;
	typedef typename BaseClass::TimePointType TimePointType;
	typedef typename BaseClass::MatrixType MatrixType;
	typedef typename BaseClass::VectorType VectorType;
	typedef typename BaseClass::ScalarType ScalarType;
	typedef typename BaseClass::size_type size_type;
	typedef typename BaseClass::DataType DataType;
	typedef typename BaseClass::ValueStorageType ValueStorageType;
	typedef typename BaseClass::VariableStorageType VariableStorageType;
	typedef typename BaseClass::JacobianStorageType JacobianStorageType;
	/** storage type for discrete number of delays */
	typedef std::vector< TimePointType > DelayStorageType;
	/** iterator over this map will iterate over delays! */
	typedef typename DelayStorageType::const_iterator const_iterator;
	/** iterator over this map will iterate over delays! */
	typedef typename DelayStorageType::iterator iterator;

    template <typename OtherSolutionSpec, typename OtherJetSpec=typename OtherSolutionSpec::JetType>
    struct rebind { typedef DiscreteDelaysFunctionalMap<FinDimMapSpec, OtherSolutionSpec, OtherJetSpec> other; };

	template<typename OtherDiscreteDelaysFunctionalMap>
	DiscreteDelaysFunctionalMap(const OtherDiscreteDelaysFunctionalMap& orig): m_map(orig.m_map), m_delays(orig.m_delays){}
	DiscreteDelaysFunctionalMap(const DiscreteDelaysFunctionalMap& orig): m_map(orig.m_map), m_delays(orig.m_delays){}
	DiscreteDelaysFunctionalMap(){ this->checkMapSignature(); }
	DiscreteDelaysFunctionalMap(const TimePointType& tau) { m_delays.push_back(tau); this->checkMapSignature(); }
	DiscreteDelaysFunctionalMap(const MapType& f): m_map(f) { this->checkMapSignature(); }
	DiscreteDelaysFunctionalMap(const MapType& f, const TimePointType& tau): m_map(f) { m_delays.push_back(tau); this->checkMapSignature(); }
	DiscreteDelaysFunctionalMap(const MapType& f, std::vector<TimePointType> delays): m_map(f), m_delays(delays){ this->checkMapSignature(); }
	DiscreteDelaysFunctionalMap(const MapType& f, int n, TimePointType delays[]): m_map(f) { for (size_type i = 0; i < n; ++i) m_delays.push_back(delays[i]); this->checkMapSignature(); }
	/** standard */
	virtual ~DiscreteDelaysFunctionalMap() {}

	/** iterator over delays */
	iterator begin() { return m_delays.begin(); }
	/** iterator over delays */
	iterator end() { return m_delays.end(); }
	/** iterator over delays */
	const_iterator begin() const { return m_delays.begin(); }
	/** iterator over delays */
	const_iterator end() const { return m_delays.end(); }

	/** output dimension of the internal map */
	size_type imageDimension() const { return m_map.imageDimension(); }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	size_type dimension() const { return (m_delays.size() + 1) * m_map.imageDimension(); }
	/** number of delays */
	size_type delaysCount() const { return m_delays.size(); }

	/** Implementation of virual func. See Interface docs. */
	VectorType operator()(const TimePointType& t, const CurveType& x) const {
		// checkCurveDimension(x, "operator()"); // TODO: (NOT URGENT) use from Interface later
		VectorType args(this->dimension());
		VectorType x0 = x.eval(t);
		auto iargs = args.begin();
		for (auto ix0 = x0.begin(); ix0 != x0.end(); ++ix0, ++iargs) *iargs = *ix0;
		for (auto idelay = m_delays.begin(); idelay != m_delays.end(); idelay++){
			VectorType xtau = x.eval(t - *idelay);
			for (auto ixtau = xtau.begin(); ixtau != xtau.end(); ++ixtau, ++iargs) *iargs = *ixtau;
		}
		VectorType result(imageDimension());
		m_map(RealType(t), args, result);
		return result;
	}

	/**
	 * Implementation of virual func. See Interface docs.
	 *
	 * TODO: (NOT URGENT) explain what is stored in output for DiscreteDelaysMap...
	 *
	 * Warning: dt should be enclosure for the step [0, h]
	 */
	virtual void collectComputationData(
					const TimePointType& t0, const TimePointType& th,		//input
					const RealType& dt, const CurveType& x, 	 			// input
					VariableStorageType& out_u, ValueStorageType& out_encl, // output
					size_type& out_admissible_order) const {				// output
		out_u.clear();
		out_encl.clear();
		size_type d = this->imageDimension();
		// we collect jets and compute admissible order
		// TODO: (NOT URGENT) rewrite using pointers to Jets (will be faster)!
		out_admissible_order = 40;
		std::vector<JetType> jets;
		for (auto tau = begin(); tau != end(); ++tau){
			jets.push_back(x.jet(t0 - *tau));
			size_type k = jets.back().order();
			if (k + 1 < out_admissible_order) out_admissible_order = k + 1;
		}
		// Then we fill in at u[0] the value at current time (always present in u, by assumption)
		out_u.push_back(x.getValueAtCurrent());
		out_encl.push_back(VectorType(d)); // just make room for future.
		// and then all jets, order by order (in this way it is easier to fill AutoDiff later)
		ValueStorageType pastData; ///< this is used in computation of rough enclosure, it only need to be the enclosures of value of solution in the past.
		for (auto jit = jets.begin(); jit != jets.end(); ++jit){
			out_u.push_back((*jit)[0]);
			VectorType enc = jit->evalAtDelta(dt);
			out_encl.push_back(enc);
			pastData.push_back(enc);
		}
		for (size_type k = 1; k < out_admissible_order; ++k){
			for (auto jit = jets.begin(); jit != jets.end(); ++jit){
				out_u.push_back((*jit)[k]);
				out_encl.push_back(jit->evalCoeffAtDelta(k, dt));
			}
		}
		// add extra data necessary in current Solvers, that is, global estimate of one higher order.
		for (auto ijet = jets.begin(); ijet != jets.end(); ++ijet)
			out_encl.push_back(ijet->evalCoeffAtDelta(out_admissible_order, dt));

		// run find enclosure, we have everything we need.
		VectorType dump_Z(d); // it is not used
		// in a rare event, that we are computing for 0 length interval...
		if (t0 == th){
			out_encl[0] = out_u[0];
		}else{
			findRoughEnclosure(t0, dt, pastData, out_u[0], dump_Z, out_encl[0]);
		}
	}

	/**
	 * Implementation of virual func. See Interface docs.
	 *
	 * TODO: (NOT URGENT) explain what is stored in output for DiscreteDelaysMap...
	 *
	 * Warning: dt should be enclosure for the step [0, h]
	 */
	virtual void collectComputationData(
					const TimePointType& t0, const TimePointType& th,	//input
					const RealType& dt, const CurveType& x, 	 		// input
					VariableStorageType& out_u, 						// output
					size_type& out_admissible_order) const {			// output
		// TODO: (NOT URGENT) DRY with the other  collectComputationData() method.
		out_u.clear();
		// we collect jets and compute admissible order
		// TODO: (NOT URGENT) rewrite using pointers to Jets (will be faster)!
		out_admissible_order = 40;
		std::vector<JetType> jets;
		for (auto tau = begin(); tau != end(); ++tau){
			jets.push_back(x.jet(t0 - *tau));
			size_type k = jets.back().order();
			if (k + 1 < out_admissible_order) out_admissible_order = k + 1;
		}
		// Then we fill in at u[0] the value at current time (always present in u, by assumption)
		out_u.push_back(x.getValueAtCurrent());
		// and then all jets, order by order (in this way it is easier to fill AutoDiff later)
		ValueStorageType pastData; ///< this is used in computation of rough enclosure, it only need to be the enclosures of value of solution in the past.
		for (auto jit = jets.begin(); jit != jets.end(); ++jit){
			out_u.push_back((*jit)[0]);
			VectorType enc = jit->evalAtDelta(dt);
			pastData.push_back(enc);
		}
		for (size_type k = 1; k < out_admissible_order; ++k)
			for (auto jit = jets.begin(); jit != jets.end(); ++jit)
				out_u.push_back((*jit)[k]);
	}

	/** Implementation of virual func. See Interface docs. */
	virtual void computeDDECoefficients(
				const RealType& t0, const ValueStorageType& u, 	// input
				ValueStorageType& coeffs) const {					// output
		// we assume u[0] is always present, and all other u[i] has the same
		checkDimension(u[0]);
		// we assume that coeffs has a proper size for the computation.
		// for example, set when user took data u with collectComputationData()
		size_type resultOrder = coeffs.size() - 1;

		typedef fadbad::T<ScalarType> AutoDiffTScalar; 		// Auto Diff T (Taylor) scalar type
		typedef typename VectorType::template rebind<AutoDiffTScalar>::other TADVectorType;
		TADVectorType args(dimension(), false);
		// propagate arguments in a given order - first value at t0
		// then jet of \tau_1, then \tau_2, etc...
		// this is change from the initial paper!, where I have done delay terms first
		// this is coherent with what is in the general DDE literature though...
		// we assume this order in u (because we produce itd this way in collectComputationData()!).
		auto current_u = u.begin();
		auto iarg = args.begin();
		for (size_type idelay = 0; idelay <= delaysCount(); ++idelay, ++current_u)
			for (auto iu = current_u->begin(); iu != current_u->end(); ++iu, ++iarg){
				(*iarg)[0] = *iu;
			}

		// we are ready to go for the recurrent relation on coefficients.
		coeffs[0] = u[0]; // by definition, see paper
		for (size_type k = 1 ; k <= resultOrder; ++k){
			// eval Autodiff routine so that we would be able to extract v[k-1] = F^{[k-1]}(...)
			TADVectorType v(imageDimension(), false);
			m_map(t0, args, v);
			// propagate AD jets and output.
			// this will be F_k from the latest paper on DDE integration.
			VectorType Fk(imageDimension());
			auto iarg = args.begin();
			for (size_type i = 0; i < imageDimension(); ++i, ++iarg){
				v[i].eval(k-1);
				Fk[i] = v[i][k-1];
				// propagate computed recurrent formula into args with proper scaling!
				(*iarg)[k] = Fk[i] / ScalarType(k);
			}
			// this is the right formula from the last paper!
			coeffs[k] = Fk / ScalarType(k);
			// do not propagate anything else if we are past the requestedOrder (it might cause errors)
			if (k == resultOrder) break;
			// propagate higher order coefficients in the past to the args Vector.
			// the recurrent Taylor formula is already there.
			for (size_type idelay = 1; idelay <= delaysCount(); ++idelay, ++current_u){
				for (auto iu = current_u->begin(); iu != current_u->end(); ++iu, ++iarg)
					(*iarg)[k] = *iu;
			}
		} /// for k loop
	} /// computeDDECoefficients

	/**
	 * see Interface docs.
	 */
	virtual void computeDDECoefficients(
				const RealType& t0, const ValueStorageType& u, 			// input
				ValueStorageType& coeffs, JacobianStorageType& Du) const {	// output
		// we assume u[0] is always present, and all other u[i] has the same
		checkDimension(u[0]);
		// we assume that coeffs has a proper size for the computation.
		// for example, set when user took data u with collectComputationData()
		size_type resultOrder = coeffs.size() - 1;
		// d = space dimension of a single variable in the map / input data in ValueStorageType (must be the same).
		size_type d = this->imageDimension();
		size_type numDelays = delaysCount();
		// those are common matrices of dimension M(d, d),
		MatrixType Id(d, d); Id.setToIdentity();
		MatrixType Zero(d, d);

		// This is obvious
		size_type partialsCount = u.size();
		// in partialsDimension we keep how many actual coefficients from all of the Jets used to
		// compute rhs of the equation. As we are computing jets up to coeffs.order() (incl.).
		// then we are going to use only coefficients of all jets up to coeffs.order() - 1 (incl).
		// Thus, the formula. Please note in all those computations we have 0,..,resultOrder-1
		// Taylor coefficients in any jet used. Therefore the number of Taylor coeffs at each
		// jet used is resultOrder. This is trivial, but I want this to be reiterated for readability.
		size_type partialsDimension = d * u.size();
		// we make room to store Jacobian matrices. Each matrix will be M(d,d).
		// Du is shorthand for $D_u coeffs$
		// in Du[i][j] we will have \partial coeffs[i]
		Du.clear(); Du.resize(resultOrder + 1);
		for (auto iDu = Du.begin(); iDu != Du.end(); ++iDu)
			iDu->resize(partialsCount);

		typedef fadbad::F<ScalarType> FADScalar; 		// Auto Diff Forward scalar type (i.e. w.r.t. initial coefficients)
		typedef fadbad::T<FADScalar> TFADScalar;		// Auto Diff Taylor type with option to compute derivatives wrt initials
		typedef typename VectorType::template rebind<TFADScalar>::other TFADVectorType;
		// this is Automatic Differentiation both w.r.t t but also w.r.t. coefficients.
		// we will use it to automatically eval the r.h.s. of the equation with all
		// partial derivatives when recursively compute Taylor coefficients at current t0
		TFADVectorType args(dimension(), false);

		// we will successively fill-in args as the recursion will go on
		auto current_u = u.begin();
		auto iarg = args.begin();
		// this will hold current number of parameter (among 0,.., partialsDimension-1)
		// we need those for AutoDiff types (they use this kind of numbering)
		size_type partialAD = 0;
		// this will hold current number of variable in u (among 0,..,partialsCount-1, each variable is d dimensional)
		// we need those for our purpose (to match scalar variables to appropriate entries in vector/matrices)
		size_type partialu = 0;
		// we fill-in initial values for the recurrence, i.e. the value at t0 stored
		// in u[0] and values at values in the past, at given delays.
		// We mark each of those values as dependent only on self. This will give us
		// possibility to automatically compute dependence of the formulas on them
		// i.e. the value \frac{\partial coeffs}{\partial u[partialu]}
		// Note: here k = 0, and we iterate idelay from 0. Later, it will be from 1. See docs there.
		for (size_type idelay = 0; idelay <= numDelays; ++idelay, ++current_u, ++partialu){
			for (auto iu = current_u->begin(); iu != current_u->end(); ++iarg, ++iu, ++partialAD){
				// (*iarg)[0] is [0]-th order Taylor coefficient
				(*iarg)[0] = *iu;
				// and we mark it as dependent only on self among all 0,..,partialsDimension-1 variables
				// (i.e. \frac{\partial{(*iarg)[0]}}{\partial{(*iarg)[0]}} = 1, and 0 otherwise)
				(*iarg)[0].diff(partialAD, partialsDimension);
			}
		}

		// we compute coeffs[0] and derivatives...
		coeffs[0] = u[0];		// coeffs[0] = x(t0) by definition of (see paper)
		Du[0][0] = Id; 			// \frac{Dcoeffs_{[0]}}{Dx(t0)} = Id (as a matrix M(d, d)) (since coeffs[0] = x(t0))
		// Next... \frac{Dcoeffs_{[0]}}{{Djet_m}_{[mk]}} = 0 (as a matrix M(d, d)) (since coeffs[0] = x(t0))
		for (size_type par = 1; par < partialsCount; ++par)
			Du[0][par] = Zero;

		// then we compute higher order coeffs, with recurrent Taylor formula (see paper)
		for (size_type k = 1; k <= resultOrder; ++k){
			// eval AutoDiff routine so we would be able to extract v[.][k-1] = F^{[k-1]}(...)
			// together with partial derivatives with respect to coefficients in the past.
			TFADVectorType v(imageDimension(), false);
			m_map(t0, args, v);

			// propagate AD jets and output.
			VectorType Fk(imageDimension());
			auto iarg = args.begin();
			for (size_type i = 0; i < imageDimension(); ++i, ++iarg){
				v[i].eval(k-1);
				Fk[i] = v[i][k-1].val(); 					// here we need to .value(), to convert from FAD to VectorType
				(*iarg)[k] = v[i][k-1] / ScalarType(k);		// here we need to propagate things with all partial derivatives! So we copy from v (appropriate AutoDiff Type)
			}

			// we output the value, properly scaled
			coeffs[k] = Fk / ScalarType(k);
			// and partial derivatives w.r.t all previous coeffs.
			MatrixType parDu(d, d);
			for (size_type par = 0; par < partialu; ++par){
				// we copy the appropriate partial derivatives into a local Matrix, which will be later
				// Du[k][par] = \frac{\partial coeffs[k]}{\partial u[par]}
				for (size_type i = 0; i < d; ++i)								// we use formula par*d+j
					for (size_type j = 0; j < d; ++j)							// to obtain the local index
						parDu[i][j] = v[i][k-1].d(par * d + j) / ScalarType(k);	// among all partialsCount variables
				Du[k][par] = parDu;
			}
			// mark that this order coeff does not depend on higher order u's
			for (size_type par = partialu; par < partialsCount; ++par)
				Du[k][par] = Zero;


			// this if is to ensure that we do not copy, if there is no more values left in u.
			if (k == resultOrder) break;

			// We have now propagated dependence of the solution on the value at
			// present time t0 to Du. Now, finally, we propagate higher order
			// coefficients in the past to args for the next level of recurrence.
			for (size_type idelay = 1; idelay <= numDelays; ++idelay, ++current_u, ++partialu){
				for (auto iu = current_u->begin(); iu != current_u->end(); ++iarg, ++iu, ++partialAD){
					// (*iarg)[0] is [0]-th order Taylor coefficient
					(*iarg)[k] = *iu;
					// and we mark it as dependent only on self among all 0,..,partialsDimension-1 variables
					// (i.e. \frac{\partial{(*iarg)[0]}}{\partial{(*iarg)[0]}} = 1, and 0 otherwise)
					(*iarg)[k].diff(partialAD, partialsDimension);
				}
			}
			// TODO: (FUTURE) rewrite on iterators only, it is possible for sure!
			// TODO: (FUTURE) in such a case probably, I would not need many of the *Count variables.

			// now we are ready for the next level of recursion.
		} // for k loop
	} // function computeDDECoefficients


	using BaseClass::operator();
	using BaseClass::computeDDECoefficients;
protected:
	MapType m_map;
	DelayStorageType m_delays;

	static const double ROUGH_MULFACT;
	static const double ROUGH_REFINE_FACTOR;
	static const double ROUGH_TRIALSTEP;
	static const size_type ROUGH_LIMIT;

	void checkDimension(size_type d) const {
		if (d != this->imageDimension()){
			std::ostringstream info;
			info << "DiscreteDelaysFunctionalMap::checkDimension(size_type): incompatible dimension: ";
			info << "is " << d << ", should be " << this->imageDimension() << ".";
			throw std::logic_error(info.str());
		}
	}

	template<typename AnyVector>
	void checkDimension(AnyVector const & v) const {
		try { checkDimension(v.dimension()); }
		catch (std::logic_error& e) { throw rethrow("DiscreteDelaysFunctionalMap::checkDimension(AnyVector)", e); }
	}

	void checkMapSignature() const {
		if (m_map.dimension() != this->dimension()){
			std::ostringstream info;
			info << "DiscreteDelaysFunctionalMap: An internal map has incompatible dimensions: ";
			info << "is " << m_map.dimension() << ", should be " << this->dimension() << ".";
			info << "Number of delays is " << this->m_delays.size() << ", output dimension is" << this->imageDimension() << ".";
			throw std::logic_error(info.str());
		}
	}

	/**
	 * Finds a roughEnclosure Z and Y such that:
	 * Y = z + h * f(t0+h, Z, enclosure_of_past_data(x)) \subset \interior Z
	 * h = [t0, t0+dt] and z \in Z, where set Z is chosen heuristically.
	 * The theorem guarantees that solution to DDE x(t)
	 * exists for initial data z (at t0) for any t \in h and x(t) \in Y.
	 * The enclosure_of_past_data(x) is the enclosure of the value of the solution in the past
	 * times over intervals of length h. See current paper for better explanation.
	 *
	 * Note: the code is direct reimplementation of CAPD code for ODEs
	 */
	void findRoughEnclosure(
				RealType const & t0, RealType const & dt,					// input
				ValueStorageType const & pastData, VectorType const & z, 	// input
				VectorType& out_Z, VectorType& out_Y) const {				// output

		// some helper constants
		size_type d = this->imageDimension();
		RealType trialStep = ScalarType(-ROUGH_TRIALSTEP, 1.0 + ROUGH_TRIALSTEP) * dt;
		RealType h = RealType(0.0, 1.0) * dt;
		VectorType Small(d);
		for (size_type i = 0; i < d; ++i)
			Small[i] = ScalarType(-1.0, 1.0) * 1e-21;

		// we will use x to evaluate the function F on it
		VectorType x(this->dimension()); auto ix = x.begin();
		// we will need to update current time t0 data each time later...
		for (auto ival = z.begin(); ival != z.end(); ++ival, ++ix)
			(*ix) = (*ival);
		// but we only need to copy past data once.
		for (auto ipast = pastData.begin(); ipast != pastData.end(); ++ipast)
			for (auto ival = ipast->begin(); ival != ipast->end(); ++ival, ++ix)
				(*ix) = (*ival);

		VectorType Fx(d);
		m_map(t0 + h, x, Fx);					// eval F(x) to Fx
		out_Z = z + trialStep * Fx + Small;  	// initial very rough guess

		bool found = false;
		size_type counter = 0; // counts until max number of heuristic trials is done (see ROUGH_LIMIT)
		while ((!found) && (counter++ < ROUGH_LIMIT)){
			// copy Z as the first variable in x
			auto ix = x.begin();
			for (auto ival = out_Z.begin(); ival != out_Z.end(); ++ival, ++ix)
				(*ix) = (*ival);
			m_map(t0 + h, x, Fx);
			out_Y = z + h * Fx;
			found = true;
			for (size_type i = 0; i < d; ++i){
				if(!(out_Y[i].subsetInterior(out_Z[i]))){
					found = false;
					// Z to small in the i-th coordinate, grow Y and use as Z for another try
					out_Z[i] = out_Y[i];
					ScalarType s;
					out_Z[i].split(s);
					s = ROUGH_MULFACT * s;
					out_Z[i] += s;
				}
			}
		}

		if(!found)
			throw std::logic_error("DiscreteDelaysFunctionalMap::findRoughEnclosure(): Rough enclosure not found!");

		// try to refine estimates
		VectorType save_Y(d);
		VectorType save_Z(d);

		// first (new) strategy to refine
		found = true;
		counter = 0;
		do{
			save_Y = out_Y;
			save_Z = out_Z;
			out_Z = out_Y; // Try current Y as Z and compute new Y
			// this is very old version
			auto ix = x.begin();
			for (auto ival = out_Z.begin(); ival != out_Z.end(); ++ival, ++ix)
				(*ix) = (*ival);
			m_map(t0 + h, x, Fx);
			out_Y = z + h * Fx;

			for (size_type i = 0; i < d; ++i)
				if(!(out_Y[i].subsetInterior(out_Z[i])))
					found = false;

		}while(found && counter++ < ROUGH_LIMIT);
		// save_y must be good candidate for Y, because
		// at the first entrance to do it was{ save_Y = Y
		// and we exited }while() either if
		// counter >= limit and then (z+h*f(z)).subsetInterior(z) -> satisfies condition for Rough Enclosure
		// or save_y = y for y = (z+h*f(z)) which satisfied y.subsetInterior(z) in previous iteration
		out_Y = save_Y;
		out_Z = save_Z;

		// second (older) strategy to refine TODO: (FUTURE) rethink: run both strategies "in parallel" then compare / intersect, rather than run in  consecutive
		found = true;
		counter = 0;
		do{
			save_Y = out_Y;
			save_Z = out_Z;
			for (size_type i = 0; i < d; ++i){
				out_Z[i] = (ROUGH_REFINE_FACTOR * out_Z[i] + (1-ROUGH_REFINE_FACTOR) * out_Y[i]);
			}
			auto ix = x.begin();
			for (auto ival = out_Z.begin(); ival != out_Z.end(); ++ival, ++ix)
				(*ix) = (*ival);
			m_map(t0 + h, x, Fx);
			out_Y = z + h * Fx;

			for (size_type i = 0; i < d; ++i)
				if(!(out_Y[i].subsetInterior(out_Z[i])))
					found = false;

		}while(found && counter++ < ROUGH_LIMIT);

		// see comment above why this gives good out_Y and out_Z.
		out_Y = save_Y;
		out_Z = save_Z;

	} // function findRoughEnclosure()

};


} // namespace ddes
} // namespace capd


#endif /* _CAPD_DDES_DISCRETEDELAYSFUNCTIONALMAP_H_ */
