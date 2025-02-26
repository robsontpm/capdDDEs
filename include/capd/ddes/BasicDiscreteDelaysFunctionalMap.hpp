/*
 * This file constitutes part of DDEs rigorous integration framework developed
 * in PhD Thesis under supervision of prof. Piotr Zgliczynski:
 *
 * 		"Rigorous Integration of Delay Differential Equations", Jagiellonian University, 2015
 *
 * When using in scientific work please consider citing my articles (preferred),
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

#ifndef _CAPD_DDES_BASICDISCRETEDELAYSFUNCTIONALMAP_HPP_
#define _CAPD_DDES_BASICDISCRETEDELAYSFUNCTIONALMAP_HPP_

#include <capd/ddes/BasicDiscreteDelaysFunctionalMap.h>


namespace capd{
namespace ddes{

template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>						// template spec
typename BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::VectorType	// return type
BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::operator()(			// method name
			const TimePointType& t,																	// args
			const CurveType& x) const {																// args...
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


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>								// template spec
void BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::collectComputationData( 	// method name
				const TimePointType& t0, const TimePointType& th,	//input
				const RealType& dt, const CurveType& x, 	 		// input
				VariableStorageType& out_u, 						// output
				size_type& out_admissible_order) const {			// output
	// TODO: (NOT URGENT) DRY with the other  collectComputationData() method.
	// TODO: (NOT URGENT) I'm not sure it will be possible, but think about it!
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
	for (size_type k = 0; k < out_admissible_order; ++k)
		for (auto jit = jets.begin(); jit != jets.end(); ++jit)
			out_u.push_back((*jit)[k]);
}


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>								// template spec
void BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::computeDDECoefficients(
			const RealType& t0, const ValueStorageType& u, 		// input
			ValueStorageType& coeffs) const {					// output
	// we assume u[0] is always present, and all other u[i] has the same dimension
	// we only check once. If other have different, then the CAPD will throw a generic exception
	checkDimension(u[0].dimension());
	// we assume that coeffs has a proper size for the computation.
	// for example, set when user took data u with collectComputationData()
	// if not, then some other generic exceptions can happen.
	size_type resultOrder = coeffs.size() - 1;

	typedef fadbad::T<ScalarType> AutoDiffTScalar; 		// Auto Diff T (Taylor) scalar type
	typedef typename VectorType::template rebind<AutoDiffTScalar>::other TADVectorType;
	TADVectorType args(dimension(), false);
	AutoDiffTScalar targ(t0); targ[1] = 1; // we have targ(t) = t, targ'(t) = 1, targ^{(k > 1)} = 0;
	// propagate arguments in a given order - first value at t0
	// then jet of x at t0-\tau_1, then t0-\tau_2, etc...
	// please note that u is a Value storage (not variable storage), that is, it
	// contains by assumption only vectors (interval hulls, not e.g. Lohner/Affine sets with a structure!)
	// moreover, we have u[0] = x(t0) \subset \R^d, u[1] = x(t-\tau_1), ..., x[m] = x(t-\tau_m).
	// we assume this order in u (because we produce it this way in collectComputationData()!).
	// We will successively fill-in TFAD args below as the recursion will go on
	// using this iterator to u variable.
	// NOTE: important fact!
	//		this iterator only goes forward, never it should get reset!
	// 		because the older data was inserted into 'args', and we only need to fill new data,
	// This is why we expect the u to be organized in this way!
	auto current_u = u.begin();
	auto iarg = args.begin();
	// here we start from 0, as we are extracting the value at t0
	// later in the code you might notice we start from 1.
	for (size_type idelay = 0; idelay <= delaysCount(); ++idelay, ++current_u)
		for (auto iu = current_u->begin(); iu != current_u->end(); ++iu, ++iarg)
			(*iarg)[0] = *iu;

	// now we are ready to go for the recurrent relation on coefficients
	// to produce first derivative. Later, after each loop, we will have
	// all derivatives up to order n and we are ready to compute n+1'st.
	coeffs[0] = u[0]; // by definition, see paper
	for (size_type k = 1 ; k <= resultOrder; ++k){
		// eval Autodiff routine so that we would be able to
		// extract v[k-1] = F^{[k-1]}(...), as in the paper.
		// The constructor with false here is important to not initialize v data
		// with anything, each position in v will be just a fresh AD variable.
		// in (2024) there was an error introduced by use of other constructor, when CAPD changed.
		TADVectorType v(imageDimension(), false);
		m_map(targ, args, v);
		// propagate AD jets and output.
		// this will be F_k from the latest paper on DDE integration (FoCM 2024).
		VectorType Fk(imageDimension());
		auto iarg = args.begin();
		for (size_type i = 0; i < imageDimension(); ++i, ++iarg){
			v[i].eval(k-1);
			Fk[i] = v[i][k-1];
			// propagate computed recurrent formula into args with proper scaling!
			(*iarg)[k] = Fk[i] / ScalarType(k);
		}
		// this is the right formula from the last paper!
		// it is the value of x^{[k]}(t0) here.
		coeffs[k] = Fk / ScalarType(k);
		// do not propagate anything else if we are past the
		// requestedOrder (it might cause errors if we do!)
		if (k == resultOrder) break;
		// propagate higher order coefficients in the past to the args Vector.
		// note that we have already updated the vector with value at t0,
		// therefore we are starting from idelay = 1.
		for (size_type idelay = 1; idelay <= delaysCount(); ++idelay, ++current_u){
			for (auto iu = current_u->begin(); iu != current_u->end(); ++iu, ++iarg)
				(*iarg)[k] = *iu;
		}
	} /// for k loop
} /// computeDDECoefficients


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>								// template spec
void BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::computeDDECoefficients(
			const RealType& t0, const ValueStorageType& u, 				// input
			ValueStorageType& coeffs, JacobianStorageType& Du) const {	// output
	 // TODO: (NOT URGENT, FUTURE, RETHINK) there is some redundancy in the code (with the other version),
	 // TODO: (NOT URGENT, FUTURE, RETHINK) but I don't believe we can make it more DRY without sacrificing speed/RAM usage...

	// we assume u[0] is always present, and all other u[i] has the same dimension
	// we only check once. If other have different, then the CAPD will throw a generic exception
	checkDimension(u[0].dimension());
	// we assume that coeffs has a proper size for the computation.
	// for example, set when user took data u with collectComputationData()
	// if not, then some other generic exceptions can happen.
	size_type resultOrder = coeffs.size() - 1;
	// d = space dimension of a single variable in the
	// map / input data in ValueStorageType (must be the same).
	size_type d = this->imageDimension();
	size_type numDelays = delaysCount();
	// those are common matrices of dimension M(d, d),
	MatrixType Id(d, d); Id.setToIdentity();
	MatrixType Zero(d, d);

	// This is obviously the number of d-dimensional variables
	// that the coefficients at t0 depend on.
	// Those will be basically 1 + num-delays * resultOrder
	size_type partialsCount = u.size();
	// in partialsDimension we keep how many actual coefficients from all of the Jets used to
	// compute rhs of the equation. As we are computing jets up to coeffs.order() (incl.).
	// then we are going to use only coefficients of all jets up to coeffs.order() - 1 (incl).
	// Thus, the formula. Please note in all those computations we have 0,..,resultOrder-1
	// Taylor coefficients in any jet used. Therefore the number of Taylor coeffs at each
	// jet used is resultOrder. This is trivial, but I want this to be reiterated for readability.
	size_type partialsDimension = d * u.size();
	// we make room to store Jacobian matrices. Each matrix will be of shape d x d.
	// Du is shorthand for $D_u coeffs = \frac{\partial coeffs_j}{\partial u_i}$
	// (the last notation means we are iterating over all possibilities and we put it into a matrix)
	// in Du[i][j] we will have $\frac{\partial coeffs[i]}{\partial u[j]}$ (a d x d Matrix)
	Du.clear(); Du.resize(resultOrder + 1);
	for (auto iDu = Du.begin(); iDu != Du.end(); ++iDu)
		iDu->resize(partialsCount);

	typedef fadbad::F<ScalarType> FADScalar; 		// Auto Diff Forward scalar type (i.e. w.r.t. initial coefficients)
	typedef fadbad::T<FADScalar> TFADScalar;		// Auto Diff Taylor type with option to compute derivatives wrt initials
	typedef typename VectorType::template rebind<TFADScalar>::other TFADVectorType;
	// this is Automatic Differentiation both w.r.t t but also w.r.t. coefficients.
	// we will use it to automatically eval the r.h.s. of the equation with all
	// partial derivatives when recursively compute Taylor coefficients at current t0
	// So TAD computes derivatives of the jet composition, but with FAD inside we can see how the
	// composition depends on the initial variables!
	// again, we need the constructor with false here to combat a
	// 'bug' in FADBAD/CAPD cooperation. The bug was that, using another constructor,
	// all AD variables were assigned to the same space and shared internals
	// form AD perspective this was ok, but this was not what we meant. We wanted to have n fresh
	// variables inside the vector, not dependent on each other, not a vector (x, x, x, x), that
	// depends on x. I leave this comment for future if we encounter similar problem.
	TFADVectorType args(dimension(), false);
	TFADScalar targ(t0); targ[1] = 1; // we have targ(t) = t, targ'(t) = 1, targ^{(k > 1)} = 0;

	// we will successively fill-in TFAD args below as the recursion will go on
	// using this iterator to u variable.
	// NOTE: important fact!
	//		this iterator only goes forward, never it should get reset!
	// 		because the older data was inserted into 'args', and we only need to fill new data,
	// This is why we expect the u to be organized in this way!
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
	// we start from 0to fill the data corresponding to the value u[0] = x(t0) (see other version of the computeDDECoefficients)
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
		TFADVectorType v(imageDimension());
		m_map(targ, args, v);

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
		// we fill from idelay = 1 because we have filled the respective value from computation above.
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


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>						// template spec
typename BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::TimePointType	// return type
BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::getMaxDelay(){
	// TODO: (NOT URGENT, RETHINK) cache max delay? Now the problem is that TimePoint has a reference to grid, and in case of 0 there is incompatibility... maybe I need to refactor m_grid in Time point to pointer and have an option to change it?
	if (m_delays.size() == 0)
		return TimePointType();
	TimePointType max_tau = *(this->begin());
	for (auto tau = this->begin() + 1; tau < this->end(); ++tau)
		if (tau > max_tau) max_tau = tau;
	return max_tau;
}


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>						// template spec
typename BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::CurveType		// return type
BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::makeCompatibleSegment(VectorType v){
	// if (v == {}) v = VectorType();
	// TODO: (SOMEHOW URGENT): implement!
	throw std::logic_error("BasicDiscreteDelaysFunctionalMap::makeCompatibleSegment(): Not Implemented Yet");
//	auto max_tau = getMaxDelay();
//	CurveType segment(-max_tau, -max_tau.gird()(0));
//	return segment;
}


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>								// template spec
void BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::checkDimension(size_type d) const {
	if (d != this->imageDimension()){
		std::ostringstream info;
		info << "DiscreteDelaysFunctionalMap::checkDimension(size_type): incompatible dimension: ";
		info << "is " << d << ", should be " << this->imageDimension() << ".";
		throw std::logic_error(info.str());
	}
}


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>								// template spec
void BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::checkMapSignature() const {
	if (m_map.dimension() != this->dimension()){
		std::ostringstream info;
		info << "DiscreteDelaysFunctionalMap: An internal map has incompatible dimensions: ";
		info << "is " << m_map.dimension() << ", should be " << this->dimension() << ".";
		info << "Number of delays is " << this->m_delays.size() << ", output dimension is" << this->imageDimension() << ".";
		throw std::logic_error(info.str());
	}
}


} // namespace ddes
} // namespace capd


#endif /* _CAPD_DDES_BASICDISCRETEDELAYSFUNCTIONALMAP_HPP_ */
