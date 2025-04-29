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
				ValueStorageType& out_v, 							// output
				VariableStorageType& out_u, 						// output
				JacobianStorageType& out_dvdu, 						// output
				size_type& out_admissible_order) const {			// output
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

	// empty the container. We exploit the fact that in case it is empty
	// then the FunctionalMap assume it is identity, i.e. v == u
	// and can optimize some computations.
	out_dvdu.clear();

	// and here we say that v = u (but we forget the variable DataClass and convert to pure Vector
	BaseClass::convert(out_u, out_v);
}


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>		// template spec
void BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::computeDDECoefficients(
			const RealType& t0, const ValueStorageType& v, 							// input
			ValueStorageType& coeffs) const {										// output
	// we assume u[0] is always present, and all other u[i] has the same dimension
	// we only check once. If other have different, then the CAPD will throw a generic exception
	checkDimension(v[0].dimension());
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
	auto current_v = v.begin();
	auto iarg = args.begin();
	// here we start from 0, as we are extracting the value at t0
	// later in the code you might notice we start from 1.
	for (size_type idelay = 0; idelay <= delaysCount(); ++idelay, ++current_v)
		for (auto iu = current_v->begin(); iu != current_v->end(); ++iu, ++iarg)
			(*iarg)[0] = *iu;

	// now we are ready to go for the recurrent relation on coefficients
	// to produce first derivative. Later, after each loop, we will have
	// all derivatives up to order n and we are ready to compute n+1'st.
	coeffs[0] = v[0]; // by definition, see paper
	for (size_type k = 1 ; k <= resultOrder; ++k){
		// eval Autodiff routine so that we would be able to
		// extract v[k-1] = F^{[k-1]}(...), as in the paper.
		// The constructor with false here is important to not initialize v data
		// with anything, each position in v will be just a fresh AD variable.
		// in (2024) there was an error introduced by use of other constructor, when CAPD changed.
		TADVectorType fv(imageDimension(), false);
		m_map(targ, args, fv);
		// propagate AD jets and output.
		// this will be F_k from the latest paper on DDE integration (FoCM 2024).
		VectorType Fk(imageDimension());
		auto iarg = args.begin();
		for (size_type i = 0; i < imageDimension(); ++i, ++iarg){
			fv[i].eval(k-1);
			Fk[i] = fv[i][k-1];
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
		for (size_type idelay = 1; idelay <= delaysCount(); ++idelay, ++current_v){
			for (auto iu = current_v->begin(); iu != current_v->end(); ++iu, ++iarg)
				(*iarg)[k] = *iu;
		}
	} /// for k loop
} /// computeDDECoefficients


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>	// template spec
void BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::computeDDECoefficients(
			const RealType& t0, const ValueStorageType& v, 						// input
			ValueStorageType& coeffs, JacobianStorageType& Dv) const {			// output
	 // TODO: (NOT URGENT, FUTURE, RETHINK) there is some redundancy in the code (with the other version),
	 // TODO: (NOT URGENT, FUTURE, RETHINK) but I don't believe we can make it more DRY without sacrificing speed/RAM usage...

	// we assume u[0] is always present, and all other u[i] has the same dimension
	// we only check once. If other have different, then the CAPD will throw a generic exception
	checkDimension(v[0].dimension());
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
	size_type partialsCount = v.size();
	// in partialsDimension we keep how many actual coefficients from all of the Jets used to
	// compute rhs of the equation. As we are computing jets up to coeffs.order() (incl.).
	// then we are going to use only coefficients of all jets up to coeffs.order() - 1 (incl).
	// Thus, the formula. Please note in all those computations we have 0,..,resultOrder-1
	// Taylor coefficients in any jet used. Therefore the number of Taylor coeffs at each
	// jet used is resultOrder. This is trivial, but I want this to be reiterated for readability.
	size_type partialsDimension = d * v.size();
	// we make room to store Jacobian matrices. Each matrix will be of shape d x d.
	// Du is shorthand for $D_u coeffs = \frac{\partial coeffs_j}{\partial u_i}$
	// (the last notation means we are iterating over all possibilities and we put it into a matrix)
	// in Du[i][j] we will have $\frac{\partial coeffs[i]}{\partial u[j]}$ (a d x d Matrix)
	Dv.clear(); Dv.resize(resultOrder + 1);
	for (auto iDv = Dv.begin(); iDv != Dv.end(); ++iDv)
		iDv->resize(partialsCount);

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
	auto current_v = v.begin();
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
	for (size_type idelay = 0; idelay <= numDelays; ++idelay, ++current_v, ++partialu){
		for (auto iv = current_v->begin(); iv != current_v->end(); ++iarg, ++iv, ++partialAD){
			// (*iarg)[0] is [0]-th order Taylor coefficient
			(*iarg)[0] = *iv;
			// and we mark it as dependent only on self among all 0,..,partialsDimension-1 variables
			// (i.e. \frac{\partial{(*iarg)[0]}}{\partial{(*iarg)[0]}} = 1, and 0 otherwise)
			(*iarg)[0].diff(partialAD, partialsDimension);
		}
	}

	// we compute coeffs[0] and derivatives...
	coeffs[0] = v[0];		// coeffs[0] = x(t0) by definition of (see paper)
	Dv[0][0] = Id; 			// \frac{Dcoeffs_{[0]}}{Dx(t0)} = Id (as a matrix M(d, d)) (since coeffs[0] = x(t0))
	// Next... \frac{Dcoeffs_{[0]}}{{Djet_m}_{[mk]}} = 0 (as a matrix M(d, d)) (since coeffs[0] = x(t0))
	for (size_type par = 1; par < partialsCount; ++par)
		Dv[0][par] = Zero;

	// then we compute higher order coeffs, with recurrent Taylor formula (see paper)
	for (size_type k = 1; k <= resultOrder; ++k){
		// eval AutoDiff routine so we would be able to extract v[.][k-1] = F^{[k-1]}(...)
		// together with partial derivatives with respect to coefficients in the past.
		TFADVectorType fv(imageDimension());
		m_map(targ, args, fv);

		// propagate AD jets and output.
		VectorType Fk(imageDimension());
		auto iarg = args.begin();
		for (size_type i = 0; i < imageDimension(); ++i, ++iarg){
			fv[i].eval(k-1);
			Fk[i] = fv[i][k-1].val(); 					// here we need to .value(), to convert from FAD to VectorType
			(*iarg)[k] = fv[i][k-1] / ScalarType(k);	// here we need to propagate things with all partial derivatives! So we copy from v (appropriate AutoDiff Type)
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
					parDu[i][j] = fv[i][k-1].d(par * d + j) / ScalarType(k);// among all partialsCount variables
			Dv[k][par] = parDu;
		}
		// mark that this order coeff does not depend on higher order u's
		for (size_type par = partialu; par < partialsCount; ++par)
			Dv[k][par] = Zero;


		// this if is to ensure that we do not copy, if there is no more values left in u.
		if (k == resultOrder) break;

		// We have now propagated dependence of the solution on the value at
		// present time t0 to Du. Now, finally, we propagate higher order
		// coefficients in the past to args for the next level of recurrence.
		// we fill from idelay = 1 because we have filled the respective value from computation above.
		for (size_type idelay = 1; idelay <= numDelays; ++idelay, ++current_v, ++partialu){
			for (auto iu = current_v->begin(); iu != current_v->end(); ++iarg, ++iu, ++partialAD){
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









//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// TODO: BELOW I WORK ON SD-DDES - I need to move this out, after I finish the experiments
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>							// template spec
typename BasicStateDependentDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::VectorType	// return type
BasicStateDependentDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::operator()(			// method name
			const TimePointType& t,																		// args
			const CurveType& x) const {																	// args...
	// checkCurveDimension(x, "operator()"); // TODO: (NOT URGENT) use from Interface later
	VectorType args(this->dimension());
	VectorType x0 = x.eval(t);
	auto iargs = args.begin();
	for (auto ix0 = x0.begin(); ix0 != x0.end(); ++ix0, ++iargs) *iargs = *ix0;
	for (auto idelay = m_delays.begin(); idelay != m_delays.end(); idelay++){
		VectorType xtau = x.eval((*idelay)->operator()(t, x));
		for (auto ixtau = xtau.begin(); ixtau != xtau.end(); ++ixtau, ++iargs) *iargs = *ixtau;
	}
	VectorType result(imageDimension());
	m_map(RealType(t), args, result);
	return result;
}


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>									// template spec
void BasicStateDependentDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::collectComputationData( // method name
				const TimePointType& t0, const TimePointType& th,	//input
				const RealType& dt, const CurveType& x, 	 		// input
				ValueStorageType& out_v,							// output
				VariableStorageType& out_u,							// output
				JacobianStorageType& out_dvdu,						// output
				size_type& out_admissible_order) const {			// output

	// TODO: remove after testing everything below
//	// TODO: (NOT URGENT) DRY with the other  collectComputationData() method.
//	// TODO: (NOT URGENT) I'm not sure it will be possible, but think about it!
//	out_u.clear();
//	// we collect jets and compute admissible order
//	// TODO: (NOT URGENT) rewrite using pointers to Jets (will be faster)!
//	out_admissible_order = 40; // TODO: (NOT URGENT) Magic Constant....
//	std::vector<JetType> jets;
//	for (auto tau = begin(); tau != end(); ++tau){
//		jets.push_back(x.jet((*tau)->inf(t0, x)));   // TODO: here I just get the lower bound on the delay. The next quest is to make it truly evaluation of the delay and the jet at the delay
//		//jets.push_back((*tau)->jet(t0, x)); // return the jet at the retarded argument
//		size_type k = jets.back().order();
//		if (k + 1 < out_admissible_order) out_admissible_order = k + 1;
//	}
//	// Then we fill in at u[0] the value at current time (always present in u, by assumption)
//	out_u.push_back(x.getValueAtCurrent());
//	// and then all jets, order by order (in this way it is easier to fill AutoDiff later)
//	for (size_type k = 0; k < out_admissible_order; ++k)
//		for (auto jit = jets.begin(); jit != jets.end(); ++jit)
//			out_u.push_back((*jit)[k]);
//
//	// and here we say that v = u (but we forget the variable DataClass and convert to pure Vector
//	BaseClass::convert(out_u, out_v);
//
//	// for test, I put a true Id here
//	// TODO: later, I need to get this from the delays...
//	auto d = imageDimension();
//	auto M = out_u.size(); // == out_v.size();
//	MatrixType Zero(d, d), Id = MatrixType::Identity(d);
//	out_dvdu.resize(M);
//	for (size_type i = 0; i < M; ++i){
//		out_dvdu[i].resize(M, Zero);
//		out_dvdu[i][i] = Id;
//	}

	// NEW VERSION FOR STATE DEPENDENT DELAYS

	// First, we extract data at t0 - which is always present in the data
	// (I always assume that x'(t) = f(x(t), ...) - that way I always have at least one variable to do Lohner stuff later)
	// TODO: it should be eval at t0! (RETHINK!) (KIND OF IMPORTANT!)
	auto x_at_t0 = x.getValueAtCurrent();
	out_v.push_back(x_at_t0);
	out_u.push_back(x_at_t0);
	auto d = out_v.back().dimension();
	MatrixType Id = MatrixType::Identity(d), Zero(d, d);
	out_dvdu.push_back({Id});
	auto& dvdu_0 = out_dvdu.back();

	// TODO: (RETHINK) what I am implementing here is in fact some form automatic differentiation of v(t, u) w.r.t. to both t and u
	// TODO: (RETHINK) but with variables u which are d-dimensional by default
	// TODO: (RETHINK) also, if one delay returns a Variable it might be the same Variable
	// TODO: (RETHINK) as in the other delay. So I store this twice. It's ok, because I will use those same representations of variables for both, but the redundancy might be large... (but for usual r.h.s. of DDEs this should happen very very rarely)
	// TODO: (RETHINK) maybe rethink to store v, u, dvdu as a common collection?
	// TODO: (RETHINK, FUTURE) maybe rethink to store u as pointers to datastructures inside the solution curve, so that I can identify them, and then remove redundancy in variables?
	out_admissible_order = 40; // TODO: (NOT URGENT) Magic Constant....
//	std::vector<ValueStorageType> values;			// represents v(u)
//	std::vector<VariableStorageType> variables;		// represents u
//	std::vector<JacobianStorageType> derivatives;	// represents \frac{\partial v}{\partial u}
	std::vector<std::tuple<ValueStorageType, VariableStorageType, JacobianStorageType>> data;
	for (auto tau = begin(); tau != end(); ++tau){
		ValueStorageType tau_v;
		VariableStorageType tau_u;
		JacobianStorageType tau_dvdu;
		size_type k;
		(*tau)->collectComputationData(t0, th, dt, x, tau_v, tau_u, tau_dvdu, k);
		data.push_back(std::make_tuple(tau_v, tau_u, tau_dvdu));
		// we put all variables to the common structure to be used later
		// see comments later why we put all the variables, but we will
		// use admissible order to filter some values.
		for (auto u_i = tau_u.begin(); u_i != tau_u.end(); ++u_i){
			MatrixType locZero(d, u_i->dimension());
			out_u.push_back(*u_i);
			dvdu_0.push_back(locZero);
		}

		// we assume that v represents values of some variable v, then v',
		// then v'', etc. (actually v^{[k]} = v^{[k]}/k!, i.e. Taylor coefficients)
		// therefore the order is 0, ..., v.size() - 1, and v.size() + 1 is the
		// admissible order, because of the smoothing of solutions
		if (k < out_admissible_order) out_admissible_order = k;
	}

	// reorganize data to match the algorithm
	// TODO: maybe reorganize algorithm to match the data?
	// NOTE: here we cut out the variables that we cannot use in our computations
	//       because of the admissible order.

	// TODO: (FUTURE) I definitely need to put variables, values and derivatives together in a more robust structure
	// TODO: (FUTURE) because now I have to put a lot of zeros. I think, I now understands why Daniel did the Node-based AutoDiff
	// TODO: (FUTURE) I need to ask him if this is indeed the reason. Then, I think I might be able to rewrite this using Node algebra
	// TODO: (FUTURE) and maybe then I will be able to rewrite this code with Nodes.
	// TODO: (RETHINK, FUTURE): after looking at Node it seems that Daniel is indexing Nodes with integers too (the dag is just a vector of Nodes).
	// TODO: (RETHINK, FUTURE): on the other hand, maybe I can create DAGs for Delays and for Function itself, and then combine them? Rethink... Maybe meet with Daniel nevertheless.

	// and then all jets, order by order (in this way it is easier to fill AutoDiff later)
	// i - "index" of the delay, this just denotes that this is index of delay, i use iterators here
	// k - index of the Taylor coefficient used (the derivative w.r.t. t of the jet of v treated as a function v = v(t)
	// j - index of the variable within the given delay
	for (size_type k = 0; k < out_admissible_order; ++k){
		size_type base_u_i = 1; //because we have x_at_t0 in out_u[0]
		for (size_type i = 0; i < data.size(); ++i){
			auto& item = data[i];
			// unpack the data at i-th tau
			auto& v_i = std::get<0>(item);
			auto& u_i = std::get<1>(item);
			auto& dvdu_i = std::get<2>(item);
			auto& v_ik = v_i[k]; 		// k-th value at i-th delay
			auto& dvdu_ik = dvdu_i[k]; 	// the derivative of k-th value at i-th delay wrt all variables at i-th delay
			out_v.push_back(v_ik);
			out_dvdu.push_back({Zero}); // dependence of v_ik on x_at_t0 is theoretically 0
			auto& out_dv_ik__du = out_dvdu.back();
			auto u_i_count = u_i.size();
			auto dvdu_ik_j = dvdu_ik.begin();
			for (int j = 1; j < out_u.size(); ++j){ // we skip the first variable = x_at_t0, therefore we start at 1
				MatrixType D(v_ik.dimension(), out_u[j].dimension()); // by default Zero
				if (base_u_i <= j && j < base_u_i + u_i_count){
					// extract data dependence, instead of assuming Zero
					// we can increase linearly iterator to dvdu_ik, as we added u_i's in this order,
					// and in this loop we are processing them in this order
					D = *dvdu_ik_j;
					++dvdu_ik_j;
				}
				out_dv_ik__du.push_back(D);
			}
			base_u_i += u_i_count; // increase the base index where u_i's starts in the out_u collection
		}
	}
}


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>		// template spec
void BasicStateDependentDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::computeDDECoefficients(
			const RealType& t0, const ValueStorageType& v, 							// input
			ValueStorageType& coeffs) const {										// output
	// we assume u[0] is always present, and all other u[i] has the same dimension
	// we only check once. If other have different, then the CAPD will throw a generic exception
	checkDimension(v[0].dimension());
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
	auto current_v = v.begin();
	auto iarg = args.begin();
	// here we start from 0, as we are extracting the value at t0
	// later in the code you might notice we start from 1.
	for (size_type idelay = 0; idelay <= delaysCount(); ++idelay, ++current_v)
		for (auto iu = current_v->begin(); iu != current_v->end(); ++iu, ++iarg)
			(*iarg)[0] = *iu;

	// now we are ready to go for the recurrent relation on coefficients
	// to produce first derivative. Later, after each loop, we will have
	// all derivatives up to order n and we are ready to compute n+1'st.
	coeffs[0] = v[0]; // by definition, see paper
	for (size_type k = 1 ; k <= resultOrder; ++k){
		// eval Autodiff routine so that we would be able to
		// extract v[k-1] = F^{[k-1]}(...), as in the paper.
		// The constructor with false here is important to not initialize v data
		// with anything, each position in v will be just a fresh AD variable.
		// in (2024) there was an error introduced by use of other constructor, when CAPD changed.
		TADVectorType fv(imageDimension(), false);
		m_map(targ, args, fv);
		// propagate AD jets and output.
		// this will be F_k from the latest paper on DDE integration (FoCM 2024).
		VectorType Fk(imageDimension());
		auto iarg = args.begin();
		for (size_type i = 0; i < imageDimension(); ++i, ++iarg){
			fv[i].eval(k-1);
			Fk[i] = fv[i][k-1];
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
		for (size_type idelay = 1; idelay <= delaysCount(); ++idelay, ++current_v){
			for (auto iu = current_v->begin(); iu != current_v->end(); ++iu, ++iarg)
				(*iarg)[k] = *iu;
		}
	} /// for k loop
} /// computeDDECoefficients


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>	// template spec
void BasicStateDependentDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::computeDDECoefficients(
			const RealType& t0, const ValueStorageType& v, 						// input
			ValueStorageType& coeffs, JacobianStorageType& Dv) const {			// output
	 // TODO: (NOT URGENT, FUTURE, RETHINK) there is some redundancy in the code (with the other version),
	 // TODO: (NOT URGENT, FUTURE, RETHINK) but I don't believe we can make it more DRY without sacrificing speed/RAM usage...

	// we assume u[0] is always present, and all other u[i] has the same dimension
	// we only check once. If other have different, then the CAPD will throw a generic exception
	checkDimension(v[0].dimension());
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
	size_type partialsCount = v.size();
	// in partialsDimension we keep how many actual coefficients from all of the Jets used to
	// compute rhs of the equation. As we are computing jets up to coeffs.order() (incl.).
	// then we are going to use only coefficients of all jets up to coeffs.order() - 1 (incl).
	// Thus, the formula. Please note in all those computations we have 0,..,resultOrder-1
	// Taylor coefficients in any jet used. Therefore the number of Taylor coeffs at each
	// jet used is resultOrder. This is trivial, but I want this to be reiterated for readability.
	size_type partialsDimension = d * v.size();
	// we make room to store Jacobian matrices. Each matrix will be of shape d x d.
	// Du is shorthand for $D_u coeffs = \frac{\partial coeffs_j}{\partial u_i}$
	// (the last notation means we are iterating over all possibilities and we put it into a matrix)
	// in Du[i][j] we will have $\frac{\partial coeffs[i]}{\partial u[j]}$ (a d x d Matrix)
	Dv.clear(); Dv.resize(resultOrder + 1);
	for (auto iDv = Dv.begin(); iDv != Dv.end(); ++iDv)
		iDv->resize(partialsCount);

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
	auto current_v = v.begin();
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
	for (size_type idelay = 0; idelay <= numDelays; ++idelay, ++current_v, ++partialu){
		for (auto iv = current_v->begin(); iv != current_v->end(); ++iarg, ++iv, ++partialAD){
			// (*iarg)[0] is [0]-th order Taylor coefficient
			(*iarg)[0] = *iv;
			// and we mark it as dependent only on self among all 0,..,partialsDimension-1 variables
			// (i.e. \frac{\partial{(*iarg)[0]}}{\partial{(*iarg)[0]}} = 1, and 0 otherwise)
			(*iarg)[0].diff(partialAD, partialsDimension);
		}
	}

	// we compute coeffs[0] and derivatives...
	coeffs[0] = v[0];		// coeffs[0] = x(t0) by definition of (see paper)
	Dv[0][0] = Id; 			// \frac{Dcoeffs_{[0]}}{Dx(t0)} = Id (as a matrix M(d, d)) (since coeffs[0] = x(t0))
	// Next... \frac{Dcoeffs_{[0]}}{{Djet_m}_{[mk]}} = 0 (as a matrix M(d, d)) (since coeffs[0] = x(t0))
	for (size_type par = 1; par < partialsCount; ++par)
		Dv[0][par] = Zero;

	// then we compute higher order coeffs, with recurrent Taylor formula (see paper)
	for (size_type k = 1; k <= resultOrder; ++k){
		// eval AutoDiff routine so we would be able to extract v[.][k-1] = F^{[k-1]}(...)
		// together with partial derivatives with respect to coefficients in the past.
		TFADVectorType fv(imageDimension());
		m_map(targ, args, fv);

		// propagate AD jets and output.
		VectorType Fk(imageDimension());
		auto iarg = args.begin();
		for (size_type i = 0; i < imageDimension(); ++i, ++iarg){
			fv[i].eval(k-1);
			Fk[i] = fv[i][k-1].val(); 					// here we need to .value(), to convert from FAD to VectorType
			(*iarg)[k] = fv[i][k-1] / ScalarType(k);	// here we need to propagate things with all partial derivatives! So we copy from v (appropriate AutoDiff Type)
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
					parDu[i][j] = fv[i][k-1].d(par * d + j) / ScalarType(k);// among all partialsCount variables
			Dv[k][par] = parDu;
		}
		// mark that this order coeff does not depend on higher order u's
		for (size_type par = partialu; par < partialsCount; ++par)
			Dv[k][par] = Zero;


		// this if is to ensure that we do not copy, if there is no more values left in u.
		if (k == resultOrder) break;

		// We have now propagated dependence of the solution on the value at
		// present time t0 to Du. Now, finally, we propagate higher order
		// coefficients in the past to args for the next level of recurrence.
		// we fill from idelay = 1 because we have filled the respective value from computation above.
		for (size_type idelay = 1; idelay <= numDelays; ++idelay, ++current_v, ++partialu){
			for (auto iu = current_v->begin(); iu != current_v->end(); ++iarg, ++iu, ++partialAD){
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


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>								// template spec
typename BasicStateDependentDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::TimePointType	// return type
BasicStateDependentDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::getMaxDelay(){
	// TODO: (NOT URGENT, RETHINK) cache max delay? Now the problem is that TimePoint has a reference to grid, and in case of 0 there is incompatibility... maybe I need to refactor m_grid in Time point to pointer and have an option to change it?
	if (m_delays.size() == 0)
		return TimePointType();
	TimePointType max_tau = *(this->begin())->inf();
	for (auto tau = this->begin() + 1; tau < this->end(); ++tau){
		auto this_max = (*tau)->inf();
		if (this_max > max_tau) max_tau = this_max;
	}
	return max_tau;
}


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>								// template spec
typename BasicStateDependentDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::CurveType		// return type
BasicStateDependentDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::makeCompatibleSegment(VectorType v){
	// if (v == {}) v = VectorType();
	// TODO: (SOMEHOW URGENT): implement!
	throw std::logic_error("BasicDiscreteDelaysFunctionalMap::makeCompatibleSegment(): Not Implemented Yet");
//	auto max_tau = getMaxDelay();
//	CurveType segment(-max_tau, -max_tau.gird()(0));
//	return segment;
}


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>								// template spec
void BasicStateDependentDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::checkDimension(size_type d) const {
	if (d != this->imageDimension()){
		std::ostringstream info;
		info << "DiscreteDelaysFunctionalMap::checkDimension(size_type): incompatible dimension: ";
		info << "is " << d << ", should be " << this->imageDimension() << ".";
		throw std::logic_error(info.str());
	}
}


template<typename FinDimMapSpec, typename SolutionCurveSpec, typename JetSpec>								// template spec
void BasicStateDependentDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec>::checkMapSignature() const {
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
