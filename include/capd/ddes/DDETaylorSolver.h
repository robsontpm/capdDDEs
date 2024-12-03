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

#ifndef _CAPD_DDES_TAYLORSOLVER_H_
#define _CAPD_DDES_TAYLORSOLVER_H_

#include <stdexcept>
#include <string>
#include <capd/ddes/DDECommon.h>
#include <capd/ddes/FunctionalMap.h>

namespace capd{
namespace ddes{

template<typename FunctionalMapSpec>
class DDETaylorSolver {
public:
	typedef FunctionalMapSpec FunctionalMapType;
	typedef typename FunctionalMapType::CurveType SolutionCurveType;
	typedef typename FunctionalMapType::CurveType CurveType;
	typedef typename FunctionalMapType::JetType JetType;
	typedef typename SolutionCurveType::SetType SetType;  	// here, we deliberately use SolutionCurveType::SetType to prevent using with nonrigorous code

	typedef typename FunctionalMapType::DataType DataType; 	// it is in fact SetType
	typedef typename FunctionalMapType::VectorType VectorType;
	typedef typename FunctionalMapType::MatrixType MatrixType;
	typedef typename FunctionalMapType::TimePointType TimePointType;
	typedef typename FunctionalMapType::RealType RealType;
	typedef typename FunctionalMapType::size_type size_type;
	typedef typename FunctionalMapType::VariableStorageType VariableStorageType;
	typedef typename FunctionalMapType::JacobianStorageType JacobianStorageType;
	typedef typename FunctionalMapType::ValueStorageType ValueStorageType;

	DDETaylorSolver(DDETaylorSolver const& solver): m_map(solver.m_map), m_maxOrder(solver.m_maxOrder) {}
	DDETaylorSolver(FunctionalMapType const& map, size_type maxOrder): m_map(map), m_maxOrder(maxOrder) {}

	/**
	 * this function produces estimates that could be used to produce solution
	 * after time h = distance between grid points in the input curve.
	 * The solution will be valid on [t0, t0+h], its enclosure will be as interval set
	 * in out_CoeffsEnclosure up to order one higher than the out_Phi_coeffs_t0
	 * The var out_Phi_coeffs_t0 will hold the Jet coefficients at point t0 computed
	 * for the middle point of the curve, out_JacPhi_coeffs_t0 will be the Jacobian of
	 * coeffs w.r.t. variables in out_u,  over the whole set of the curve.
	 * The value out_Rem_coeffs_t0 is the remainder, but it will be 0 always in this
	 * implementation, as Jets are computed exactly for grid points.
	 * Then out_Phi_x[0] is the value at t = t0+h computed by the Taylor method,
	 * i.e. evaluation of out_Phi_coeffs_t0(t0+h). out_JacPhi_x[0] and out_Rem_x[0]
	 * are like in Taylor method for ODEs.
	 *
	 * Curve.move() procedure then should be able to extract information needed
	 * to construct a succesive CurvePiece of the solution over [t0, t0+h) and
	 * to set current value by use of out_{Phi|JacPhi|Rem}_x as it is customary in ODEs
	 * code.
	 *
	 * TODO: (FAR FUTURE) rethink that all data in and out encapsulated in some other structure?
	 */
	void encloseSolution(
			TimePointType const& 		in_t0,
			SolutionCurveType const& 	in_curve,
			TimePointType&				out_th,
			RealType&					out_HH,
			VariableStorageType&		out_u,
			ValueStorageType&			out_u_encl,
			ValueStorageType& 			out_Phi_coeffs_t0,
			JacobianStorageType& 		out_JacPhi_coeffs_t0,
			ValueStorageType& 			out_Rem_coeffs_t0,
			ValueStorageType& 			out_Y, 				// this will be the enclosure of the jet of the solution over time [t0, t0+h]
			ValueStorageType& 			out_Phi_z,			// for possible future use, I assume that there might be more than only value of x(t+h) !
			JacobianStorageType& 		out_JacPhi_z,		// I need storage for as many as there will be entries in out_Phi_z
			ValueStorageType& 			out_Rem_z			// I need storage for as many as there will be entries in out_Phi_z
			)
	{
		size_type d = m_map.imageDimension();
		VectorType zero(d); MatrixType Zero(d, d);

		RealType h = in_curve.getStep();
		out_th = in_t0 + in_curve.getStep();
		out_HH = ecloseStep(h);

		size_type order;
		// collecting computation data will make things faster and allow to collect u_encl - needed in rigorous computations.
		m_map.collectComputationData(in_t0, out_th, out_HH, in_curve, out_u, out_u_encl, order);
		// order is the final result order of the new representation at the current time
		// i.e. it will be one higher than the smallest orders at each delay.
		// then we take into account the possibility that order is too high for the computation to be feasible
		// order = order > m_maxOrder ? m_maxOrder : order;
		size_type coeffs_order = order;
		bool extraOrder = false;
		if (order > m_maxOrder){
			// extra rank computation as in old method
			extraOrder = true;
			order = m_maxOrder;
			coeffs_order = m_maxOrder + 1;
		}

		if (out_Phi_coeffs_t0.size() != coeffs_order + 1)
			out_Phi_coeffs_t0.resize(coeffs_order + 1);

		// make sure that JacPhi can store enough data
		out_JacPhi_coeffs_t0.resize(coeffs_order + 1);
		for (auto iJac = out_JacPhi_coeffs_t0.begin(); iJac != out_JacPhi_coeffs_t0.end(); ++iJac)
			iJac->resize(out_u.size(), Zero);

		// those will be used for computation of DDECeoficients
		// uMid - middle point of variables in u
		// uHull - interval hull of u
		// uEncl - enclosure of variables in u on whole intervals of length out_HH
		ValueStorageType uMid, uHull, uEncl;

		variableToValueMiddle(out_u, uMid);
		variableToValueHull(out_u, uHull);

		// compute Jacobian of the numerical method on the hull,
		// out_Phi_coeffs_t0 will be discarded later (computed only for midpoint)
		m_map.computeDDECoefficients(in_t0, uHull, out_Phi_coeffs_t0, out_JacPhi_coeffs_t0);

		// this is to test older codes.
		// this is version where I compute by interval method (not mean value form)
		// mean value form is below.
		VectorType extraOrderAtT0J1(d);
		if (extraOrder)
			extraOrderAtT0J1 = out_Phi_coeffs_t0[coeffs_order];

		// compute coefficients for the middle of the set, overwriting previous out_Phi_coeffs_t0
		m_map.computeDDECoefficients(in_t0, uMid, out_Phi_coeffs_t0);

		// now we have Taylor coefficients at t0, we now simply use standard
		// Taylor method, as in case of ODEs (with the except, that we must
		// compute roughEnclosure and enclosure of the solution over [t0, t0+h] differently

		out_Phi_z.clear(); out_Phi_z.resize(1, zero);	// there will be only one dependent value, i.e. x(t0+h)
		out_JacPhi_z.clear(); out_JacPhi_z.resize(1); 	// so, there will be only one Jacobian entry
		out_JacPhi_z[0].resize(out_u.size(), Zero); 	// but it will be dependent on  the same variables in u

		// we explicitely eval jet at t0 as a Taylor series to assure that the
		// jacobian computed along is valid.
		// TODO: (NOT URGENT) rewrite as iterators?
		// TODO: change it to: do {} while(k--) loop, instead of a hack for unsigned. Test.
		size_type k = order;
		while(true){
			out_Phi_z[0] = out_Phi_coeffs_t0[k] + h * out_Phi_z[0];
			for (size_type j = 0; j < out_u.size(); ++j){
				out_JacPhi_z[0][j] = out_JacPhi_coeffs_t0[k][j] + h * out_JacPhi_z[0][j];
			}
			if (k == 0) break; else --k; // this prevents from infinite loop in case size_type is unsigned
		}

		// there is no remainder for coeffs - the values does not depend
		// on the rough enclosure on the Xi part of the set. Therefore, we set to zero.
		out_Rem_coeffs_t0.clear(); out_Rem_coeffs_t0.resize(coeffs_order + 1, zero);
		out_Rem_z.clear(); out_Rem_z.resize(1, zero);
		// in out_CoeffsEnclosure we can compute up to order + 1 (incl.)
		// because we have enclosures on the coefficients up to order order
		// for all delayed points (Xi is the enclosure on the highest possible order)
		// Y is the name of the variable in the current paper/notes on extension of the solver.
		// please note, that rough enclosure is computed inside m_map.collectComputationData(...).
		out_Y.clear(); out_Y.resize(coeffs_order + 2, zero);

		// then w propagate it with the map
		m_map.computeDDECoefficients(in_t0, out_u_encl, out_Y);

		if (extraOrder){
			// TODO: (NOT URGENT, RETHINK) this is just to be sure the code is the same as old one, but probably not needed, since expansion of the representation is more natural...
			// implement extra expansion when order = maxOrder (as in old codes)

			// THIS IS BAD HACK JUST FOR TEST PURPOSES (move logic in solver to get better results)
			// HACK STARTS HERE
			MatrixType extra_C(d, in_curve.storageN0());
			// compute value of the extra coefficient used later in Xi
			for (size_type iu = 0; iu < out_u.size(); ++iu){
				MatrixType C_u = out_u[iu].get_C();  // TODO: (NOT URGENT) get_C copies matrix, now it is safe interface, but for speed rethink it should return reference / pointer?
				extra_C += out_JacPhi_coeffs_t0[order + 1][iu] * C_u;
			}
			MatrixType S = extra_C; S *= 0.;
			split(extra_C, S);
			// in this I do not care about extra_B - i will cast it to Vector anyway.
			VectorType extra_x(d), extra_r(d);
			extra_x  = out_Phi_coeffs_t0[coeffs_order] + out_Rem_coeffs_t0[coeffs_order];
			split(extra_x, extra_r);
			// split(extra_C, S); // TODO: (IMPORTANT) this was redundant! We lost content of S! TEST then remove this line.
			// TODO: also, rething removing this part as stated elsewhere!
			extra_r += S * in_curve.get_r0();
			for (size_type iu = 0; iu < out_u.size(); ++iu)
				extra_r += (out_JacPhi_coeffs_t0[coeffs_order][iu] * out_u[iu].get_B()) * out_u[iu].get_r();

			// HACK ENDS HERE

			// now modify new_Xi to reflect this. We have two versions.
			// first version of alternate out_Y[order + 1] computed with Interval arythm.
			VectorType extra_YJ1 = extraOrderAtT0J1 + (order + 2) * out_HH * out_Y[order + 2];
			// first version of alternate out_Y[order + 1] computed with Interval arythm.
			VectorType extra_YJ2 = extra_x + (extra_C * in_curve.get_r0()) + extra_r + (order + 2) * out_HH * out_Y[order + 2];

			try {
				out_Y[order + 1] = capd::vectalg::intersection(extra_YJ1, out_Y[order + 1]);
			} catch (...) {
				throw std::logic_error("CRITICAL: empty intersection of sets that THEORETICALLY ARE PROVED TO BE NON-EMPTY. Please contact authors (robert.szczelina@uj.edu.pl), as the source code could contain some error.");
			}

			// we need to clean the extra coefficients from data, so that set.move() is not affected
			out_Y.pop_back();
			out_Phi_coeffs_t0.pop_back();
			out_JacPhi_coeffs_t0.pop_back();
			out_Rem_coeffs_t0.pop_back();
		}

		// we extract the remainder of x(t) over [t0, th]:
		out_Rem_z[0] = power(h, order + 1) * out_Y[order + 1];
	}

	/** see docs for encloseSolution(TimePointType const& t0, ...). It does the same but for t0 = in_curve.t0() */
	void encloseSolution(
			SolutionCurveType const& 	in_curve,
			TimePointType&				out_th,
			RealType&					out_HH,
			VariableStorageType&		out_u,
			ValueStorageType&			out_u_encl,
			ValueStorageType& 			out_Phi_coeffs_t0,
			JacobianStorageType& 		out_JacPhi_coeffs_t0,
			ValueStorageType& 			out_Rem_coeffs_t0,
			ValueStorageType& 			out_CoeffsEnclosure,
			ValueStorageType& 			out_Phi_z,
			JacobianStorageType& 		out_JacPhi_z,
			ValueStorageType& 			out_Rem_z)
	{
		encloseSolution(
				in_curve.rightDomain(), in_curve,
				out_th, out_HH,		// th - time at which Phi_z ic computed. out_HH - difference between t0 and th
				out_u, out_u_encl, 	// variables used in computation, for further use
				out_Phi_coeffs_t0, out_JacPhi_coeffs_t0, out_Rem_coeffs_t0, out_CoeffsEnclosure, // Phi_j
				out_Phi_z, out_JacPhi_z, out_Rem_z  // Phi_z
		);
	}

	const FunctionalMapType& getMap() const { return m_map; }
	FunctionalMapType& getMap() { return m_map; }

private:
	FunctionalMapType m_map;
	size_type m_maxOrder;

	/** helper function. TODO: (NOT URGENT, FUTURE) move it to the definition of VariableSTorageType, this will be more work, as I use std::vector for those */
	void variableToValueMiddle(VariableStorageType const & u, ValueStorageType& uMid){
		uMid.clear();
		for (auto iu = u.begin(); iu != u.end(); ++iu)
			uMid.push_back(iu->midPoint());
	}
	void variableToValueHull(VariableStorageType const & u, ValueStorageType& uHull){
		uHull.clear();
		for (auto iu = u.begin(); iu != u.end(); ++iu)
			uHull.push_back(iu->hull());

	}

	// TODO: (NOT URGENT, FUTURE) define class Variable/Value/Jacobian Storage Binding and there put helpful conversions
	// TODO: (NOT URGENT, FUTURE) Idea is that the Variable type represents some abstract representation of sets
	// TODO: (NOT URGENT, FUTURE) while Value and Jacobian represent concrete representations of Interval Arithmetics.
	// TODO: (NOT URGENT, FUTURE) the binding is to glue them together for simpler usage.
};



} // namespace ddes
} // namespace capd


#endif /* _CAPD_DDES_TAYLORSOLVER_H_ */
