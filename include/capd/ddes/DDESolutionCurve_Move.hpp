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

#ifndef _CAPD_DDES_DDESOLUTIONCURVE_MOVE_HPP_
#define _CAPD_DDES_DDESOLUTIONCURVE_MOVE_HPP_

#include <capd/ddes/DDESolutionCurve.h>
#include <capd/ddes/DDECommon.hpp>

////////////////////////////////////////////////////////////////////////
// Below are implementations of more involved functions               //
// (to be separated as usual from the declarations in the .h file     //
// For more complete documentation/explanation please consult .h file //
////////////////////////////////////////////////////////////////////////

namespace capd{
namespace ddes{

template<typename SetSpec>
template<typename DynSysSpec>
DDESolutionCurve<SetSpec>& DDESolutionCurve<SetSpec>::move(DynSysSpec& solver){
	typedef typename DynSysSpec::VariableStorageType Variables;
	typedef typename DynSysSpec::JacobianStorageType Jacobians;
	typedef typename DynSysSpec::ValueStorageType Values;
	Variables u; Values u_encl;
	Jacobians D_uPhi_j0, D_uPhi_z;
	Values Phi_z, Rem_z;
	Values Phi_j0, Rem_j0, Y;
	TimePointType t_h = m_grid.point(0);
	RealType HH;
	// we call solver to enclose all the data needed in a common data structures
	// for ODEs it was a lot easier, as it was exactly known, what are the
	// values used in the computations were, it was just the value of the solution
	// at time t0. Now, we do not know what values the rhs of the equation uses
	// because it might have a lot of delays, or integrals, or whatever
	// we trust solver to extract all the relevant data for us.
	// but we need to pass him a mutable storages that can hold a lot
	// of Vectors / Matrices. The meaning of u, D_uPhi_j, etc. is
	// described in the docs of solver and in the new notes / paper.
	solver.encloseSolution(
			*this,							// this is the only input. rest of variables is output.
			t_h, HH,						// t_h is the new TimePoint of new jet j0, HH is [0, h]
			u, u_encl,						// u is value at current time of variables used in computation, u_encl is over current time + HH
			Phi_j0, D_uPhi_j0, Rem_j0, Y,  	// corresponds to the expansion at current time
			Phi_z, D_uPhi_z, Rem_z			// corresponds to the new value at new time (after step), corresponds to result set in ODE setting
	);
	m_lastEnclosure = u_encl[0];
	// now, implement Lohners algorithm with the data provided...

	// this is the order of the resulting Jet at t0
	size_type order = Phi_j0.size() - 1;
	size_type uCount = u.size();
	size_type d = dimension();
	size_type N0 = storageN0();
	// for clarity of writing
	VectorType const& r0 = *m_r0;

	// Rem_j0 == 0 here (no remainder from numerical method, as the j0
	// is computed from rhs of solution with AD),
	// but we keep track of it in computations nevertheless, if there will be
	// in future some other algorithm that returns some overestimates
	MatrixType Zero(d, N0);
	MatrixType Id(d, d); Id.setToIdentity();
	// new_ will indicate part of the new Jet to be added to solution
	// new_C_j0 is a collection of C matrices for sets defining jet at current time
	// (therefore we have order + 1 of them).
	// _z is the value at current time (i.e. 0-th coefficient of a new jet,
	// it corresponds to what is the result set in ODE).
	std::vector<MatrixType> new_C_j0(order + 1, Zero);
	MatrixType new_C_z = Zero;
	MatrixType S = Zero;
	// WARNING: in the loop below we assume that r0 is common for all
	// WARNING: sets in u! Otherwise it would be incorrect!
	for (size_type iu = 0; iu < uCount; ++iu){
		MatrixType C_u = u[iu].get_C();  // TODO: (FUTURE) get_C copies matrix, now it is safe interface, but for speed rethink it should return reference / pointer?
		new_C_z += D_uPhi_z[0][iu] * C_u;
		for (size_type k = 0; k <= order; ++k)
			new_C_j0[k] += D_uPhi_j0[k][iu] * C_u;
	}

	// this will be commonly used by all QR computations, so we store it in tmp variable
	VectorType u0r = u[0].get_r();

	// center the value at the current time
	VectorType new_z = Phi_z[0] + Rem_z[0];
	VectorType new_r(d); split(new_z, new_r);
	// we do the QR decomposition with B chosen as the matrix of 0-th variable in u (current time value)
	// this way we reproduce the algorithm from ODEs (in case there are no delays)
	MatrixType new_B_z = (D_uPhi_z[0][0] * u[0].get_B());
	MatrixType new_invB_z(d, d);
	this->Policy::computeBinvB(new_B_z, new_invB_z, u0r);
	new_r = new_invB_z * new_r;
	for (size_type iu = 0; iu < uCount; ++iu)
		new_r += (new_invB_z * D_uPhi_z[0][iu] * u[iu].get_B()) * u[iu].get_r();
	// next, incorporate the remains of C * r0
	// new_C_z is point matrix, while we can put
	// ((Jac * old_C) - m(Jac * old_C)) * r0 = SC * r0 into remainder!
	MatrixType SC = Zero;
	split(new_C_z, SC);
	new_r += (new_invB_z * SC) * r0;

	// then we need to redo this for j0 (new Jet at current time)
	std::vector<MatrixType> new_B_j0(order + 1, Id);
	std::vector<MatrixType> new_invB_j0(order + 1, Id);
	// here I cannot change so simply the order of loops, so I leave it outside for Phi_z
	// it should not be so important, as C matrix in a set might be large
	// (big number of columns), whereas B matrix is simply M(d, d) (small).
	for (size_type k = 0; k <= order; ++k){
		// following two lines centers the set, important!
		// Phi_j0 will be point set, and Rem_j0 an interval set centered at 0
		Phi_j0[k] += Rem_j0[k];
		split(Phi_j0[k], Rem_j0[k]);

		new_B_j0[k] = (D_uPhi_j0[k][0] * u[0].get_B());
		this->Policy::computeBinvB(new_B_j0[k], new_invB_j0[k], u0r);

		Rem_j0[k] = new_invB_j0[k] * Rem_j0[k];
		for (size_type iu = 0; iu < uCount; ++iu)
			Rem_j0[k] += (new_invB_j0[k] * D_uPhi_j0[k][iu] * u[iu].get_B()) * u[iu].get_r();

		// next, incorporate the remains of C * r0
		// new_C[k] is point matrix, while we can put
		// ((Jac * old_C) - m(Jac * old_C)) * r0 = S * r0 into remainder!
		split(new_C_j0[k], S);
		Rem_j0[k] += (new_invB_j0[k] * S) * r0;
	}

	// now compute Xi (i.e. copy value form Jet enclosure at given order)
	VectorType* new_Xi = new VectorType(Y[order + 1]);

	// now, produce the new representation
	CurvePieceType* j0 = new CurvePieceType(this->t0(), Phi_j0, new_C_j0, m_r0, false, new_B_j0, new_invB_j0, Rem_j0, new_Xi, true);
	addPiece(j0, true);
	setValueAtCurrent(SetType(new_z, new_C_z, m_r0, new_B_z, new_invB_z, new_r));

	// this should never happen:
	if (!(t_h == this->t0())) throw std::logic_error("DDESolutionCurve::move(): method has incompatible step size.");

//	/// TODO: THIS IS JUST FOR TEST
//	{
//		Values encData;
//		CurvePieceType* jit = j0;
//		auto dt = HH;
//		VectorType enc = jit->evalAtDelta(dt);
//		encData.push_back(enc);
//		for (size_type k = 1; k < Y.size(); ++k)
//			encData.push_back(jit->evalCoeffAtDelta(k, dt));
//
//		// Y holds the computed value with roughEnclosure
//		std::cout << "INTERSECTION TEST: dimensions Y: " << Y.size() << " encData: " << encData.size() << std::endl;
//		for (int k = 0; k < Y.size(); ++k){
//			std::cout << "k: " << k << " dh = " <<  dt << std::endl;
//			for (int j = 0; j < Y[k].dimension(); ++j){
//				std::cout << Y[k][j] << " vs " <<  encData[k][j] << " (diam: " << (Y[k][j].rightBound() - Y[k][j].leftBound()) << " vs " <<  encData[k][j].rightBound() - encData[k][j].leftBound() << ")"  << std::endl;
//			}
//			std::cout << std::endl;
//		}
//		std::cout << "INTERSECTION TEST END" << std::endl;
//	}
//	// TODO: dotad



	return *this;
}

template<typename SetSpec>
template<typename DynSysSpec>
void DDESolutionCurve<SetSpec>::epsilonShift(
			DynSysSpec const& solver,
			TimePointType const& at_t0, RealType const& epsilon,
			DDESolutionCurve& out_result) const {
	// this will test if the set is of good shape, more or less
	try { out_result.set_r0(*(this->m_r0)); }
	catch (std::logic_error &e) { throw rethrow("DDESolutionCurve::epsilonShift: bad r0 shape in result curve", e); }

	try {
		// TODO: (FUTURE) this is only forward epsilon step. Implement backward step later
		// TODO: (FUTURE) use solver and history to recompute Xi over shorter intervals t_i + [0, \epsi], when computing coefficients.
		// TODO: (NOT URGENT) add checks for sanity of the input before proceeding.
		TimePointType how_far = out_result.t0() - out_result.pastTime();
		RealType HH = ecloseStep(RealType(m_grid(1)));
		RealType EPSI = ecloseStep(epsilon);
		auto ijet_out = out_result.begin();
		auto ijet_this = this->at(at_t0 - how_far);
		for (; ijet_out != out_result.end(); ++ijet_out){
			(**ijet_this).evalAtDelta(epsilon, (**ijet_out)[0]);
			size_type order = (**ijet_out).order();
			for (size_type k = 1; k <= order; ++k){
				(**ijet_this).evalCoeffAtDelta(k, epsilon, (**ijet_out)[k]);
			}
			VectorType xi1 = (**ijet_this).evalCoeffAtDelta(order + 1, HH);
			++ijet_this;
			VectorType xi2 = (**ijet_this).evalCoeffAtDelta(order + 1, EPSI);
			VectorType xi(xi1.dimension());
			capd::vectalg::intervalHull(xi1, xi2, xi);
			(**ijet_out).set_Xi(xi);
		}
		// we assume this extra jet is already computed,
		// otherwise, we would need to call solver to extend solution - we do not want this.
		(**ijet_this).evalAtDelta(epsilon, out_result.m_valueAtCurrent);
		out_result.updateCommonR0();
	} catch (std::logic_error &e) {
		throw rethrow("DDESolutionCurve::epsilonShift: in shift algorithm algorithm", e);
	}
}

template<typename SetSpec>
void DDESolutionCurve<SetSpec>::reduce(
			DDESolutionCurve& out_result) const {
	try {
		if (out_result.length() != this->length()){
			std::ostringstream info;
			info << "DDESolutionCurve::reduce(): different length at output: " << out_result.length();
			info << ", expected: " << this->length();
			throw std::logic_error(info.str());
		}
		auto ijet_out = out_result.begin();
		auto ijet_this = this->begin();
		size_type jetIndex = 0;
		RealType HH = ecloseStep(RealType(m_grid(1)));
		for (; ijet_out != out_result.end(); ++ijet_out, ++ijet_this, ++jetIndex){
			size_type order = (*ijet_out)->order();
			if (order > (*ijet_this)->order()){
				std::ostringstream info;
				info << "DDESolutionCurve::reduce(): output order " << order;
				info << " higher than input " << (*ijet_this)->order();
				info << " at jet index " << jetIndex << ".";
				throw std::logic_error(info.str());
			}
			auto icoeff_out = (*ijet_out)->beginJet();
			auto icoeff_this = (*ijet_this)->beginJet();
			for (; icoeff_out != (*ijet_out)->endJet(); ++icoeff_out, ++icoeff_this){
				(*icoeff_out) = (*icoeff_this);
			}
			VectorType xi = (**ijet_this).evalCoeffAtDelta(order + 1, HH);
			(**ijet_out).set_Xi(xi);
		}
		// we assume this extra jet is already computed,
		// otherwise, we would need to call solver to extend solution - we do not want this.
		out_result.m_valueAtCurrent = this->m_valueAtCurrent;
		out_result.set_r0(this->get_r0());
		out_result.updateCommonR0();
	} catch (std::logic_error &e) {
		throw rethrow("DDESolutionCurve::reduce(): in shift algorithm algorithm", e);
	}
}

} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_DDESOLUTIONCURVE_HPP_ */
