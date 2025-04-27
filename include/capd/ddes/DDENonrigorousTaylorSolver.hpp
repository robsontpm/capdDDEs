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

#ifndef _CAPD_DDES_NONRIGOROUSTAYLORSOLVER_HPP_
#define _CAPD_DDES_NONRIGOROUSTAYLORSOLVER_HPP_

#include <capd/ddes/DDENonrigorousTaylorSolver.h>
#include <capd/ddes/DDECommon.hpp>
#include <capd/ddes/FunctionalMap.hpp>

////////////////////////////////////////////////////////////////////////
// Below are implementations of more involved functions               //
// (to be separated as usual from the declarations in the .h file     //
// For more complete documentation/explanation please consult .h file //
////////////////////////////////////////////////////////////////////////

namespace capd{
namespace ddes{


template<typename FunctionalMapSpec>
void DDENonrigorousTaylorSolver<FunctionalMapSpec>::oneStep(
		TimePointType const& 		in_t0,
		TimePointType const&		in_th,
		RealType const&				in_h,
		ValueStorageType const&		in_v,
		ValueStorageType& 			out_Phi_coeffs_t0,
		ValueStorageType& 			out_Phi_z)			// for possible future use, I assume that there might be more than only value of x(t+h) !
{
	auto& h = in_h;
	auto& Phi_coeffs_t0 = out_Phi_coeffs_t0;
	auto& Phi_z = out_Phi_z;
	auto& v = in_v;

	// compute both value and Jacobian (in nonrigorous version u is a point vector!)
	m_map.computeDDECoefficients(in_t0, v, Phi_coeffs_t0);
	// now we have Taylor coefficients at t0, we now simply use standard
	// Taylor method, as in case of ODEs
	// TODO: (NOT URGENT) rewrite as iterators?
	size_type k = Phi_coeffs_t0.size() - 1;
	while(true){
		out_Phi_z[0] = Phi_coeffs_t0[k] + h * Phi_z[0];
		if (k == 0) break; else --k; // this prevents from infinite loop in case size_type is unsigned
	}
	// and we are done!
}


template<typename FunctionalMapSpec>
void DDENonrigorousTaylorSolver<FunctionalMapSpec>::oneStep(
		TimePointType const& 		in_t0,
		TimePointType const&		in_th,
		RealType const&				in_h,
		ValueStorageType const&		in_v,
		VariableStorageType const&	in_u,
		JacobianStorageType const&	in_dvdu,
		ValueStorageType& 			out_Phi_coeffs_t0,
		JacobianStorageType& 		out_JacPhi_coeffs_t0,
		ValueStorageType& 			out_Phi_z,			// for possible future use, I assume that there might be more than only value of x(t+h) !
		JacobianStorageType& 		out_JacPhi_z)		// I need storage for as many as there will be entries in out_Phi_z
{
	auto& h = in_h;
	auto& v = in_v;
	auto& u = in_u;
	auto& dvdu = in_dvdu;

	// compute both value and Jacobian (in nonrigorous version u is a point vector!)
	m_map.computeDDECoefficients(in_t0, v, u, dvdu, out_Phi_coeffs_t0, out_JacPhi_coeffs_t0);
	// now we have Taylor coefficients at t0, we now simply use standard
	// Taylor method, as in case of ODEs
	// we explicitely eval jet at t0 as a Taylor series to assure that the
	// jacobian computed along is valid.
	// TODO: (NOT URGENT) rewrite as iterators?
	size_type k = out_Phi_coeffs_t0.size() - 1;
	while(true){
		out_Phi_z[0] = out_Phi_coeffs_t0[k] + h * out_Phi_z[0];
		for (size_type j = 0; j < u.size(); ++j){
			out_JacPhi_z[0][j] = out_JacPhi_coeffs_t0[k][j] + h * out_JacPhi_z[0][j];
		}
		if (k == 0) break; // this prevents from infinite loop in case size_type is unsigned
		--k;
	}

	// and we are done!
}


template<typename FunctionalMapSpec>
void DDENonrigorousTaylorSolver<FunctionalMapSpec>::extractVariationalMatrix(
            SolutionCurve const& curve,
            MatrixType& fullV, MatrixType& reducedV,
            std::vector<size_type> reducedShape, int p_howFar) const {
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


template<typename FunctionalMapSpec>
void DDENonrigorousTaylorSolver<FunctionalMapSpec>::extractVariationalMatrix(
			SolutionCurve const& curve,
			MatrixType& fullV) const {
    try {
    	std::vector<size_type> noShape;
    	extractVariationalMatrix(curve, fullV, fullV, noShape);
    } catch (std::logic_error & e){
    	throw rethrow("DDENonrigorousTaylorSolver::extractVariationalMatrix(SolutionSpec, Matrix&)", e);
    }
}


template<typename FunctionalMapSpec>
void DDENonrigorousTaylorSolver<FunctionalMapSpec>::operator()(
			CurveType& in_out_curve,
			size_type steps) {
	// NOTE: this is not so DRY compared to Jacobian version, but I feel it might be
	//       difficult to make it more DRY. Please remember to update the other code
	//       in case you need to update something  here.
	auto& curve = in_out_curve;

	TimePointType t0 = curve.t0();
	TimePointType th = t0 + curve.getStep();
	RealType h = curve.getStep();

	size_type d = m_map.imageDimension();
	VectorType zero(d);

	// collecting computation data
	size_type coeffs_order = m_maxOrder;
	VariableStorageType u;

	for (size_type i = 0; i < steps; ++i){
		VariableStorageType u; ValueStorageType v; JacobianStorageType dvdu;
		m_map.collectComputationData(t0, th, h, curve, v, u, dvdu, coeffs_order);
		if (coeffs_order > m_maxOrder) coeffs_order = m_maxOrder;

		ValueStorageType Phi_coeffs_t0(coeffs_order + 1);
		ValueStorageType Phi_z(1); Phi_z[0] = zero;

		this->oneStep(t0, th, h, v, Phi_coeffs_t0, Phi_z);

		VariableStorageType new_jet;
		FunctionalMapType::deconvert(Phi_coeffs_t0, new_jet);
		curve.addPiece(new JetType(th, new_jet));
		curve.setValueAtCurrent(Phi_z[0]);

		t0 = th; th = t0 + curve.getStep();
	}
}


template<typename FunctionalMapSpec>
void DDENonrigorousTaylorSolver<FunctionalMapSpec>::operator()(
			CurveType& in_out_curve,
			JacobianStorageType& D){
	// NOTE: this is not so DRY compared to Jacobian version, but I feel it might be
	//       difficult to make it more DRY. Please remember to update the other code
	//       in case you need to update something  here.
	auto& curve = in_out_curve;

	TimePointType t0 = curve.t0();
	TimePointType th = t0 + curve.getStep();
	RealType h = curve.getStep();

	size_type d = m_map.imageDimension();
	VectorType zero(d);
	MatrixType Zero(d, d);

	// collecting computation data
	size_type coeffs_order = m_maxOrder;
	VariableStorageType u; ValueStorageType v; JacobianStorageType dvdu;
	m_map.collectComputationData(t0, th, h, curve, v, u, dvdu, coeffs_order);
	if (coeffs_order > m_maxOrder) coeffs_order = m_maxOrder;

	ValueStorageType Phi_coeffs_t0(coeffs_order + 1);
	ValueStorageType Phi_z(1); Phi_z[0] = zero;
	JacobianStorageType JacPhi_z(1); JacPhi_z[0].resize(u.size(), Zero);
	JacobianStorageType JacPhi_coeffs_t0(coeffs_order + 1);
	for (auto jaci = JacPhi_coeffs_t0.begin(); jaci != JacPhi_coeffs_t0.end(); ++jaci)
		jaci->resize(u.size(), Zero);

	this->oneStep(t0, th, h, v, u, dvdu, Phi_coeffs_t0, JacPhi_coeffs_t0, Phi_z, JacPhi_z);

	typedef typename DataType::MatrixStorageType MatrixStorageType;
	VariableStorageType new_jet;
	// NOTE: important!
	//		we assume that other coeffs (PhiCoeffs)
	//      are dependent on the same list of variables as the z(t)!
	size_type varCount = u[0].getMatrixData().size();
	for (size_type k = 0; k <= coeffs_order; ++k){
		MatrixStorageType DJet_k(varCount, Zero);
		for (size_type i = 0; i < varCount; ++i)
			for (size_type j = 0; j < u.size(); ++j)
				DJet_k[i] += JacPhi_coeffs_t0[k][j] * u[j].getMatrixData()[i];

		new_jet.push_back(DataType(Phi_coeffs_t0[k], DJet_k));
		if (k > 0) D.push_back(DJet_k);
	}
	curve.addPiece(new JetType(th, new_jet));

	MatrixStorageType Dz(varCount, Zero);
	for (size_type i = 0; i < varCount; ++i)
		for (size_type j = 0; j < u.size(); ++j)
			Dz[i] += JacPhi_z[0][j] * u[j].getMatrixData()[i];

	curve.setValueAtCurrent(DataType(Phi_z[0], Dz));
	D.push_back(Dz);
}


} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_NONRIGOROUSTAYLORSOLVER_HPP_ */
