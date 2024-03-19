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

#ifndef _CAPD_DDES_NONRIGOROUSTAYLORSOLVER_H_
#define _CAPD_DDES_NONRIGOROUSTAYLORSOLVER_H_

#include <stdexcept>
#include <string>
#include <capd/ddes/DDECommon.h>
#include <capd/ddes/FunctionalMap.h>

namespace capd{
namespace ddes{

/**
 * a class to hold additional data to compute Jacobian of the
 * (nonrigorous / approximate) flow w.r.t. initial data.
 */
template<typename VectorSpec, typename MatrixSpec>
class VectorWithJacData : public VectorSpec {
public:
	typedef MatrixSpec MatrixType;
	typedef typename MatrixSpec::ScalarType ScalarType;
	typedef VectorWithJacData<VectorSpec, MatrixSpec> Class;
	typedef VectorSpec BaseClass;
	typedef typename BaseClass::size_type size_type;
	typedef std::vector<MatrixSpec> MatrixStorageType;

	VectorWithJacData(size_type d = 0): BaseClass(d) {}
	VectorWithJacData(const VectorWithJacData& orig): BaseClass(orig), m_Jac(orig.m_Jac) {}
	VectorWithJacData(const VectorSpec& orig): BaseClass(orig) {}
	VectorWithJacData(const VectorSpec& v, const std::vector<MatrixSpec>& D): BaseClass(v), m_Jac(D) { }
	VectorWithJacData(const VectorSpec& v, const MatrixSpec& D): BaseClass(v) {
		try {
			setMatrix(D);
		} catch (std::logic_error& e){
			throw rethrow("VectorWithJacData::__construct__(): ", e);
		}
	}

	VectorWithJacData& setMatrix(MatrixType const& D){
		size_type d = this->dimension();
		if (D.numberOfRows() != d)
			throw std::logic_error("VectorWithJacData::setMatrix(): D has bad number of rows");
		size_type I = 0;
		m_Jac.clear();
		while(I < D.numberOfColumns()){
			MatrixSpec M(d, d);
			// first go along columns (i, I),
			// see that I changes same number as i, so after loop I is increased by d
			// then we go with j along each row for each column.
			// this way we partitions matrix D into matrices M of dim (d,d)
			for (size_type i = 0; i < d && I < D.numberOfColumns(); ++i, ++I)
				for (size_type j = 0; j < d; j++)
					M[j][i] = D[j][I];
			m_Jac.push_back(M);
		}
		return *this;
	}

	VectorWithJacData& operator=(const VectorWithJacData& orig){
		BaseClass::operator=(BaseClass(orig));
		m_Jac = orig.m_Jac;
		return *this;
	}

	VectorWithJacData& operator=(const BaseClass& orig){
		BaseClass::operator=(orig);
		m_Jac.clear();
		return *this;
	}

	VectorWithJacData& operator*=(ScalarType const& c){
		this->VectorSpec::operator*=(c);
		for (auto i = m_Jac.begin(); i != m_Jac.end(); ++i) (*i) *= c;
		return *this;
	}

	VectorWithJacData& operator+=(VectorWithJacData const& other){
		this->VectorSpec::operator+=(other);
		auto i = m_Jac.begin();
		auto o = other.m_Jac.begin();
		for (; i != m_Jac.end() && o != other.m_Jac.end(); ++i, ++o) (*i) += (*o);
		return *this;
	}

	operator VectorSpec() { return *this; }
	operator MatrixSpec() {
		size_type d = this->dimension();
		size_type N0 = d * m_Jac.size();
		size_type I = 0;
		MatrixSpec M(d, N0);
		for (auto jaci = m_Jac.begin(); jaci != m_Jac.end(); ++jaci)
			for (size_type i = 0; i < d; ++i, ++I)
				for (size_type j = 0; j < d; ++j)
					M[j][I] = (*jaci)[j][i];
		return M;
	}

	const MatrixStorageType& getMatrixData() const { return m_Jac; }
	MatrixStorageType& getMatrixData() { return m_Jac; }

protected:
	std::vector<MatrixSpec> m_Jac;
};

template<typename FunctionalMapSpec>
class DDENonrigorousTaylorSolver {
public:
	typedef FunctionalMapSpec FunctionalMapType;
	typedef typename FunctionalMapType::CurveType SolutionCurveType;
	typedef typename FunctionalMapType::CurveType CurveType;
	typedef typename FunctionalMapType::JetType JetType;
	typedef typename FunctionalMapType::DataType DataType; 	// it is in fact SetType
	typedef typename FunctionalMapType::VectorType VectorType;
	typedef typename FunctionalMapType::MatrixType MatrixType;
	typedef typename FunctionalMapType::TimePointType TimePointType;
	typedef typename FunctionalMapType::RealType RealType;
	typedef typename FunctionalMapType::size_type size_type;
	typedef typename FunctionalMapType::VariableStorageType VariableStorageType;
	typedef typename FunctionalMapType::JacobianStorageType JacobianStorageType;
	typedef typename FunctionalMapType::ValueStorageType ValueStorageType;

	DDENonrigorousTaylorSolver(DDENonrigorousTaylorSolver const& solver): m_map(solver.m_map), m_maxOrder(solver.m_maxOrder) {}
	DDENonrigorousTaylorSolver(FunctionalMapType const& map, size_type maxOrder = 20): m_map(map), m_maxOrder(maxOrder) {}

	/** a nice operator form of the one step solver */
	void operator()(CurveType&	 in_out_curve){
		// TODO: DRY (see other)
		auto& curve = in_out_curve;

		TimePointType t0 = curve.t0();
		TimePointType th = t0 + curve.getStep();
		RealType h = curve.getStep();

		size_type d = m_map.imageDimension();
		VectorType zero(d);

		// collecting computation data
		size_type coeffs_order = m_maxOrder;
		VariableStorageType u;
		m_map.collectComputationData(t0, th, h, curve, u, coeffs_order);
		if (coeffs_order > m_maxOrder) coeffs_order = m_maxOrder;
		ValueStorageType v; FunctionalMapType::convert(u, v);

		ValueStorageType Phi_coeffs_t0(coeffs_order + 1);
		ValueStorageType Phi_z(1); Phi_z[0] = zero;

		this->oneStep(t0, th, h, v, Phi_coeffs_t0, Phi_z);

		VariableStorageType new_jet;
		FunctionalMapType::deconvert(Phi_coeffs_t0, new_jet);
		curve.addPiece(new JetType(th, new_jet));
		curve.setValueAtCurrent(Phi_z[0]);
	}

	/**
	 * a nice operator form of the one step solver
	 * NOTE IMPORTANT: we assume that DataType of Jet can hold value and Matrix
	 *                 we provide candidate for this: VectorWithJacobianData
	 *                 you cannot use this function with the basic data structure for SolutionCurve
	 *                 this is also motivated with the optimization in computation of rhs (if we do not use all past data, but only e.g. single discrete delay)
	 *                 the VectorWithJacobianData is to bind value of a variable with the derivative of coefficients w.r.t. initial data.
	 *                 Please check examples on how to use this in your code.
	 * TODO: (FUTURE): make this more user friendly...
	 */
	void operator()(CurveType&	 in_out_curve, JacobianStorageType& D){
		// TODO: DRY (see other)
		auto& curve = in_out_curve;

		TimePointType t0 = curve.t0();
		TimePointType th = t0 + curve.getStep();
		RealType h = curve.getStep();

		size_type d = m_map.imageDimension();
		VectorType zero(d);
		MatrixType Zero(d, d);

		// collecting computation data
		size_type coeffs_order = m_maxOrder;
		VariableStorageType u;
		m_map.collectComputationData(t0, th, h, curve, u, coeffs_order);
		if (coeffs_order > m_maxOrder) coeffs_order = m_maxOrder;
		ValueStorageType v; FunctionalMapType::convert(u, v);

		ValueStorageType Phi_coeffs_t0(coeffs_order + 1);
		ValueStorageType Phi_z(1); Phi_z[0] = zero;
		JacobianStorageType JacPhi_z(1); JacPhi_z[0].resize(u.size(), Zero);
		JacobianStorageType JacPhi_coeffs_t0(coeffs_order + 1);
		for (auto jaci = JacPhi_coeffs_t0.begin(); jaci != JacPhi_coeffs_t0.end(); ++jaci)
			jaci->resize(u.size(), Zero);

		this->oneStep(t0, th, h, v, Phi_coeffs_t0, JacPhi_coeffs_t0, Phi_z, JacPhi_z);

		// FunctionalMapType::deconvert(Phi_coeffs_t0, new_jet);
		typedef typename DataType::MatrixStorageType MatrixStorageType;
		VariableStorageType new_jet;
		// IMPORTANT: we assume that other coeffs (PhiCoeffs) are dependent on the same list of variables as the z(t)
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

	/**
	 * this function produces approximation to the solution
	 * after time h = distance between grid points in the input curve.
	 * The var out_Phi_coeffs_t0 will hold the Jet coefficients at point
	 * t0, then out_JacPhi_coeffs_t0 will be the Jacobian of
	 * coeffs w.r.t. variables in out_u.
	 * Then out_Phi_x[0] is the value at t = t0+h computed by the Taylor method,
	 * i.e. evaluation of out_Phi_coeffs_t0(t0+h). out_JacPhi_x[0] is
	 * like in Taylor method for ODEs.
	 *
	 * It is designed to be used internally. Operators() and epsilonShift() use this.
	 */
	void oneStep(
			TimePointType const& 		in_t0,
			TimePointType const&		in_th,
			RealType const&				in_h,
			ValueStorageType const&		in_u,
			ValueStorageType& 			out_Phi_coeffs_t0,
			JacobianStorageType& 		out_JacPhi_coeffs_t0,
			ValueStorageType& 			out_Phi_z,			// for possible future use, I assume that there might be more than only value of x(t+h) !
			JacobianStorageType& 		out_JacPhi_z)		// I need storage for as many as there will be entries in out_Phi_z
	{
		auto& h = in_h;
		auto& u = in_u;

		// compute both value and Jacobian (in nonrigorous version u is a point vector!)
		m_map.computeDDECoefficients(in_t0, u, out_Phi_coeffs_t0, out_JacPhi_coeffs_t0);
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

	/** same as the other oneStep, but without computing JacPhi */
	void oneStep(
			TimePointType const& 		in_t0,
			TimePointType const&		in_th,
			RealType const&				in_h,
			ValueStorageType const&		in_u,
			ValueStorageType& 			out_Phi_coeffs_t0,
			ValueStorageType& 			out_Phi_z)			// for possible future use, I assume that there might be more than only value of x(t+h) !
	{
		auto& h = in_h;
		auto& Phi_coeffs_t0 = out_Phi_coeffs_t0;
		auto& Phi_z = out_Phi_z;
		auto& u = in_u;

		// compute both value and Jacobian (in nonrigorous version u is a point vector!)
		m_map.computeDDECoefficients(in_t0, u, Phi_coeffs_t0);
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

	FunctionalMapType const& getMap() const { return m_map; }
	FunctionalMapType & getMap() { return m_map; }

private:
	/** used in basic computations */
	FunctionalMapType m_map;
	/** used in computation with Jacobian / Monodromy matrix */
//	JacobianNodeFunctionalMapType m_JacMap;
	/** maximal order of jet produced at each step. This might be reduced by the max allowed order due to the Curve continuity. */
	size_type m_maxOrder;
};



} // namespace ddes
} // namespace capd


#endif /* _CAPD_DDES_NONRIGOROUSTAYLORSOLVER_H_ */
