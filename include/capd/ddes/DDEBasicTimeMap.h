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

#ifndef _CAPD_DDES_DDE_BASICTIME_MAP_H_
#define _CAPD_DDES_DDE_BASICTIME_MAP_H_

#include <capd/ddes/DDECommon.h>
#include <capd/poincare/TimeMap.h>

namespace capd{
namespace ddes{

using namespace poincare;

/**
 * Implementation of Time Map for use with nonrigorous DDESolvers.
 *
 * In theory should work with any compatible DDEDynSys.
 */
template<typename DynSysSpec>
class DDEBasicTimeMap : public capd::poincare::TimeMap<DynSysSpec> {
public:
	typedef capd::poincare::TimeMap<DynSysSpec> BaseClass; 	///< I like to have those like in Java
	typedef DDEBasicTimeMap<DynSysSpec> Class;				///< I like to have those like in Java
	typedef DynSysSpec DynSysType;							///< ddesCAPD compatible name, as in DDEBasicPoincareMap, TODO: (NOT URGENT, REFATOR, RETHINK) -> make use of CAPD compatible?
	typedef typename DynSysSpec::CurveType CurveType; 		///< CAPD compatible name, as in TimeMap
	typedef typename CurveType::VectorType VectorType;		///< CAPD compatible name, as in TimeMap
	typedef typename CurveType::MatrixType MatrixType;		///< CAPD compatible name, as in TimeMap
	typedef typename CurveType::ScalarType ScalarType;		///< CAPD compatible name, as in TimeMap
	typedef typename CurveType::size_type size_type;		///< CAPD compatible name, as in TimeMap
	typedef typename CurveType::RealType RealType;							///< ddesCAPD need this
	typedef typename DynSysType::JacobianStorageType JacobianStorageType;	///< ddesCAPD need this
	typedef typename CurveType::TimePointType TimePointType;				///< ddesCAPD need this

	typedef DynSysType Solver;		///< CAPD compatible name, as in TimeMap
	typedef typename Solver::VectorFieldType VectorFieldType;	///< CAPD compatible name, as in TimeMap
	typedef CurveType SolutionCurve;							///< CAPD compatible name, as in TimeMap
	typedef typename CurveType::HessianType HessianType;        ///< CAPD compatible name, as in TimeMap
	typedef typename CurveType::JetType JetType;                ///< CAPD compatible name, as in TimeMap

	DDEBasicTimeMap(DynSysType& solver): BaseClass(solver) {

	}

//	/** compute a constant function */
//	VectorType operator()(ScalarType time, VectorType& v);
//	VectorType operator()(ScalarType time, VectorType& v, ScalarType& in_out_time);
//	VectorType operator()(ScalarType time, VectorType& v, SolutionCurve& solution);
};

} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_DDE_BASICPOINCARE_MAP_H_ */
