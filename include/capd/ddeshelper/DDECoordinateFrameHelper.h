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

#ifndef _CAPD_DDECOORDINATEFRAMEHELPER_H_
#define _CAPD_DDECOORDINATEFRAMEHELPER_H_

#include <capd/ddes/ddeslib.h>
#include <vector>

namespace capd {
namespace ddeshelper {

/**
 * Internal class to hold a affine coordinate change
 * on a (p,n)-fset.
 *
 * TODO: more docs...
 */
template<
	typename MatrixSpec,
	typename VectorSpec = typename MatrixSpec::VectorType,
	typename ScalarSpec = typename VectorSpec::ScalarType
>
class CoordinateFrame {
public:
	typedef MatrixSpec MatrixType;
	typedef VectorSpec VectorType;
	typedef ScalarSpec ScalarType;
	int M;
	VectorType reference;
	VectorType vsection;
	ScalarType secvalue;
	capd::poincare::CrossingDirection crossingDirection;
	MatrixType coords;
	MatrixType inverseCoords;
	std::string filepath;

	CoordinateFrame(VectorType x, VectorType sec, ScalarType secval, capd::poincare::CrossingDirection d, MatrixType C, MatrixType invC):
		M(x.dimension()), reference(x), vsection(sec), secvalue(secval), crossingDirection(d), coords(C), inverseCoords(invC)
	{}

	CoordinateFrame(VectorType x, MatrixType C, MatrixType invC, capd::poincare::CrossingDirection d = capd::poincare::MinusPlus):
		M(x.dimension()), reference(x), vsection(C.column(0)), secvalue(C.column(0) * x), crossingDirection(d), coords(C), inverseCoords(invC)
	{}

	explicit CoordinateFrame(int M):
		M(M), reference(M), vsection(M), crossingDirection(capd::poincare::MinusPlus), coords(M, M), inverseCoords(M, M)
	{
		coords.setToIdentity();
		inverseCoords.setToIdentity();
		vsection[0] = 1.0;
	}

	void setCrossingDirection(capd::poincare::CrossingDirection d){ crossingDirection = d; }

	/** TODO: add out<< operator... TODO: docs... */
	friend std::istream& operator>>(std::istream& in, CoordinateFrame& coords){
		int dir;
		in >> dir;
		if (dir < 0) coords.crossingDirection = capd::poincare::MinusPlus;
		else if (dir > 0) coords.crossingDirection = capd::poincare::PlusMinus;
		else if (dir == 0) coords.crossingDirection = capd::poincare::Both;
		in >> coords.secvalue;
		relativeReadData(in, coords.vsection);
		relativeReadData(in, coords.reference);
		relativeReadData(in, coords.coords);
		relativeReadData(in, coords.inverseCoords);
		return in;
	}

	friend bool operator==(CoordinateFrame& first, CoordinateFrame& second){
		return (first.reference == second.reference) &&
				(first.coords == second.coords);
	}
	friend bool operator!=(CoordinateFrame& first, CoordinateFrame& second){
		return !(first == second);
	}

	VectorType inCoords(VectorType& x){ return inverseCoords * (x - reference); }
	VectorType inGlobal(VectorType& x){ return reference + coords * x; }

	VectorType x0(){ return reference; }
	MatrixType C(){ return coords; }
	MatrixType invC(){ return inverseCoords; }
	VectorType sectionVector(){ return vsection; }
	ScalarType sectionValue(){ return secvalue; }
};

} // namespace ddeshelper
} // namespace capd

#endif /* _CAPD_DDECOORDINATEFRAMEHELPER_H_ */
