/*
 * DDEHelperRigorous.hpp
 *
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

#ifndef _CAPD_DDEHELPERRIGOROUS_HPP_
#define _CAPD_DDEHELPERRIGOROUS_HPP_

#include <capd/ddeshelper/DDEHelperRigorous.h>

namespace capd {
namespace ddeshelper {

template<typename EqSpec, int delaysSpec, typename MatrixSpec, typename VectorSpec, typename PoliciesSpec>
void RigorousHelper<EqSpec, delaysSpec, MatrixSpec, VectorSpec, PoliciesSpec>::loadSetup(std::string filepath){
	std::ifstream in(filepath);
	try {
		rawLoadSetup(in, m_p, m_n, m_params, m_coords);

		// TODO: because I have to implement binary save of the value on section i hereby dump the read value and compute it rigorously
		// TODO: change it in the final version
		m_coords.secvalue = m_coords.vsection * m_coords.reference;

		updateGrid();
		updateSteps();
	} catch (std::logic_error &e) {
		throw capd::ddes::rethrow("RigorousHelper::loadSetup(filepath)", e);
	}
}

template<typename EqSpec, int delaysSpec, typename MatrixSpec, typename VectorSpec, typename PoliciesSpec>
void RigorousHelper<EqSpec, delaysSpec, MatrixSpec, VectorSpec, PoliciesSpec>::rawLoadSetup(
			std::istream& in,
			int &p, int &n,
			Vector &params,
			CoordinateSystem& coords
			){
	int numpars, d;
	try {
		in >> numpars; if (numpars != PARAMS_COUNT) {
			std::ostringstream info;
			info << "RigorousHelper::rawLoadSetup(): incompatible number of parameters, is " << numpars << ", should be " << PARAMS_COUNT << ".";
			throw std::logic_error(info.str());
		}
		params = ParamsVector(PARAMS_COUNT);
		in >> params;
		in >> d; if (d != DIMENSION) {
			std::ostringstream info;
			info << "RigorousHelper::rawLoadSetup(): incompatible dimension, is " << d << ", should be " << DIMENSION << ".";
			throw std::logic_error(info.str());
		}
		in >> p >> n;
		coords = CoordinateSystem(d * (1 + p * (n+1)));
		in >> coords;
	} catch (std::logic_error& e) {
		throw capd::ddes::rethrow("RigorousHelper::rawLoadSetup(): ", e);
	} catch (...) {
		throw std::logic_error("RigorousHelper::rawLoadSetup(): unspecified input error.");
	}
}

} // namespace ddeshelper
} // namespace capd


#endif /* _CAPD_DDEHELPERRIGOROUS_HPP_ */

