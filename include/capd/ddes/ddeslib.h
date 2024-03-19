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

#ifndef _CAPD_DDES_DDESSLIB_H_
#define _CAPD_DDES_DDESSLIB_H_

#include <capd/ddes/DDECommon.hpp>
#include <capd/ddes/storage/storage.h>
#include <capd/ddes/DDEForwardTaylorCurvePiece.hpp>
#include <capd/ddes/DDEPiecewisePolynomialCurve.hpp>
#include <capd/ddes/DDESolutionCurve.hpp>
#include <capd/ddes/SampleEqns.hpp>
#include <capd/ddes/FunctionalMap.hpp>
#include <capd/ddes/BasicDiscreteDelaysFunctionalMap.h> // TODO: (!!!URGENT) PRAWDOPODOBNIE USUNAC JAK SKONCZE REFACTOR!
#include <capd/ddes/DiscreteDelaysFunctionalMap.hpp>
#include <capd/ddes/DDENonrigorousTaylorSolver.hpp>
#include <capd/ddes/DDETaylorSolver.hpp>
#include <capd/ddes/DDEJetSection.hpp>
#include <capd/ddes/DDEBasicPoincareMap.hpp> // TODO: (!!!URGENT) PRAWDOPODOBNIE USUNAC JAK SKONCZE REFACTOR! (PoincareMap bedzie sciagac dane z tego)
#include <capd/ddes/DDEPoincareMap.hpp>

#endif /* _CAPD_DDES_DDESSLIB_H_ */
