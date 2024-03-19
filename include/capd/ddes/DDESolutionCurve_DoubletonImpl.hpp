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

#ifndef _CAPD_DDES_DDESOLUTIONCURVE_DOUBLETONIMPL_HPP_
#define _CAPD_DDES_DDESOLUTIONCURVE_DOUBLETONIMPL_HPP_

#include <capd/ddes/DDESolutionCurve.h>
#include <capd/ddes/DDECommon.hpp>
#include <capd/ddes/DDEForwardTaylorCurvePiece.hpp>

////////////////////////////////////////////////////////////////////////
// Below are implementations of more involved functions               //
// (to be separated as usual from the declarations in the .h file     //
// For more complete documentation/explanation please consult .h file //
////////////////////////////////////////////////////////////////////////

namespace capd{
namespace ddes{

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::affineTransform(MatrixType const &M, VectorType const &v) {
	throw std::logic_error("DDESolutionCurve::affineTransform: Not implemented yet");
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::translate(typename DDESolutionCurve<SetSpec>::VectorType const &v) {
	throw std::logic_error("DDESolutionCurve::translate: Not implemented yet");
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::VectorType
DDESolutionCurve<SetSpec>::get_x() const {
	VectorType result = makeStorage_x();
	size_type I = 0;
	VectorType v = this->m_valueAtCurrent.get_x();
	for (size_type i = 0; i < v.dimension(); ++I, ++i)
		result[I] = v[i];
	for (auto j = rbegin(); j != rend(); ++j){
		VectorType v = (*j)->get_x();
		for (size_type i = 0; i < v.dimension(); ++I, ++i)
			result[I] = v[i];
	}
	return result;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::MatrixType
DDESolutionCurve<SetSpec>::get_C() const {
	MatrixType result = makeStorage_C();
	size_type I = 0;
	MatrixType C = this->m_valueAtCurrent.get_C();
	for (size_type i = 0; i < C.numberOfRows(); ++I, ++i)
		for (size_type j = 0; j < C.numberOfColumns(); ++j)
			result[I][j] = C[i][j];
	for (auto jet = rbegin(); jet != rend(); ++jet){
		MatrixType C = (*jet)->get_C();
		for (size_type i = 0; i < C.numberOfRows(); ++I, ++i)
			for (size_type j = 0; j < C.numberOfColumns(); ++j)
				result[I][j] = C[i][j];
	}
	return result;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::MatrixType
DDESolutionCurve<SetSpec>::get_B() const {
	MatrixType result = makeStorage_B();
	size_type I = 0;
	MatrixType B = this->m_valueAtCurrent.get_B();
	for (size_type i = 0; i < B.numberOfRows(); ++i)
		for (size_type j = 0; j < B.numberOfColumns(); ++j)
			result[I+i][I+j] = B[i][j];
	I += B.numberOfRows();
	for (auto jet = rbegin(); jet != rend(); ++jet){
		MatrixType B = (*jet)->get_B();
		for (size_type i = 0; i < B.numberOfRows(); ++i)
			for (size_type j = 0; j < B.numberOfRows(); ++j)
				result[I+i][I+j] = B[i][j];
		I += B.numberOfRows();
	}
	return result;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::VectorType
DDESolutionCurve<SetSpec>::get_r() const {
	VectorType result = makeStorage_r();
	size_type I = 0;
	VectorType v = this->m_valueAtCurrent.get_r();
	for (size_type i = 0; i < v.dimension(); ++I, ++i)
		result[I] = v[i];
	for (auto j = rbegin(); j != rend(); ++j){
		VectorType v = (*j)->get_r();
		for (size_type i = 0; i < v.dimension(); ++I, ++i)
			result[I] = v[i];
	}
	return result;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::VectorType
DDESolutionCurve<SetSpec>::get_r0() const {
	return *m_r0;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::MatrixType
DDESolutionCurve<SetSpec>::get_Binv() const {
	MatrixType result = makeStorage_B();
	size_type I = 0;
	MatrixType Binv = this->m_valueAtCurrent.get_Binv();
	for (size_type i = 0; i < Binv.numberOfRows(); ++i)
		for (size_type j = 0; j < Binv.numberOfColumns(); ++j)
			result[I+i][I+j] = Binv[i][j];
	I += Binv.numberOfRows();
	for (auto jet = rbegin(); jet != rend(); ++jet){
		MatrixType Binv = (*jet)->get_Binv();
		for (size_type i = 0; i < Binv.numberOfRows(); ++i)
			for (size_type j = 0; j < Binv.numberOfRows(); ++j)
				result[I+i][I+j] = Binv[i][j];
		I += Binv.numberOfRows();
	}
	return result;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_x(VectorType const &x) {
	if (x.dimension() != storageDimension()){
		std::ostringstream info;
		info << "DDESolutionCurve::set_x: incompatible x dimension, ";
		info<< "is: " << x.dimension() << " should be: " << storageDimension() << ".";
		throw std::logic_error(info.str());
	}
	size_type I = 0;
	size_type d = dimension();
	VectorType* v = new VectorType(d);
	for (size_type i = 0; i < d; ++I, ++i)
		(*v)[i] = x[I];
	m_valueAtCurrent.set_x(v, true);
	for (auto jet = rbegin(); jet != rend(); ++jet){
		for (auto coeff = (*jet)->beginJet(); coeff != (*jet)->endJet(); ++coeff){
			v = new VectorType(d);
			for (size_type i = 0; i < d; ++I, ++i)
				(*v)[i] = x[I];
			coeff->set_x(v, true);
		}
	}
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_C(MatrixType const &C) {
	if (storageN0() != C.numberOfColumns()){
		std::ostringstream info;
		info << "DDESolutionCurve::set_C: incompatible dimension, ";
		info << "number of columns in C: " << C.numberOfColumns() << ", ";
		info << "should be dimension of r0 (N0): " << storageN0() << ".";
		throw std::logic_error(info.str());
	}
	if (C.numberOfRows() != this->storageDimension()){
		std::ostringstream info;
		info << "DDESolutionCurve::set_C: incompatible dimension. ";
		info << "Number of rows in C: " << C.numberOfRows() << ", ";
		info << "should be of dimension: " << this->storageDimension() << ".";
		throw std::logic_error(info.str());
	}
	size_type d = this->dimension();
	size_type N0 = C.numberOfColumns();
	// then we split C into blocks and update Jets and valueAtCurrent
	MatrixType* new_C_block = new MatrixType(d, N0);
	size_type iC = 0;
	for (size_type i = 0; i < d; ++i, ++iC)
		new_C_block->row(i) = C.row(iC);
	m_valueAtCurrent.set_C(new_C_block, true);
	// we traverse jets in reverse order, i.e. from first past currentTime to the pastTime.
	for (auto ijet = this->rbegin(); ijet != this->rend(); ++ijet){
		for (auto icoeff = (**ijet).beginJet(); icoeff < (**ijet).endJet(); ++icoeff){
			new_C_block = new MatrixType(d, N0);
			for (size_type i = 0; i < d; ++i, ++iC)
				new_C_block->row(i) = C.row(iC);
			icoeff->set_C(new_C_block, true);
		}
	}
	// no need to update r0
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_r0(VectorType const &r0) {
	if (r0.dimension() != storageN0())
		throw std::logic_error("DDESolutionCurve::set_r0: new r0 has incompatible dimension to the whole structure");
	*m_r0 = r0; updateCommonR0(); // we update for the situation if the SetType (SetSpec) does not guarantee shared pointers to r0!
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_Cr0(MatrixType const &C, VectorType const &r0) {
	size_type new_N0 = r0.dimension();
	if (new_N0 != C.numberOfColumns()){
		std::ostringstream info;
		info << "DDESolutionCurve::set_Cr0: C and r0 dimensions do not match. ";
		info << "Number of columns in C: " << C.numberOfColumns() << ", ";
		info << "should be dimension of r0 (N0): " << new_N0 << ".";
		throw std::logic_error(info.str());
	}
	if (C.numberOfRows() != this->storageDimension()){
		std::ostringstream info;
		info << "DDESolutionCurve::set_Cr0: C do not match set storage dimension. ";
		info << "Number of rows in C: " << C.numberOfRows() << ", ";
		info << "should be of dimension: " << this->storageDimension() << ".";
		throw std::logic_error(info.str());
	}
	size_type d = this->dimension();
	// we copy by value
	*m_r0 = r0;
	// then we split C into blocks and update Jets and valueAtCurrent
	MatrixType* new_C_block = new MatrixType(d, new_N0);
	size_type iC = 0;
	for (size_type i = 0; i < d; ++i, ++iC)
		new_C_block->row(i) = C.row(iC);
	m_valueAtCurrent.set_Cr0(new_C_block, m_r0, true, false);
	// we traverse jets in reverse order, i.e. from first past currentTime to the pastTime.
	for (auto ijet = this->rbegin(); ijet != this->rend(); ++ijet){
		for (auto icoeff = (**ijet).beginJet(); icoeff < (**ijet).endJet(); ++icoeff){
			new_C_block = new MatrixType(d, new_N0);
			for (size_type i = 0; i < d; ++i, ++iC)
				new_C_block->row(i) = C.row(iC);
			icoeff->set_Cr0(new_C_block, m_r0, true, false);
		}
	}
	// no need to update r0, setCr0() for SetType took care of it in the loops.
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_B(MatrixType const &B) {
	if (B.numberOfColumns() != B.numberOfRows())
		throw std::logic_error("DDESolutionCurve::set_B(): B must be square");
	if (B.numberOfColumns() != storageDimension())
		throw std::logic_error("DDESolutionCurve::set_B(): B has dimensions not compatible with the set");
	size_type offCount = 0;
	std::vector<MatrixType> Bs = extractDiagonalBlocks(B, dimension(), offCount);
	if (offCount)
		throw std::logic_error("DDESolutionCurve::set_B(): off-block-diagonal element of B is nonzero. Check documentation why this is an error.");
	auto b = Bs.begin();
	m_valueAtCurrent.set_B(*b++);
	auto j = this->rbegin();
	for (; j != rend(); ++j)
		for (auto c = (*j)->beginJet(); c != (*j)->endJet(); ++c){
			c->set_B(*b);
	}
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_r(VectorType const &r) {
	if (r.dimension() != storageDimension()){
		std::ostringstream info;
		info << "DDESolutionCurve::set_r(): incompatible x dimension, ";
		info<< "is: " << r.dimension() << " should be: " << storageDimension() << ".";
		throw std::logic_error(info.str());
	}
	if (!capd::vectalg::containsZero(r))
		throw std::logic_error("DDESolutionCurve::set_r(): r does not contain 0 vector.");
	size_type I = 0;
	size_type d = dimension();
	VectorType* v = new VectorType(d);
	for (size_type i = 0; i < d; ++I, ++i)
		(*v)[i] = r[I];
	m_valueAtCurrent.set_r(v, true);
	for (auto jet = rbegin(); jet != rend(); ++jet){
		for (auto coeff = (*jet)->beginJet(); coeff != (*jet)->endJet(); ++coeff){
			VectorType* v = new VectorType(d);
			for (size_type i = 0; i < d; ++I, ++i)
				(*v)[i] = r[I];
			coeff->set_r(v, true);
		}
	}
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_Binv(MatrixType const &invB) {
	if (invB.numberOfColumns() != invB.numberOfRows())
		throw std::logic_error("DDESolutionCurve::set_invB(): invB must be square");
	if (invB.numberOfColumns() != storageDimension())
		throw std::logic_error("DDESolutionCurve::set_invB(): invB has dimensions not compatible with the set");
	size_type offCount = 0;
	std::vector<MatrixType> Bs = extractDiagonalBlocks(invB, dimension(), offCount);
	if (offCount)
		throw std::logic_error("DDESolutionCurve::set_invB(): off-block-diagonal element of invB is nonzero. Check documentation why this is an error.");
	auto b = Bs.begin();
	m_valueAtCurrent.set_Binv(*b++);
	auto j = this->rbegin();
	for (; j != rend(); ++j)
		for (auto c = (*j)->beginJet(); c != (*j)->endJet(); ++c)
			c->set_Binv(*b++);
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_x(VectorType* x, bool passOwnership) {
	if (x->dimension() != storageDimension()){
		std::ostringstream info;
		info << "DDESolutionCurve::set_x(ptr, flag): incompatible x dimension, ";
		info<< "is: " << x->dimension() << " should be: " << storageDimension() << ".";
		throw std::logic_error(info.str());
	}
	if (!x)
		throw std::logic_error("DDESolutionCurve::set_x(ptr, flag): x is NULL");
	size_type I = 0;
	VectorType v(dimension());
	for (size_type i = 0; i < v.dimension(); ++I, ++i)
		v[i] = (*x)[I];
	m_valueAtCurrent.set_x(v);
	for (auto jet = rbegin(); jet != rend(); ++jet){
		VectorType v( (*jet)->storageDimension());
		for (size_type i = 0; i < v.dimension(); ++I, ++i)
			v[i] = (*x)[I];
		(*jet)->set_x(v);
	}
	helper_safe_delete(x, passOwnership);
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_C(MatrixType* C, bool passOwnership) {
	if (!C) throw std::logic_error("DDESolutionCurve::set_C(ptr, flag): C is NULL");
	try { set_C(*C); } catch (std::logic_error& e) { rethrow("DDESolutionCurve::set_C(ptr, flag):", e); }
	helper_safe_delete(C, passOwnership);
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_r0(VectorType* r0, bool passOwnership) {
	if (r0->dimension() != storageN0())
		throw std::logic_error("DDESolutionCurve::set_r0(ptr, flag): new r0 has incompatible dimension to the whole structure");
	if (!r0)
		throw std::logic_error("DDESolutionCurve::set_r0(ptr, flag): r0 is NULL");
	if (passOwnership){
		deallocateR0();
		m_r0 = r0;
	}else{
		reallocateR0();
		*m_r0 = *r0;
	}
	updateCommonR0(); // for underlying sets that do not share explicitely the r0
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_Cr0(MatrixType* C, VectorType* r0, bool passCOwnership, bool passR0Ownership) {
	//
	// TODO: (!!!URGENT) this does not take into account passOwnership and also forgets previous ptr! FIX FIX FIX! FIX ALSO OTHERS!
	throw std::logic_error("DDESolutionCurve::set_Cr0(ptr, ptr, flag, flag): Not implemented yet");
	size_type new_N0 = r0->dimension();
	if (new_N0 != C->numberOfColumns()){
		std::ostringstream info;
		info << "DDESolutionCurve::set_Cr0(ptr, ptr, flag, flag): C and r0 dimensions do not match. ";
		info << "Number of columns in C: " << C->numberOfColumns() << ", ";
		info << "should be dimension of r0 (N0): " << new_N0 << ".";
		throw std::logic_error(info.str());
	}
	if (C->numberOfRows() != this->storageDimension()){
		std::ostringstream info;
		info << "DDESolutionCurve::set_Cr0(ptr, ptr, flag, flag): C do not match set storage dimension. ";
		info << "Number of rows in C: " << C->numberOfRows() << ", ";
		info << "should be of dimension: " << this->storageDimension() << ".";
		throw std::logic_error(info.str());
	}
	size_type d = this->dimension();
	// we copy by value, as SolutionCurve do not implement ownership logic
	(*m_r0) = (*r0);
	// then we split C into blocks and update Jets and valueAtCurrent
	MatrixType* new_C_block = new MatrixType(d, new_N0);
	size_type iC = 0;
	for (size_type i = 0; i < d; ++i, ++iC)
		new_C_block->row(i) = C->row(iC);
	m_valueAtCurrent.set_Cr0(new_C_block, m_r0, true, false);
	// we traverse jets in reverse order, i.e. from first past currentTime to the pastTime.
	for (auto ijet = this->rbegin(); ijet != this->rend(); ++ijet){
		for (auto icoeff = (**ijet).beginJet(); icoeff < (**ijet).endJet(); ++icoeff){
			new_C_block = new MatrixType(d, new_N0);
			for (size_type i = 0; i < d; ++i, ++iC)
				new_C_block->row(i) = C->row(iC);
			icoeff->set_Cr0(new_C_block, m_r0, true, false);
		}
	}
	helper_safe_delete(C, passCOwnership);
	helper_safe_delete(r0, passR0Ownership);
	// no need to update r0, setCr0() for SetType took care of it in the loops.
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_B(MatrixType* B, bool passOwnership) {
	if (!B) throw std::logic_error("DDESolutionCurve::set_B(ptr, flag): B is NULL");
	try { set_B(*B); } catch (std::logic_error& e) { rethrow("DDESolutionCurve::set_B(ptr, flag):", e); }
	helper_safe_delete(B, passOwnership);
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_r(VectorType* r, bool passOwnership) {
	if (!r) throw std::logic_error("DDESolutionCurve::set_r(ptr, flag): r is NULL");
	try { set_r(*r); } catch (std::logic_error& e) { rethrow("DDESolutionCurve::set_r(ptr, flag):", e); }
	helper_safe_delete(r, passOwnership);
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::set_Binv(MatrixType* invB, bool passOwnership) {
	if (!invB) throw std::logic_error("DDESolutionCurve::set_Binv(ptr, flag): invB is NULL");
	try { set_Binv(*invB); } catch (std::logic_error& e) { rethrow("DDESolutionCurve::set_Binv(ptr, flag):", e); }
	helper_safe_delete(invB, passOwnership);
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::VectorType* DDESolutionCurve<SetSpec>::take_x(){
	return new VectorType(this->get_x());
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::MatrixType* DDESolutionCurve<SetSpec>::take_C(){
	return new MatrixType(this->get_C());
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::VectorType* DDESolutionCurve<SetSpec>::take_r0(){
	return new VectorType(this->get_r0());
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::MatrixType* DDESolutionCurve<SetSpec>::take_B(){
	return new MatrixType(this->get_B());
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::VectorType* DDESolutionCurve<SetSpec>::take_r(){
	return new VectorType(this->get_r());
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::MatrixType* DDESolutionCurve<SetSpec>::take_Binv(){
	return new MatrixType(this->get_Binv());
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::mul(ScalarType const &c) {
	m_valueAtCurrent.mul(c);
	for (auto j = begin(); j != end(); ++j) (*j)->mul(c);
	return *this;
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::add(VectorType const &set) {
	// TODO: (NOT URGENT) implement
	throw std::logic_error("DDESolutionCurve::add(vector): Not implemented yet");
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::add(BaseClass const &set) {
	// TODO: (NOT URGENT) implement
	throw std::logic_error("DDESolutionCurve::add(set): Not implemented yet");
}

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::BaseClass&
DDESolutionCurve<SetSpec>::mulThenAdd(ScalarType const &c, BaseClass const &set) {
	// TODO: (NOT URGENT) implement
	throw std::logic_error("DDESolutionCurve::mulThenAdd: Not implemented yet");
}

} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_DDESOLUTIONCURVE_DOUBLETONIMPL_HPP_ */
