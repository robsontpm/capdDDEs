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

#ifndef _CAPD_DDES_SHAREDDOUBLETON_H_
#define _CAPD_DDES_SHAREDDOUBLETON_H_

#include <capd/ddes/DDECommon.h>
#include <capd/ddes/storage/DoubletonInterface.h>
#include <bitset>

namespace capd{
namespace ddes{

/**
 * This class stores a set X \in \R^d of the form:
 *
 *    X = x + C * r0 + B * r
 *
 * with $x \in \R^d$, $r0, r$ - closed intervals centered at 0,
 * $r0 \subset \R^{N_0}$, $C \in Lin(\R^{N_0}, \R^d)$, $B \in Lin(\R^d, \R^d)$, $r \in \R^d$,
 * $B$ being a matrix easy to invert, e.g. B = ID (Interval form of the remainder) or $B$ orthogonal
 * (Doubleton set with QR decomposition).
 *
 * NOTE: compare to the classic Lohner Doubleton Set in CAPD, where authors assume $N_0 = d$)
 * NOTE: the set allows for its components to be set outside of the object and controlled there.
 *       We use pointers for this purpose. However, we always work with the set to assure that
 *       no pointer is NULL. If user does not supply us with a pointer, then we are creating a
 *       default one instead. We use the fact that CAPD can have Vectors of dimension 0 and Matrices
 *       of dimension (d, 0) and (0, d) and algebra is well defined for them.
 *
 * NOTE: This is quite long and tedious code in C++ because of C++
 *       This code is heavily tested for correctness and memory leaks.
 *       Reader (e.g. reviewer of the manuscript for publication) should not worry to check that code
 *
 * TODO: (NOT URGENT) move big implementation into .hpp part.
 *
 */
template<typename MatrixSpec, typename PolicySpec = capd::dynset::IdQRPolicy>
class SharedDoubleton :
		public DoubletonInterface<MatrixSpec, PolicySpec>{
public:
	typedef MatrixSpec MatrixType;
	typedef typename MatrixType::RowVectorType VectorType;
	typedef typename MatrixType::ScalarType ScalarType;
	typedef typename MatrixType::size_type size_type;
	typedef SharedDoubleton Class;
	typedef DoubletonInterface<MatrixSpec, PolicySpec> BaseClass;
	typedef VectorType* VectorTypePtr;
	typedef MatrixType* MatrixTypePtr;
	typedef std::bitset<6> OwnershipType;
	typedef PolicySpec Policy;
	typedef PolicySpec QRPolicy;

	/** assign operator - it needs to deallocate memory if necessary, then setup set anew */
	Class& operator=(Class const & orig){
		deallocate();
		m_owner = orig.m_owner;
		rawSetup(
			orig.dimension(), orig.storageN0(),
			orig.m_x, orig.m_C, orig.m_r0, orig.m_B, orig.m_r, orig.m_Binv,
			m_owner
		);
		return *this;
	}
	/** default constructor makes a 1D point-set at 0 */
	SharedDoubleton():
			m_x(0), m_C(0), m_r0(0), m_B(0), m_r(0), m_Binv(0),
			m_owner(OWN_ALL)
	{
		allocate(0, 0);
	}

	/**
	 * this is to be the conversion the same as in CAPD, that is, the vector is split into mid(x) + C * r0,
	 *
	 * where C = Id, r0 = x - mid(x). The set retains the ownership of the C * r0 part.
	 */
	SharedDoubleton(VectorType const & x):
			m_x(0), m_C(0), m_r0(0), m_B(0), m_r(0), m_Binv(0),
			m_owner(OWN_ALL)
	{
		allocate(x.dimension(), x.dimension(), OWN_ALL);
		m_C->setToIdentity();
		split(x, *m_x, *m_r0);
		// not necessary to sanityCheck(), everything must be sane
	}
	/**
	 * setup set with a given vector and initialize r0.
	 * You can pass NULL as r0, then vector of dim=0 will be used (the only point of \R^0 vector space).
	 *
	 * If it is an interval vector, then it be split into midPoint and the error part, which goes to B*r part (error part).
	 *
	 * If you want some structure (i.e. Lohner part C * r0) then use other constructors.
	 */
	SharedDoubleton(VectorType const & x, VectorType* set_external_r0, bool passOwnership = false):
			m_x(0), m_C(0), m_r0(set_external_r0), m_B(0), m_r(0), m_Binv(0),
			m_owner(OWN_ALL)
	{
		// by default (see initializers list), we own r0,
		// but if user wants otherwise, then we act accordingly
		if (set_external_r0 && !passOwnership)
			m_owner = ~OWN_r0;
		else if (!set_external_r0)
			m_r0 = new VectorType();
		// if m_r0 is already set (either by set_external_r0 or by us in the last statement), then we do not allocate
		// Note: allocate will set also B and invB to Id.
		allocate(x.dimension(), m_r0 ? m_r0->dimension() : 0, m_r0 ? ~OWN_r0 : OWN_ALL);
		split(x, *m_x, *m_r);
		// not necessary to sanityCheck(), everything must be sane
	}
	/** copies the original set (ownership is preserved) */
	SharedDoubleton(Class const & orig):
			m_x(0), m_C(0), m_r0(0), m_B(0), m_r(0), m_Binv(0),
			m_owner(orig.m_owner)
	{
		rawSetup(
			orig.dimension(), orig.storageN0(),
			orig.m_x, orig.m_C, orig.m_r0, orig.m_B, orig.m_r, orig.m_Binv,
			orig.m_owner
		);
		// not necessary to sanityCheck(), everything must be sane
	}
	/** setup this set with a given data, set owns everything */
	SharedDoubleton(VectorType const & x, MatrixType const & C, VectorType const & r0, MatrixType const & B, VectorType const & r):
			m_x(0), m_C(0), m_r0(0), m_B(0), m_r(0), m_Binv(0),
			m_owner(OWN_ALL)
	{
		rawSetup(x.dimension(), r0.dimension(), &x, &C, &r0, &B, &r, NULL, m_owner);
		this->sanityCheck("::__construct__");
	}
	/** setup this set with a given data, set owns everything */
	SharedDoubleton(VectorType const & x, MatrixType const & C, VectorType const & r0, MatrixType const & B, MatrixType const & Binv, VectorType const & r):
			m_x(0), m_C(0), m_r0(0), m_B(0), m_r(0), m_Binv(0),
			m_owner(OWN_ALL)
	{
		rawSetup(x.dimension(), r0.dimension(), &x, &C, &r0, &B, &r, &Binv, m_owner);
		this->sanityCheck("::__construct__");
	}
	/** setup this set with a given data but the set does not own the r0 (user is responsible for deleting it) */
	SharedDoubleton(VectorType const & x, MatrixType const & C, VectorType* r0, MatrixType const & B, VectorType const & r):
			m_x(0), m_C(0), m_r0(0), m_B(0), m_r(r0), m_Binv(0),
			m_owner(OWN_ALL & ~OWN_r0)
	{
		rawSetup(x.dimension(), r0 ? r0->dimension() : 0, &x, &C, r0, &B, &r, NULL, m_owner);
		this->sanityCheck("::__construct__");
	}
	/** setup this set with a given data but the set does not own the r0 (user is responsible for deleting it) */
	SharedDoubleton(VectorType const & x, MatrixType const & C, VectorType* r0, MatrixType const & B, MatrixType const & Binv, VectorType const & r):
			m_x(0), m_C(0), m_r0(0), m_B(0), m_r(r0), m_Binv(0),
			m_owner(OWN_ALL & ~OWN_r0)
	{
		rawSetup(x.dimension(), r0 ? r0->dimension() : 0, &x, &C, r0, &B, &r, &Binv, m_owner);
		this->sanityCheck("::__construct__");
	}
	/** setup this set with a given data, set is owner of everything */
	SharedDoubleton(VectorType const & x, MatrixType const & C, VectorType const & r0):
			m_x(0), m_C(0), m_r0(0), m_B(0), m_r(0), m_Binv(0),
			m_owner(OWN_ALL)
	{
		rawSetup(x.dimension(), r0.dimension(), &x, &C, &r0, NULL, NULL, NULL, m_owner);
		this->sanityCheck("::__construct__");
	}
	/** setup this set with a given data but the set does not own the r0 (user is responsible for deleting it) */
	SharedDoubleton(VectorType const & x, MatrixType const & C, VectorType* r0):
			m_x(0), m_C(0), m_r0(0), m_B(0), m_r(r0), m_Binv(0),
			m_owner(OWN_ALL & ~OWN_r0)
	{
		rawSetup(x.dimension(), r0 ? r0->dimension() : 0, &x, &C, r0, NULL, NULL, NULL, m_owner);
		this->sanityCheck("::__construct__");
	}
	/** setup this set with a given data but the set does not own the data (user is responsible for deleting) */
	SharedDoubleton(VectorType* x, MatrixType* C, VectorType* r0, MatrixType* B, VectorType* r):
			m_x(0), m_C(0), m_r0(0), m_B(0), m_r(r0), m_Binv(0),
			m_owner(OWN_NONE)
	{
		rawSetup(x->dimension(), r0 ? r0->dimension() : 0, x, C, r0, B, r, NULL, m_owner);
		this->sanityCheck("::__construct__");
	}
	/** setup this set as zero vector, but the structure of given dimensions. If second arg is < 0 then d is used instead. */
	SharedDoubleton(size_type d, size_type N0 = -1):
			m_x(0), m_C(0), m_r0(0), m_B(0), m_r(0), m_Binv(0),
			m_owner(OWN_ALL){
		if (N0 < 0) N0 = d;
		rawSetup(d, N0, NULL, NULL, NULL, NULL, NULL, NULL, m_owner);
		// no need for sanity check
	}
	/** standard thing */
	virtual ~SharedDoubleton(){ deallocate(); }

	using BaseClass::dimension; ///< import dimension from interface
	/** see interface description */
	size_type storageN0() const { return m_r0->dimension(); }
	/** see interface description */
	size_type storageDimension() const { return m_x->dimension(); }

	/** unmark itself as owner and returns the pointer. From now on user is responsible for deleting. Use with caution. */
	VectorType* take_x() { m_owner.set(OWNERBIT_x, false); return this->m_x; }
	/** unmark itself as owner and returns the pointer. From now on user is responsible for deleting. Use with caution. */
	MatrixType* take_C() { m_owner.set(OWNERBIT_C, false); return this->m_C; }
	/** unmark itself as owner and returns the pointer. From now on user is responsible for deleting. Use with caution. */
	VectorType* take_r0() { m_owner.set(OWNERBIT_r0, false); return this->m_r0; }
	/** unmark itself as owner and returns the pointer. From now on user is responsible for deleting. Use with caution. */
	MatrixType* take_B() { m_owner.set(OWNERBIT_B, false); return this->m_B; }
	/** unmark itself as owner and returns the pointer. From now on user is responsible for deleting. Use with caution. */
	VectorType* take_r() { m_owner.set(OWNERBIT_r, false); return this->m_r; }
	/** unmark itself as owner and returns the pointer. From now on user is responsible for deleting. Use with caution. */
	MatrixType* take_Binv() { m_owner.set(OWNERBIT_Binv, false); return this->m_Binv; }

	/** see interface docs */
	VectorType get_x() const { return *this->m_x; }
	/** see interface docs */
	MatrixType get_C() const { return *this->m_C; }
	/** see interface docs */
	VectorType get_r0() const { return *this->m_r0; }
	/** see interface docs */
	MatrixType get_B() const { return *this->m_B; }
	/** see interface docs */
	VectorType get_r() const { return *this->m_r; }
	/** see interface docs */
	MatrixType get_Binv() const { return *this->m_Binv; }
	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
	BaseClass& set_x(VectorType const &x) { assureOwner(OWN_x); *m_x = x; sanityCheck("::set_x"); return *this; }
	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
	BaseClass& set_C(MatrixType const &C) { assureOwner(OWN_C); *m_C = C; sanityCheck("::set_C"); return *this; }
	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
	BaseClass& set_r0(VectorType const &r0) { assureOwner(OWN_r0); *m_r0 = r0; sanityCheck("::set_r0"); return *this; }
	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
	BaseClass& set_B(MatrixType const &B) { assureOwner(OWN_B); *m_B = B; updateBinv(); sanityCheck("::set_B"); return *this; }
	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
	BaseClass& set_r(VectorType const &r) { assureOwner(OWN_r); *m_r = r; sanityCheck("::set_r"); return *this; }
	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
	BaseClass& set_Cr0(MatrixType const &C, VectorType const &r0) {
		assureOwner(OWN_C | OWN_r0);
		*m_C = C; *m_r0 = r0;
		sanityCheck("::set_Cr0"); return *this;
	}
	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
	BaseClass& set_Binv(MatrixType const &Binv) { assureOwner(OWN_Binv); *m_Binv = Binv; sanityCheck("::set_Binv"); return *this; }
	/* Reimplemented: is faster than original in interface */
	VectorType midPoint() const { return *m_x; }
	/* Reimplemented: is faster than original in interface */
	VectorType hull() const { return (*m_x) + (*m_C) * (*m_r0) + (*m_B) * (*m_r); }

	/** Note: the element must be compatible with other elements (dimensions!) */
	BaseClass& set_x(VectorType* x, bool passOwnership = false) {
		deallocate(OWN_x); m_x = x; m_owner.set(OWNERBIT_x, passOwnership); sanityCheck("::set_x(ptr)"); return *this;
	}
	/** Note: the element must be compatible with other elements (dimensions!) */
	BaseClass& set_C(MatrixType* C, bool passOwnership = false) {
		deallocate(OWN_C); m_C = C; m_owner.set(OWNERBIT_C, passOwnership); sanityCheck("::set_C(ptr)"); return *this;
	}
	/** Note: the element must be compatible with other elements (dimensions!) */
	BaseClass& set_r0(VectorType* r0, bool passOwnership = false) {
		deallocate(OWN_r0); m_r0 = r0; m_owner.set(OWNERBIT_r0, passOwnership); sanityCheck("::set_r0(ptr)"); return *this;
	}
	/** Note: the element must be compatible with other elements (dimensions!) */
	BaseClass& set_Cr0(MatrixType* C, VectorType* r0, bool passCOwnership = false, bool passR0Ownership = false) {
		deallocate(OWN_C | OWN_r0);
		m_C = C; m_r0 = r0;
		m_owner.set(OWNERBIT_C, passCOwnership);
		m_owner.set(OWNERBIT_r0, passR0Ownership);
		sanityCheck("::set_Cr0(ptr)");
		return *this;
	}
	/** Note: the element must be compatible with other elements (dimensions!) */
	BaseClass& set_B(MatrixType* B, bool passOwnership = false) {
		deallocate(OWN_B); m_B = B; m_owner.set(OWNERBIT_B, passOwnership); updateBinv(); sanityCheck("::set_B(ptr)"); return *this;
	}
	/** Note: the element must be compatible with other elements (dimensions!) */
	BaseClass& set_r(VectorType* r, bool passOwnership = false) {
		deallocate(OWN_r); m_r = r; m_owner.set(OWNERBIT_r, passOwnership); sanityCheck("::set_r(ptr)"); return *this;
	}
	/** Note: the element must be compatible with other elements (dimensions!) */
	BaseClass& set_Binv(MatrixType* Binv, bool passOwnership = false) {
		deallocate(OWN_Binv); m_Binv = Binv; m_owner.set(OWNERBIT_Binv, passOwnership); sanityCheck("::set_Binv(ptr)"); return *this;
	}

	/** checked by pointer equality */
	bool common_x(VectorType const* x) const { return x == m_x; }
	/** checked by pointer equality */
	bool common_C(MatrixType const* C) const { return C == m_C; }
	/** checked by pointer equality */
	bool common_r0(VectorType const* r0) const { return r0 == m_r0; }
	/** checked by pointer equality */
	bool common_B(MatrixType const* B) const { return B == m_B; }
	/** checked by pointer equality */
	bool common_r(VectorType const* r) const { return r == m_r; }
	/** checked by pointer equality */
	bool common_Binv(MatrixType const* Binv) const { return Binv == m_Binv; }

	/** reimplemented for better perfoprmance */
	BaseClass& add(VectorType const & v){
		size_type d = this->dimension();
		if (d != v.dimension()) throw std::logic_error("SharedDoubleton::add(vector): Input vector has incompatible dimension.");
		assureOwner(OWN_x | OWN_r);
		(*m_x) += v; VectorType s(v.dimension());
		capd::vectalg::split(*m_x, s);
		(*m_r) += (*m_Binv) * s;
		return *this;
	}
	/** reimplemented for better perfoprmance */
	BaseClass& add(BaseClass const & set){
		size_type N0 = this->storageN0();
		size_type d = this->dimension();
		if (set.storageN0() != N0) throw std::logic_error("SharedDoubleton::add(set): Sets has incompatible N0 dimensions.");
		if (set.dimension() != d) throw std::logic_error("SharedDoubleton::add(set): Sets has incompatible dimensions.");
		VectorType other_r0 = set.get_r0();
		if (set.common_r0(m_r0)){
			// testing by set.common_r() virtual gives advantage, that we do not need to know underlyings of set.
			// if set is of the good class type, i.e. SharedDoubleton, then we have optimal situation.
			assureOwner(OWN_x | OWN_r | OWN_C);
		}else{
			assureOwner(OWN_x | OWN_r | OWN_C | OWN_r0);
			capd::vectalg::intervalHull(*m_r0, other_r0, *m_r0);
		}
		(*m_C) += set.get_C(); MatrixType S(d, N0);
		capd::vectalg::split(*m_C, S);

		(*m_x) += set.get_x(); VectorType s(d);
		capd::vectalg::split(*m_x, s);
		(*m_r) += ((*m_Binv) * set.get_B()) * set.get_r() + (*m_Binv) * S * (*m_r0) + (*m_Binv) * s;
		return *this;
	}
	/** reimplemented for better performance */
	BaseClass& mul(ScalarType const & c){
		assureOwner(OWN_x | OWN_C | OWN_r);
		size_type d = this->dimension();
		size_type N0 = this->storageN0();
		VectorType s(d);
		MatrixType S(d, N0);
		(*m_x) *= c;
		capd::vectalg::split(*m_x, s);
		(*m_C) *= c;
		capd::vectalg::split(*m_C, S);
		(*m_r) *= c;
		(*m_r) += ((*m_Binv) * S) * (*m_r0) + (*m_Binv) * s;
		return *this;
	}
	/** reimplemented for better perfoprmance */
	BaseClass& mulThenAdd(ScalarType const & c, BaseClass const & set){
		// TODO: (!!!URGENT) reimplement for better performance!
		//throw std::logic_error("SharedDoubleton::mulThenAdd(scalar, set): Not reimplemented yet!");
		this->mul(c); this->add(set); return *this;
	}

	/** applies in a smart way to this set X the affine transform f(y) = M * (X - v) */
	BaseClass& affineTransform(MatrixType const &M, VectorType const &v){
		assureOwner(OWN_x | OWN_C | OWN_r | OWN_B);
		if (M.numberOfRows() != dimension() && M.numberOfColumns() != dimension())
			throw std::logic_error("SharedDoubleton::affineTransform: requires square matrix of same dimension as the set");
		if (v.dimension() != dimension())
			throw std::logic_error("SharedDoubleton::affineTransform: requires v to be of the same dimension as the set");
		*m_x = M * (*m_x  - v);
		if (storageN0()){
			*m_B = M * (*m_B);
			*m_C = M * (*m_C);
		}
		updateBinv();
		// TODO: (!!!URGENT) split m_x and incorporate rest to m_B * m_r?
		// TODO: (!!!URGENT) split m_C and incorporate rest * m_r0 to m_B * m_r?
		// TODO: (!!!URGENT) return m_B to a good shape? (QR strategy)
		return *this;
	}
	/** applies in a smart way to this set X the transform f(y) = X - v */
	BaseClass& translate(VectorType const &v){
		assureOwner(OWN_x);
		if (v.dimension() != dimension())
			throw std::logic_error("SharedDoubleton::translate: requires v to be of the same dimension as the set");
		*m_x += v;
		return *this;
	}

	// TODO: (???) better out handling (add read (cin>>) ability?)
	std::string show() const {
		std::ostringstream oss;
		oss << "SharedDoubleton(x + C * r0 + B * r):\nDimension / N_0: ";
		oss << (*this); return oss.str();
	}
//	// TODO: (???) better out handling (add read (cin>>) ability?)
//	friend std::ostream& operator<<(std::ostream & out, Class const & c) {
//		out << c.dimension() << " " << c.storageN0() << " " << c.m_owner << "\n";
//		out << (*c.m_x);
//		if (c.storageN0()){
//			out << "\n";
//			out << (*c.m_C) << "\n";
//			out << (*c.m_r0) << "\n";
//			out << (*c.m_B) << "\n";
//			out << (*c.m_r) << "\n";
//			out << (*c.m_Binv);
//		}
//		return out;
//	}
	/** see base class */
	virtual void reinitialize(size_type d, size_type N0){
		reallocate(d, N0);
	}

	void reinit(const VectorType& x, const MatrixType& C, const VectorType& r0, const MatrixType& B, const VectorType& r){
		throw std::logic_error("Implementation is incompatible with rawSetup. Rewrite.");
//		deallocate();
//		rawSetup(true, &x, &C, &B, &r, true, &r0);
	}
	void reinit(const VectorType* x, const MatrixType* C, const VectorType* r0, const MatrixType* B, const VectorType* r	){
		throw std::logic_error("Implementation is incompatible with rawSetup. Rewrite.");
//		deallocate();
//		rawSetup(false, x, C, B, r, false, r0);
	}

protected:
	VectorTypePtr m_x;
	MatrixTypePtr m_C;
	VectorTypePtr m_r0;
	MatrixTypePtr m_B;
	VectorTypePtr m_r;
	MatrixTypePtr m_Binv;
	OwnershipType m_owner; 	// ownership
	enum {					// ownership ENUM for readability,
		OWNERBIT_x = 5,		// sequence is reversed to match bitset order
		OWNERBIT_C = 4,		// i.e. least significant on the right and is 0
		OWNERBIT_r0 = 3,
		OWNERBIT_B = 2,
		OWNERBIT_r = 1,
		OWNERBIT_Binv = 0
	};
	static OwnershipType OWN_x; 	/// bitset "I own only x" - used for binary operations
	static OwnershipType OWN_C; 	/// similarly to above - used for binary operations
	static OwnershipType OWN_r0; 	/// similarly to above - used for binary operations
	static OwnershipType OWN_B; 	/// similarly to above - used for binary operations
	static OwnershipType OWN_r; 	/// similarly to above - used for binary operations
	static OwnershipType OWN_Binv; 	/// similarly to above - used for binary operations
	static OwnershipType OWN_ALL; 	/// bitset "I own everything" - used for binary operations
	static OwnershipType OWN_NONE; 	/// bitset "I own nothing" - used for binary operations

	/** assures m_Binv holds real inverse of B */
	void updateBinv(){
		assureOwner(OWN_Binv);
		if (this->m_B){
			MatrixType old_B = *(this->m_B);
			this->Policy::computeBinvB(*(this->m_B), *(this->m_Binv), (*this->m_r));
			// we had part old_B * m_r, but we want to have (new) m_B * new_r enclose the old part
			// note that Policy::computeBinvB can alter both m_B and m_invB
			// if m_B was in a good shape before computing invB then we do not need to correct m_r
			if (this->m_r){
				if (old_B != *this->m_B)
					(*this->m_r) = ((*this->m_Binv) * old_B) * (*this->m_r);
			}
		} else {
			throw std::logic_error("SharedDoubleton::updateBinv(): called, but matrix B not set!");
		}
	}
	/**
	 * Low level setup function to be called in constructors for DRY.
	 * This does not deallocate memory! For safer version see public setup().
	 */
	void rawSetup(
			size_type dim, size_type N0,
			const VectorType* x,
			const MatrixType* C,
			const VectorType* r0,
			const MatrixType* B,
			const VectorType* r,
			const MatrixType* Binv,
			OwnershipType ownership){
		if (!x) ownership |= OWN_x;
		if (!C) ownership |= OWN_C;
		if (!r0) ownership |= OWN_r0;
		if (!B) ownership |= OWN_B;
		if (!r) ownership |= OWN_r;
		if (!Binv) ownership |= OWN_Binv;
		allocate(dim, N0, ownership);
		if (x){ if (ownership[OWNERBIT_x]) *m_x = *x; else m_x = const_cast<VectorType*>(x); }
		if (C){ if (ownership[OWNERBIT_C]) *m_C = *C; else m_C = const_cast<MatrixType*>(C); }
		if (r0){ if (ownership[OWNERBIT_r0]) *m_r0 = *r0; else m_r0 = const_cast<VectorType*>(r0); }
		if (B){ if (ownership[OWNERBIT_B]) *m_B = *B; else m_B = const_cast<MatrixType*>(B); }
		if (r){ if (ownership[OWNERBIT_r]) *m_r = *r; else m_r = const_cast<VectorType*>(r); }
		// r must be set before Binv for updateBinv() to work properly
		if (Binv){
			if (ownership[OWNERBIT_Binv]) *m_Binv = *Binv; else m_Binv = const_cast<MatrixType*>(Binv);
		} else {
			updateBinv();
		}
		m_owner = ownership;
	}
	/** does nothing if is owner of others, otherwise it reasigns the memory and copies values form the external data */
	void assureOwner(OwnershipType const & what = OWN_ALL){
		// we assign to tmp all elements for simplicity, this des not cost much.
		VectorType* tmp_x = m_x;
		MatrixType* tmp_C = m_C;
		VectorType* tmp_r0 = m_r0;
		MatrixType* tmp_B = m_B;
		VectorType* tmp_r = m_r;
		MatrixType* tmp_Binv = m_Binv;
		// Below: to_allocate will give the bits where what is 1 and m_owner is 0
		// please check this! Those are elements that needs to be allocated
		// because they were not owned by the object and we want them to be owned
		// (wee need to copy them by value to newly allocated space)
		// other items are either not required to be owned (0 in what)
		// or are already owned (1 in m_owner).
		// ~(~what | m_owner) = 1 <=> what = 1, owner = 0
		OwnershipType to_allocate = ~(~what | m_owner);
		// below: allocate this is safe, as it will only perform new on the
		// elements that are already owned by someone else (0 in m_owner)
		allocate(dimension(), storageN0(), to_allocate);
		// now we copy the data from those vectors that needs to be copied
		if (to_allocate[OWNERBIT_x]) *m_x = *tmp_x;
		if (to_allocate[OWNERBIT_C]) *m_C = *tmp_C;
		if (to_allocate[OWNERBIT_r0]) *m_r0 = *tmp_r0;
		if (to_allocate[OWNERBIT_B]) *m_B = *tmp_B;
		if (to_allocate[OWNERBIT_r]) *m_r = *tmp_r;
		if (to_allocate[OWNERBIT_Binv]) *m_Binv = *tmp_Binv;
		m_owner |= to_allocate;
	}
	/**
	 * helper function to check if the object represents sane data
	 * (all dimensions compatible and all pointers set if neccesary )
	 * should be called after constructing or manipulating object
	 */
	void sanityCheck(std::string const& what = "") const {
		try {
			if (!m_x) throw std::domain_error("SharedDoubleton"+what+": x is NULL.");
			if (!m_C) throw std::domain_error("SharedDoubleton"+what+": C is NULL.");
			if (!m_B) throw std::domain_error("SharedDoubleton"+what+": B is NULL.");
			if (!m_Binv) throw std::domain_error("SharedDoubleton"+what+": B inverse is NULL.");
			if (!m_r) throw std::domain_error("SharedDoubleton"+what+": r is NULL.");
			if (!m_r0) throw std::domain_error("SharedDoubleton"+what+": r0 is NULL.");
			if (storageN0() != m_C->numberOfColumns()) { throw std::domain_error("SharedDoubleton"+what+": C and r0 have incompatible dimensions."); }
			if (dimension() != m_B->numberOfRows() || dimension() != m_B->numberOfColumns()) {
				std::ostringstream info;
				info << "SharedDoubleton" << what << ": B has incompatible dimensions.";
				info << "was " << dimension() << ", expected " << m_B->numberOfRows() << ".";
				throw std::domain_error(info.str());
			}
			if (dimension() != m_Binv->numberOfRows() || dimension() != m_Binv->numberOfColumns()) {
				std::ostringstream info;
				info << "SharedDoubleton" << what << ": B inverse has incompatible dimensions.";
				info << "was " << dimension() << ", expected " << m_Binv->numberOfRows() << ".";
				throw std::domain_error(info.str());
			}
			if (dimension() != m_C->numberOfRows()) { throw std::domain_error("SharedDoubleton"+what+": C and x has incompatible dimensions."); }
			if (dimension() != m_r->dimension()) { throw std::domain_error("SharedDoubleton"+what+": r has incompatible dimensions."); }
		} catch (...){
			// TODO: (NOT URGENT) maybe add more info to exception?
//			std::cerr << (m_C ? *m_C : MatrixType()) << "|";
//			std::cerr << (m_r0 ? *m_r0 : VectorType()) << "|";
			throw;
		}
	}
	/** safely removes all specified data (default: all) from the object, checking the ownership if neccesary. */
	void deallocate(OwnershipType const & what = OWN_ALL){
		if (what[OWNERBIT_x]) 	capd::ddes::helper_safe_delete(m_x, m_owner[OWNERBIT_x]);
		if (what[OWNERBIT_C]) 	capd::ddes::helper_safe_delete(m_C, m_owner[OWNERBIT_C]);
		if (what[OWNERBIT_r0]) 	capd::ddes::helper_safe_delete(m_r0, m_owner[OWNERBIT_r0]);
		if (what[OWNERBIT_B]) 	capd::ddes::helper_safe_delete(m_B, m_owner[OWNERBIT_B]);
		if (what[OWNERBIT_r]) 	capd::ddes::helper_safe_delete(m_r, m_owner[OWNERBIT_r]);
		if (what[OWNERBIT_Binv]) capd::ddes::helper_safe_delete(m_Binv, m_owner[OWNERBIT_Binv]);
	}
	/** allocate specified elements of the object (default: all) so it is d- dimensional, and has N0 dimensional r0 part. */
	void allocate(size_type d, size_type N0 = 0, OwnershipType const & what = OWN_ALL){
		if (what[OWNERBIT_x])		{ m_x = new VectorType(d); }
		if (what[OWNERBIT_C]) 		{ m_C = new MatrixType(d, N0); }
		if (what[OWNERBIT_r0]) 		{ m_r0 = new VectorType(N0); }
		if (what[OWNERBIT_B]) 		{ m_B = new MatrixType(d, d); m_B->setToIdentity(); }
		if (what[OWNERBIT_r]) 		{ m_r = new VectorType(d); }
		if (what[OWNERBIT_Binv])	{ m_Binv = new MatrixType(d, d); m_Binv->setToIdentity(); }
	}
	/** dealocates and allocates new memory for specified elements (default: all) */
	void reallocate(size_type d, size_type N0 = 0, OwnershipType const & what = OWN_ALL){
		deallocate(what);
		allocate(d, N0, what);
	}
};

// Here we initialize static items. Template language is a mess, read carefully :P
template<typename MatrixSpec, typename PolicySpec>
typename SharedDoubleton<MatrixSpec, PolicySpec>::OwnershipType
SharedDoubleton<MatrixSpec, PolicySpec>::OWN_x("100000");
template<typename MatrixSpec, typename PolicySpec>
typename SharedDoubleton<MatrixSpec, PolicySpec>::OwnershipType
SharedDoubleton<MatrixSpec, PolicySpec>::OWN_C("010000");
template<typename MatrixSpec, typename PolicySpec>
typename SharedDoubleton<MatrixSpec, PolicySpec>::OwnershipType
SharedDoubleton<MatrixSpec, PolicySpec>::OWN_r0("001000");
template<typename MatrixSpec, typename PolicySpec>
typename SharedDoubleton<MatrixSpec, PolicySpec>::OwnershipType
SharedDoubleton<MatrixSpec, PolicySpec>::OWN_B("000100");
template<typename MatrixSpec, typename PolicySpec>
typename SharedDoubleton<MatrixSpec, PolicySpec>::OwnershipType
SharedDoubleton<MatrixSpec, PolicySpec>::OWN_r("000010");
template<typename MatrixSpec, typename PolicySpec>
typename SharedDoubleton<MatrixSpec, PolicySpec>::OwnershipType
SharedDoubleton<MatrixSpec, PolicySpec>::OWN_Binv("000001");
template<typename MatrixSpec, typename PolicySpec>
typename SharedDoubleton<MatrixSpec, PolicySpec>::OwnershipType
SharedDoubleton<MatrixSpec, PolicySpec>::OWN_ALL("111111");
template<typename MatrixSpec, typename PolicySpec>
typename SharedDoubleton<MatrixSpec, PolicySpec>::OwnershipType
SharedDoubleton<MatrixSpec, PolicySpec>::OWN_NONE("000000");

} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_SHAREDDOUBLETON_H_ */
