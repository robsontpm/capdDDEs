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

#ifndef _CAPD_DDES_BASICDOUBLETON_H_
#define _CAPD_DDES_BASICDOUBLETON_H_

#include <capd/ddes/storage/DoubletonInterface.h>

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
 * NOTE: This is the most basic implementation, I use bare Vectors and Matrices
 *       and all the data is stored inside the set. For DDEs as described in notes it is
 *       somehow innefficient, as each component will contain its copy of r0 vector
 *       which in the computation in fact is common to all the components.
 *       Therefore we will store O(size of repr.) elements of size dim(r0)
 *       instead of O(1). SharedDoubleton is designed for this purpose, but is is
 *       highly complicated as of now (and I have some problems with memory leak - o be corrected later).
 *       I have decided to implement simple solution to test things without
 *       need to worry about the problems with dynamical memory allocation.
 *       The code is also simpler to read.
 *
 * NOTE: This is quite long and tedious code in C++ because of C++
 *       This should be heavily tested for correctedness and memory leaks
 *       Reader (e.g. reviewer of the manuscript for publication) should not worry to check that code
 *
 * NOTE: THIS PARTICULAR IMPLEMENTATION does not use shared memory, everything is static
 *       (from the point of view of implkemnentation, afaik CAPD uses shared ptr inside the library)
 *       Its main purpose is testing, as it will be (probably) slower and (for sure) more memory consuming (2x at least)
 *       than a SharedDoubleton version. SharedDoubleton will be used in applications.
 *
 * TODO: (NOT URGENT) implement QR policy in BinvB (now IdPolicy is hardcoded). This is not urgent, as I use SharedDoubleton in my computations.
 */
template<typename MatrixSpec, typename PoliciesSpec = capd::dynset::IdQRPolicy>
class BasicDoubleton : public DoubletonInterface<MatrixSpec, PoliciesSpec>{
public:
	typedef MatrixSpec MatrixType;
	typedef typename MatrixType::RowVectorType VectorType;
	typedef typename MatrixType::ScalarType ScalarType;
	typedef typename MatrixType::size_type size_type;
	typedef BasicDoubleton Class;
	typedef DoubletonInterface<MatrixSpec> BaseClass;
	typedef VectorType* VectorTypePtr;
	typedef MatrixType* MatrixTypePtr;
	typedef std::bitset<5> OwnershipType;
	typedef PoliciesSpec Policy;
	typedef PoliciesSpec QRPolicy;

	/** assign operator - it needs to deallocate memory if neccessery, then setup set anew */
	Class& operator=(Class const & orig){ setupFromOther(orig); return *this; }
	/** default constructor makes a 1D point-set at 0 */
	BasicDoubleton() {}
	/**
	 * setup set with a given vector. You can pass common_r0 and passOwnership = false (default)
	 * if you want to set externally controlled r0 element. If passOwnership = true
	 * then the set will take responsibility for deallocating item.
	 */
	BasicDoubleton(VectorType const & x, VectorType* common_r0 = 0, bool passOwnership = false){
		size_type N0 = common_r0 ? common_r0->dimension() : 0;
		setupDimension(x.dimension(), N0); // sets m_B to Id
		split(x, m_x, m_r);
		if (common_r0) m_r0 = *common_r0;
		helper_safe_delete(common_r0, passOwnership);
	}
	/** copies the original set (ownership is preserved) */
	BasicDoubleton(Class const & orig){
		setupFromOther(orig);
		// no need to sanity check. Orig must be sane.
	}
	/** setup this set with a given data, set owns everything */
	BasicDoubleton(VectorType const & x, MatrixType const & C, VectorType const & r0, MatrixType const & B, VectorType const & r){
		setupFromData(x, C, r0, B, r);
		this->sanityCheck("::__construct__");
	}
	/** setup this set with a given data, set owns everything */
	BasicDoubleton(VectorType const & x, MatrixType const & C, VectorType const & r0, MatrixType const & B, MatrixType const & Binv, VectorType const & r){
		setupFromData(x, C, r0, B, r);
		this->sanityCheck("::__construct__");
	}
	/** setup this set with a given data but the set does not own the r0 (user is responsible for deleting it) */
	BasicDoubleton(VectorType const & x, MatrixType const & C, VectorType* r0, MatrixType const & B, VectorType const & r){
		if (r0)
			setupFromData(x, C, *r0, B, r);
		else
			setupFromData(x, C, VectorType(), B, r);
		this->sanityCheck("::__construct__");
	}
	/** setup this set with a given data but the set does not own the r0 (user is responsible for deleting it) */
	BasicDoubleton(VectorType const & x, MatrixType const & C, VectorType* r0, MatrixType const & B, MatrixType const & Binv, VectorType const & r){
		if (r0)
			setupFromData(x, C, *r0, B, r);
		else
			setupFromData(x, C, VectorType(), B, r);
		this->sanityCheck("::__construct__");
	}
	/** setup this set with a given data, set is owner of everything */
	BasicDoubleton(VectorType const & x, MatrixType const & C, VectorType const & r0){
		MatrixType B(x.dimension(), x.dimension()); B.setToIdentity();
		setupFromData(x, C, r0, B, VectorType(x.dimension()));
		this->sanityCheck("::__construct__");
	}
	/** setup this set with a given data but the set does not own the r0 (user is responsible for deleting it) */
	BasicDoubleton(VectorType const & x, MatrixType const & C, VectorType* r0){
		MatrixType B(x.dimension(), x.dimension()); B.setToIdentity();
		if (r0)
			setupFromData(x, C, *r0, B, VectorType(x.dimension()));
		else
			setupFromData(x, C, VectorType(0), B, VectorType(x.dimension()));
		this->sanityCheck("::__construct__");
	}
	/** setup this set with a given data but the set does not own the data (user is responsible for deleting) */
	BasicDoubleton(VectorType* x, MatrixType* C, VectorType* r0, MatrixType* B, VectorType* r){
		if (!x) throw std::logic_error("BasicDoubleton::__construct__: x cannot be NULL");
		if (!C) throw std::logic_error("BasicDoubleton::__construct__: C cannot be NULL");
		MatrixType BB;
		if (B){
			BB = *B;
		} else {
			BB = MatrixType(x->dimension(), x->dimension());
			BB.setToIdentity();
		}
		if (r0)
			setupFromData(x, C, *r0, BB, r ? *r : VectorType(x->dimension()));
		else
			setupFromData(x, C, VectorType(0), BB, r ? *r : VectorType(x->dimension()));
		this->sanityCheck("::__construct__");
	}
	/** setup this set as zero vector, but the structure of given dimensions. If second arg is < 0 then d is used instead. */
	BasicDoubleton(size_type d, size_type N0 = -1){
		if (N0 < 0) N0 = d;
		setupDimension(d, N0); // sets m_B to Id
	}
	/** standard thing */
	virtual ~BasicDoubleton(){ }

	using BaseClass::dimension; ///< import dimension from interface
	/** see interface */
	size_type storageN0() const { return m_r0.dimension(); }
	/** see interface */
	size_type storageDimension() const { return m_x.dimension(); }

	/** see interface */
	VectorType get_x() const { return this->m_x; }
	/** see interface */
	MatrixType get_C() const { return this->m_C; }
	/** see interface */
	VectorType get_r0() const { return this->m_r0; }
	/** see interface */
	MatrixType get_B() const { return this->m_B; }
	/** see interface */
	VectorType get_r() const { return this->m_r; }
	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
	BaseClass& set_x(VectorType const &x) { m_x = x; sanityCheck("::set_x"); return *this; }
	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
	BaseClass& set_C(MatrixType const &C) { m_C = C; sanityCheck("::set_C"); return *this; }
	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
	BaseClass& set_r0(VectorType const &r0) { m_r0 = r0; sanityCheck("::set_r0"); return *this; }
	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
	BaseClass& set_B(MatrixType const &B) { m_B = B; sanityCheck("::set_B"); return *this; }
	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
	BaseClass& set_r(VectorType const &r) { m_r = r; sanityCheck("::set_r"); return *this; }
	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
	BaseClass& set_Cr0(MatrixType const &C, VectorType const &r0) {
		m_C = C; m_r0 = r0;
		sanityCheck("::set_Cr0");
		return *this;
	}
	// TODO: (NOT URGENT) implement Binv as a member...
//	/** Note: the element must be compatible with other elements (dimensions!) or we throw exception! */
//	BaseClass& set_Binv(MatrixType const &Binv) { this->m_Binv = Binv; return *this; }
	/** see interface */
	VectorType midPoint() const { return m_x; }
	/** see interface */
	VectorType hull() const { return m_x + m_C * m_r0 + m_B * m_r; }
	using BaseClass::set_x; ///< using the base class set_*(pointer, ownershhip) method...
	using BaseClass::set_C; ///< using the base class set_*(pointer, ownershhip) method...
	using BaseClass::set_r0; ///< using the base class set_*(pointer, ownershhip) method...
	using BaseClass::set_B; ///< using the base class set_*(pointer, ownershhip) method...
	using BaseClass::set_r; ///< using the base class set_*(pointer, ownershhip) method...
	using BaseClass::set_Cr0; ///< using the base class set_*(pointer, ownershhip) method...
	using BaseClass::set_Binv; ///< using the base class set_*(pointer, ownershhip) method...
	using BaseClass::take_x; ///< using the base class take_* method...
	using BaseClass::take_C; ///< using the base class take_* method...
	using BaseClass::take_r0; ///< using the base class take_* method...
	using BaseClass::take_B; ///< using the base class take_* method...
	using BaseClass::take_r; ///< using the base class take_* method...
	using BaseClass::take_Binv; ///< using the base class take_* method...

	/** checked by object equality */
	virtual bool common_x(VectorType const* x) const { return this->m_x == *x; }
	/** checked by object equality */
	virtual bool common_C(MatrixType const* C) const { return this->m_C == *C; }
	/** checked by object equality */
	virtual bool common_r0(VectorType const* r0) const { return this->m_r0 == *r0; }
	/** checked by object equality */
	virtual bool common_B(MatrixType const* B) const { return this->m_B == *B; }
	/** checked by object equality */
	virtual bool common_r(VectorType const* r) const { return this->m_r == *r; }

	/** add a Set or Vector to this representation. It takes the representation into consideration. */
	using BaseClass::add;
	/** multiply set by a scalar. Should take set structure into consideration. */
	using BaseClass::mul;
	/** multiply set by a scalar, then add another set. Very simple basic implementation by first mul then add. */
	using BaseClass::mulThenAdd;

	/** applies in a smart way to this set X the affine transform f(y) = M * (X - v) */
	BaseClass& affineTransform(MatrixType const &M, VectorType const &v){
		if (M.numberOfRows() != dimension() && M.numberOfColumns() != dimension())
			throw std::logic_error("BasicDoubleton::affineTransform: requires square matrix of same dimension as the set");
		if (v.dimension() != dimension())
			throw std::logic_error("BasicDoubleton::affineTransform: requires v to be of the same dimension as the set");
		m_x = M * (m_x  - v);
		m_B = M * m_B;
		if (storageN0()) m_C = M * m_C;
		// TODO: (!!!URGENT) split m_x and incorporate rest to m_B * m_r?
		// TODO: (!!!URGENT) split m_C and incorporate rest * m_r0 to m_B * m_r?
		// TODO: (!!!URGENT) return m_B to a good shape? (QR strategy)
		return *this;
	}
	/** applies in a smart way to this set X the transform f(y) = X - v */
	BaseClass& translate(VectorType const &v){
		if (v.dimension() != dimension())
			throw std::logic_error("BasicDoubleton::translate: requires v to be of the same dimension as the set");
		m_x += v; return *this;
	}

	/** show info on this set */
	std::string show() const {
		std::ostringstream oss;
		oss << "BasicDoubleton(x + C * r0 + B * r):\nDimension / N_0: ";
		oss << (*this); return oss.str();
	}

	/** see base class */
	virtual void reinitialize(size_type d, size_type N0){
		setupDimension(d, N0);
	}

protected:
	VectorType m_x;
	MatrixType m_C;
	VectorType m_r0;
	MatrixType m_B;
	VectorType m_r;

	/**
	 * helper function to check if the object represents sane data
	 * (all dimensions compatible and all pointers set if neccesary )
	 * should be called after constructing or manipulating object
	 */
	void sanityCheck(std::string const& what = "") const {
		if (storageN0() != m_C.numberOfColumns()) {
			std::ostringstream extraInfo; extraInfo << storageN0() << " vs " << m_C.numberOfColumns();
			throw std::domain_error("BasicDoubleton"+what+": C and r0 have incompatible dimensions: " + extraInfo.str());
		}
		if (dimension() != m_B.numberOfRows() || dimension() != m_B.numberOfColumns()) {
			std::ostringstream info;
			info << "BasicDoubleton" << what << ": B has incompatible dimensions.";
			info << "was " << dimension() << ", expected " << m_B.numberOfRows() << ".";
			throw std::domain_error(info.str());
		}
		if (dimension() != m_C.numberOfRows()) {
			std::ostringstream extraInfo; extraInfo << dimension() << " vs " << m_C.numberOfRows();
			throw std::domain_error("BasicDoubleton"+what+": C and x has incompatible dimensions: " + extraInfo.str());
		}
		if (dimension() != m_r.dimension()) {
			std::ostringstream extraInfo; extraInfo << dimension() << " vs " << m_r.dimension();
			throw std::domain_error("BasicDoubleton"+what+": r has incompatible dimensions." + extraInfo.str());
		}
	}
	/** makes all vector and matrices to be of appropriate dimension */
	void setupDimension(size_type d, size_type N0 = 0){
		m_x = VectorType(d);
		m_C = MatrixType(d, N0);
		m_r0 = VectorType(N0);
		m_B = MatrixType(d, d);
		m_B.setToIdentity();
		m_r = VectorType(d);
	}
	/** helper in setup */
	void setupFromData(
			VectorType  const& x,
			MatrixType  const& C,
			VectorType  const& r0,
			MatrixType  const& B,
			VectorType  const& r)
	{
		m_x = x;
		m_C = C;
		m_r0 = r0;
		m_B = B;
		m_r = r;
	}
	/** helper in copying */
	void setupFromOther(Class const& other){
		setupFromData(other.m_x, other.m_C, other.m_r0, other.m_B, other.m_r);
	}
	/** helper in copying from base class (interface) */
	void setupFromOther(BaseClass const& other){
		setupFromData(other.get_x(), other.get_C(), other.get_r0(), other.get_B(), other.get_r());
	}
};

} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_BASICDOUBLETON_H_ */
