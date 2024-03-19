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

#ifndef _CAPD_DDES_DOUBLETONINTERFACE_H_
#define _CAPD_DDES_DOUBLETONINTERFACE_H_

#include <capd/ddes/DDECommon.h>
#include <bitset>

namespace capd{
namespace ddes{

/**
 * interface I use in my code to store doubleton data of the form:
 *
 *    X = x + C * r0 + B * r
 *
 * I use storage* to name methods, as for example dimension() is usually related
 * to the algebraic dimension of the space in CAPD. In my case I will have
 * doubletons that represent data by a huge number of coefficients, but in small dimensional
 * space, i.e. piecewise Taylor curves of order n in d -dimensional space.
 * so the storageDimension will be (n+1) * d;
 *
 * NOTE: This is quite long and tedious code in C++ because of C++
 *       This should be heavily tested for correctedness and memory leaks
 *       Reader (e.g. reviewer of the manuscript for publication) should not worry to check that code
 *
 * NOTE: THIS IS JUST AN INTERFACE for a doubleton set. We provide two specific versions:
 *       Basic* and Shared* - first for testing, second for applications.
 *
 * TODO: (NOT URGENT) move big implementation into .hpp part.
 * TODO: (FAR FUTURE, RETHINK) make set_*(ptr*) methods return bool if either they truly accepted to hold the ownership and the ptr and it is really shared. Then user would be able to take actions accordingly. Now, it is not allways so obvious.
 */
template<typename MatrixSpec, typename PoliciesSpec = capd::dynset::IdQRPolicy>
class DoubletonInterface : public PoliciesSpec {
public:
	typedef MatrixSpec MatrixType;
	typedef typename MatrixType::RowVectorType VectorType;
	typedef typename MatrixType::ScalarType ScalarType;
	typedef typename MatrixType::size_type size_type;
	typedef DoubletonInterface Class;
	typedef DoubletonInterface BaseClass;
	typedef PoliciesSpec Policy;
	typedef PoliciesSpec QRPolicy;

	static std::string badge() { return "DoubletonInterface"; }

	/** reinitializes the Doubleton to a gibven structure (all should be 0 or Id for B, Binv)  */
	virtual void reinitialize(size_type d, size_type N0) = 0;

	/** returns a dimension of the stored data. Needs to be implemented.  */
	virtual size_type storageDimension() const = 0;
	/** by default dimension is equal to storage dimension */
	virtual size_type dimension() const { return storageDimension();  };
	/** returns a dimension of the r0 vector. Needs to be implemented.  */
	virtual size_type storageN0() const = 0;
	/** makes a vector that can store x part */
	virtual VectorType makeStorage_x() const { return VectorType(storageDimension()); }
	/** makes a vector that can store C part */
	virtual MatrixType makeStorage_C() const { return MatrixType(storageDimension(), storageN0()); }
	/** makes a vector that can store r0 part */
	virtual VectorType makeStorage_r0() const { return VectorType(storageN0()); }
	/** makes a vector that can store B part */
	virtual MatrixType makeStorage_B() const { return MatrixType(storageDimension(), storageDimension()); }
	/** makes a vector that can store r part */
	virtual VectorType makeStorage_r() const { return VectorType(storageDimension()); }
	/** returns mid point of the set, default: returns get_X(); */
	virtual VectorType midPoint() const { return get_x(); };
	/** returns mid point of the set, default: returns value obtained using get_*() methods. You might want to reimplement for better performance(?) */
	virtual VectorType hull() const { return get_x() + get_C() * get_r0() + get_B() * get_r(); };
	/** returns Vector representation of the set, default: returns hull() */
	virtual operator VectorType() const { return hull(); };

	/** returns x value of the doubleton. Needs to be implemented. */
	virtual VectorType get_x() const = 0;
	/** returns B value of the doubleton. Needs to be implemented. */
	virtual MatrixType get_C() const = 0;
	/** returns B value of the doubleton. Needs to be implemented. */
	virtual MatrixType get_B() const = 0;
	/** returns inverse of B value of the doubleton. Default implementation by capd inverse. */
	virtual MatrixType get_Binv() const { return capd::matrixAlgorithms::inverseMatrix(get_B()); }
	/** returns r value of the doubleton. Needs to be implemented. */
	virtual VectorType get_r() const = 0;
	/** returns r0 value of the doubleton. Needs to be implemented. */
	virtual VectorType get_r0() const = 0;
	/** Must assure that arg is compatible with other parts (dimensions!) or throw exception! Should return *this for a cascade. */
	virtual Class& set_x(VectorType const &x) = 0;
	/** Must assure that arg is compatible with other parts (dimensions!) or throw exception! Should return *this for a cascade. */
	virtual Class& set_C(MatrixType const &C) = 0;
	/** Must assure that arg is compatible with other parts (dimensions!) or throw exception! Should return *this for a cascade. */
	virtual Class& set_r0(VectorType const &r0) = 0;
	/** Must assure that arg is compatible with other parts (dimensions!) or throw exception! Should return *this for a cascade. */
	virtual Class& set_B(MatrixType const &B) = 0;
	/** Must assure that arg is compatible with other parts (dimensions!) or throw exception! Should return *this for a cascade. */
	virtual Class& set_r(VectorType const &r) = 0;
	/**
	 * This is to allow to safely change the structure of the set. If changing r0 dimension,
	 * then C matrix must be updated to match the dimension, and this must be an atomic operation.
	 * Must assure that arg is compatible with other parts (dimensions!) or throw exception! Should return *this for a cascade.
	 */
	virtual Class& set_Cr0(MatrixType const &C, VectorType const &r0) = 0;
	/**
	 * A version of the interface that allows for sharing of the pointers from external source.
	 * passOwnership is a flag to indicate if the class should take resposibility of freeing
	 * the memory allocation of the precceding pointer. A default implementation provided for
	 * convenience uses the corresponding set_* method and then frees the memory of
	 * the argument if passOwnership = true.
	 * Must assure that args are compatible with other parts (dimensions!) or throw exception!
	 * Should return *this for a cascade.
	 */
	virtual Class& set_x(VectorType* x, bool passOwnership = false){
		try{
			callSetMethodForPtr(x, passOwnership, &Class::set_x);
		} catch (std::logic_error& e){
			throw std::logic_error(std::string("DoubletonInterface::set_x(ptr, owner): ") + e.what());
		}
		return *this;
	}
	/**
	 * A version of the interface that allows for sharing of the pointers from external source.
	 * passOwnership is a flag to indicate if the class should take resposibility of freeing
	 * the memory allocation of the precceding pointer. A default implementation provided for
	 * convenience uses the corresponding set_* method and then frees the memory of
	 * the argument if passOwnership = true.
	 * Must assure that args are compatible with other parts (dimensions!) or throw exception!
	 * Should return *this for a cascade.
	 */
	virtual Class& set_C(MatrixType* C, bool passOwnership = false){
		try{
			callSetMethodForPtr(C, passOwnership, &Class::set_C);
		} catch (std::logic_error& e){
			throw std::logic_error(std::string("DoubletonInterface::set_C(ptr, owner): ") + e.what());
		}
		return *this;
	}
	/**
	 * A version of the interface that allows for sharing of the pointers from external source.
	 * passOwnership is a flag to indicate if the class should take resposibility of freeing
	 * the memory allocation of the precceding pointer. A default implementation provided for
	 * convenience uses the corresponding set_* method and then frees the memory of
	 * the argument if passOwnership = true.
	 * Must assure that args are compatible with other parts (dimensions!) or throw exception!
	 * Should return *this for a cascade.
	 */
	virtual Class& set_r0(VectorType* r0, bool passOwnership = false){
		try{
			callSetMethodForPtr(r0, passOwnership, &Class::set_r0);
		} catch (std::logic_error& e){
			throw std::logic_error(std::string("DoubletonInterface::set_r0(ptr, owner): ") + e.what());
		}
		return *this;
	}
	/**
	 * A version of the interface that allows for sharing of the pointers from external source.
	 * passOwnership is a flag to indicate if the class should take resposibility of freeing
	 * the memory allocation of the precceding pointer. A default implementation provided for
	 * convenience uses the corresponding set_* method and then frees the memory of
	 * the argument if passOwnership = true.
	 * Must assure that args are compatible with other parts (dimensions!) or throw exception!
	 * Should return *this for a cascade.
	 */
	virtual Class& set_B(MatrixType* B, bool passOwnership = false){
		try{
			callSetMethodForPtr(B, passOwnership, &Class::set_B);
		} catch (std::logic_error& e){
			throw std::logic_error(std::string("DoubletonInterface::set_B(ptr, owner): ") + e.what());
		}
		return *this;
	}
	/**
	 * A version of the interface that allows for sharing of the pointers from external source.
	 * passOwnership is a flag to indicate if the class should take resposibility of freeing
	 * the memory allocation of the precceding pointer. A default implementation provided for
	 * convenience uses the corresponding set_* method and then frees the memory of
	 * the argument if passOwnership = true.
	 * Must assure that args are compatible with other parts (dimensions!) or throw exception!
	 * Should return *this for a cascade.
	 */
	virtual Class& set_r(VectorType* r, bool passOwnership = false){
		try{
			callSetMethodForPtr(r, passOwnership, &Class::set_r);
		} catch (std::logic_error& e){
			throw std::logic_error(std::string("DoubletonInterface::set_r(ptr, owner): ") + e.what());
		}
		return *this;
	}
	/**
	 * A version of the interface that allows for sharing of the pointers from external source.
	 * passOwnership is a flag to indicate if the class should take resposibility of freeing
	 * the memory allocation of the precceding pointer. A default implementation provided for
	 * convenience uses the corresponding set_* method and then frees the memory of
	 * the argument if passOwnership = true.
	 * Must assure that args are compatible with other parts (dimensions!) or throw exception!
	 * Should return *this for a cascade.
	 */
	virtual Class& set_Cr0(MatrixType* C, VectorType* r0, bool passOwnership1 = false, bool passOwnership2 = false){
		try{
			if (!C) throw std::logic_error("C argument is NULL");
			if (!r0) throw std::logic_error("r0 argument is NULL");
			try { set_Cr0(*C, *r0); } catch (...){ throw std::logic_error("incompatible dimensions"); }
		} catch (std::logic_error& e){
			helper_safe_delete(C, passOwnership1);
			helper_safe_delete(r0, passOwnership2);
			throw e;
		}
		helper_safe_delete(C, passOwnership1);
		helper_safe_delete(r0, passOwnership2);
		return *this;
	}
	/**
	 * A version of the interface that allows for sharing of the pointers from external source.
	 * passOwnership is a flag to indicate if the class should take resposibility of freeing
	 * the memory allocation of the precceding pointer. A default implementation provided for
	 * convenience uses the corresponding set_* method and then frees the memory of
	 * the argument if passOwnership = true.
	 * Must assure that args are compatible with other parts (dimensions!) or throw exception!
	 * Should return *this for a cascade.
	 */
	virtual Class& set_Binv(MatrixType* Binv, bool passOwnership = false){
		try{
			callSetMethodForPtr(Binv, passOwnership, &Class::set_Binv);
		} catch (std::logic_error& e){
			throw std::logic_error(std::string("DoubletonInterface::set_Binv(ptr, owner): ") + e.what());
		}
		return *this;
	}

	/**
	 * for future implementations. By default it does nothing - Binv is computed always by inverse matrix.
	 * NOTE: user of this function is responsible for providing TRUE inverse. Code does not check anything in this regard.
	 */
	virtual Class& set_Binv(MatrixType const &Binv) { /* purposefully do nothing */ return *this; };

	/** tests if the corresponding part of the set is the same entity as the parameter */
	virtual bool common_x(VectorType const* x) const { return false; }
	/** tests if the corresponding part of the set is the same entity as the parameter */
	virtual bool common_C(MatrixType const* C) const { return false; }
	/** tests if the corresponding part of the set is the same entity as the parameter */
	virtual bool common_r0(VectorType const* r0) const { return false; }
	/** tests if the corresponding part of the set is the same entity as the parameter */
	virtual bool common_B(MatrixType const* B) const { return false; }
	/** tests if the corresponding part of the set is the same entity as the parameter */
	virtual bool common_r(VectorType const* r) const { return false; }

	/** computes v . thisSet, it is virtual to allow optimization. Now uses get_* in computation (might be slow) */
	virtual ScalarType dot(VectorType const& v) const {
		MatrixType CT = transpose(get_C());
		MatrixType BT = transpose(get_B());
		return (v * get_x()) + ((CT * v) * get_r0()) + ((BT * v) * get_r());
	}

	/** add vector to this representation. It is distributed between x and B*r part by default */
	virtual Class& add(VectorType const & v){
		size_type d = this->dimension();
		if (d != v.dimension()) throw std::logic_error("DoubletonInterface::add(vector): Input vector has incompatible dimension.");
		VectorType new_x = get_x() + v, s(d);
		capd::vectalg::split(new_x, s);
		this->set_x(new_x);
		this->set_r(get_r() + get_Binv() * s);
		return *this;
	}
	/** add a Set to this representation. It takes the representation into consideration. */
	virtual Class& add(Class const & set){
		size_type N0 = this->storageN0();
		size_type d = this->dimension();
		if (set.storageN0() != N0) throw std::logic_error("DoubletonInterface::add(set): Sets has incompatible N0 dimensions.");
		if (set.dimension() != d) throw std::logic_error("DoubletonInterface::add(set): Sets has incompatible dimensions.");
		VectorType other_r0 = set.get_r0();
		VectorType this_r0 = this->get_r0();
		if (!this->common_r0(&other_r0))
			capd::vectalg::intervalHull(this_r0, other_r0, this_r0);
		this->set_r0(this_r0);

		MatrixType new_C = this->get_C() + set.get_C(), S(d, N0);
		capd::vectalg::split(new_C, S);
		this->set_C(new_C);

		MatrixType invB = this->get_Binv();
		VectorType new_x = get_x() + set.get_x(), s(d);
		this->set_x(new_x);
		capd::vectalg::split(new_x, s);
		VectorType new_r = this->get_r() + (invB * set.get_B()) * set.get_r() + invB * (s + S * this_r0);
		this->set_r(new_r);

		return *this;
	}
	/** multiply set by a scalar. Should take set structure into consideration. */
	virtual Class& mul(ScalarType const & c){
		size_type d = this->dimension();
		size_type N0 = this->storageN0();
		MatrixType invB = this->get_Binv();
		VectorType new_x = c * get_x(), s(d);
		capd::vectalg::split(new_x, s);
		MatrixType new_C = c * get_C(), S(d, N0);
		capd::vectalg::split(new_C, S);
		VectorType new_r = c * get_r() + invB * (S * this->get_r0() + s);
		this->set_x(new_x);
		this->set_C(new_C);
		this->set_r(new_r);
		// TODO: if c == 0. it would be good to reset B and invB for better performance later
		return *this;
	}
	/** multiply set by a scalar, then add another set. Very simple basic implementation by first mul then add. */
	virtual Class& mulThenAdd(ScalarType const & c, Class const & set){ this->mul(c); this->add(set); return *this; }
	/** applies in a smart way to this set X the affine transform f(y) = M * (X - v). Needs to be implemented.  */
	virtual Class& affineTransform(MatrixType const &M, VectorType const &v) = 0;
	/** applies in a smart way to this set X the transform f(y) = X - v. Needs to be implemented.  */
	virtual Class& translate(VectorType const &v) = 0;

	/** unmark itself as owner and returns the pointer. From now on user is responsible for deleting. It provides basic implementation as returning copy of the data by get_* */
	virtual VectorType* take_x() { return new VectorType(get_x()); }
	/** unmark itself as owner and returns the pointer. From now on user is responsible for deleting. It provides basic implementation as returning copy of the data by get_* */
	virtual MatrixType* take_C() { return new MatrixType(get_C()); }
	/** unmark itself as owner and returns the pointer. From now on user is responsible for deleting. It provides basic implementation as returning copy of the data by get_* */
	virtual VectorType* take_r0() { return new VectorType(get_r0()); }
	/** unmark itself as owner and returns the pointer. From now on user is responsible for deleting. It provides basic implementation as returning copy of the data by get_* */
	virtual MatrixType* take_B() { return new MatrixType(get_B()); }
	/** unmark itself as owner and returns the pointer. From now on user is responsible for deleting. It provides basic implementation as returning copy of the data by get_* */
	virtual VectorType* take_r() { return new VectorType(get_r()); }
	/** unmark itself as owner and returns the pointer. From now on user is responsible for deleting. It provides basic implementation as returning copy of the data by get_* */
	virtual MatrixType* take_Binv() { return new MatrixType(get_Binv()); }

	/** virtual destructor for strict warning compliance */
	virtual ~DoubletonInterface() {};

	/** Should be good for any Doubleton that inherits */
	friend std::ostream& operator<<(std::ostream & out, Class const & c) {
		out << Class::badge() << std::endl;
		out << c.get_x() << std::endl;
		out << c.get_C() << std::endl;
		out << c.get_r0() << std::endl;
		out << c.get_B() << std::endl;
		out << c.get_Binv() << std::endl;
		out << c.get_r();
		return out;
	}

	/** Should be good for any Doubleton that inherits */
	friend std::istream& operator>>(std::istream & in, Class & c) {
		helper_dump_badge(in);
		MatrixType B, Binv, C;
		VectorType x, r0, r;
		in >> x;
		in >> C;
		in >> r0;
		in >> B;
		in >> Binv;
		in >> r;
		c.reinitialize(x.dimension(), r0.dimension());
		c.set_x(x);
		c.set_Cr0(C, r0);
		c.set_B(B);
		c.set_Binv(Binv);
		c.set_r(r);
		return in;
	}
	/** alternative to mul(); uses mul to gain polymorphism */
	DoubletonInterface& operator*=(ScalarType const& c){ return this->mul(c); }
private:
	/**
	 * This is very private function to have a common implementation
	 * of almost the same code. It could be achieved with #define, but I have
	 * decided to go with the templates and the pointers to members.
	 */
	template<typename TypeSpec>
	void callSetMethodForPtr(TypeSpec* item, bool passOwnership, Class& (Class::* call)(TypeSpec const &)){
		try{
			if (!item) throw std::logic_error("argument is NULL");
			try { (this->*call)(*item); } catch (...){ throw std::logic_error("incompatible dimensions"); }
		} catch (std::logic_error& e){
			helper_safe_delete(item, passOwnership);
			throw e;
		}
		helper_safe_delete(item, passOwnership);
	}
};

} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_DOUBLETONINTERFACE_H_ */
