/*
 * DDESet.h
 *
 * This file constitutes part of DDEs rigorous integration framework developed
 * in PhD Thesis:
 *
 * 		"Rigorous Integration of Delay Differential Equations", Jagiellonian University, 2015
 *
 * When using in scientific work please consider citing my articles (preffered),
 * PhD thesis and/or my webpage. For the publications describing the source code
 * please refer to http://scirsc.org/p/papers.
 *
 * This work would not be possible without aid and expertise of people involved in
 * CAPD developing library (Computer Assisted Proofs in Dynamics).
 * Please refer to http://capd.ii.uj.edu.pl and consider citing also this library when
 * using those codes in any scientific work.
 *
 * Author: Robert Szczelina, PhD
 * Małopolska Center of Biotechnology, Jagiellonian University AND
 * (former) Faculty of Mathematics and Computer Science, Jagiellonian University
 * email: 	robert.szczelina@uj.edu.pl
 * www: 	scirsc.org
 *
 * This source code is provided under GNU GPL license (v.2 or whatever)
 */

#ifndef _CAPD_HOOKEDSET_H_
#define _CAPD_HOOKEDSET_H_

#include <iostream>
#include <vector>

namespace capd{
namespace dynset{


/**
 * remembers sets in a public
 * std::vector history. You can access
 * this vector after iteration.
 */
template<typename SetSpec>
class RememberHook : public SetSpec::HookType{
public:
	std::vector<SetSpec> history;
	RememberHook() {}
	virtual void operator()(SetSpec& the_set){
		history.push_back(the_set);
	}
	virtual ~RememberHook(){}

};  // class RememberHook ends here


/**
 * adapter for the standard sets from CAPD to the representation
 */
template<typename BaseSetSpec>
class HookedSet : public BaseSetSpec{
public:
	typedef BaseSetSpec BaseSetType;
	typedef BaseSetType super;
	typedef typename BaseSetType::DynSysType DynSysType;

	typedef typename super::MatrixType MatrixType;
	typedef typename MatrixType::RowVectorType VectorType;
	typedef typename MatrixType::ScalarType ScalarType;
	typedef typename super::NormType NormType;

	class HookType{public: typedef HookedSet HookedSetType; virtual void operator()(HookedSetType& item) = 0;} *after_move_hook, *before_move_hook;
	void setBeforeMoveHook(HookType* hook) { before_move_hook = hook;}
	void setAfterMoveHook(HookType* hook) {	after_move_hook = hook; }

	explicit HookedSet(const BaseSetType& orig):
			super(orig),
			after_move_hook(0),
			before_move_hook(0) {
	}

	/** konstruktor kopiujacy */
	HookedSet(const HookedSet& original):
			super(original),
			after_move_hook(original.after_move_hook),
			before_move_hook(original.before_move_hook){
	}

	HookedSet& operator=(const HookedSet& that){
		this->super::operator=((BaseSetType)that);
		this->after_move_hook = that.after_move_hook;
		this->before_move_hook = that.before_move_hook;
		return *this;
	}

	HookedSet& operator=(const BaseSetType& that){
		this->super::operator=(that);
		this->after_move_hook = 0;
		this->before_move_hook = 0;
		return *this;
	}

	/** we can simply delegate it to the base class */
	operator VectorType() const{
		return (VectorType)((super)*this);
	}

	/** we can delegate it to the base class, additionaly before and after we execute hooks */
	void move(DynSysType & dynsys){
		if (before_move_hook) (*before_move_hook)(*this);
		super::move(dynsys);
		if (after_move_hook) (*after_move_hook)(*this);
	}
	/** we can delegate it to the base class, additionaly before and after we execute hooks */
	void move(DynSysType & dynsys, BaseSetType& result) const {
		if (before_move_hook) (*before_move_hook)(*this);
		super::move(dynsys, result);
		if (after_move_hook) (*after_move_hook)(result);
	}
	/** we can delegate it to the base class, additionaly before and after we execute hooks */
	void move(DynSysType & dynsys, HookedSet& result) const {
		if (before_move_hook) (*before_move_hook)(*this);
		super::move(dynsys, result);
		if (after_move_hook) (*after_move_hook)(result);
	}

	/** prints the representation as a Matrix */
	friend std::ostream& operator<<(std::ostream& out, HookedSet& item){
		HookedSet::VectorType v_data(item);
		for (int k = 0; k <= item.m_index.N + 1; k++){
			for (int i = item.m_index.P; i > 0 ; i--){
				out << v_data(item.m_index.index(i, k)) << " ";
			}
			out << "\n";
		}
		return out;
	}

}; // class DDESet ends here


} //namespace dynset
} //namespace capd


#endif /* _CAPD_DYNSET_DDESET_H_ */
