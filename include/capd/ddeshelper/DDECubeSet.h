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

#ifndef _CAPD_DDECUBESET_H_
#define _CAPD_DDECUBESET_H_

#include <capd/ddes/ddeslib.h>
#include <set>
#include <vector>
#include <algorithm>

namespace capd {
namespace ddeshelper {

template<typename IType>
class GeneratorRange {
	void commonSetup(){
		index = 0;
		if (start > end){ std::swap(start, end); }
	}
public:
	GeneratorRange(IType lo=0, IType up=0, IType dt=1):
			start(lo),
			end(up),
			inc(dt) {
		commonSetup();
	}

	IType start;
	IType end;
	IType inc;
	int index;

	int i(){ return index; }
	IType v(){ return start + inc * index; }
	IType lo(){ return start; }
	IType up(){ return end; }

	void next(){
		if (!isFinished()) ++index;
	}
	bool isFinished(){
		return v() >= end;
	}
	void reset(){
		index = 0;
	}
	bool isTrivial(){ return start == end; }

	template<typename iterator>
	void fill(iterator from, iterator to){
		for (iterator p = from; p < to; p++)
			*p = *this;
	}

	/**
	 * returns true if the incrementation was successful and the next element is available
	 * returns false if the incrementation reached last element and reseted all elements (we are back to all indices 0)
	 * */
	template<typename iterator>
	static bool increment(iterator from, iterator to){
		iterator p = from;
		while (p < to && p->isFinished()){
			p->reset();
			++p;
		}
		if (p == to) return false;
		p->next();
		return true;
	}

	friend std::ostream& operator<<(std::ostream& out, GeneratorRange const& range){
		out << "GeneratiorRange(" << range.start << ", " << range.end << "," << range.inc << ")";
		return out;
	}
};

template<typename VectorSpec>
class CubeSet {
public:
	typedef VectorSpec VectorType;
	typedef typename VectorType::ScalarType ScalarType;
	typedef typename VectorType::size_type size_type;
	typedef std::vector<int> CubeSignatureType;
	typedef std::set<CubeSignatureType> CSContainerType;
	typedef CSContainerType::const_iterator CubeSignatureConstIterator;
	typedef CSContainerType::iterator CubeSignatureIterator;

	class CubeSetIterator{
	public:
		CubeSetIterator(CubeSetIterator const& orig): m_ptr(orig.m_ptr), m_owner(orig.m_owner) {}
		operator VectorType() const { return m_owner.get_dx(*m_ptr); }
		CubeSignatureType const& operator*() const { return *m_ptr; }
		CubeSetIterator& operator++(){ ++m_ptr; return *this; }
		CubeSetIterator operator++(int){ auto cpy = *this; ++m_ptr; return cpy; }
		friend bool operator==(CubeSetIterator const& lhs, CubeSetIterator const& rhs) {
			return (&lhs.m_owner == &rhs.m_owner) && (lhs.m_ptr == rhs.m_ptr);
		}
		friend bool operator!=(CubeSetIterator const& lhs, CubeSetIterator const& rhs) {
			return !(lhs == rhs);
		}
	private:
		CubeSignatureConstIterator m_ptr;
		CubeSet const& m_owner;
		CubeSetIterator(CubeSet& owner, CubeSignatureConstIterator&& ptr): m_ptr(ptr), m_owner(owner) {}
	friend class CubeSet;
	};

	CubeSet() {}
	CubeSet(VectorType const& x0, VectorType const& r0): m_x0(x0), m_r0(r0) {
		assertDimension(m_r0, "In constructor(x0, r0)");
		symetrize();
	}
	CubeSet(VectorType const& r0): m_x0(r0.dimension()), m_r0(r0) { symetrize(); }

	CubeSet& insert(CubeSignatureType const& cube){
		// assertDimension(std::max<size_type>(int(cube.size(), m_r0.dimension()), "In insert(cube), cube dimension is too big;");
		m_cubes.insert(cube); return *this;
	}
	/**
	 * Max_cuts == -1 means we mean cut in all neccessary dimensions. This is default behaviour.
	 * returns 	the last index at what the cut was needed. If it is smaller than max_cuts requirement
	 * 			then we are ok with the set. If it is greater on equal, then we need to probably
	 * 			change the m_r0, as the box cannot be covered with the current small boxes.
	 */
	int insert_cover(VectorType box, int max_cuts=-1){
		assertDimension(box, "In insert_cover(box), box dimension is wrong;");
		std::vector<GeneratorRange<int>> series;
		int last_cut_needed = -1;
		auto done = [max_cuts, &series](){ return (max_cuts > -1 && series.size() >= max_cuts); };
		for (int i = 0; i < m_r0.dimension(); ++i){
			auto r = m_r0[i].rightBound(); // we assume m_r0 is symmetric!
			if (r == 0. && !done()){
				series.push_back(GeneratorRange<int>(0, 0, 1));
			}else{
				auto diff = box[i] - m_x0[i];
				// let diff = [left, right]
				// we looking for i \in \Z such that -r -i * r < left,
				// so we get i < d, with d = (left/r+1)/2, analogously
				// for right we get i > d = (right/r-1)/2
				// in case d \in \Z, we need to include extra piece, just in case
				auto lod = (((diff / r) + 1) / 2.0).leftBound();
				auto upd = (((diff / r) - 1) / 2.0).rightBound();
				int loi = int(floor(lod)); if (lod == loi) --loi;
				int upi = int(ceil(upd));  if (upd == upi) ++upi;

				if (loi != upi)
					last_cut_needed = i;
				if (!done())
					series.push_back(GeneratorRange<int>(loi, upi, 1));

				// TODO: remove debug after tests!
//				if (i == 1){
//					std::cout << "box[i]  = " << box[i] << std::endl;
//					std::cout << "x0[i]   = " << m_x0[i] << std::endl;
//					std::cout << "diff    = " << diff << std::endl;
//					std::cout << "m_r0[i] = " << m_r0[i] << std::endl;
//					std::cout << "loi     = " << loi << " " << loi * 2 * m_r0[i].rightBound() + m_x0[i] + m_r0[i].leftBound() << " vs box " << box[i].leftBound() << std::endl;
//					std::cout << "upi     = " << upi << " " << upi * 2 * m_r0[i].rightBound() + m_x0[i] + m_r0[i].rightBound() << " vs box " << box[i].rightBound() << std::endl;
//				}
			}
		}
		std::cout << "All cut" << std::endl;

		do {
			CubeSignatureType cube;
			for (int i = 0; i < series.size(); ++i)
				cube.push_back(series[i].v());
			this->insert(cube);
		} while(GeneratorRange<int>::increment(series.begin(), series.end()));
		return last_cut_needed;
	}

	CubeSetIterator begin() { return CubeSetIterator(*this, csbegin()); }
	CubeSetIterator end() { return CubeSetIterator(*this, csend()); }

	CubeSignatureIterator csbegin() { return m_cubes.begin(); }
	CubeSignatureIterator csend() { return m_cubes.end(); }
	CubeSignatureConstIterator csbegin() const { return m_cubes.begin(); }
	CubeSignatureConstIterator csend() const { return m_cubes.end(); }

	VectorType get_x0() const { return m_x0; }
	VectorType get_r0() const { return m_r0; }
	VectorType get_dx(CubeSignatureType const& signature) const {
		int i = 0;
		VectorType dx = m_x0;
		for (auto& di: signature){
			dx[i] += m_r0[i].rightBound() * 2. * di;
			++i;
		}
		return dx;
	}

	bool operator==(CubeSet const& other) const {
		return this->subcubeset(other) && other.subcubeset(*this);
//		bool result = m_r0 == other.m_r0 && m_x0 == other.m_x0;
//		for (auto& p1: m_cubes)
//			result &= (other.m_cubes.find(p1) != other.m_cubes.end());
//		for (auto& p2: other.m_cubes)
//			result &= (m_cubes.find(p2) != m_cubes.end());
//		return result;
	}

	/** this only outputs cube signatures and add to the set */
	friend std::ostream& operator<<(std::ostream &out, CubeSet const& items){
		auto print_vector = [&](CubeSignatureType const& signature){
			out << "[";
			if (signature.size()){
				auto elem = signature.begin();
				out << (*elem++);
				for (; elem != signature.end(); ++elem)
					out << "," << (*elem);
			}
			out << "]";
		};
		out << "{";
		if (items.m_cubes.size()){
			out << "\n";
			auto signature = items.csbegin();
			print_vector(*signature++);
			for (; signature != items.csend(); ++signature){
				out << ",\n";
				print_vector(*signature);
			}
			out << "\n";
		}
		out << "}";
		return out;
	}
	/** this only reads cube signatures and adds to the set */
	friend std::istream& operator>>(std::istream &in, CubeSet &items){
		auto scan_for_token = [&](char expected){
			char token; in >> token;
			if (token != expected)
				in.putback(token);
			return token == expected;
		};
		if (!scan_for_token('{')) { return in; } // nothing to read
		bool done = false;
		while(!done){
			if (!scan_for_token('[')) { break; } // nothing to read, break free
			CubeSignatureType signature;
			while (true){
				if (scan_for_token(']')) {
					// it can be empty list, so we test first
					// add the current items
					items.insert(signature);
					std::cout << "Adding " << signature << std::endl;
					if (!scan_for_token(','))
						done = true; // we require lists to be separated by ,
					break; // go, look for the next signature
				}
				int item; in >> item;
				signature.push_back(item);
				scan_for_token(',');
			}
		}
		scan_for_token('}'); // remove eventual closing bracket
		return in;
	}

	/** this only check if x0's, r0's are equal (same cube structure), and that all cubes are in the other */
	bool subcubeset(CubeSet const& other) const {
		bool result = (m_x0 == other.m_x0 && m_r0 == other.m_r0);
		for (auto& c: m_cubes)
			result &= (other.m_cubes.find(c) != other.m_cubes.end());
		return result;
	}

	// TODO: ...
//	bool subset(CubeSet const& other) const {
//		for ()
//		return false; // TODO: implemenmt
//	}

	template<typename Iterator>
	CubeSet& insert(Iterator from, Iterator to){
		m_cubes.insert(from, to);
		return *this;
	}

	int size() const { return m_cubes.size(); }
private:
	VectorType m_x0;
	VectorType m_r0;
	CSContainerType m_cubes;

	void symetrize(){ m_r0 *= ScalarType(-1.,1.); } // TODO: make it better: hull = x0 + r0; x0 = mid(hull), r0 = hull - x0; ???

	void assertDimension(size_type const& dim, std::string extra_msg = "") const {
		if (dim != m_r0.dimension()){
			std::ostringstream msg;
			msg << "CubeSet::assertDimension(): " << extra_msg;
			msg << ", wrong dimension, is " << dim;
			msg << ", should be " << m_r0.dimension();
			throw std::logic_error(msg.str());
		}
	}
	void assertDimension(VectorType const& v, std::string extra_msg = "") const {
		assertDimension(v.dimension(), extra_msg);
	}
};

} // namespace ddeshelper
} // namespace capd

#endif /* _CAPD_DDECUBESET_H_ */
