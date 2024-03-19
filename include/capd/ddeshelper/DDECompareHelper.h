/*
 * DDECompareHelper.h
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

#ifndef _CAPD_DDECOMPAREHELPER_H_
#define _CAPD_DDECOMPAREHELPER_H_

#include <capd/ddes/ddeslib.h>

namespace capd {
namespace ddeshelper {

template<typename T>
std::string strfix(T const& a, typename std::string::size_type len){
	std::ostringstream oss; oss.precision(15);
	oss << a;
	std::string result = oss.str();
	while (result.length() < len) result += " ";
	if (result.length() > len) result = result.substr(0, len);
	return result;
}

enum SetRelation {
	suPset = 0, suBset = 1, miss_L = 2, miss_R = 3, int_L = 4, int_R = 5, unkn = 6,
};

template<typename M>
std::string matrixInfo(M const& A){
	std::ostringstream oss;
	oss << "row x cols: " << A.numberOfRows() << " x " << A.numberOfColumns();
	return oss.str();
}

template<typename V>
std::string vectorInfo(V const& v){
	std::ostringstream oss;
	oss << "dimension:  " << v.dimension();
	return oss.str();
}

/** checks what is relation of intervals a to b */
template<typename ScalarType>
std::pair<SetRelation, std::string> checkRelation(ScalarType const& a, ScalarType const& b){
	std::string op = "";
	SetRelation rel = SetRelation::unkn;
	if (a.contains(b)){
		op = "suPset";
		rel = SetRelation::suPset;
	} else if (b.contains(a)){
		op = "suBset";
		rel = SetRelation::suBset;
	} else if (a.leftBound() > b.rightBound()) {
		op = "miss R";
		rel = SetRelation::miss_R;
	} else if (b.leftBound() > a.rightBound()) {
		op = "miss L";
		rel = SetRelation::miss_L;
	} else if (a.contains(b.leftBound())) {
		op = "int L ";
		rel = SetRelation::int_L;
	} else if (a.contains(b.rightBound())) {
		op = "int R ";
		rel = SetRelation::int_R;
	} else {
		op = "unknown (should not happen)";
		rel = SetRelation::unkn;
	}
	return std::make_pair(rel, op);
}

template<typename RealT>
std::string printInterval(RealT v){
	std::ostringstream oss; oss.precision(15);
	oss << "[" << strfix(v.leftBound(), 25) << "," << strfix(v.rightBound(), 25) << "] (diam: " << v.rightBound() - v.leftBound() << ")";
	return oss.str();
}

/**
 * allows to easily check how the other vector relates to the original.
 */
template<typename VectorSpec>
class DDECompareHelper {
public:
	typedef VectorSpec VectorType;
	typedef typename VectorType::size_type size_type;
	typedef typename VectorType::ScalarType ScalarType;

	DDECompareHelper(VectorType const& orig, VectorType const& other): m_orig(orig), m_other(other), m_counter(7, 0) {
		for (size_type i = 0; i < m_orig.dimension(); ++i){ // TODO: use iterators
			auto relation = checkRelation(m_other[i], m_orig[i]);
			++m_counter[relation.first];
			m_relations.push_back(relation);
		}
	}

	void printRelation(int i, std::ostream &out) const {
		out << m_relations[i].second << " at " << i << std::endl;
		out << "    " << printInterval(m_other[i]) << std::endl;
		out << "    " << m_relations[i].second << std::endl;
		out << "    " << printInterval(m_orig[i]) << std::endl;
		{
			ScalarType a1, a2, r1, r2;
			m_other[i].split(a1, r1);
			m_orig[i].split(a2, r2);
			if (r1.rightBound() > r2.rightBound()){ auto tmp = r2; r2 = r1; r1 = tmp; }
			std::ostringstream margin; margin.precision(10);
			margin << (1.0 - (r1.rightBound() / r2.rightBound()));
			if (r2.rightBound() != 0.){
				out << "    (diameter margin: " << margin.str() << ")" << std::endl;
			}
		}
		out  << std::endl;
	}


	void printRelations(std::ostream &out) const {
		for (size_type i = 0; i < m_relations.size(); ++i)
			printRelation(i, out);
	}

	void printSummary(std::ostream &out) const {
		out << "suPsets  " << m_counter[SetRelation::suPset] << std::endl;
		out << "suBsets  " << m_counter[SetRelation::suBset] << std::endl;
		out << "miss L   " << m_counter[SetRelation::miss_L] << std::endl;
		out << "miss R   " << m_counter[SetRelation::miss_R] << std::endl;
		out << "int L    " << m_counter[SetRelation::int_L] << std::endl;
		out << "int R    " << m_counter[SetRelation::int_R] << std::endl;
		out << "unknowns " << m_counter[SetRelation::unkn] << std::endl;
	}

	std::string inlineSummary() const {
		std::ostringstream out;
		out << "suPs=" << m_counter[SetRelation::suPset] << ", ";
		out << "suBs=" << m_counter[SetRelation::suBset] << ", ";
		out << "miss=" << m_counter[SetRelation::miss_L] + m_counter[SetRelation::miss_R] << ", ";
		out << "ints=" << m_counter[SetRelation::int_L] + m_counter[SetRelation::int_R] << ", ";
		out << "unks=" << m_counter[SetRelation::unkn] << "; ";
		return out.str();
	}

	operator std::string() const {
		std::ostringstream oss;
		printSummary(oss);
		printRelations(oss);
		printSummary(oss);
		return oss.str();
	}

	friend std::ostream& operator<<(std::ostream& out, DDECompareHelper<VectorSpec> const& comp){
		out << std::string(comp);
		return out;
	}

	int getSubsetsCount() const { return m_counter[SetRelation::suBset]; }
	int getSupsetsCount() const { return m_counter[SetRelation::suPset]; }
	int getMissLeftCount() const { return m_counter[SetRelation::miss_L]; }
	int getMissRightCount() const { return m_counter[SetRelation::miss_R]; }
	int getMissCount() const { return getMissLeftCount() + getMissRightCount(); }
	int getIntersectionLeftCount() const { return m_counter[SetRelation::int_L]; }
	int getIntersectionRightCount() const { return m_counter[SetRelation::int_R]; }
	int getIntersectionsCount() const { return getIntersectionLeftCount() + getIntersectionRightCount(); }
	/** this should always return 0, otherwise some serious errors happen */
	int getUnknownsCount() const { return m_counter[SetRelation::unkn]; }

	bool isSubset(size_type i) const { return m_relations[i].first == SetRelation::suBset; }
	bool isSupset(size_type i) const { return m_relations[i].first == SetRelation::suPset; }
	bool isMissLeft(size_type i) const { return m_relations[i].first == SetRelation::miss_L; }
	bool isMissRight(size_type i) const { return m_relations[i].first == SetRelation::miss_R; }
	bool isMiss(size_type i) const { return isMissLeft(i) || isMissRight(i); }
	bool isIntersectLeft(size_type i) const { return m_relations[i].first == SetRelation::int_L; }
	bool isIntersectRight(size_type i) const { return m_relations[i].first == SetRelation::int_R; }
	bool isIntersection(size_type i) const { return isIntersectLeft(i) || isIntersectRight(i); }
	bool isUnknown(size_type i) const { return m_relations[i].first == SetRelation::unkn; }

private:
	VectorType m_orig;
	VectorType m_other;
	std::vector<int> m_counter;
	std::vector<std::pair<SetRelation, std::string>> m_relations;

}; // DDECompareHelper

} // namespace ddeshelper
} // namespace capd

#endif /* _CAPD_DDECOMPAREHELPER_H_ */
