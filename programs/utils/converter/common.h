/*
 * common.h
 *
 *  Created on: Sep 3, 2020
 *      Author: robson
 */

#ifndef EXAMPLES_CONVERTER_COMMON_H_
#define EXAMPLES_CONVERTER_COMMON_H_

#include <capd/capdlib.h>
#include <capd/mpcapdlib.h>
#include <capd/ddeshelper/ddeshelperlib.h>

using namespace std;
using namespace capd;
using namespace multiPrec;

typedef capd::vectalg::Matrix<capd::Interval, 0, 0> MIMatrix;
typedef capd::vectalg::Matrix<double, 0, 0>  MDMatrix;
typedef capd::vectalg::Vector<capd::Interval, 0> MIVector;
typedef capd::vectalg::Vector<double, 0>  MDVector;

template<typename T>
std::string strfix(T const& a, std::string::size_type len){
	std::ostringstream oss; oss.precision(15);
	oss << a;
	std::string result = oss.str();
	while (result.length() < len) result += " ";
	if (result.length() > len) result = result.substr(0, len);
	return result;
}

capd::Interval mpi_to_interval(const capd::MpInterval& a){
	return capd::Interval(toDouble(a.leftBound(), capd::multiPrec::MpReal::RoundDown), toDouble(a.rightBound(), capd::multiPrec::MpReal::RoundUp));
	//return Interval(toDouble(a.leftBound(), MpReal::RoundDown), toDouble(a.rightBound(), MpReal::RoundUp));
}

MpInterval interval_to_mpi(const Interval& a){
	return MpInterval(a.leftBound(), a.rightBound());
}

MpInterval interval_to_mpi(const Interval& a, int precission){
	int old_precision = MpReal::getDefaultPrecision();
	MpReal::setDefaultPrecision(precission);
	return MpInterval(a.leftBound(), a.rightBound());
	MpReal::setDefaultPrecision(old_precision);
}

template<typename MatrixSpec>
void to_mpi_matrix(const MatrixSpec& inA, MpIMatrix& outA, int precission = 512){
	int old_precision = MpReal::getDefaultPrecision();
	MpReal::setDefaultPrecision(precission);
	int size = inA.numberOfColumns(); // assume square matrix
	outA = MpIMatrix(size, size);
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			outA[i][j] = MpInterval(leftBound(inA[i][j]), rightBound(inA[i][j]));
	MpReal::setDefaultPrecision(old_precision);
}

template<typename MatrixSpec>
void from_mpi_matrix(const MpIMatrix& inA, MatrixSpec& outA){
	int size = inA.numberOfColumns(); // assume square matrix
	outA = MatrixSpec(size, size);
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			outA[i][j] = mpi_to_interval(inA[i][j]);
}

template<typename MatrixSpec>
void to_mpi_matrix(const MatrixSpec& inA, MpMatrix& outA, int precission = 512){
	int old_precision = MpReal::getDefaultPrecision();
	MpReal::setDefaultPrecision(precission);
	int size = inA.numberOfColumns(); // assume square matrix
	outA = MpIMatrix(size, size);
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			outA[i][j] = MpFloat(inA[i][j]);
	MpReal::setDefaultPrecision(old_precision);
}

template<typename MatrixSpec>
void from_mpi_matrix(const MpMatrix& inA, MatrixSpec& outA){
	int size = inA.numberOfColumns(); // assume square matrix
	outA = MatrixSpec(size, size);
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			outA[i][j] = toDouble(inA[i][j], capd::multiPrec::MpReal::RoundNearest);
}

template<typename VectorSpec>
void to_mpi_vector(const VectorSpec& in_v, MpIVector& out_v, int precission = 512){
	int old_precision = MpReal::getDefaultPrecision();
	MpReal::setDefaultPrecision(precission);
	int size = in_v.dimension();
	out_v = MpIVector(size);
	for (int i = 0; i < size; ++i)
		out_v[i] = MpInterval(leftBound(in_v[i]), rightBound(in_v[i]));
	MpReal::setDefaultPrecision(old_precision);
}

template<typename VectorSpec>
void from_mpi_vector(const MpIVector& in_v, VectorSpec& out_v){
	int size = in_v.dimension();
	out_v = VectorSpec(size);
	for (int i = 0; i < size; ++i)
		out_v[i] = mpi_to_interval(in_v[i]);
}

template<typename VectorSpec>
void to_mpi_vector(const VectorSpec& in_v, MpVector& out_v, int precission = 512){
	int old_precision = MpReal::getDefaultPrecision();
	MpReal::setDefaultPrecision(precission);
	int size = in_v.dimension();
	out_v = MpVector(size);
	for (int i = 0; i < size; ++i)
		out_v[i] = MpFloat(in_v[i]);
	MpReal::setDefaultPrecision(old_precision);
}

template<typename VectorSpec>
void from_mpi_vector(const MpVector& in_v, VectorSpec& out_v){
	int size = in_v.dimension();
	out_v = VectorSpec(size);
	for (int i = 0; i < size; ++i)
		out_v[i] = toDouble(in_v[i], capd::multiPrec::MpReal::RoundNearest);
}

MpIMatrix mpi_inverse(MpIMatrix& A, int precission = 512){
	int old_precision = MpReal::getDefaultPrecision();
	MpReal::setDefaultPrecision(precission);
	MpIMatrix mpA = A;
//	mpA.Transpose();
	MpIMatrix mpAinv = capd::matrixAlgorithms::inverseMatrix(A);
//	mpAinv.Transpose();
	MpReal::setDefaultPrecision(old_precision);
	return mpAinv;
}

template<typename M>
void testInverse(M const& A, M const& Ainv){
	int size = A.numberOfRows(); // assume square matrix
	M Id(size, size); Id.setToIdentity();
	M AAinv = A * Ainv - Id;
	M AinvA = Ainv * A - Id;

	if (capd::vectalg::containsZero(AAinv)){
		cout << "AA^{-1} - Id contains zero matrix, OK" << endl;
	} else {
		cout << "AA^{-1} - Id DOES NOT contains zero matrix, ERROR!!!" << endl;
	}

	if (capd::vectalg::containsZero(AinvA)){
		cout << "A^{-1}A - Id contains zero matrix, OK" << endl;
	} else {
		cout << "A^{-1}A - Id DOES NOT contains zero matrix, ERROR!!!" << endl;
	}

	auto maxAAinvDiam = AAinv[0][0].rightBound() - AAinv[0][0].leftBound();
	auto maxAinvADiam = AinvA[0][0].rightBound() - AinvA[0][0].leftBound();
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j){
			auto x = AAinv[i][j].rightBound() - AAinv[i][j].leftBound();
			if (maxAAinvDiam < x) maxAAinvDiam = x;
			x = AinvA[i][j].rightBound() - AinvA[i][j].leftBound();
			if (maxAinvADiam < x) maxAinvADiam = x;
		}

	cout << "Max diam AAinv:                  " << maxAAinvDiam << endl;
	cout << "Max diam AinvA (more important): " << maxAinvADiam << endl;
}


enum SetRelation {
	suPset = 0, suBset = 1, miss_L = 2, miss_R = 3, int_L = 4, int_R = 5, unkn = 6, equal = 7,
};

template<typename Scalar>
std::pair<SetRelation, std::string> check_inclusion(Scalar a, Scalar b, std::vector<int> &counter){
	if (counter.size() != 8) counter.resize(8, 0);
	std::string op = "";
	SetRelation rel = SetRelation::unkn;
	if (a == b){
		op = "equal ";
		rel = SetRelation::equal;
	}else if (a.contains(b)){
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
	++counter[rel];
	return std::make_pair(rel, op);
}

void print_counter(std::ostream &out, std::vector<int> const& counter){
	out << "equals   " << counter[SetRelation::equal] << endl;
	out << "suPsets  " << counter[SetRelation::suPset] << endl;
	out << "suBsets  " << counter[SetRelation::suBset] << endl;
	out << "miss L   " << counter[SetRelation::miss_L] << endl;
	out << "miss R   " << counter[SetRelation::miss_R] << endl;
	out << "int L    " << counter[SetRelation::int_L] << endl;
	out << "int R    " << counter[SetRelation::int_R] << endl;
	out << "unknowns " << counter[SetRelation::unkn] << endl;
}

template<typename Real>
std::string printInterval(Real v){
	std::ostringstream oss; oss.precision(15);
	oss << "[" << strfix(v.leftBound(), 25) << "," << strfix(v.rightBound(), 25) << "]" << " (diam: " << strfix(v.rightBound() - v.leftBound(), 25) << ")";
	return oss.str();
}


#endif /* EXAMPLES_CONVERTER_COMMON_H_ */
