/*
 * utils.h
 *
 *  Created on: Nov 30, 2015
 *      Author: robson
 */

#ifndef _CAPD_CUSTOM_UTILS_H_
#define _CAPD_CUSTOM_UTILS_H_

#include <scirsc/pretty/strhlp.h>
#include <capd/capdlib.h>
#include <fstream>

using namespace std;
using namespace capd;
using namespace intervals;
using namespace scirsc; // for now, TODO: namespaces

//TODO: namespaces???
//TODO: some of those are already implemented in CAPD

double upper(const double& smth){ return smth; }
template<typename T>
double upper(const T& smth){ return smth.rightBound(); }

double lower(const double& smth){ return smth; }
template<typename T>
double lower(const T& smth){ return smth.leftBound(); }

double middle(const double& smth){ return smth; }
template<typename T>
double middle(const T& smth){
	//typedef typename T::BoundType BT;
	//T m; BT d;
	// capd::intervals::split(smth, m, d);
	T m = (smth.rightBound() + smth.leftBound()) / 2.0;
	return m.leftBound();
}

double diameter(const double& smth){ return 0.0; } // assume 0 for non-intervals...

double diameter(const MpInterval& smth){
	typename MpInterval::BoundType d;
	d = smth.rightBound() - smth.leftBound();
	return capd::multiPrec::toDouble(d, MpReal::RoundUp);
}

template<typename T>
double diameter(const T& smth){
	//T m; typename T::BoundType d;
	//capd::intervals::split(smth, m, d);
	typename T::BoundType d;
	d = smth.rightBound() - smth.leftBound();
	return d;
}

/**
 *
 */
template<typename T>
std::string pm(const T& d){
	std::ostringstream oss;
	#ifdef GLOBAL_PRECISION
	oss.precision(GLOBAL_PRECISION);
	#else
	oss.precision(15);
	#endif
	oss << middle(d) << " +/- " << (diameter(d)/2.0);
	return oss.str();
}


template<typename VectorSpec>
void vectorMaxDiam(const VectorSpec& v, double& max_diam, int& max_i){
	//std::cout << "v: "<< v.dimension() << "\n" << std::flush;
	max_diam = v[0].rightBound() - v[0].leftBound();
    max_i = 0;
    for (int i = 1; i < v.dimension(); i++){
  	  double diam = v[i].rightBound() - v[i].leftBound();
  	  if (diam > max_diam) {
  		  max_diam = diam;
  		  max_i = i;
  	  }
    }
}

template<typename MatrixSpec>
void matrixMaxDiam(const MatrixSpec& A, double& max_diam, int& max_i, int& max_j){
	//std::cout << "A: "<< A.numberOfRows() << ", "<< A.numberOfColumns() << "\n" << std::flush;
	max_diam = A[0][0].rightBound() - A[0][0].leftBound();
    max_i = 0;
    max_j = 0;
    for (int i = 0; i < A.numberOfRows(); i++){
    	for (int j = 0; j < A.numberOfColumns(); j++){
			double diam = A[i][j].rightBound() - A[i][j].leftBound();
			if (diam > max_diam) {
				max_diam = diam;
				max_i = i;
				max_j = j;
			}
    	}
    }
}

template<typename VectorType>
std::vector<double> vlo(VectorType v){
	std::vector<double> r(v.dimension());
	for (int i = 0; i < v.dimension(); i++){
		r[i] = v[i].leftBound();
	}
	return r;
}

template<typename VectorType>
std::vector<double> vup(VectorType v){
	std::vector<double> r(v.dimension());
	for (int i = 0; i < v.dimension(); i++){
		r[i] = v[i].rightBound();
	}
	return r;
}

template<typename VectorType>
std::vector<double> vdiam(VectorType v){
	std::vector<double> r(v.dimension());
	for (int i = 0; i < v.dimension(); i++){
		r[i] = v[i].rightBound() - v[i].leftBound();
	}
	return r;
}

template<typename VectorType>
std::vector<double> vmid(VectorType v){
	std::vector<double> r(v.dimension());
	for (int i = 0; i < v.dimension(); i++){
		r[i] = v[i].mid().leftBound();
	}
	return r;
}

/**
 * Makes the following matrix:
 * A 0
 * 0 B
 */
template<typename MatrixSpec>
MatrixSpec make2BlockMatrix(const MatrixSpec& A, const MatrixSpec& B){
	int nrows = A.numberOfRows() + B.numberOfRows();
	int ncols = A.numberOfColumns() + B.numberOfColumns();
	MatrixSpec C(nrows, ncols);
	for (int i = 1; i <= A.numberOfRows(); i++)
		for (int j = 1; j <= A.numberOfColumns(); j++)
			C(i, j) = A(i, j);
	int rowshift = A.numberOfRows();
	int colshift = A.numberOfColumns();
	for (int i = 1; i <= B.numberOfRows(); i++)
			for (int j = 1; j <= B.numberOfColumns(); j++)
				C(rowshift + i, colshift + j) = B(i, j);
	return C;
}

/**
 * adds n x n identity Id as a block to the matrix.
 * Returns matrix of the form:
 * A 0
 * 0 Id
 */
template<typename MatrixSpec>
MatrixSpec addIdentityBlock(const MatrixSpec& A, int n){
	MatrixSpec Id(n, n); Id.setToIdentity();
	return make2BlockMatrix(A, Id);
}

/**
 * If C is
 * A ?
 * ? B
 * then this function extracts a nrows x ncols block
 * from matrix C to A and the second block goes to B
 * It ignores blocks marked as ?
 */
template<typename MatrixSpec>
void split2BlockMatrix(const MatrixSpec& C, MatrixSpec& A, MatrixSpec &B, int nrows, int ncols){
	int nrowsB = C.numberOfRows() - nrows;
	int ncolsB = C.numberOfColumns() - ncols;
	A = MatrixSpec(nrows, ncols);
	B = MatrixSpec(nrowsB, ncolsB);
	for (int i = 1; i <= nrows; i++)
			for (int j = 1; j <= ncols; j++)
				A(i, j) = C(i, j);
	for (int i = 1; i <= nrowsB; i++)
			for (int j = 1; j <= ncolsB; j++)
				B(i, j) = C(nrows + i, ncols + j);
}

template<typename VectorSpec>
typename VectorSpec::ScalarType sum(const VectorSpec& v){
	typename VectorSpec::ScalarType res = 0.0;
	for (int i = 0; i < v.dimension(); i++){
		res += v[i];
	}
	return res;
}

template<typename ScalarSpec, typename VectorSpec>
ScalarSpec linFit(VectorSpec x, VectorSpec fx, ScalarSpec& a, ScalarSpec& b){
	ScalarSpec sumX = sum(x);
	ScalarSpec sumY = sum(fx);
	ScalarSpec xy = x * fx;
	ScalarSpec xx = x * x;
	ScalarSpec yy = fx * fx;
	ScalarSpec n = x.dimension();
	ScalarSpec denom = (n * xx - sumX * sumX);
	b = (sumY * xx - sumX * xy) / denom;
	a = (n * xy - sumX * sumY) / denom;
	ScalarSpec err = 0.0;
	for (int i = 0; i < x.dimension(); i++){
		ScalarSpec diff = fx[i] - (a * x[i] + b);
		err += diff * diff;
	}
	return sqrt(err);
}

void eigenSort(DVector& real, DVector& imag, DMatrix& real_vectors, DMatrix& imag_vectors){
	DVector modulus = real; for (int i = 0; i < modulus.dimension(); i++) modulus[i] = modulus[i] * modulus[i] + imag[i] * imag[i];
	capd::vectalg::Vector<int, 0> perm(modulus.dimension()); modulus.sorting_permutation(perm);
	DMatrix tmp_real_v = real_vectors;
	DMatrix tmp_imag_v = imag_vectors;
	DVector tmp_real = real;
	DVector tmp_imag = imag;
	for (int i = 0; i < modulus.dimension(); i++){
		real_vectors.column(i) = tmp_real_v.column(perm[i]);
		imag_vectors.column(i) = tmp_imag_v.column(perm[i]);
		real[i] = tmp_real[perm[i]];
		imag[i] = tmp_imag[perm[i]];
	}
}
void eigenRealMatrix(DVector& real, DVector& imag, DMatrix& real_vectors, DMatrix& imag_vectors){
	int size = real.dimension();
	for (int i = 0; i < size; i++){
		if (imag[i] != 0.){
			real_vectors.column(i+1) = imag_vectors.column(i+1);
			i++;
		}
	}
}

void eigenMatrixTest(DMatrix DP){
	int size = DP.numberOfColumns();
	DVector real(size);
	DVector imag(size);
	DVector norm(size);
	DMatrix real_vectors(size, size);
	DMatrix imag_vectors(size, size);
	capd::alglib::computeEigenvaluesAndEigenvectors(DP, real, imag, real_vectors, imag_vectors);
	std::cout << "computeEigenvaluesAndEigenvectors(DP)\n";
	strhlp::prettyPrintEigenvaluesAndEigenvectors(cout, real, imag, real_vectors, imag_vectors);
	// to be sure - sort eigenvalues depending on modulus (sometimes 0 were put before other items)
	eigenRealMatrix(real, imag, real_vectors, imag_vectors);
	eigenSort(real, imag, real_vectors, imag_vectors);
	DMatrix C = real_vectors;
	for (int i = 0; i < size; i++) C.column(i).normalize();


	std::cout << "computeEigenvaluesAndEigenvectors(DP^T)\n";
	DP.Transpose();
	capd::alglib::computeEigenvaluesAndEigenvectors(DP, real, imag, real_vectors, imag_vectors);
	strhlp::prettyPrintEigenvaluesAndEigenvectors(cout, real, imag, real_vectors, imag_vectors);
	eigenRealMatrix(real, imag, real_vectors, imag_vectors);
	eigenSort(real, imag, real_vectors, imag_vectors);
	DMatrix invC = real_vectors;
	for (int i = 0; i < size; i++) invC.column(i).normalize();
	invC.Transpose();

	cout << "eigenMatrixTest C^-1 * C: \n" << flush;
	strhlp::prettyPrint(cout, invC * C) << "\n\n";
	cout << "eigenMatrixTest: C * C^-1\n" << flush;
	strhlp::prettyPrint(cout, C * invC) << "\n\n";
}

void computeEigenvaluesAndEigenvectorsByLeftInverse(DMatrix DP, DVector& real, DVector& imag, DMatrix& real_vectors, DMatrix& imag_vectors){
	DP.Transpose();
	capd::alglib::computeEigenvaluesAndEigenvectors(DP, real, imag, real_vectors, imag_vectors);
	strhlp::prettyPrintEigenvaluesAndEigenvectors(cout, real, imag, real_vectors, imag_vectors);
	eigenRealMatrix(real, imag, real_vectors, imag_vectors);
	eigenSort(real, imag, real_vectors, imag_vectors);
	real_vectors.Transpose();
	real_vectors = capd::matrixAlgorithms::inverseMatrix(real_vectors);
}

void computeEigenvaluesAndEigenvectorsDirectly(DMatrix DP, DVector& real, DVector& imag, DMatrix& real_vectors, DMatrix& imag_vectors){
	capd::alglib::computeEigenvaluesAndEigenvectors(DP, real, imag, real_vectors, imag_vectors);
	strhlp::prettyPrintEigenvaluesAndEigenvectors(cout, real, imag, real_vectors, imag_vectors);
	eigenRealMatrix(real, imag, real_vectors, imag_vectors);
	eigenSort(real, imag, real_vectors, imag_vectors);
}

void eigenMatrixForDP(DMatrix DP, DMatrix& oreal, DMatrix& oimag){
	int size = DP.numberOfColumns();
	DVector real(size);
	DVector imag(size);
	DVector norm(size);
	DMatrix real_vectors(size, size);
	DMatrix imag_vectors(size, size);
	capd::alglib::computeEigenvaluesAndEigenvectors(DP, real, imag, real_vectors, imag_vectors);
	// to be sure - sort eigenvalues depending on modulus (sometimes 0 were put before other items)
	DVector modulus = real; for (int i = 0; i < modulus.dimension(); i++) modulus[i] = modulus[i] * modulus[i] + imag[i] * imag[i];
	capd::vectalg::Vector<int, 0> perm(modulus.dimension()); modulus.sorting_permutation(perm);
	DMatrix tmp_real_v = real_vectors;
	DMatrix tmp_imag_v = imag_vectors;
	DVector tmp_real = real;
	DVector tmp_imag = imag;
	for (int i = 0; i < modulus.dimension(); i++){
		real_vectors.column(i) = tmp_real_v.column(perm[i]);
		imag_vectors.column(i) = tmp_imag_v.column(perm[i]);
		real[i] = tmp_real[perm[i]];
		imag[i] = tmp_imag[perm[i]];
	}
	oreal = real_vectors;
	oimag = imag_vectors;
}

DMatrix eigenRealMatrixForDP(DMatrix DP){
	int size = DP.numberOfColumns();
	DVector real(size);
	DVector imag(size);
	DVector norm(size);
	DMatrix real_vectors(size, size);
	DMatrix imag_vectors(size, size);
	computeEigenvaluesAndEigenvectorsDirectly(DP, real, imag, real_vectors, imag_vectors);
	return real_vectors;
}

///**
// * helps to cope with computing, sorting and selecting various eigenvalues;
// */
//class EigenvectorsHelper {
//public:
//	// typedefs for possible future templatization...
//	typedef DVector VectorType;
//	typedef DMatrix MatrixType;
//	typedef double ScalarType;
//
//	EigenvectorsHelper(MatrixType matrix): m_A(matrix), size(matrix.numberOfColumns()) {
//		real = DVector(size);
//		imag = DVector(size);
//		norm = DVector(size);
//		realEigenvectors = MatrixType(size, size);
//		imagEigenvectors = MatrixType(size, size);
//		computeEigenvaluesAndEigenvectorsDirectly(m_A, real, imag, realEigenvectors, imagEigenvectors);
//	}
//private:
//	MatrixType m_A;
//
//	int size;
//
//	VectorType real;
//	VectorType imag;
//	VectorType norm;
//	MatrixType realEigenvectors;	// in fact: real parts of the eigenvectors
//	MatrixType imagEigenvectors; 	// in fact: imag parts of the eigenvectors
//};

/*** returns interval containing number dec.frac */
interval rigval(int dec, int frac, int p = 10){
	interval sign(1.0);
	if (frac < 0) { frac = -frac; }
	if (dec < 0) { dec = -dec; sign = -1; }
	interval res;
	for (; p <= 1000000000; p *= 10){
		res = interval(frac) / interval(p); if (res < 1.0) break;
	}
	if (!(res < 1.0)) throw -1; // TODO: better exception handling
	return sign * (interval(dec) + res);
}
/** simplest name for common function */
interval RV(int dec, int frac, int p = 10){
	return rigval(dec, frac, p);
}
/** returns interval containing [0, 1/p] */
interval rigsmall(int p = 1000000000){
	return rigval(0, 1, p);
}
/** simplest name for common function */
interval RI(int p){
	return rigsmall(p);
}


/**
 * class to represent rigorous values of paramaters
 * which are pretty in print (fixed precision) but
 * non-rigorous when inserted directly into programs
 *
 * we assume that no more than 8 digits for decimal part and no more than 8 digits for fractional part
 *
 * Sample usages:
 *
 * 		interval a = RigParam("1.1");
 * 		interval b = RigParam("-1.3");
 * 		interval c = RigParam("-0.000013");
 */
class RigParam{
public:
	int dec, frac, fracn;
	short sign;
	interval rigval;
	RigParam(std::string value):
			sign(1),
			dec(0),
			frac(0),
			rigval(0),
			fracn(10){
		std::istringstream ss(value);
		int* cur = &dec;
		int c;
		do c = ss.get(); while(ss && (c == ' ' || c == '+')); 									// discard any + or ' ' character at the begining...
		while(ss){																				// begin reading
			if      (c == '.' && cur == &dec)	{ cur = &frac; fracn = 1; }						// we have . for the first time -> change to frac part
			else if (c == '-' && sign > 0)		{ sign = -1; }									// we have '-' sign -> first time
			else if ('0' <= c && c <= '9') 		{ (*cur) = 10 * (*cur) + (c - '0'); fracn*=10; }// we have digit -> extend current
			else 								{ break; }										// no more valid data -> brak free
			c = ss.get(); 																		// get next char
		}
		if (cur == &dec) fracn = 10; 	// if we had no point . then fracn is 10;
		if (fracn < 10) fracn = 10;		// if we had no digit after .
		rigval = interval(sign) * (interval(dec) + interval(frac) / interval(fracn));	// n is allways != 0
	}
	operator interval() const{
		return this->rigval;
	}
	std::string str() const{
		std::ostringstream sf; sf << (fracn + frac);
		std::ostringstream ss;
		ss << (sign < 0 ? "-" : "+") << dec << "." << sf.str().substr(1);
		return ss.str();
	}
	operator std::string() const{
		return this->str();
	}
	friend std::ostream& operator<<(std::ostream& out, const RigParam& rp){
		out << rp.str();
		return out;
	}
};
// short name
typedef RigParam RP;


template<typename T>
union BinaryData {
	T value;
	char repr[sizeof(T)];

	BinaryData(){};
	BinaryData(T v): value(v) {};
};
template<typename T>
streamsize binsize(const BinaryData<T>& d){
	return sizeof(T);
}

/**
 * User must supply out which supports write and is set to binary!
 * It does not store Vector size, it might be used for both vectors and matrices
 */
template<typename VectorType>
void saveBinary(std::ostream& out, const VectorType& vector){
	typedef typename VectorType::ScalarType ScalarType;
	for (auto it = vector.begin(); it != vector.end(); ++it){
		BinaryData<ScalarType> value(*it);
		out.write(value.repr, binsize(value));
	}
}
/**
 * For simplicity...
 */
template<typename VectorType>
void saveBinary(const std::string& filepath, const VectorType& vector){
	std::ofstream out(filepath, ios::binary);
	saveBinary(out, vector);
	out.close();
}

/**
 * User must supply in which supports read and is set to binary!
 * Vector must be of appropriate size
 * It works for both matrix and vector! (any Container!)
 */
template<typename VectorType>
void readBinary(std::istream& in, VectorType& vector){
	typedef typename VectorType::ScalarType ScalarType;
	for (auto it = vector.begin(); it != vector.end(); ++it){
		BinaryData<ScalarType> value;
		in.read(value.repr, binsize(value));
		*it = value.value;
	}
}
/**
 * For simplicity...
 */
template<typename VectorType>
void readBinary(const std::string& filepath, VectorType& vector){
	std::ifstream in(filepath, ios::binary);
	readBinary(in, vector);
	in.close();
}


/////////////////////////////////////////////////////////////////
// TESTS ////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//{
//	DVector v(3); v[0] = 1.0; v[1] = 2.0; v[2] = 3.0;
//	std::ofstream out("test1.bin", ios::binary);
//	saveBinary(out, v);
//	out.close();
//
//	DVector u(3); u[0] = 15.0; u[1] = 15.0; u[2] = 15.0;
//	std::ifstream in("test1.bin", ios::binary);
//	readBinary(in, u);
//	in.close();
//
//	std::cout << u << " should be {1.0,2.0,3.0} [" << (u == v ? "OK": "BAD") << "]\n";
//}
//
//{
//	IVector v(2); v[0] = interval(1.0, 1.5); v[1] = interval(2.0, 3.0);;
//	std::ofstream out("test2.bin", ios::binary);
//	saveBinary(out, v);
//	out.close();
//
//	IVector u(2); u[0] = 15.0; u[1] = 15.0;
//	std::ifstream in("test2.bin", ios::binary);
//	readBinary(in, u);
//	in.close();
//
//	std::cout << u << " should be {[1.0,1.5],[2.0,3.0]} [" << (u == v ? "OK": "BAD") << "]\n";
//}
//
//{
//	DVector v(2); v[0] = 1.0; v[1] = 2.0;
//	DVector w(3); w[0] = 3.0; w[1] = 4.0; w[2] = 5.0;
//	std::ofstream out("test3.bin", ios::binary);
//	saveBinary(out, v);
//	saveBinary(out, w);
//	out.close();
//
//	DVector u(2); u[0] = 15.0; u[1] = 15.0;
//	std::ifstream in("test3.bin", ios::binary);
//	readBinary(in, u);
//	std::cout << u << " should be {1.0,2.0} [" << (u == v ? "OK": "BAD") << "]\n";
//  u = DVector(3);
//	readBinary(in, u);
//	std::cout << u << " should be {3.0,4.0,5.0} [" << (u == w ? "OK": "BAD") << "]\n";
//
//	in.close();
//}
//
//{
//	DMatrix A(2,2);
//	A[0][0] = 1.0; A[0][1] = 2.0;
//	A[1][0] = 3.0; A[1][1] = 4.0;
//	std::ofstream out("test4.bin", ios::binary);
//	saveBinary(out, A);
//	out.close();
//
//	DMatrix B(2, 2);
//	B[0][0] = 10.0; B[0][1] = 10.0;
//	B[1][0] = 10.0; B[1][1] = 10.0;
//	std::ifstream in("test4.bin", ios::binary);
//	readBinary(in, B);
//	in.close();
//
//	std::cout << B << " should be {{1.0,2.0},{3.0,4.0}} [" << (A == B ? "OK": "BAD") << "]\n";
//}
/////////////////////////////////////////////////////////////////////////////////////////////


#endif /* _CAPD_CUSTOM_UTILS_H_ */
