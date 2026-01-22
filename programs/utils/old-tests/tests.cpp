#define MYDEBUG(item) std::cout << #item << ": " << item << "    <END>" << std::endl;

#include <iostream>
#include <iomanip>
#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>
#include <capd/map/Map.hpp>

using namespace std;

//#define WITH_GNUPLOT
#define FOR_VALGRIND

// TODO: better tests!

typedef capd::intervals::Interval<double, capd::rounding::DoubleRounding>  Interval;
typedef Interval Scalar;
typedef Interval Real;
typedef capd::vectalg::Vector<Interval, 0> Vector;
typedef capd::vectalg::Matrix<Interval, 0, 0>  Matrix;
typedef capd::ddes::BasicDoubleton<Matrix> BasicSetType;
typedef capd::ddes::SharedDoubleton<Matrix> SharedSetType;

typedef capd::ddes::DDESolutionCurve<BasicSetType> BasicSolution;
typedef capd::ddes::DDESolutionCurve<SharedSetType> SharedSolution;
typedef typename BasicSolution::GridType Grid;
typedef typename BasicSolution::TimePointType TimePoint;
typedef typename BasicSolution::CurvePieceType BasicCurvePiece;
typedef typename SharedSolution::CurvePieceType SharedCurvePiece;

void test_Grid();
template<typename SetT>
void test_Doubleton(std::string);
template<typename DataType, typename VectorType = typename DataType::VectorType, typename MatrixType = typename DataType::MatrixType>
void test_GenericJet(std::string info, DataType& set);
template<typename CurvePiece>
void test_CurvePieces(std::string);
template<typename CurvePiece>
void test_CurvePiecesCommonR0(std::string);
template<typename CurvePiece>
void test_CurvePiecesEvals(std::string);
template<typename Solution>
void test_SolutionCurve(std::string);
template<typename Solution>
void test_SolutionCurveSubcurve(std::string);
template<typename Solution>
void test_FunctionalMap(std::string);
template<typename Solution>
void test_Solver(std::string, int numIters = 1);
template<typename Solution>
void test_SolverEpsilon(std::string, int numIters, double epsi);
template<typename Solution>
void test_ODETaylor(std::string, Real h = 1./32., int order = 4, int numIters = 128);
template<typename Solution>
void test_JetSection(std::string info);
template<typename Solution>
void test_PoincareMap(std::string info);
// TODO: exhaustive testiong of all functions / constructors
// TODO: make it Boost.Tests or something similar.
// TODO: Test addPiece(ptr, false), with valgrind, there is some memory leak there (in Xi probably)
int main(int, char**){
	test_Grid();
	test_Doubleton<BasicSetType>("Basic");
	test_Doubleton<SharedSetType>("Shared");
	Vector x(2), r0(1); Matrix C(2,1); C[0][0] = 2.0; C[1][0] = 3.0; x[0] = 2.0; x[1] = 5.0; r0[0] = 0.0;
	capd::DVector dv(2); 		test_GenericJet<capd::DVector, capd::DVector, capd::DMatrix>("DVector", dv);
	capd::IVector iv(2); 		test_GenericJet<capd::IVector, capd::IVector, capd::IMatrix>("IVector", iv);
	BasicSetType bs(x, C, r0); 	test_GenericJet<BasicSetType>("BasicSetType", bs);
	SharedSetType ss(x, C, r0); test_GenericJet<SharedSetType>("SharedSetType", ss);
	test_CurvePieces<BasicCurvePiece>("Basic");
	test_CurvePieces<SharedCurvePiece>("Shared");
	test_CurvePiecesCommonR0<BasicCurvePiece>("Basic");
	test_CurvePiecesCommonR0<SharedCurvePiece>("Shared");
	test_CurvePiecesEvals<BasicCurvePiece>("Basic");
	test_CurvePiecesEvals<SharedCurvePiece>("Shared");
	test_SolutionCurve<BasicSolution>("Basic");
	test_SolutionCurve<SharedSolution>("Shared");
	test_SolutionCurveSubcurve<BasicSolution>("Basic");
	test_SolutionCurveSubcurve<SharedSolution>("Shared");
	test_FunctionalMap<BasicSolution>("Basic");
	test_FunctionalMap<SharedSolution>("Shared");
	test_Solver<BasicSolution>("Basic", 1);
	test_Solver<SharedSolution>("Shared", 1);
	#ifdef FOR_VALGRIND
	test_Solver<BasicSolution>("Basic", 5);
	test_Solver<SharedSolution>("Shared", 5);
	#else
	test_Solver<BasicSolution>("Basic", -10);
	test_Solver<SharedSolution>("Shared", -10);
	#endif

	test_ODETaylor<BasicSolution>("Basic", 1./64., 4, 512);
	test_ODETaylor<SharedSolution>("Shared", 1./64., 20, 512);

	test_SolverEpsilon<BasicSolution>("Basic", -10, 1.0/64.);
	test_SolverEpsilon<SharedSolution>("Shared", -10, 1.0/64.);

	test_JetSection<BasicSolution>("Basic");
	test_JetSection<SharedSolution>("Shared");

	test_PoincareMap<BasicSolution>("Basic");
	test_PoincareMap<BasicSolution>("Shared");
	return 0;}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// TEST ROUTINES ARE DEFINED HERE //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void test_Grid(){
	cout << "Testing grid" << endl;
	Real h = Real(2.0);
	Grid grid(h);
	auto t0 = grid.point(0);
	t0 += grid.point(1); cout << t0 << endl;
	t0 += grid.point(1); cout << t0 << endl;
	t0 += grid.point(1); cout << t0 << endl;
	t0 -= grid.point(5); cout << t0 << endl;
	t0 = grid.point(5) + grid.point(-2); cout << t0 << endl;

	auto t1 = grid.point(5); cout << t1 << endl;

	ostringstream oss;
	oss << t0 << " test" << std::endl << "ala ma kota" << std::endl;

	string ossstr = oss.str();
	cout << ossstr << endl;

	istringstream iss(ossstr);
	iss >> t1;
	cout << t0 << " " << t1 << endl;

	Grid grid2(h*2);
	auto t2 = grid2(0);
	istringstream iss2(ossstr);
	try { iss2 >> t2; throw -1; } catch (std::logic_error& e) { cout << "Expected exception OK: " << e.what() << endl; }
	cout << "Finished testing grid" << endl;
}

template<typename SetT>
void test_Doubleton(std::string info){
	cout << "Testing Doubleton : " << info << endl;
	Vector v = Vector(3); v[0] = 0.5; v[1] = -0.5; v[2] = 1.0;
	Vector r = Vector(2); r[0] = Interval(-0.5,0.5); r[1] = Interval(-0.5,0.5);
	Matrix Cgood(3, 2);
	Cgood[0][0] = 1.0; Cgood[0][1] = 0.0;
	Cgood[1][0] = 0.0; Cgood[1][1] = 1.0;
	Cgood[2][0] = 1.0; Cgood[2][1] = 1.0;

	Matrix Bgood(3, 3); Bgood.setToIdentity();
	Matrix M(3, 3); M[0][0] = 2.0; M[1][1] = 3.0;  M[2][2] = 4.0;
	Vector rr = Vector(3); rr[0] = Interval(-0.05,0.05); rr[1] = Interval(-0.05,0.05); rr[2] = Interval(-0.05,0.05);
	try{
		cout << "Testing basics..." << endl;
		SetT data0;
		cout << "0: " << flush << data0.get_x() << endl;
		cout << "0: " << flush << data0.get_r0() << endl;
		cout << "0: " << flush << data0.hull() << endl;
		cout << "0: " << flush << data0 << endl;
		SetT data1(Vector(2));
		cout << "1: " << flush << data1.get_x() << endl;
		cout << "1: " << flush << data1.get_r0() << endl;
		cout << "1: " << flush << data1.hull() << endl;
		cout << "1: " << flush << data1 << endl;
		SetT data2(v + Cgood * r);
		cout << "2: " << flush  << (v + Cgood * r) << endl;
		cout << "2: " << flush  << data2 << endl;
		SetT data3(v, Cgood, r);
		cout << "3: " << flush  << data3 << endl;
		SetT data4(v, Cgood, r, Bgood, rr);
		cout << "4: " << flush  << data4 << endl;
		SetT data5a(v, Cgood, &r, Bgood, rr);
		SetT data5b(v, Cgood, &r, Bgood, rr);
		cout << "5: " << flush  << r << endl;
		cout << "5: " << flush  << data5a << endl;
		cout << "5: " << flush  << data5b << endl;
		r[0] *= 2.0; r[1] *= 0.5;
		cout << "5: " << flush  << r << endl;
		cout << "5: " << flush  << data5a << endl;
		cout << "5: " << flush  << data5b << endl;
		SetT data6(v, Cgood, r);
		cout << "6: " << flush  << data6.hull() << endl;
		data6.translate(v);
		cout << "6: " << flush  << data6.hull() << endl;
		data6.translate(-v);
		cout << "6: " << flush  << data6.hull() << endl;
		SetT data7(v, Cgood, r);
		cout << "7: " << flush  << data7.hull() << endl;
		data7.affineTransform(M, v);
		cout << "7: " << flush  << data7.hull() << endl;
		SetT data8(data1);
		cout << "8: " << flush  << data1.hull() << endl;
		cout << "8: " << flush  << data8.hull() << endl;
		SetT data9;
		data9 = data2;
		cout << "9: " << flush  << data2.hull() << endl;
		cout << "9: " << flush  << data9.hull() << endl;
		Vector* external_r0 = new Vector(2);
		Matrix* external_C = new Matrix(v.dimension(), external_r0->dimension());
		(*external_r0)[0] = Scalar(-1, 1);
		(*external_r0)[1] = Scalar(-2, 2);
		SetT* data10 = new SetT(v, external_r0);
		cout << "10:" << flush  << data10->hull() << endl;
		delete data10; // should not mess up with external_r0
		cout << "10:" << flush  << *external_r0 << endl;
		delete external_r0;
		external_r0 = new Vector(2);
		(*external_r0)[0] = Scalar(-1, 1);
		(*external_r0)[1] = Scalar(-2, 2);
		SetT* data11 = new SetT(v);
		data11->set_Cr0(external_C, external_r0);
		cout << "11:" << flush  << data11->hull() << endl;
		delete data11; // should not mess up with external_r0
		cout << "11:" << flush  << *external_r0 << endl;
		delete external_r0;
		delete external_C;
		cout << "Testing basics finished..." << endl;
	} catch (std::logic_error& exception){
		cout << "exception in basic tests: " << exception.what() << endl;
	}
	cout << "r-after basic test: " << flush  << r << endl;

	try {
		cout << "Operations test started" << endl;
		int d = 2, N0 = 2;
		Vector x1(d); 		x1[0] = 1.0; 				x1[1] = 2.0;
		Vector x2(d); 		x2[0] = -1.0; 				x2[1] = -1.0;
		Matrix C1(d, d);	C1[0][0] = 1.0;				C1[0][1] = 1.0;
							C1[1][0] = 0.5;				C1[1][1] = 2.0;
		Matrix C2(d, d);	C2[0][0] = -1.0;			C2[0][1] = 0.0;
							C2[1][0] = 0.0;				C2[1][1] = -1.0;
		Vector r0(N0); 		r0[0] = Real(-1.0, 1.0); 	r0[1] = Real(-2.0, 2.0);
		Matrix B(d, d);		B.setToIdentity();
		Vector r1(d); 		r1[0] = Real(-0.5, 0.5); 	r1[1] = Real(-1.5, 1.5);
		Vector r2(d); 		r2[0] = Real(-1.5, 1.5); 	r2[1] = Real(-0.0, 0.0);

		SetT set1(x1, C1, &r0, B, r1);
		SetT set2(x2, C2, &r0, B, r2);

		cout << "set1: " << set1 << endl;
		cout << "set2: " << set2 << endl;

		SetT result = set1;
		cout << "mul before " << result << endl;
		cout << "mul before " << result.hull() << endl;
		result.mul(10.0);
		cout << "mul after  " << result << endl;
		cout << "mul after  " << result.hull() << endl;

		result = set1;
		cout << "add vector before " << result << endl;
		cout << "add vector before " << result.hull() << endl;
		result.add(r2);
		cout << "add vector after  " << result << endl;
		cout << "add vector after  " << result.hull() << endl;

		result = set1;
		cout << "add set before " << result << endl;
		cout << "add set before " << result.hull() << endl;
		result.add(set2);
		cout << "add set after  " << result << endl;
		cout << "add set after  " << result.hull() << endl;

		result = set1;
		cout << "mulThenAdd set before " << result << endl;
		cout << "mulThenAdd set before " << result.hull() << endl;
		result.mulThenAdd(0.5, set2);
		cout << "mulThenAdd set after  " << result << endl;
		cout << "mulThenAdd set after  " << result.hull() << endl;

		cout << "Operations test finished" << endl;
	} catch (std::logic_error& exception){
		cout << "exception in operations test: " << exception.what() << endl;
	}

	try {
		cout << "Operations (2) test started" << endl;
		int d = 2, N0 = 2;
		Vector x1(d); 		x1[0] = 1.0; 				x1[1] = 2.0;
		Vector x2(d); 		x2[0] = -1.0; 				x2[1] = -1.0;
		Matrix C1(d, d);	C1.setToIdentity();
		Matrix C2(d, d);	C2.setToIdentity();
		Vector r0(N0); 		r0[0] = Real(-1.0, 1.0); 	r0[1] = Real(-2.0, 2.0);
		Matrix B(d, d);		B.setToIdentity();
		Vector r1(d); 		r1[0] = Real(-1.0, 1.0); 	r1[1] = Real(-2.0, 2.0);
		Vector r2(d); 		r2[0] = Real(-2.0, 2.0); 	r2[1] = Real(-0.0, 0.0);

		SetT set1(x1, C1, &r0, B, r1);
		SetT set2(x2, C2, &r0, B, r2);

		cout << "set1: " << set1 << endl;
		cout << "set2: " << set2 << endl;

		SetT result = set1;
		cout << "mul before " << result << endl;
		cout << "mul before " << result.hull() << endl;
		result.mul(10.0);
		cout << "mul after  " << result << endl;
		cout << "mul after  " << result.hull() << endl;

		result = set1;
		cout << "add vector before " << result << endl;
		cout << "add vector before " << result.hull() << endl;
		result.add(r2);
		cout << "add vector after  " << result << endl;
		cout << "add vector after  " << result.hull() << endl;

		result = set1;
		cout << "add set before " << result << endl;
		cout << "add set before " << result.hull() << endl;
		result.add(set2);
		cout << "add set after  " << result << endl;
		cout << "add set after  " << result.hull() << endl;

		result = set1;
		cout << "mulThenAdd set before " << result << endl;
		cout << "mulThenAdd set before " << result.hull() << endl;
		result.mulThenAdd(2.0, set2);
		cout << "mulThenAdd set after  " << result << endl;
		cout << "mulThenAdd set after  " << result.hull() << endl;

		cout << "Operations test finished" << endl;
	} catch (std::logic_error& exception){
		cout << "exception in operations (2) test: " << exception.what() << endl;
	}

	{
		int d = 2, N0 = 1;
		Vector v(d); for (int i = 0; i < d; i++) v[i] = i+1;
		Vector r0(N0); for (int i = 0; i < N0; i++) r0[i] = Real(-1., 1.0) * 0.5 * (i+1);
		Matrix C(d, N0); for (int i = 0; i < d; i++) for (int j = 0; j < N0; j++) C[i][j] = (i+1)*(j+1);
		SetT set1(v);
		set1.set_Cr0(C, r0);
		SetT set2(Vector(2));

		cout << set1 << endl;
		cout << set2 << endl;

		ostringstream oss; oss << set1 << " test" << endl << "Ala ma kota" << endl;
		string ossstr = oss.str();
		cout << ossstr << endl;
		istringstream iss(ossstr);
		iss >> set2;
		cout << set2 << endl << endl;
	}

	cout << "Finished Testing Doubleton : " << info << endl;
}

template<typename DataType, typename VectorType, typename MatrixType>
void test_GenericJet(std::string info, DataType& set){
	// TODO: finish tests and implementation of the GenericJets, then refactor TaylorForwardCUrvePiece.
	typedef typename MatrixType::ScalarType ScalarType;
	typedef typename MatrixType::ScalarType RealType;
	typedef typename capd::ddes::DiscreteTimeGrid<RealType> GridType;
	typedef typename GridType::TimePointType TimePointType;
	typedef capd::ddes::GenericJet<TimePointType, DataType, VectorType, MatrixType> SomeJet;
	cout << "Testing GenericJet: " << info << endl;
	RealType h = 2.0;
	GridType grid(h);
	VectorType v(set);
	int d = v.dimension();
	for (int i = 0; i < d; ++i) v[i] = 1.0 * (i+1);
	cout << " - constructors " << endl;
	{ SomeJet test; 					cout << "   default:     " << test.show() << endl; }
	{ SomeJet test(grid(2)); 			cout << "   t:           " << test.show() << endl; }
	{ SomeJet test(grid(2), 2); 		cout << "   t dim:       " << test.show() << endl; }
	{ SomeJet test(grid(2), 2, 3); 		cout << "   t dim order: " << test.show() << endl; }
	{ SomeJet test(grid(2), 3, v); 		cout << "   t order v:   " << test.show() << endl; }
	{
	SomeJet test(grid(0), 3, v);		cout << "   t order v:   " << test.show() << endl;
	SomeJet testcpy1(test);				cout << "   copy:        " << testcpy1.show() << endl;
	}
	{
	DataType* c = new DataType[3]; for (int i = 0; i < 3; ++i) c[i] = DataType((1.0 + i) * v);
	SomeJet test(grid(0), c, c+3);		cout << "   t, ptr, ptr: " << test.show() << endl;
	delete[] c;
	}
	{
		cout << "PTR NULL" << endl;
		DataType* ptr = NULL;
		try { SomeJet test(grid(1), ptr, ptr); throw -1; }
		catch(std::logic_error &e) { cout << "PTR NULL exception thrown, OK, e:\n" << e.what() << endl; }
		catch(...) { cout << "ERROR: bad exception caught! " << endl; }
		cout << "DIM DIFF" << endl;
		DataType* c = new DataType[3];
		for (int i = 0; i < 2; ++i)
			c[i] = DataType((1.0 + i) * v);
		c[2] = DataType(VectorType(2*d));
		try { SomeJet test(grid(1), c, c+3); throw -1; }
		catch(std::logic_error &e) { cout << "DIM DIFF exception thrown, OK, e:\n" << e.what() << endl; }
		catch(...) { cout << "ERROR: bad exception caught! " << endl; }
		cout << "BAD DIR " << endl;
		try { SomeJet test(grid(1), c+3, c); throw -1; }
		catch(std::logic_error &e) { cout << "BAD DIR  exception thrown, OK, e:\n" << e.what() << endl; }
		catch(...) { cout << "ERROR: bad exception caught! " << endl; }
		delete[] c;
	}
	{
	std::vector<DataType> coeffs; for (int i = 0; i < 3; ++i) coeffs.push_back(DataType((1.0 + i) * v));
	SomeJet test(grid(0), coeffs);		cout << "   t, vector:   " << test.show() << endl;
	}
	{
		std::vector<DataType> coeffs;
		try { SomeJet test(grid(0), coeffs); throw -1; }
		catch(std::logic_error &e) { cout << "exception thrown, OK, e:\n" << e.what() << endl; }
		catch(...) { cout << "ERROR: bad exception caught! " << endl; }
		coeffs.push_back(v); coeffs.push_back(DataType(VectorType(2*d)));
		try { SomeJet test(grid(0), coeffs); throw -1; }
		catch(std::logic_error &e) { cout << "exception thrown, OK, e:\n" << e.what() << endl; }
		catch(...) { cout << "ERROR: bad exception caught! " << endl; }
	}
	{
	cout << " - copy operator" << endl;
	SomeJet test(grid(1), 3, v);
	SomeJet testcpy1(grid(1));
	cout << "   test:        " << test << endl;
	cout << "   testcpy1:    " << testcpy1 << endl;
	cout << "   copying...   " << endl;
	testcpy1 = test;
	cout << "   test:        " << test << endl;
	cout << "   testcpy1:    " << testcpy1 << endl;
	}
	{
	cout << " - vector conversion" << endl;
	SomeJet test(grid(1), 3, v);
	cout << "   test:        " << test << endl;
	cout << "   hull:        " << test.hull() << endl;
	cout << "   Vector(test):" << VectorType(test) << endl;
	}

	// TODO: other operations.

	{
	// UWAGA: ten test jest zrobiony troche naokolo, bo CAPD zle sobie roadzi z pustymi wektorami i macierzami
	//        np. zapisuje kazda macierz M(d, 0) jako {}, z czego nie da sie odzyskac wymiarow oryginalnych.
	//        wydaje mi sie tez, ze wtedy przeskakuje nad danymi (powinna chyba wypisywac cos w stylu
	//        {
	//        {}
    //		  }
	cout << " - input / output" << endl;
	SomeJet test(grid(1), 3, set);
	std::vector<DataType> coeffs; coeffs.push_back(set); coeffs.push_back(set); coeffs.push_back(set); coeffs.push_back(set);
	test.setupCoeffs(coeffs);
	SomeJet test2(grid(2), 5, 2*v);
	ostringstream oss; oss << test << " test" << endl << "ala ma kote" << endl;
	string ossstr = oss.str();
	cout << ossstr << endl;
	istringstream iss(ossstr);
	iss >> test2;
	cout << test << endl;
	cout << test2 << endl;
	}

	cout << "Finished Testing GenericJet: " << info << endl;
}

template<typename CurvePiece>
void test_CurvePieces(std::string info){
	cout << "Testing curve pieces: " << info << endl;
	Real h = Real(2.0);
	Grid grid(h);
	Vector value(2);
	value[0] = 1.0; value[1] = 2.0;
	cout << " - basic constructor (constant function)" << endl;
	CurvePiece test(grid(0), 3, value);
	cout << " - basic set structure" << endl;
	cout << "a" << test.get_x() << "a" << endl;
	cout << "b" << test.get_r0() << "b" << endl;
	cout << "c" << test.get_r() << "c" << endl;
	cout << "d" << test.get_C() << "d" << endl;
	cout << "e" << test.get_B() << "e" << endl;

	cout << " - basic output" << endl;
	cout << test.midPoint() << " FINISEHD" << endl;
	cout << test.hull() << " FINISEHD" <<endl;
	cout << test.show() << " FINISEHD" << endl;

	cout << " - copy operator" << endl;
	CurvePiece testcpy1;
	cout << testcpy1.midPoint() << " testcpy1 a FINISEHD" << endl;
	cout << testcpy1.hull() << " testcpy1 a FINISEHD" <<endl;
	cout << testcpy1.show() << " testcpy1 a FINISEHD" << endl;
	testcpy1 = test;
	cout << testcpy1.midPoint() << " testcpy1 b FINISEHD" << endl;
	cout << testcpy1.hull() << " testcpy1 b FINISEHD" <<endl;
	cout << testcpy1.show() << " testcpy1 b FINISEHD" << endl;

	cout << " - copy constructor" << endl;
	CurvePiece testcpy2(test);
	cout << testcpy2.midPoint() << " testcpy1 b FINISEHD" << endl;
	cout << testcpy2.hull() << " testcpy1 b FINISEHD" <<endl;
	cout << testcpy2.show() << " testcpy1 b FINISEHD" << endl;

	cout << " - midCurve" << endl;
	cout << test.midCurve() << " FINISEHD" << endl;
	cout << testcpy1.midCurve() << " FINISEHD" << endl;
	cout << testcpy2.midCurve() << " FINISEHD" << endl;

	cout << " - exceptions" << endl;
	try { cout << "A" << test.jetAt(grid.point(0)) << endl; } catch (std::logic_error& err){ cout << "A" << err.what() << endl; }
	try { cout << "B" << test.jetAt(grid.point(2)) << endl; } catch (std::logic_error& err){ cout << "B" << err.what() << endl; }
	try { cout << "C" << test.jetAt(0.5) << endl; } catch (std::logic_error& err){ cout << "C" << err.what() << endl; }

	cout << " - eval, taylor, summa" << endl;
	for (double t = -1.5; t <= 1.5; t += 0.25){
		cout << "T(" << t << ") " << test.taylor(t) << endl;
		cout << "S(" << t << ") " << test.summa(t) << endl;
		cout << "E(" << t << ") " << test.eval(t) << endl;
	}

	cout << " - midCurve(), further operations" << endl;
	CurvePiece testMid = test.midCurve();
	cout << testMid.midPoint() << endl;
	cout << testMid.hull() << endl;
	cout << testMid.show() << endl;
	try { cout << "A" << testMid.jetAt(grid.point(0)) << endl; } catch (std::logic_error& err){ cout << "A" << err.what() << endl; }
	try { cout << "B" << testMid.jetAt(grid.point(2)) << endl; } catch (std::logic_error& err){ cout << "B" << err.what() << endl; }
	try { cout << "C" << testMid.jetAt(0.5) << endl; } catch (std::logic_error& err){ cout << "C" << err.what() << endl; }
	for (double t = -1.5; t <= 1.5; t += 0.25){
		cout << "T(" << t << ") " << testMid.taylor(t) << endl;
		cout << "S(" << t << ") " << testMid.summa(t) << endl;
		cout << "E(" << t << ") " << testMid.eval(t) << endl;
	}

	cout << " - cube1D" << endl;
	CurvePiece cube1D(grid(0), 1, 3);
	Vector x(4); x[0] = 2.0; x[1] = 0.0; x[2] = 0.0; x[3] = 1.0;
	Vector r0(1); r0[0] = Interval(-0.01, 0.01);
	Matrix C(4, 1); C[0][0] = 1.0; C[1][0] = 1.0; C[2][0] = 2.0; C[3][0] = 6.0;
	cout << " - cube1D (2)" << endl;
	cube1D.set_x(x);
	cube1D.set_Cr0(C, r0);
	cout << " - cube1D (3)" << endl;
	cout << cube1D.hull() << "\n";
	cout << cube1D.midPoint() << "\n";
	cout << cube1D.show() << "\n";

	#ifdef WITH_GNUPLOT
		std::ofstream outf("test.dat");
		capd::ddeshelper::value_to_gnuplot(outf, grid.point(0), grid.point(200), 1000, cube1D);
		outf.close();
		std::ofstream outg("test.gp");
		outg << "plot 'test.dat' using 1:($3-$4):($3+$4) with filledcu, 'test.dat' using 1:3 with lines, x**3+2 with points" << endl;
		outg.close();
		system("gnuplot -p test.gp");
	#endif

	cout << "I/O Test START" << endl;
	CurvePiece testIO(grid(0), 1, 2);
	cout << cube1D << endl;
	cout << testIO << endl;

	ostringstream oss;
	oss << cube1D << " test" << endl << "ala ma kota" << endl;
	string ossstr = oss.str();
	cout << ossstr << endl;
	istringstream iss(ossstr);
	iss >> testIO;
	cout << testIO << endl;
	cout << "I/O Test END" << endl;

	cout << "Finished Testing curve pieces: " << info << endl;
}

template<typename CurvePiece>
void test_CurvePiecesCommonR0(std::string info){
	cout << "Testing curve pieces common r0: " << info << endl;
	Real h = Real(2.0);
	Grid grid(h);
	Vector value(2);
	value[0] = 1.0; value[1] = 2.0;
	CurvePiece test(grid.point(0), 3, value);
	Vector *common_r0 = new Vector(test.get_r0());
	test.set_r0(common_r0);

	cout << test.midPoint() << " FINISEHD" << endl;
	cout << test.hull() << " FINISEHD" <<endl;
	cout << test.show() << " FINISEHD" << endl;

	cout << "Testing copy constructor" << endl;
	CurvePiece testcpy2(test);
	cout << testcpy2.midPoint() << " testcpy1 b FINISEHD" << endl;
	cout << testcpy2.hull() << " testcpy1 b FINISEHD" <<endl;
	cout << testcpy2.show() << " testcpy1 b FINISEHD" << endl;

	cout << "Testing copy operator" << endl;
	CurvePiece testcpy1;
	cout << testcpy1.midPoint() << " testcpy1 a FINISEHD" << endl;
	cout << testcpy1.hull() << " testcpy1 a FINISEHD" <<endl;
	cout << testcpy1.show() << " testcpy1 a FINISEHD" << endl;
	testcpy1 = test;
	cout << testcpy1.midPoint() << " testcpy1 b FINISEHD" << endl;
	cout << testcpy1.hull() << " testcpy1 b FINISEHD" <<endl;
	cout << testcpy1.show() << " testcpy1 b FINISEHD" << endl;

	cout << test.midCurve() << " FINISEHD" << endl;
	cout << testcpy1.midCurve() << " FINISEHD" << endl;
	cout << testcpy2.midCurve() << " FINISEHD" << endl;

	try { cout << "A" << test.jetAt(grid.point(0)) << endl; } catch (std::logic_error& err){ cout << "A" << err.what() << endl; }
	try { cout << "B" << test.jetAt(grid.point(2)) << endl; } catch (std::logic_error& err){ cout << "B" << err.what() << endl; }
	try { cout << "C" << test.jetAt(0.5) << endl; } catch (std::logic_error& err){ cout << "C" << err.what() << endl; }
	for (double t = -1.5; t <= 1.5; t += 0.25){
		cout << "T(" << t << ") " << test.taylor(t) << endl;
		cout << "S(" << t << ") " << test.summa(t) << endl;
		cout << "E(" << t << ") " << test.eval(t) << endl;
	}

	CurvePiece testMid = test.midCurve();
	cout << testMid.midPoint() << endl;
	cout << testMid.hull() << endl;
	cout << testMid.show() << endl;
	try { cout << "A" << testMid.jetAt(grid.point(0)) << endl; } catch (std::logic_error& err){ cout << "A" << err.what() << endl; }
	try { cout << "B" << testMid.jetAt(grid.point(2)) << endl; } catch (std::logic_error& err){ cout << "B" << err.what() << endl; }
	try { cout << "C" << testMid.jetAt(0.5) << endl; } catch (std::logic_error& err){ cout << "C" << err.what() << endl; }
	for (double t = -1.5; t <= 1.5; t += 0.25){
		cout << "T(" << t << ") " << testMid.taylor(t) << endl;
		cout << "S(" << t << ") " << testMid.summa(t) << endl;
		cout << "E(" << t << ") " << testMid.eval(t) << endl;
	}

	delete common_r0;
	cout << "Finished testing curve pieces common r0: " << info << endl;
}

template<typename CurvePiece>
void test_CurvePiecesEvals(std::string info){

	typedef typename CurvePiece::SetType SetType;

	int d = 1, n0 = 3;
	cout << "Testing curve pieces EVALS: " << info << endl;
	Real h = Real(2.0);
	Grid grid(h);


	CurvePiece cube1D(grid(0), d, 3);
	Vector x(4); x[0] = 2.0; x[1] = 0.0; x[2] = 0.0; x[3] = 1.0;
	Vector r0(n0);
	r0[0] = Interval(-0.01, 0.01);
	r0[1] = Interval(-0.01, 0.01);
	r0[2] = Interval(-0.01, 0.01);
	Matrix C(4, n0);
	C[0][0] = 0.5; C[1][0] = 0.25; C[2][0] = 1.0; C[3][0] = 3.0;
	C[0][1] = 0.25; C[1][1] = 0.50; C[2][1] = 1.0; C[3][1] = 2.0;
	C[0][2] = 0.25; C[1][2] = 0.25; C[2][2] = 1.0; C[3][2] = 1.0;
	cube1D.set_x(x);
	cube1D.set_Cr0(C, r0);
	cout << cube1D.hull() << "\n";
	cout << cube1D.midPoint() << "\n";
	cout << cube1D.show() << "\n";

	#ifdef WITH_GNUPLOT
		std::ofstream outf("cube_evals_x.dat");
		capd::ddeshelper::value_to_gnuplot(outf, grid.point(0), grid.point(200), 1000, cube1D);
		outf.close();
		std::ofstream outg("cube_evals_x.gp");
		outg << "plot 'cube_evals_x.dat' using 1:($3-$4):($3+$4) with filledcu, 'cube_evals_x.dat' using 1:3 with lines, x**3+2 with points" << endl;
		outg.close();
		system("gnuplot -p cube_evals_x.gp");
	#endif

	Real dt = 0.25;
	SetType out1(d, n0);
	cube1D.evalAtDelta(dt, out1);

	cout << "evalAtDelta(set):       " << Vector(out1) << endl;
	cout << "vector evalDtAtDelta(): " << cube1D.evalAtDelta(dt) << endl;

	SetType out2(d, n0);
	cube1D.evalCoeffAtDelta(2, dt, out2);

	cout << "evalCoeffAtDelta(set):     " << Vector(out2) << endl;
	cout << "vector evalCoeffAtDelta(): " << cube1D.evalCoeffAtDelta(2, dt) << endl;

	SetType out3(d, n0);
	cube1D.taylorAtDelta(dt, out3);

	cout << "taylorAtDelta(set):     " << Vector(out3) << endl;
	cout << "vector taylorAtDelta(): " << cube1D.taylorAtDelta(dt) << endl;

	cout << "Finished Testing curve pieces EVALS: " << info << endl;
}

template<typename Solution>
void test_SolutionCurve(std::string info){
	typedef typename Solution::CurvePieceType CurvePiece;
	cout << "Testing solution curve: " << info << endl;
	Real h = Real(2.0);
	Grid grid(h);

	Scalar v0 = 2.0; Vector v0v(1); v0v[0] = v0;
	cout << "- making solution" << endl;
	Solution sol(grid, 1, 1);

	cout << "- making cube" << endl;
	CurvePiece cube1D(grid.point(0), 1, 3);
	Vector x(4); x[0] = v0; x[1] = 0.0; x[2] = 0.0; x[3] = 1.0;
	Vector r0(1); r0[0] = Interval(-0.01, 0.01);
	Matrix C(4, 1); C[0][0] = 1.0; C[1][0] = 1.0; C[2][0] = 2.0; C[3][0] = 6.0;
	cube1D.set_x(x);
	cube1D.set_Cr0(C, r0);

	cout << " - extending solution" << endl;
	cout << sol.pastTime() << " " << sol.currentTime() << endl;
	sol.addPiece(cube1D);	cout << "Added cube" << endl;
	sol.set_r0(cube1D.get_r0());
	cout << sol.pastTime() << " " << sol.currentTime() << endl;
	sol.addPiece(cube1D);	cout << "Added cube" << endl;
	cout << sol.pastTime() << " " << sol.currentTime() << endl;
	sol.addPiece(cube1D).addPiece(cube1D).addPiece(cube1D).addPiece(cube1D);	cout << "Added some cube" << endl;
	cout << sol.pastTime() << " " << sol.currentTime() << endl;

	cout << " - eval" << endl;
	cout << sol.eval(0.00) << endl;
	cout << sol.eval(0.01) << endl;
	cout << sol.eval(0.02) << endl;
	cout << sol.eval(0.03) << endl;
	try{
		cout << sol.eval(100) << endl;
	} catch (std::exception& e){
		cout << e.what() << endl;
	}

	#ifdef WITH_GNUPLOT
	{
	std::ofstream outf("test.dat"); //outf.precision(15);
	capd::ddeshelper::value_to_gnuplot(outf, sol.pastTime(), sol.currentTime(), Interval(0.1 * h), sol);
	outf.close();
	std::ofstream outg("test.gp"); //outf.precision(15);
	outg << "plot 'test.dat' using 1:($3-$4):($3+$4) with filledcu, 'test.dat' using 1:3 with lines" << endl;
	outg.close();
	system("gnuplot -p test.gp");
	}
	#endif
	cout << " - set value at current" << endl;
	sol.setValueAtCurrent(v0v);
	#ifdef WITH_GNUPLOT
	{
	std::ofstream outf("test.dat"); //outf.precision(15);
	capd::ddeshelper::value_to_gnuplot(outf, sol.pastTime(), sol.currentTime(), Interval(0.1 * h), sol);
	outf.close();
	std::ofstream outg("test.gp"); //outf.precision(15);
	outg << "plot 'test.dat' using 1:($3-$4):($3+$4) with filledcu, 'test.dat' using 1:3 with lines" << endl;
	outg.close();
	system("gnuplot -p test.gp");
	}
	#endif

	cout << " - j's" << endl;
	cout << sol.rightDomain() << endl;
	cout << sol.j(grid.point(0)) << endl;
	cout << sol.j(grid.point(1)) << endl;
	cout << sol.j(grid.point(2)) << endl;
	cout << sol.j(grid.point(0), 0) << endl;
	cout << sol.j(grid.point(0), 1) << endl;
	cout << sol.j(grid.point(0), 2) << endl;

	cout << " - copy constructor" << endl;
	Solution testcpy2(sol);
	cout << testcpy2.show() << endl;

	cout << " - copy operator" << endl;
	Solution testcpy1(grid);
	testcpy1 = sol;
	cout << testcpy1.show() << endl;

	cout << " - mid curve" << endl;
	cout << testcpy2.midCurve().show() << endl;

	typedef typename Solution::SetType SetType;
	cout << " - constructor constant over interval (vector)" << endl; {
		Interval epsi(-0.0001, 0.0001);
		Vector v(2); for (int j = 0; j < 2; j++) v[j] = 1.1 + epsi;
		Solution X(grid, grid(0), grid(10), 10, v);
		cout << X.show() << endl;
		cout << "DONE" << endl;
	}

	cout << " - constructor constant over interval (set)" << endl; {
		Interval epsi(-0.0001, 0.0001);
		Vector r0(2); for (int j = 0; j < 2; j++) r0[j] = epsi;
		Matrix C(2,2); C.setToIdentity();
		Vector v(2); for (int j = 0; j < 2; j++) v[j] = 1.1;
		SetType set(v, C, r0);
		Solution X(grid, grid(0), grid(10), 10, set);
		cout << X.show() << endl;
		cout << "DONE" << endl;
	}

	cout << " - constructor constant over interval (const set)" << endl; {
		Interval epsi(-0.0001, 0.0001);
		Vector r0(2); for (int j = 0; j < 2; j++) r0[j] = epsi;
		Matrix C(2,2); C.setToIdentity();
		Vector v(2); for (int j = 0; j < 2; j++) v[j] = 1.1;
		const SetType set(v, C, r0);
		Solution X(grid, grid(0), grid(10), 10, set);
		cout << X.show() << endl;
		cout << "DONE" << endl;
	}

	cout << " - constructor from data" << endl;
	cout << "TODO" << endl;


	cout << "Finished testing solution curve: " << info << endl;
}

template<typename Solution>
void test_SolutionCurveSubcurve(std::string info){
	typedef typename Solution::CurvePieceType CurvePiece;
	cout << "Testing solution curve subcurve: " << info << endl;
	Real h = Real(2.0);
	Grid grid(h);

	Scalar v0 = 2.0; Vector v0v(1); v0v[0] = v0;
	cout << "- making solution" << endl;
	Solution sol(grid, 1, 1);

	cout << "- making cube" << endl;
	CurvePiece cube1D(grid.point(0), 1, 3);
	Vector x(4); x[0] = v0; x[1] = 0.0; x[2] = 0.0; x[3] = 1.0;
	Vector r0(1); r0[0] = Interval(-0.01, 0.01);
	Matrix C(4, 1); C[0][0] = 1.0; C[1][0] = 1.0; C[2][0] = 2.0; C[3][0] = 6.0;
	cube1D.set_x(x);
	cube1D.set_Cr0(C, r0);

	cout << " - extending solution" << endl;
	cout << sol.pastTime() << " " << sol.currentTime() << endl;
	sol.addPiece(cube1D);	cout << "Added cube" << endl;
	sol.set_r0(cube1D.get_r0());
	cout << sol.pastTime() << " " << sol.currentTime() << endl;
	sol.addPiece(cube1D);	cout << "Added cube" << endl;
	cout << sol.pastTime() << " " << sol.currentTime() << endl;
	sol.addPiece(cube1D).addPiece(cube1D).addPiece(cube1D).addPiece(cube1D);	cout << "Added some cubes" << endl;
	cout << sol.pastTime() << " " << sol.currentTime() << endl;

	auto t0 = grid(1);
	auto t1 = grid(3);

	cout << "SOL" <<  endl;
	cout << sol << endl;

	cout << "subcurve 1 between " << t0 << " and " << sol.currentTime() <<  endl;
	Solution sub1 = sol.subcurve(t0);
	cout << sub1 << endl;

	cout << "subcurve 1 between " << t0 << " and " << t1 <<  endl;
	Solution sub2 = sol.subcurve(t0, t1);
	cout << sub2 << endl;

	#ifdef WITH_GNUPLOT
	{
	{ std::ofstream outf("suball.dat"); capd::ddeshelper::value_to_gnuplot(outf, sol.pastTime(), sol.currentTime(), Interval(0.1 * h), sol); outf.close(); }
	{ std::ofstream outf("sub1.dat"); capd::ddeshelper::value_to_gnuplot(outf, sub1.pastTime(), sub1.currentTime(), Interval(0.1 * h), sub1); outf.close(); }
	{ std::ofstream outf("sub2.dat"); capd::ddeshelper::value_to_gnuplot(outf, sub2.pastTime(), sub2.currentTime(), Interval(0.1 * h), sub2); outf.close(); }
	std::ofstream outg("testsubcurve.gp"); //outf.precision(15);
	outg << "set multiplot layout 3,1" << endl;
	outg << "set xrange [" <<  Real(sol.pastTime()).leftBound() << ":" <<  Real(sol.currentTime()).rightBound() << "]" << endl;
	outg << "plot 'suball.dat' using 1:($3-$4):($3+$4) with filledcu, 'suball.dat' using 1:3 with lines" << endl;
	outg << "plot 'sub1.dat' using 1:($3-$4):($3+$4) with filledcu, 'sub1.dat' using 1:3 with lines" << endl;
	outg << "plot 'sub2.dat' using 1:($3-$4):($3+$4) with filledcu, 'sub2.dat' using 1:3 with lines" << endl;
	outg.close();
	system("gnuplot -p testsubcurve.gp");
	}
	#endif


	cout << "Finished testing solution curve subcurve: " << info << endl;
}

template<typename Solution>
void test_FunctionalMap(std::string info){
	typedef typename Solution::CurvePieceType CurvePiece;
	typedef typename Solution::SetType SetType;
	cout << "Testing functional map: " << info << endl;

	Real h = Real(1.0 / 32.0);
	Grid grid(h);
	CurvePiece cube1D(grid.point(0), 1, 3);
	Vector x(4); x[0] = 2.0; x[1] = 0.0; x[2] = 0.0; x[3] = 1.0;
	Vector r0(1); r0[0] = Interval(-0.01, 0.01);
	Matrix C(4, 1); C[0][0] = 1.0; C[1][0] = 1.0; C[2][0] = 2.0; C[3][0] = 6.0;
	cube1D.set_x(x);
	cube1D.set_Cr0(C, r0);

	Solution sol(grid, 1, 1); sol.set_r0(r0); // TODO: this is important to set common r0, otherwise, addPiece will forgot about r0 from added Piece! Add EXCEPTION or some handling!
	cout << sol.pastTime() << " " << sol.currentTime() << endl;
	for (int j = 0; j < 10; j++)
		sol.addPiece(cube1D);

	Vector value0(1); value0[0] = 2.0 + Interval(-0.01, 0.01);
	sol.setValueAtCurrent(value0);

	typedef capd::ddes::ToyModel Eq;
	Eq f;
	Vector v(2);
	v[0] = 1.1; v[1] = 1.1;
	Vector w(1);
	f(grid.point(0), v, w);
	cout << "f(x) = " << w << endl;

	typedef capd::ddes::DiscreteDelaysFunctionalMap<Eq, Solution> DDEq;
	typedef typename DDEq::VariableStorageType Variables;
	typedef typename DDEq::JacobianStorageType Jacobians;
	typedef typename DDEq::ValueStorageType Values;
	auto tau = grid.point(4);
	DDEq rhs(f, tau);
	cout << sol.eval(sol.t0()) << " " << sol.eval(sol.t0() - tau) << endl;
	v[0] = sol.eval(sol.t0())[0];
	v[1] = sol.eval(sol.t0() - tau)[0];
	cout << "curve(0), curve(-tau) = " << v << endl;
	f(sol.t0(), v, w);
	cout << "f(curve(0), curve(-tau)) = " << w << endl;
	cout << "f(curve) = " << rhs(sol) << endl;

	Values coeffs, coeffs2, coeffsM;
	cout << "computeDDECoefficients-simple started" << endl;
	rhs.computeDDECoefficients(sol.t0(), sol, coeffs);
	rhs.computeDDECoefficients(sol.t0(), sol, coeffs2);
	rhs.computeDDECoefficients(sol.t0(), sol.midCurve(), coeffsM);
	cout << "computeDDECoefficients-simple finished" << endl;
	CurvePiece jet1(sol.t0(), coeffs);
	CurvePiece jet2(sol.t0(), coeffs2);
	CurvePiece jetM(sol.t0(), coeffsM);
	cout << "Without AD" << endl;
	cout << "A1 " << jet1 << endl;
	cout << "A2 " << jet2 << endl;
	cout << "B1 " << jet1.hull() << endl;
	cout << "B2 " << jet2.hull() << endl;
	cout << "C1 " << jet1.midPoint() << endl;
	cout << "C2 " << jet2.midPoint() << endl;

	std::cout << "should be: \n";
	std::cout << "B1 {[1.99,2.01],[3.9601,4.0401],[3.93025,4.07035],[2.5802,2.754],[1.72092,1.94691]} \n";
	std::cout << "B2 {[1.99,2.01],[3.9601,4.0401],[3.93025,4.07035],[2.5802,2.754],[1.72092,1.94691]} \n";
	std::cout << "C1 {[2,2],[4.0001,4.0001],[4.0003,4.0003],[2.6671,2.6671],[1.83392,1.83392]} \n";
	std::cout << "C2 {[2,2],[4.0001,4.0001],[4.0003,4.0003],[2.6671,2.6671],[1.83392,1.83392]} \n";

	cout << "AM " << jetM << endl;
	cout << "BM " << jetM.hull() << endl;
	cout << "CM " << jetM.midPoint() << endl;

	Values coeffs3;
	Variables u;
	Jacobians Du;
	cout << "computeDDECoefficients-uDu started" << endl;
	rhs.computeDDECoefficients(sol.t0(), sol, coeffs2, u, Du);
	rhs.computeDDECoefficients(sol.t0(), sol, coeffs3, u, Du);
	cout << "computeDDECoefficients-uDu finished" << endl;
	jet2 = CurvePiece(sol.t0(), coeffs2);
	CurvePiece jet3(sol.t0(), coeffs3);

	cout << "With AD" << endl;
	cout << "A1 " << jet1 << endl;
	cout << "A2 " << jet2 << endl;
	cout << "A3 " << jet3 << endl;
	cout << "B1 " << jet1.hull() << endl;
	cout << "B2 " << jet2.hull() << endl;
	cout << "B3 " << jet3.hull() << endl;
	cout << "C1 " << jet1.midPoint() << endl;
	cout << "C2 " << jet2.midPoint() << endl;
	cout << "C3 " << jet3.midPoint() << endl;
	cout << sol.eval(sol.t0()) << endl;
	cout << sol.eval(sol.pastTime()) << endl;

	for (int j = 0; j < u.size(); ++j)
		cout << "u[" << j << "] = " << u[j] << endl;

	for (int k = 0; k < coeffs.size(); ++k)
		for (int j = 0; j < u.size(); ++j)
			cout << "dcoeffs[" << k << "] / du[" << j << "] = " << Du[k][j] << endl;

	cout << "Finished testing functional map: " << info << endl;
}

template<typename Solution>
void test_Solver(std::string info, int numIters){
	typedef typename Solution::CurvePieceType CurvePiece;
	typedef typename Solution::SetType SetType;
	cout << "Testing SOLVER of MackeyGlass " << info << endl;
	Real par_tau 	= 2.0;
	Real par_beta 	= 2.0;
	Real par_gamma 	= 1.0;
	Real par_n 		= 6.0;
	int p = 32; int n = 4; int d = 1;
	Grid grid(par_tau / p);
	auto tau = grid.point(p);
	auto t_0 = grid.point(0);
	if (numIters < 0) numIters = -numIters * p;

	typedef capd::ddes::MackeyGlass<Real> Eq;
	typedef capd::ddes::DiscreteDelaysFunctionalMap<Eq, Solution> DDEq;
	typedef capd::ddes::DDETaylorSolver<DDEq> Solver;
	typedef typename Solver::VariableStorageType Variables;
	typedef typename Solver::JacobianStorageType Jacobians;
	typedef typename Solver::ValueStorageType Values;
	typedef typename Solver::JetType Jet;
	typedef typename Solver::size_type size_type;

	Interval epsi(-0.0001, 0.0001);
	Vector v(d); for (int j = 0; j < d; j++) v[j] = 1.1 + epsi;
	Solution X(grid, -tau, t_0, n, v);
	std::cout << X.show() << std::endl;

	DDEq dde(Eq(par_beta, par_gamma, par_n), tau);
	Solver solver(dde, 10);

	typedef typename Solver::VariableStorageType Variables;
	typedef typename Solver::JacobianStorageType Jacobians;
	typedef typename Solver::ValueStorageType Values;
	Variables u; Values u_encl;
	Jacobians D_uPhi_j0, D_uPhi_z;
	Values Phi_z, Rem_z;
	Values Phi_j0, Rem_j0, Y;
	TimePoint t_h = grid.point(0);
	Real HH;
	solver.encloseSolution(X, t_h, HH, u, u_encl, Phi_j0, D_uPhi_j0, Rem_j0, Y, Phi_z, D_uPhi_z, Rem_z);

	std::cout << X.show() << std::endl;

//	for (auto iu = u.begin(); iu != u.end(); ++iu)
//		MYDEBUG(*iu);
//		// std::cout << *iu << std::endl;
//	MYDEBUG(t_h);
//	MYDEBUG(HH);
//	for (auto ij0 = Phi_j0.begin(); ij0 != Phi_j0.end(); ++ij0)
//		MYDEBUG(*ij0);
//	for (auto iDu = D_uPhi_j0.begin(); iDu != D_uPhi_j0.end(); ++iDu)
//		for (auto jDu = iDu->begin(); jDu != iDu->end(); ++jDu)
//			MYDEBUG(*jDu);
//	for (auto irem = Rem_j0.begin(); irem != Rem_j0.end(); ++irem)
//		MYDEBUG(*irem);
//
//	MYDEBUG(Phi_z);
//	for (auto iDu_z = D_uPhi_z.begin(); iDu_z != D_uPhi_z.end(); ++iDu_z)
//		MYDEBUG(*iDu_z);
//	MYDEBUG(Rem_z);
//	MYDEBUG(Y);

	size_type itersCount = 0;
	try{
		cout << "move test" << endl;
		while (itersCount++ < numIters){
			X.move(solver);
		}
		cout << "move test end" << endl;
	} catch (...) {
		std::cout << "Error integrating at numIters = " << itersCount << " / " << numIters << std::endl;
	}

	#ifdef WITH_GNUPLOT
	std::ofstream outf("test.dat"); //outf.precision(15);
	capd::ddeshelper::value_to_gnuplot(outf, X.pastTime(), X.currentTime(), grid.point(1), X);
	outf.close();
	std::ofstream outg("test.gp"); //outf.precision(15);
	outg << "set yrange [-1.5:1.5]" << endl;
	outg << "plot 'test.dat' using 1:($3-$4):($3+$4) with filledcu, 'test.dat' using 1:3 with lines" << endl;
	outg.close();
	system("gnuplot -p test.gp");
	#endif

	cout << "Finished testing SOLVER of MackeyGlass " << info << endl;
}

template<typename Solution>
void test_SolverEpsilon(std::string info, int numIters, double epsi){
	typedef typename Solution::CurvePieceType CurvePiece;
	typedef typename Solution::SetType SetType;
	cout << "Testing SOLVER EPSILON of MackeyGlass " << info << endl;
	Real par_tau 	= 2.0;
	Real par_beta 	= 2.0;
	Real par_gamma 	= 1.0;
	Real par_n 		= 6.0;
	int p = 32; int n = 4; int d = 1;
	Grid grid(par_tau / p);
	auto h = grid.point(1);
	auto tau = grid.point(p);
	auto t_0 = grid.point(0);
	if (numIters < 0) numIters = -numIters * p;

	typedef capd::ddes::MackeyGlass<Real> Eq;
	typedef capd::ddes::DiscreteDelaysFunctionalMap<Eq, Solution> DDEq;
	typedef capd::ddes::DDETaylorSolver<DDEq> Solver;
	typedef typename Solver::VariableStorageType Variables;
	typedef typename Solver::JacobianStorageType Jacobians;
	typedef typename Solver::ValueStorageType Values;
	typedef typename Solver::JetType Jet;
	typedef typename Solver::size_type size_type;

	Real repsi(-0.001,0.001);
	Vector v(d); for (int j = 0; j < d; j++) v[j] = 1.1;
	Solution X(grid, -tau, t_0, n, v);
	size_type storage_d = X.storageDimension();
	Matrix C(storage_d, storage_d); C.setToIdentity();
	Vector r0(storage_d); for (int i = 0; i < storage_d; ++i) r0[i] = repsi;
	X.set_Cr0(C, r0);

	cout << X << endl;

	DDEq dde(Eq(par_beta, par_gamma, par_n), tau);
	Solver solver(dde, 10);

	cout << "moving into position... " << flush;
	while (--numIters) X.move(solver);
	cout << "DONE!" << endl;

	Vector zero(1);
	//Solution Y(grid, -tau, t_0, n, zero);
	Solution Y(grid, X.t0() - tau - h, X.t0() - h, n, zero, X.storageN0());
	X.epsilonShift(solver, Real(epsi), Y);
	//solver.epsilonShift(X, epsi, Y);

	#ifdef WITH_GNUPLOT
	{
	std::ofstream outX("epsitestX.dat"); //outf.precision(15);
	std::ofstream outY("epsitestY.dat"); //outf.precision(15);
	capd::ddeshelper::value_to_gnuplot(outX, X.currentTime() - tau - h, X.currentTime() - h, 0.1 * Real(h), X);
	capd::ddeshelper::value_to_gnuplot(outY, Y.pastTime(), Y.currentTime(), 0.1 * Real(h), Y);
	outX.close();
	outY.close();
	std::ofstream outg("epsitest.gp"); //outf.precision(15);
	outg << "set yrange [-0.1:1.5]" << endl;
	outg << "plot 'epsitestX.dat' using 1:($3-$4):($3+$4) with filledcu, 'epsitestX.dat' using 1:3 with lines";
	outg << ", 'epsitestY.dat' using ($1+" << epsi << "):($3-$4):($3+$4) with filledcu, 'epsitestY.dat' using ($1+" << epsi << "):3 with lines" << endl;
	outg.close();
	system("gnuplot -p epsitest.gp");
	}
	#endif

	cout << "Finished testing SOLVER EPSILON of MackeyGlass " << info << endl;
}
template<typename Solution>
void test_ODETaylor(std::string info, Real h, int order, int numIters){
	cout << "Testing Compare ODETaylor(CAPD) vs DDETaylor in ODE setting " << info << endl;
	int d = 2;
	double r = 0.1;
	Interval epsi(-r, r);
	capd::interval epsi_capd(-r, r);
	Vector initial(d);
	capd::IVector initial_capd(d);
	for (int j = 0; j < d; j++){
		initial[j] = 1.1 + epsi;
		initial_capd[j] = 1.1 + epsi_capd;
	}
	capd::interval h_capd(h.leftBound(), h.rightBound());

	typedef capd::IOdeSolver CAPDOdeSolver;
	// TODO: ponizszy nie dziala - linker error wychodzi...
	//typedef capd::dynsys::OdeSolver<capd::IMap, capd::dynsys::ILastTermsStepControl, capd::dynsys::FirstOrderEnclosure> CAPDOdeSolver;
	// I will use it to extract definition for my set
	typedef capd::dynset::C0DoubletonSet<capd::IMatrix, capd::C0Intv2Policies> RectSetIdQR;
	RectSetIdQR X_capd(initial_capd);
	capd::IMatrix C_capd = X_capd.get_C();
	capd::IVector r0_capd = X_capd.get_r0();
	capd::IVector x_capd = X_capd.get_x();

	Vector r0_dde(d);
	Vector x_dde(d);
	Matrix C_dde(d, d);
	for (int i = 0; i < d; ++i){
		r0_dde[i] = Interval(r0_capd[i].leftBound(), r0_capd[i].rightBound());
		x_dde[i] = Interval(x_capd[i].leftBound(), x_capd[i].rightBound());
		for (int j = 0; j < d; ++j)
			C_dde[i][j] = Interval(C_capd[i][j].leftBound(), C_capd[i][j].rightBound());
	}

	{
	// dde version
	typedef typename Solution::CurvePieceType CurvePiece;
	typedef typename Solution::SetType SetType;
	Grid grid(h);
	auto t_0 = grid.point(0);
	if (numIters < 0) numIters = -numIters;

	typedef capd::ddes::ODEPendulum Eq;
	typedef capd::ddes::DiscreteDelaysFunctionalMap<Eq, Solution> DDEq;
	typedef capd::ddes::DDETaylorSolver<DDEq> Solver;
	typedef typename Solver::VariableStorageType Variables;
	typedef typename Solver::JacobianStorageType Jacobians;
	typedef typename Solver::ValueStorageType Values;
	typedef typename Solver::JetType Jet;
	typedef typename Solver::size_type size_type;

	Solution X(grid, t_0, SetType(x_dde, C_dde, r0_dde));

	DDEq dde(Eq(), 0, {});
	Solver solver(dde, order);

	size_type itersCount = 0;
	std::ofstream outenc("enclosuresDDE.dat"); //outf.precision(15);
	try{
		cout << "move test" << endl;
		while (itersCount++ < numIters){
			X.move(solver);
			capd::ddeshelper::to_dat(outenc, Real(X.t0())); outenc << " "; capd::ddeshelper::to_dat(outenc, X.getLastEnclosure()); outenc << endl;
		}
		cout << "move test end" << endl;
	} catch (...) {
		std::cout << "Error integrating DDE at numIters = " << itersCount << " / " << numIters << std::endl;
	}
	outenc.close();

	#ifdef WITH_GNUPLOT
	std::ofstream outf("testDDE.dat"); outf.precision(15);
	capd::ddeshelper::value_to_gnuplot(outf, X.pastTime(), X.currentTime(), grid.point(1), X);
	outf.close();
	std::ofstream outg("testDDE.gp"); //outf.precision(15);
	outg << "set yrange [-2.0:2.0]" << endl;
	outg << "plot 'testDDE.dat' using 1:($3-$4):($3+$4) with filledcu, 'testDDE.dat' using 1:3 with lines" << endl;
	outg.close();
	system("gnuplot -p testDDE.gp");
	#endif
	}

	{
	// CAPD version, with IdQRPolicy (same as basic DDE code).
	capd::IMap pendulum("var:x,y;fun:y,-x;");
	CAPDOdeSolver solver(pendulum, order);
	solver.turnOffStepControl();
	solver.setStep(capd::interval(h.leftBound()));

	RectSetIdQR X(X_capd);

	std::ofstream outf("testODE.dat"); outf.precision(15);
	std::ofstream outenc("enclosuresODE.dat"); //outf.precision(15);
	int itersCount = 0;

	capd::interval t = 0.0;
	capd::IVector v = X;
	capd::ddeshelper::to_dat(outf, t); outf << " "; capd::ddeshelper::to_dat(outf, v); outf << endl;
	try{
		cout << "move test" << endl;
		while (itersCount++ < numIters){
			t += h_capd;
			X.move(solver);
			capd::ddeshelper::to_dat(outenc, t); outenc << " "; capd::ddeshelper::to_dat(outenc, X.getLastEnclosure()); outenc << endl;
			v = X;
			capd::ddeshelper::to_dat(outf, t); outf << " "; capd::ddeshelper::to_dat(outf, v); outf << endl;
		}
		cout << "move test end" << endl;
	} catch (...) {
		std::cout << "Error integrating ODE at numIters = " << itersCount << " / " << numIters << std::endl;
	}
	outf.close();
	outenc.close();

	#ifdef WITH_GNUPLOT
	std::ofstream outg("testODE.gp"); //outf.precision(15);
	outg << "set yrange [-2.0:2.0]" << endl;
	outg << "plot 'testODE.dat' using 1:($3-$4):($3+$4) with filledcu, 'testODE.dat' using 1:3 with lines" << endl;
	outg.close();
	system("gnuplot -p testODE.gp");
	#endif
	}


	cout << "Finished testing Compare ODETaylor(CAPD) vs DDETaylor in ODE setting " << info << endl;
}

template<typename Solution>
void test_JetSection(std::string info){
	cout << "Testing JetSection " << info << endl;

	typedef typename Solution::CurvePieceType CurvePiece;
	typedef typename Solution::SetType SetType;
	typedef typename SetType::ScalarType ScalarType;
	Real par_tau 	= 2.0;
	Real par_beta 	= 2.0;
	Real par_gamma 	= 1.0;
	Real par_n 		= 6.0;
	int p = 32; int n = 4; int d = 1;
	Grid grid(par_tau / p);
	auto tau = grid.point(p);
	auto t_0 = grid.point(0);

	typedef capd::ddes::MackeyGlass<Real> Eq;
	typedef capd::ddes::DiscreteDelaysFunctionalMap<Eq, Solution> DDEq;
	typedef capd::ddes::DDETaylorSolver<DDEq> Solver;
	typedef typename Solver::VariableStorageType Variables;
	typedef typename Solver::JacobianStorageType Jacobians;
	typedef typename Solver::ValueStorageType Values;
	typedef typename Solver::JetType Jet;
	typedef typename Solver::size_type size_type;

	Vector v(d); for (int j = 0; j < d; j++) v[j] = 1.1;
	Solution X(grid, -tau, t_0, n, v);
	std::cout << X.show() << std::endl;

	DDEq dde(Eq(par_beta, par_gamma, par_n), tau);
	Solver solver(dde, 10);

	typedef capd::ddes::DDEJetSection<Solution> JetSection;
	typedef typename JetSection::JetType SecJet;
	JetSection section(d, 0);
	ScalarType SX = section(X);
	std::cout << "SX = " << SX << endl;

	Vector s(d); s(1.0) = -1.0;
	section.set_s(s);
	SX = section(X);
	std::cout << "SX = " << SX << endl;

	section.set_c(1.0);
	SX = section(X);
	std::cout << "SX = " << SX << endl;

	SecJet j0(grid(0), d, n);
	j0[0][0] = 1.0;
	j0[1][0] = 1.0;
	section.extend(j0);
	SX = section(X);
	std::cout << "SX = " << SX << endl;

	SecJet j1(grid(0), d, n-2);
	j1[0][0] = 1.0;
	j1[1][0] = 1.0;
	section.extend(j1);
	SX = section(X);
	std::cout << "SX = " << SX << endl;

	section.set_c(section(X));

	#ifdef FOR_VALGRINND
	size_type numIters = 10;
	#else
	size_type numIters = p * 5;
	#endif
	size_type itersCount = 0;
	try{
		cout << "move test with section" << endl;
		while (itersCount++ < numIters){
			X.move(solver);
			SX = section(X);
			std::cout << "SX = " << SX << endl;
		}
	} catch (...) {
		std::cout << "Error integrating at numIters = " << itersCount << " / " << numIters << std::endl;
	}

	cout << "Finished testing JetSection " << info << endl;
}

template<typename Solution>
void test_PoincareMap(std::string info){
	cout << "Testing PoincareMap" << info << endl;

	typedef typename Solution::CurvePieceType CurvePiece;
	typedef typename Solution::SetType SetType;
	typedef typename SetType::ScalarType ScalarType;
	Real par_tau 	= 2.0;
	Real par_beta 	= 2.0;
	Real par_gamma 	= 1.0;
	Real par_n 		= 6.0;
	int p = 32; int n = 4; int d = 1;
	Grid grid(par_tau / p);
	auto tau = grid.point(p);
	auto t_0 = grid.point(0);

	typedef capd::ddes::MackeyGlass<Real> Eq;
	typedef capd::ddes::DiscreteDelaysFunctionalMap<Eq, Solution> DDEq;
	typedef capd::ddes::DDETaylorSolver<DDEq> Solver;
	typedef typename Solver::VariableStorageType Variables;
	typedef typename Solver::JacobianStorageType Jacobians;
	typedef typename Solver::ValueStorageType Values;
	typedef typename Solver::JetType Jet;
	typedef typename Solver::size_type size_type;

	Vector v(d); for (int j = 0; j < d; j++) v[j] = 1.1;
	Solution X(grid, -tau, t_0, n, v);
	std::cout << X.show() << std::endl;

	std::cout << X.storageN0() << endl;
	Solution on_section(grid, -tau, t_0, n, (0. * v) );
	std::cout << on_section.show() << std::endl;

	DDEq dde(Eq(par_beta, par_gamma, par_n), tau);
	Solver solver(dde, 10);

	typedef capd::ddes::DDEJetSection<Solution> JetSection;
	typedef typename JetSection::JetType SecJet;
	typedef capd::ddes::DDEBasicPoincareMap<Solver, JetSection> PoincareMap;
	JetSection section(d, 0, 1.0);

	#ifdef FOR_VALGRIND
	const int REQ_STEPS = p;
	#else
	const int REQ_STEPS = p * (n+1);
	#endif
	PoincareMap pm(solver, section);
	pm.setDirection(pm.detectCrossingDirection(X));
	pm.setRequiredSteps(REQ_STEPS);
	pm.setMaxSteps(5 * REQ_STEPS);

	Real reachTime;
	pm(X, on_section, reachTime);
	cout << "reachTime = " << reachTime << endl;
	cout << on_section.show() << endl;

	#ifdef FOR_VALGRIND
	const int NUM_ITERS = 0;
	#else
	const int NUM_ITERS = 5;
	#endif
	for (int i = 0; i < NUM_ITERS; i++){
		X = on_section;
		pm(X, on_section, reachTime);
		cout << "reachTime = " << reachTime << endl;
	}

	#ifdef WITH_GNUPLOT
	std::ofstream outfX("testPM_X.dat"); outfX.precision(15);
	capd::ddeshelper::value_to_gnuplot(outfX, on_section.pastTime(), on_section.currentTime(), grid.point(1), X);
	std::ofstream outfPX("testPM_PX.dat"); outfPX.precision(15);
	capd::ddeshelper::value_to_gnuplot(outfPX, on_section.pastTime(), on_section.currentTime(), grid.point(1), on_section);
	outfX.close(); outfPX.close();
	std::ofstream outg("testPM.gp"); //outf.precision(15);
	outg << "set yrange [-2.0:2.0]" << endl;
	outg << "plot 'testPM_X.dat' using 1:($3-$4):($3+$4) with filledcu, 'testPM_X.dat' using 1:3 with lines, ";
	outg << "     'testPM_PX.dat' using 1:($3-$4):($3+$4) with filledcu, 'testPM_PX.dat' using 1:3 with lines" << endl;
	outg.close();
	system("gnuplot -p testPM.gp");
	#endif

	cout << "Finished testing PoincareMap " << info << endl;
}
