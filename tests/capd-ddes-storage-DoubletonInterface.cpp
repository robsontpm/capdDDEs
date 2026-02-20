#define BOOST_TEST_MODULE DoubletonInterfaceTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/ddes/storage/BasicDoubleton.h"

BOOST_AUTO_TEST_SUITE(DoubletonInterfaceTestSuite)

// Use BasicDoubleton as the concrete implementation for testing the interface
typedef capd::ddes::BasicDoubleton<capd::IMatrix> Doubleton;
typedef Doubleton::VectorType VectorType;
typedef Doubleton::MatrixType MatrixType;
typedef Doubleton::ScalarType ScalarType;

BOOST_AUTO_TEST_CASE(DefaultImplementationTest) {
    // Setup a doubleton X = x + C*r0 + B*r
    int dim = 2;
    int dim_r0 = 1;
    VectorType x(dim); x[0] = 1.0; x[1] = 2.0;
    MatrixType C(dim, dim_r0); C[0][0] = 0.5; C[1][0] = 0.5;
    VectorType r0(dim_r0); r0[0] = capd::Interval(-1, 1);
    MatrixType B(dim, dim); B.setToIdentity();
    VectorType r(dim); r[0] = capd::Interval(-0.1, 0.1); r[1] = capd::Interval(-0.1, 0.1);

    Doubleton d(x, C, r0, B, r);

    // Test midPoint (default: returns x)
    BOOST_CHECK_EQUAL(d.midPoint(), x);

    // Test hull
    // hull = x + C*r0 + B*r
    VectorType expectedHull = x + C * r0 + B * r;
    // We use inclusion check for intervals
    BOOST_CHECK(subset(d.hull(), expectedHull));
    BOOST_CHECK(subset(expectedHull, d.hull()));

    // Test conversion operator
    VectorType v = d;
    BOOST_CHECK(subset(v, expectedHull));
}

BOOST_AUTO_TEST_CASE(GetBinvDefaultTest) {
    int dim = 2;
    Doubleton d(dim, 0); // B is Identity

    // Default get_Binv returns inverse of B
    MatrixType Binv = d.get_Binv();
    MatrixType expectedBinv(dim, dim); expectedBinv.setToIdentity();

    // We check if expectedBinv (exact identity) is contained in computed Binv (interval enclosure)
    BOOST_CHECK(subset(expectedBinv, Binv));

    // Set B to something else
    MatrixType B(dim, dim); B[0][0]=2; B[1][1]=2;
    d.set_B(B);
    Binv = d.get_Binv();
    expectedBinv[0][0]=0.5; expectedBinv[1][1]=0.5;
    BOOST_CHECK(subset(expectedBinv, Binv));
}

BOOST_AUTO_TEST_CASE(PointerSettersTest) {
    int dim = 2;
    Doubleton d(dim, 0);

    // Test set_x(ptr, passOwnership=true)
    VectorType* px = new VectorType(dim);
    (*px)[0] = 10.0;
    d.set_x(px, true);
    BOOST_CHECK_EQUAL(d.get_x()[0], 10.0);
    // px should be deleted, we can't verify that easily but valgrind would

    // Test set_x(ptr, passOwnership=false)
    VectorType x(dim); x[0] = 20.0;
    d.set_x(&x, false);
    BOOST_CHECK_EQUAL(d.get_x()[0], 20.0);

    // Similarly for set_C
    MatrixType* pC = new MatrixType(dim, 0); // N0=0
    d.set_C(pC, true);
    BOOST_CHECK_EQUAL(d.get_C().numberOfRows(), dim);

    // set_r0
    VectorType* pr0 = new VectorType((VectorType::size_type)0);
    d.set_r0(pr0, true);

    // set_B
    MatrixType* pB = new MatrixType(dim, dim); pB->setToIdentity();
    d.set_B(pB, true);

    // set_r
    VectorType* pr = new VectorType(dim);
    d.set_r(pr, true);

    // set_Cr0
    MatrixType* pC2 = new MatrixType(dim, 1);
    VectorType* pr02 = new VectorType(1);
    d.set_Cr0(pC2, pr02, true, true);
    BOOST_CHECK_EQUAL(d.get_C().numberOfColumns(), 1);
    BOOST_CHECK_EQUAL(d.get_r0().dimension(), 1);

    // set_Binv
    MatrixType* pBinv = new MatrixType(dim, dim);
    d.set_Binv(pBinv, true);
    // Note: BasicDoubleton doesn't store Binv, so set_Binv does nothing effectively (implementation detail of BasicDoubleton's base?)
    // Actually DoubletonInterface::set_Binv(ptr) calls set_Binv(ref) which is empty in DoubletonInterface.
    // So this test mainly checks it compiles and runs without memory leak (delete called).
}

BOOST_AUTO_TEST_CASE(OwnershipTransferTakeTest) {
    int dim = 2;
    VectorType x(dim); x[0] = 5.0;
    Doubleton d(x);

    // take_x
    VectorType* px = d.take_x();
    BOOST_CHECK(px != nullptr);
    BOOST_CHECK_EQUAL((*px)[0], 5.0);
    delete px;

    MatrixType* pC = d.take_C();
    BOOST_CHECK(pC != nullptr);
    delete pC;

    VectorType* pr0 = d.take_r0();
    BOOST_CHECK(pr0 != nullptr);
    delete pr0;

    MatrixType* pB = d.take_B();
    BOOST_CHECK(pB != nullptr);
    delete pB;

    VectorType* pr = d.take_r();
    BOOST_CHECK(pr != nullptr);
    delete pr;

    MatrixType* pBinv = d.take_Binv();
    BOOST_CHECK(pBinv != nullptr);
    delete pBinv;
}

BOOST_AUTO_TEST_CASE(CommonCheckTest) {
    int dim = 2;
    Doubleton d(dim, 0);
    VectorType x = d.get_x();

    // Default implementation returns false
    BOOST_CHECK(d.common_x(&x) == true); // BasicDoubleton overrides it to return true if equal
    // Wait, DoubletonInterface returns false. BasicDoubleton overrides it.
    // So we are testing BasicDoubleton's override here mostly.

    VectorType y(dim); y[0] = 100;
    BOOST_CHECK(d.common_x(&y) == false);
}

BOOST_AUTO_TEST_CASE(DotProductTest) {
    int dim = 2;
    Doubleton d(dim, 0);
    // x=0, B=Id, r=0 -> set is {0}
    VectorType v(dim); v[0]=1; v[1]=1;

    ScalarType dot = d.dot(v);
    BOOST_CHECK_EQUAL(dot, 0.0);

    VectorType x(dim); x[0]=1; x[1]=2;
    d.set_x(x);
    // dot = v*x + ... = 1*1 + 1*2 = 3
    dot = d.dot(v);
    BOOST_CHECK_EQUAL(dot, 3.0);
}

BOOST_AUTO_TEST_CASE(OperationsTest) {
    int dim = 2;
    Doubleton d(dim, 0);
    VectorType x(dim); x[0]=1;
    d.set_x(x);

    // Add vector
    VectorType v(dim); v[0]=1;
    d.add(v); // x becomes 2
    BOOST_CHECK_EQUAL(d.get_x()[0], 2.0);

    // Mul scalar
    d.mul(2.0); // x becomes 4
    BOOST_CHECK_EQUAL(d.get_x()[0], 4.0);

    // operator*=
    d *= 0.5; // x becomes 2
    BOOST_CHECK_EQUAL(d.get_x()[0], 2.0);

    // Add set
    Doubleton d2(dim, 0);
    VectorType x2(dim); x2[0]=3;
    d2.set_x(x2);
    d.add(d2); // x should include 2+3=5
    BOOST_CHECK_EQUAL(d.get_x()[0], 5.0);

    // mulThenAdd
    d.mulThenAdd(2.0, d2); // (5*2) + 3 = 13
    BOOST_CHECK_EQUAL(d.get_x()[0], 13.0);
}

BOOST_AUTO_TEST_CASE(StreamOperatorTest) {
    // Use N0=1 to avoid issues with 0-dimension stream I/O if any
    int dim = 2;
    int N0 = 1;
    Doubleton d(dim, N0);
    std::ostringstream oss;
    oss << d;
    BOOST_CHECK(!oss.str().empty());

    std::istringstream iss(oss.str());
    Doubleton d2(dim, N0);
    iss >> d2;

    BOOST_CHECK_EQUAL(d.get_x(), d2.get_x());
    BOOST_CHECK_EQUAL(d.get_C(), d2.get_C());
    BOOST_CHECK_EQUAL(d.get_r0(), d2.get_r0());
    BOOST_CHECK_EQUAL(d.get_B(), d2.get_B());
    // r might differ slightly due to numeric operations or reinitialization, but should be close/same here
    BOOST_CHECK_EQUAL(d.get_r(), d2.get_r());
}

BOOST_AUTO_TEST_CASE(SetCr0ExceptionTest) {
    int dim = 2;
    Doubleton d(dim, 0);

    MatrixType* pC = new MatrixType(dim, 2); // dim_r0 = 2
    VectorType* pr0 = new VectorType(1); // dim_r0 = 1, mismatch!

    // set_Cr0 should throw logic_error and delete pointers
    BOOST_CHECK_THROW(d.set_Cr0(pC, pr0, true, true), std::logic_error);

    // Since we can't check if they were deleted, we trust valgrind.
    // If we didn't pass ownership, we would leak here if exception is thrown and we don't catch it and delete.

    MatrixType* pC2 = new MatrixType(dim, 2);
    VectorType* pr02 = new VectorType(1);
    try {
        d.set_Cr0(pC2, pr02, true, true);
    } catch (...) {
        // Expected
    }
}

BOOST_AUTO_TEST_CASE(MakeStorageTest) {
    int dim = 2;
    int N0 = 1;
    Doubleton d(dim, N0);

    BOOST_CHECK_EQUAL(d.makeStorage_x().dimension(), dim);
    BOOST_CHECK_EQUAL(d.makeStorage_C().numberOfRows(), dim);
    BOOST_CHECK_EQUAL(d.makeStorage_C().numberOfColumns(), N0);
    BOOST_CHECK_EQUAL(d.makeStorage_r0().dimension(), N0);
    BOOST_CHECK_EQUAL(d.makeStorage_B().numberOfRows(), dim);
    BOOST_CHECK_EQUAL(d.makeStorage_B().numberOfColumns(), dim);
    BOOST_CHECK_EQUAL(d.makeStorage_r().dimension(), dim);
}

BOOST_AUTO_TEST_SUITE_END()
