#define BOOST_TEST_MODULE BasicDoubletonTest
#include <boost/test/included/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/ddes/storage/BasicDoubleton.h"

BOOST_AUTO_TEST_SUITE(BasicDoubletonTestSuite)

typedef capd::ddes::BasicDoubleton<capd::IMatrix> Doubleton;
typedef Doubleton::VectorType VectorType;
typedef Doubleton::MatrixType MatrixType;

// Test the default constructor
BOOST_AUTO_TEST_CASE(DefaultConstructorTest) {
    Doubleton doubleton;
    BOOST_CHECK_EQUAL(doubleton.dimension(), 0);
}

// Test the constructor with a vector
BOOST_AUTO_TEST_CASE(ConstructorWithVectorTest) {
    VectorType x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    Doubleton doubleton(x);
    BOOST_CHECK_EQUAL(doubleton.dimension(), 3);
}

// Test the copy constructor
BOOST_AUTO_TEST_CASE(CopyConstructorTest) {
    VectorType x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    Doubleton doubleton1(x);
    Doubleton doubleton2(doubleton1);
    BOOST_CHECK_EQUAL(doubleton2.dimension(), 3);
}

// Test the assignment operator
BOOST_AUTO_TEST_CASE(AssignmentOperatorTest) {
    VectorType x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    Doubleton doubleton1(x);
    Doubleton doubleton2;
    doubleton2 = doubleton1;
    BOOST_CHECK_EQUAL(doubleton2.dimension(), 3);
}

// Test the affineTransform method
BOOST_AUTO_TEST_CASE(AffineTransformTest) {
    VectorType x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    MatrixType M(3, 3);
    M.setToIdentity();
    VectorType v(3);
    v[0] = 1.0;
    v[1] = 1.0;
    v[2] = 1.0;
    Doubleton doubleton(x);
    doubleton.affineTransform(M, v);
    BOOST_CHECK_EQUAL(doubleton.dimension(), 3);
}

// Test the translate method
BOOST_AUTO_TEST_CASE(TranslateTest) {
    VectorType x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    VectorType v(3);
    v[0] = 1.0;
    v[1] = 1.0;
    v[2] = 1.0;
    Doubleton doubleton(x);
    doubleton.translate(v);
    BOOST_CHECK_EQUAL(doubleton.dimension(), 3);
}

// Test the reinitialize method
BOOST_AUTO_TEST_CASE(ReinitializeTest) {
    VectorType x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    Doubleton doubleton(x);
    doubleton.reinitialize(2, 2);
    BOOST_CHECK_EQUAL(doubleton.dimension(), 2);
}

// Test the show method
BOOST_AUTO_TEST_CASE(ShowTest) {
    VectorType x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    Doubleton doubleton(x);
    std::string result = doubleton.show();
    BOOST_CHECK(!result.empty());
}

BOOST_AUTO_TEST_SUITE_END()
