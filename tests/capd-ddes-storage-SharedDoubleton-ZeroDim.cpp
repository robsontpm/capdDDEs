#define BOOST_TEST_MODULE SharedDoubletonZeroDimTest
#include <boost/test/included/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/ddes/storage/SharedDoubleton.h"

BOOST_AUTO_TEST_SUITE(SharedDoubletonZeroDimTestSuite)

typedef capd::ddes::SharedDoubleton<capd::IMatrix> Doubleton;
typedef Doubleton::VectorType VectorType;
typedef Doubleton::size_type size_type;

// Test Zero Dimension
BOOST_AUTO_TEST_CASE(ZeroDimTest) {
    // Constructor with dimensions
    Doubleton db(0, 0);
    BOOST_CHECK_EQUAL(db.dimension(), 0);
    BOOST_CHECK_EQUAL(db.storageN0(), 0);

    // Add zero vector
    VectorType v((size_type)0);
    db.add(v); // Should not crash

    // Mul
    db.mul(2.0); // Should not crash

    // Constructor with zero vector
    VectorType x((size_type)0);
    Doubleton db2(x);
    BOOST_CHECK_EQUAL(db2.dimension(), 0);
    BOOST_CHECK_EQUAL(db2.storageN0(), 0);
}

BOOST_AUTO_TEST_SUITE_END()
