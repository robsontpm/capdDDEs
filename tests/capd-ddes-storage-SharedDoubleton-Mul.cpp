#define BOOST_TEST_MODULE SharedDoubletonMulTest
#include <boost/test/included/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/ddes/storage/SharedDoubleton.h"

BOOST_AUTO_TEST_SUITE(SharedDoubletonMulTestSuite)

typedef capd::ddes::SharedDoubleton<capd::IMatrix> Doubleton;
typedef Doubleton::VectorType VectorType;
typedef Doubleton::ScalarType ScalarType;

// Test Mul Set Isolated
BOOST_AUTO_TEST_CASE(MulTest) {
    VectorType x(2); x[0]=1.0; x[1]=2.0;
    Doubleton db(x);

    ScalarType c = 2.0;

    // This operation likely caused bad_alloc in main suite
    db.mul(c);

    // x should be scaled
    BOOST_CHECK(db.get_x()[0] == 2.0);
    BOOST_CHECK(db.get_x()[1] == 4.0);
}

BOOST_AUTO_TEST_SUITE_END()
