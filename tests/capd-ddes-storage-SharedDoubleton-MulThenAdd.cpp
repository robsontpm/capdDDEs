#define BOOST_TEST_MODULE SharedDoubletonMulThenAddTest
#include <boost/test/included/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/ddes/storage/SharedDoubleton.h"

BOOST_AUTO_TEST_SUITE(SharedDoubletonMulThenAddTestSuite)

typedef capd::ddes::SharedDoubleton<capd::IMatrix> Doubleton;
typedef Doubleton::VectorType VectorType;
typedef Doubleton::ScalarType ScalarType;

// Test MulThenAdd Set Isolated
BOOST_AUTO_TEST_CASE(MulThenAddTest) {
    VectorType x1(2); x1[0]=1.0; x1[1]=2.0;
    Doubleton db1(x1);

    VectorType x2(2); x2[0]=0.5; x2[1]=0.5;
    Doubleton db2(x2);

    ScalarType c = 2.0;

    // This operation likely caused bad_alloc in main suite
    db1.mulThenAdd(c, db2);

    // x1 should be (x1 * 2) + x2
    // (1.0 * 2) + 0.5 = 2.5
    // (2.0 * 2) + 0.5 = 4.5
    BOOST_CHECK(db1.get_x()[0] == 2.5);
    BOOST_CHECK(db1.get_x()[1] == 4.5);
}

BOOST_AUTO_TEST_SUITE_END()
