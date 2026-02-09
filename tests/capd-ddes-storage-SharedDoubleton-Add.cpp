#define BOOST_TEST_MODULE SharedDoubletonAddTest
#include <boost/test/included/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/ddes/storage/SharedDoubleton.h"

BOOST_AUTO_TEST_SUITE(SharedDoubletonAddTestSuite)

typedef capd::ddes::SharedDoubleton<capd::IMatrix> Doubleton;
typedef Doubleton::VectorType VectorType;

// Test Add Set Isolated
BOOST_AUTO_TEST_CASE(AddSetTest) {
    VectorType x1(2); x1[0]=1.0; x1[1]=2.0;
    Doubleton db1(x1);

    VectorType x2(2); x2[0]=0.5; x2[1]=0.5;
    Doubleton db2(x2);

    // This is the operation causing bad_alloc in the main suite
    db1.add(db2);

    // x should be sum
    BOOST_CHECK(db1.get_x()[0] == 1.5);
    BOOST_CHECK(db1.get_x()[1] == 2.5);
}

BOOST_AUTO_TEST_SUITE_END()
