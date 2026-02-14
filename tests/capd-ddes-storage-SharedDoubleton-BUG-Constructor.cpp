#define BOOST_TEST_MODULE SharedDoubletonBUGConstructorTestSuite
#include <boost/test/included/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/ddes/storage/SharedDoubleton.h"
#include <exception>

BOOST_AUTO_TEST_SUITE(SharedDoubletonBUGConstructorTestSuite)

typedef capd::ddes::SharedDoubleton<capd::IMatrix> Doubleton;
typedef Doubleton::VectorType VectorType;
typedef Doubleton::MatrixType MatrixType;
typedef Doubleton::size_type size_type;
typedef Doubleton::ScalarType ScalarType;

BOOST_AUTO_TEST_CASE(KnownBugConstructorTest) {
    // This constructor causes bad_alloc due to unsigned integer underflow
    // SharedDoubleton(size_type d, size_type N0 = -1)
    try {
        Doubleton db(2); // Should trigger bad_alloc
        BOOST_WARN_MESSAGE(false, "Bug fixed? Constructor(d) did not throw bad_alloc.");
    } catch (std::bad_alloc&) {
        BOOST_TEST_MESSAGE("Caught expected bad_alloc from buggy constructor");
    } catch (...) {
        BOOST_WARN_MESSAGE(false, "Constructor(d) threw unknown exception.");
    }
}

BOOST_AUTO_TEST_SUITE_END()
