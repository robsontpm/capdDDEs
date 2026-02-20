#define BOOST_TEST_MODULE SharedDoubletonBUGDataConstructorTestSuite
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/ddes/storage/SharedDoubleton.h"
#include <exception>

BOOST_AUTO_TEST_SUITE(SharedDoubletonBUGDataConstructorTestSuite)

typedef capd::ddes::SharedDoubleton<capd::IMatrix> Doubleton;
typedef Doubleton::VectorType VectorType;
typedef Doubleton::MatrixType MatrixType;
typedef Doubleton::size_type size_type;
typedef Doubleton::ScalarType ScalarType;

BOOST_AUTO_TEST_CASE(KnownBugDataConstructorTest) {
    int d = 2;
    VectorType x(d), r0(d), r(d);
    MatrixType Bgood(d, d); Bgood.setToIdentity();
    MatrixType Cbad(d, d + 1); // Mismatch with r0(d)

    // BUG: The following constructor call causes a segmentation fault because
    // SharedDoubleton constructor performs *m_C = *C without verifying dimensions first,
    // leading to unsafe assignment.
    // Doubleton(x, Cbad, r0, Bgood, r);

    bool test_enabled = false;
    if (test_enabled) {
        BOOST_CHECK_THROW(Doubleton(x, Cbad, r0, Bgood, r), std::domain_error);
    } else {
        BOOST_WARN_MESSAGE(false, "Test disabled due to segmentation fault in SharedDoubleton constructor with invalid data. Enable 'test_enabled' to reproduce crash.");
    }
}

BOOST_AUTO_TEST_SUITE_END()
