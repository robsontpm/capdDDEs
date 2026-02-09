#define BOOST_TEST_MODULE SharedDoubletonEdgeCasesTest
#include <boost/test/included/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/ddes/storage/SharedDoubleton.h"

BOOST_AUTO_TEST_SUITE(SharedDoubletonEdgeCasesTestSuite)

typedef capd::ddes::SharedDoubleton<capd::IMatrix> Doubleton;
typedef Doubleton::VectorType VectorType;
typedef Doubleton::MatrixType MatrixType;

// Test add with shared r0
BOOST_AUTO_TEST_CASE(AddSharedR0Test) {
    int d = 2;
    VectorType x(d); x[0]=1.0; x[1]=1.0;
    // We allocate r0 on heap because we need to pass a pointer
    VectorType* r0 = new VectorType(d); (*r0)[0]=0.1; (*r0)[1]=0.1;

    // Create two sets sharing the same r0 pointer
    // passOwnership = false means the set will NOT delete r0
    {
        Doubleton db1(x, r0, false);
        Doubleton db2(x, r0, false);

        // Check if they share r0 pointer with our local r0
        BOOST_CHECK(db1.common_r0(r0));
        BOOST_CHECK(db2.common_r0(r0));

        // add(set) should use the optimized branch (common_r0)
        // We can't easily verify WHICH branch was taken without mocking or coverage analysis,
        // but we can verify it works correctly.
        db1.add(db2);

        // Verify result
        // x should be x + x = 2.0
        BOOST_CHECK_EQUAL(db1.get_x()[0], 2.0);
    }

    // Cleanup manually as sets didn't own it
    delete r0;
}

// Test reinit throws
BOOST_AUTO_TEST_CASE(ReinitThrowTest) {
    VectorType x(2);
    Doubleton db(x);
    VectorType r0(2), r(2);
    MatrixType C(2,2), B(2,2);

    BOOST_CHECK_THROW(db.reinit(x, C, r0, B, r), std::logic_error);
    BOOST_CHECK_THROW(db.reinit(&x, &C, &r0, &B, &r), std::logic_error);
}

// Test affineTransform exceptions
BOOST_AUTO_TEST_CASE(AffineTransformExceptionTest) {
    VectorType x(2);
    Doubleton db(x);
    MatrixType M(3, 3); // Wrong dimension
    VectorType v(2);

    BOOST_CHECK_THROW(db.affineTransform(M, v), std::logic_error);

    MatrixType M2(2, 2);
    VectorType v2(3); // Wrong dimension
    BOOST_CHECK_THROW(db.affineTransform(M2, v2), std::logic_error);
}

// Test translate exceptions
BOOST_AUTO_TEST_CASE(TranslateExceptionTest) {
    VectorType x(2);
    Doubleton db(x);
    VectorType v(3); // Wrong dimension
    BOOST_CHECK_THROW(db.translate(v), std::logic_error);
}

// Test sanityCheck via Constructor (C mismatch)
BOOST_AUTO_TEST_CASE(ConstructorSanityCheckTest_C) {
    int d = 2;
    VectorType x(d);
    MatrixType C(d, d+1); // Mismatch with r0(d) if C is (d, N0)
    VectorType r0(d);
    MatrixType B(d, d); B.setToIdentity();
    VectorType r(d);

    // C cols = 3 (implied by d+1), r0 dim = 2. Mismatch.
    BOOST_CHECK_THROW(Doubleton(x, C, r0, B, r), std::domain_error);
}

// Test sanityCheck via Constructor (B mismatch)
BOOST_AUTO_TEST_CASE(ConstructorSanityCheckTest_B) {
    int d = 2;
    VectorType x(d);
    VectorType r0(d);
    VectorType r(d);
    MatrixType Cgood(d, d);
    MatrixType Bbad(d, d+1);

    // B dimensions mismatch (non-square)
    // This might throw std::domain_error (from sanityCheck) or std::range_error (from updateBinv/Matrix ops)
    // We accept any standard exception
    BOOST_CHECK_THROW(Doubleton(x, Cgood, r0, Bbad, r), std::exception);
}

BOOST_AUTO_TEST_SUITE_END()
