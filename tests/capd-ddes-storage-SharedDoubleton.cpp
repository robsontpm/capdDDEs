#define BOOST_TEST_MODULE SharedDoubletonTestSuite
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/ddes/storage/SharedDoubleton.h"
#include <exception>

BOOST_AUTO_TEST_SUITE(SharedDoubletonTestSuite)

typedef capd::ddes::SharedDoubleton<capd::IMatrix> Doubleton;
typedef Doubleton::VectorType VectorType;
typedef Doubleton::MatrixType MatrixType;
typedef Doubleton::size_type size_type;
typedef Doubleton::ScalarType ScalarType;

// ==========================================
// Basic Constructor Tests
// ==========================================

BOOST_AUTO_TEST_CASE(DefaultConstructorTest) {
    Doubleton doubleton;
    BOOST_CHECK_EQUAL(doubleton.dimension(), 0);
    BOOST_CHECK_EQUAL(doubleton.storageN0(), 0);
}

BOOST_AUTO_TEST_CASE(ConstructorWithVectorTest) {
    VectorType x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;

    // Test with just x
    Doubleton doubleton(x);
    BOOST_CHECK_EQUAL(doubleton.dimension(), 3);
    BOOST_CHECK_EQUAL(doubleton.storageN0(), 3);
    BOOST_CHECK(doubleton.get_x() == x);

    // Test with x and common_r0
    VectorType r0(2);
    r0[0] = 0.1; r0[1] = 0.2;
    VectorType* r0_ptr = new VectorType(r0);

    // passOwnership = true
    Doubleton doubleton2(x, r0_ptr, true);
    BOOST_CHECK_EQUAL(doubleton2.dimension(), 3);
    BOOST_CHECK_EQUAL(doubleton2.storageN0(), 2);
    BOOST_CHECK(doubleton2.get_r0() == r0);

    // passOwnership = false
    VectorType* r0_ptr2 = new VectorType(r0);
    {
        Doubleton doubleton3(x, r0_ptr2, false);
        BOOST_CHECK_EQUAL(doubleton3.dimension(), 3);
        BOOST_CHECK_EQUAL(doubleton3.storageN0(), 2);
        BOOST_CHECK(doubleton3.get_r0() == r0);
    }
    // r0_ptr2 should still be valid
    BOOST_CHECK(*r0_ptr2 == r0);
    delete r0_ptr2;
}

BOOST_AUTO_TEST_CASE(CopyConstructorTest) {
    VectorType x(3);
    x[0] = 1.0; x[1] = 2.0; x[2] = 3.0;
    Doubleton doubleton1(x);
    Doubleton doubleton2(doubleton1);
    BOOST_CHECK_EQUAL(doubleton2.dimension(), 3);
    BOOST_CHECK(doubleton2.get_x() == x);

    Doubleton d3(x);
    Doubleton d4(d3);
    VectorType new_x(3); new_x[0] = 5.0;
    d3.set_x(new_x);
    BOOST_CHECK(d4.get_x() == x); // d4 should not change (deep copy)
}

BOOST_AUTO_TEST_CASE(AssignmentOperatorTest) {
    VectorType x(3);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    Doubleton doubleton1(x);
    Doubleton doubleton2;
    doubleton2 = doubleton1;
    BOOST_CHECK_EQUAL(doubleton2.dimension(), 3);
    BOOST_CHECK(doubleton2.get_x() == x);
}

BOOST_AUTO_TEST_CASE(DataConstructorTest) {
    int d = 2;
    int N0 = 1;
    VectorType x(d), r0(N0), r(d);
    MatrixType C(d, N0), B(d, d);
    B.setToIdentity();

    x[0] = 1.0; x[1] = 2.0;
    r0[0] = 0.5;
    C[0][0] = 0.1; C[1][0] = 0.2;
    r[0] = 0.01; r[1] = 0.02;

    Doubleton db1(x, C, r0, B, r);
    BOOST_CHECK_EQUAL(db1.dimension(), d);
    BOOST_CHECK_EQUAL(db1.storageN0(), N0);
    BOOST_CHECK(db1.get_x() == x);
    BOOST_CHECK(db1.get_C() == C);
    BOOST_CHECK(db1.get_r0() == r0);
    BOOST_CHECK(db1.get_B() == B);
    BOOST_CHECK(db1.get_r() == r);
}

BOOST_AUTO_TEST_CASE(PointerConstructorTest) {
    int d = 2;
    int N0 = 1;
    VectorType *x = new VectorType(d); (*x)[0] = 1.0;
    MatrixType *C = new MatrixType(d, N0);
    VectorType *r0 = new VectorType(N0);
    MatrixType *B = new MatrixType(d, d); B->setToIdentity();
    VectorType *r = new VectorType(d);

    Doubleton db(x, C, r0, B, r);

    BOOST_CHECK_EQUAL(db.dimension(), d);
    BOOST_CHECK(db.get_x() == *x);

    (*x)[1] = 99.0;
    BOOST_CHECK(db.get_x()[1] == 99.0);

    // Manually cleanup (ownership was OWN_NONE)
    delete x; delete C; delete r0; delete B; delete r;
}

BOOST_AUTO_TEST_CASE(DimensionConstructorTest) {
    Doubleton db1(3, 0);
    BOOST_CHECK_EQUAL(db1.dimension(), 3);
    BOOST_CHECK_EQUAL(db1.storageN0(), 0);

    Doubleton db2(3, 2);
    BOOST_CHECK_EQUAL(db2.dimension(), 3);
    BOOST_CHECK_EQUAL(db2.storageN0(), 2);
}

// ==========================================
// Advanced Constructor Tests (from Constructors.cpp)
// ==========================================

BOOST_AUTO_TEST_CASE(ExtraConstructorsTest) {
    int d = 2;
    VectorType x(d);
    MatrixType C(d, d);
    VectorType r0(d);
    MatrixType B(d, d); B.setToIdentity();
    MatrixType Binv(d, d); Binv.setToIdentity();
    VectorType r(d);

    // Constructor with Binv (const ref)
    Doubleton db6(x, C, r0, B, Binv, r);
    BOOST_CHECK_EQUAL(db6.dimension(), d);

    // Constructor with r0 pointer (no Binv)
    VectorType* r0_ptr = new VectorType(r0);
    {
        Doubleton db7(x, C, r0_ptr, B, r);
        BOOST_CHECK(db7.get_r0() == r0);
        BOOST_CHECK(db7.common_r0(r0_ptr));
    }
    delete r0_ptr;

    // Constructor with r0 pointer and Binv (const ref)
    VectorType* r0_ptr2 = new VectorType(r0);
    {
        Doubleton db8(x, C, r0_ptr2, B, Binv, r);
        BOOST_CHECK(db8.common_r0(r0_ptr2));
    }
    delete r0_ptr2;

    // Constructor with x, C, r0 (const ref)
    Doubleton db9(x, C, r0);
    BOOST_CHECK_EQUAL(db9.dimension(), d);

    // Constructor with x, C, r0 pointer
    VectorType* r0_ptr3 = new VectorType(r0);
    {
        Doubleton db10(x, C, r0_ptr3);
        BOOST_CHECK(db10.common_r0(r0_ptr3));
    }
    delete r0_ptr3;
}

// ==========================================
// Zero Dimension Tests (from ZeroDim.cpp)
// ==========================================

BOOST_AUTO_TEST_CASE(ZeroDimTest) {
    Doubleton db(0, 0);
    BOOST_CHECK_EQUAL(db.dimension(), 0);
    BOOST_CHECK_EQUAL(db.storageN0(), 0);

    VectorType v((size_type)0);
    db.add(v); // Should not crash
    db.mul(2.0); // Should not crash

    VectorType x((size_type)0);
    Doubleton db2(x);
    BOOST_CHECK_EQUAL(db2.dimension(), 0);
    BOOST_CHECK_EQUAL(db2.storageN0(), 0);
}

// ==========================================
// Getter / Setter Tests
// ==========================================

BOOST_AUTO_TEST_CASE(GetterSetterTest) {
    VectorType x(2); x[0]=1; x[1]=2;
    MatrixType C(2, 2); C.setToIdentity();
    VectorType r0(2);
    MatrixType B(2, 2); B.setToIdentity();
    VectorType r(2);

    Doubleton db(2, 2);
    db.set_x(x); BOOST_CHECK(db.get_x() == x);
    db.set_C(C); BOOST_CHECK(db.get_C() == C);
    db.set_r0(r0); BOOST_CHECK(db.get_r0() == r0);
    db.set_B(B); BOOST_CHECK(db.get_B() == B);
    db.set_r(r); BOOST_CHECK(db.get_r() == r);

    // Test set with pointer and ownership transfer
    VectorType* x_ptr = new VectorType(x);
    db.set_x(x_ptr, true);
    BOOST_CHECK(db.get_x() == x);

    // Test set with pointer without ownership
    VectorType* r_ptr = new VectorType(r);
    db.set_r(r_ptr, false);
    BOOST_CHECK(db.get_r() == r);
    db.set_r(r); // Reset to owned copy to safely delete ptr
    delete r_ptr;
}

BOOST_AUTO_TEST_CASE(SettersRefTest) {
    Doubleton db(2, 2);
    MatrixType C(2, 2); C.setToIdentity();
    VectorType r0(2); r0[0]=1.0;

    db.set_Cr0(C, r0);
    BOOST_CHECK(db.get_C() == C);
    BOOST_CHECK(db.get_r0() == r0);

    MatrixType Binv(2, 2); Binv.setToIdentity();
    db.set_Binv(Binv);
    BOOST_CHECK(db.get_Binv() == Binv);
}

BOOST_AUTO_TEST_CASE(SettersPointerTest) {
    Doubleton db(2, 2);
    VectorType* x = new VectorType(2);
    MatrixType* C = new MatrixType(2, 2);
    VectorType* r0 = new VectorType(2);
    MatrixType* B = new MatrixType(2, 2); B->setToIdentity();
    VectorType* r = new VectorType(2);
    MatrixType* Binv = new MatrixType(2, 2); Binv->setToIdentity();

    db.set_x(x, true);
    db.set_C(C, true);
    db.set_r0(r0, true);
    db.set_B(B, true);
    db.set_r(r, true);
    db.set_Binv(Binv, true);

    BOOST_CHECK(db.common_x(x));
    BOOST_CHECK(db.common_C(C));
    BOOST_CHECK(db.common_r0(r0));
    BOOST_CHECK(db.common_B(B));
    BOOST_CHECK(db.common_r(r));
    BOOST_CHECK(db.common_Binv(Binv));

    // Test set_Cr0 pointers
    Doubleton db2(2, 2);
    MatrixType* C2 = new MatrixType(2, 2);
    VectorType* r02 = new VectorType(2);
    db2.set_Cr0(C2, r02, true, true);
    BOOST_CHECK(db2.common_C(C2));
    BOOST_CHECK(db2.common_r0(r02));
}

BOOST_AUTO_TEST_CASE(TakeMethodsTest) {
    VectorType x(2); x[0]=10;
    Doubleton db(x);
    VectorType* x_ptr = db.take_x();
    BOOST_CHECK(x_ptr != nullptr);
    BOOST_CHECK(*x_ptr == x);
    BOOST_CHECK(db.get_x() == x);
    delete x_ptr;

    // Test other take methods
    Doubleton db2(2, 2);
    MatrixType* C = db2.take_C(); delete C;
    VectorType* r0 = db2.take_r0(); delete r0;
    MatrixType* B = db2.take_B(); delete B;
    VectorType* r = db2.take_r(); delete r;
    MatrixType* Binv = db2.take_Binv(); delete Binv;
    VectorType* x2 = db2.take_x(); delete x2;
}

BOOST_AUTO_TEST_CASE(HullMidPointTest) {
    int d = 2;
    VectorType x(d); x[0]=1.0; x[1]=2.0;
    VectorType r0(d); r0[0]=0.1; r0[1]=0.1;
    MatrixType C(d,d); C.setToIdentity();
    VectorType r(d); r[0]=0.01; r[1]=0.01;
    MatrixType B(d,d); B.setToIdentity();

    Doubleton db(x, C, r0, B, r);

    BOOST_CHECK(db.midPoint() == x);

    VectorType h = db.hull();
    BOOST_CHECK(h.dimension() == d);
}

BOOST_AUTO_TEST_CASE(UpdateBinvTest) {
    int d = 2;
    VectorType x(d); VectorType r0(d); VectorType r(d); r[0]=1.0; r[1]=1.0;
    MatrixType C(d,d); MatrixType B(d,d); B.setToIdentity();
    Doubleton db(x, C, r0, B, r);

    MatrixType B2(d,d); B2[0][0]=2.0; B2[1][1]=2.0;
    db.set_B(B2);

    // IdQRPolicy normalization check
    VectorType r_final = db.get_r();
    VectorType expected(d); expected[0]=2.0; expected[1]=2.0;
    BOOST_CHECK(r_final == expected);
}

// ==========================================
// Logic Tests (Add, Mul, Transform)
// ==========================================

BOOST_AUTO_TEST_CASE(AddVectorTest) {
    VectorType x(2); x[0]=1.0; x[1]=2.0;
    Doubleton db(x);
    VectorType v(2); v[0]=0.5; v[1]=0.5;
    db.add(v);
    BOOST_CHECK(db.get_x()[0] == 1.5);
    BOOST_CHECK(db.get_x()[1] == 2.5);
}

BOOST_AUTO_TEST_CASE(AddSetTest) {
    VectorType x1(2); x1[0]=1.0; x1[1]=2.0;
    Doubleton db1(x1);
    VectorType x2(2); x2[0]=0.5; x2[1]=0.5;
    Doubleton db2(x2);

    db1.add(db2);

    BOOST_CHECK(db1.get_x()[0] == 1.5);
    BOOST_CHECK(db1.get_x()[1] == 2.5);
}

BOOST_AUTO_TEST_CASE(AddSharedR0Test) {
    int d = 2;
    VectorType x(d); x[0]=1.0; x[1]=1.0;
    VectorType* r0 = new VectorType(d); (*r0)[0]=0.1; (*r0)[1]=0.1;

    {
        Doubleton db1(x, r0, false);
        Doubleton db2(x, r0, false);
        BOOST_CHECK(db1.common_r0(r0));
        BOOST_CHECK(db2.common_r0(r0));

        db1.add(db2);
        BOOST_CHECK_EQUAL(db1.get_x()[0], 2.0);
    }
    delete r0;
}

BOOST_AUTO_TEST_CASE(MulTest) {
    VectorType x(2); x[0]=1.0; x[1]=2.0;
    Doubleton db(x);
    ScalarType c = 2.0;
    db.mul(c);
    BOOST_CHECK(db.get_x()[0] == 2.0);
    BOOST_CHECK(db.get_x()[1] == 4.0);
}

BOOST_AUTO_TEST_CASE(MulThenAddTest) {
    VectorType x1(2); x1[0]=1.0; x1[1]=2.0;
    Doubleton db1(x1);
    VectorType x2(2); x2[0]=0.5; x2[1]=0.5;
    Doubleton db2(x2);
    ScalarType c = 2.0;

    db1.mulThenAdd(c, db2);

    // (1.0 * 2) + 0.5 = 2.5
    // (2.0 * 2) + 0.5 = 4.5
    BOOST_CHECK(db1.get_x()[0] == 2.5);
    BOOST_CHECK(db1.get_x()[1] == 4.5);
}

BOOST_AUTO_TEST_CASE(AffineTransformTest) {
    VectorType x(2); x[0] = 1.0; x[1] = 2.0;
    MatrixType M(2, 2); M[0][0] = 2.0; M[1][1] = 2.0;
    VectorType v(2); v[0] = 1.0; v[1] = 1.0;

    Doubleton doubleton(x);
    // x - v = (0, 1), M*(x-v) = (0, 2)
    doubleton.affineTransform(M, v);
    VectorType expected(2); expected[0] = 0.0; expected[1] = 2.0;
    BOOST_CHECK(doubleton.get_x() == expected);

    // Test Exceptions
    MatrixType M_bad(3, 3);
    BOOST_CHECK_THROW(doubleton.affineTransform(M_bad, v), std::logic_error);
    VectorType v_bad(3);
    BOOST_CHECK_THROW(doubleton.affineTransform(M, v_bad), std::logic_error);
}

BOOST_AUTO_TEST_CASE(TranslateTest) {
    VectorType x(2); x[0] = 1.0; x[1] = 2.0;
    VectorType v(2); v[0] = 1.0; v[1] = 1.0;
    Doubleton doubleton(x);
    doubleton.translate(v);
    VectorType expected(2); expected[0] = 2.0; expected[1] = 3.0;
    BOOST_CHECK(doubleton.get_x() == expected);

    VectorType bad_v(3);
    BOOST_CHECK_THROW(doubleton.translate(bad_v), std::logic_error);
}

BOOST_AUTO_TEST_CASE(ShowTest) {
    VectorType x(3);
    Doubleton doubleton(x);
    std::string result = doubleton.show();
    BOOST_CHECK(!result.empty());
}

// ==========================================
// Exception & Edge Case Tests
// ==========================================

BOOST_AUTO_TEST_CASE(ReinitializeTest) {
    VectorType x(3);
    Doubleton doubleton(x);
    doubleton.reinitialize(2, 2);
    BOOST_CHECK_EQUAL(doubleton.dimension(), 2);
    BOOST_CHECK_EQUAL(doubleton.storageN0(), 2);

    VectorType r0(2), r(2); MatrixType C(2,2), B(2,2);
    BOOST_CHECK_THROW(doubleton.reinit(x, C, r0, B, r), std::logic_error);
    BOOST_CHECK_THROW(doubleton.reinit(&x, &C, &r0, &B, &r), std::logic_error);
}

BOOST_AUTO_TEST_CASE(SanityCheckTest) {
    int d = 2; int N0 = 1;

    {
        Doubleton db(d, N0);
        BOOST_CHECK_THROW(db.set_x(nullptr, false), std::domain_error);
    }
    {
        Doubleton db(d, N0);
        BOOST_CHECK_THROW(db.set_C(nullptr, false), std::domain_error);
    }
    {
        Doubleton db(d, N0);
        BOOST_CHECK_THROW(db.set_B(nullptr, false), std::logic_error);
    }
    {
        Doubleton db(d, N0);
        BOOST_CHECK_THROW(db.set_r(nullptr, false), std::domain_error);
    }
    {
        Doubleton db(d, N0);
        MatrixType B_bad(d+1, d);
        BOOST_CHECK_THROW(db.set_B(B_bad), std::exception);
    }
    {
        Doubleton db(d, N0);
        MatrixType Binv_bad(d+1, d);
        BOOST_CHECK_THROW(db.set_Binv(Binv_bad), std::domain_error);
    }
    {
        Doubleton db(d, N0);
        MatrixType C_bad(d+1, N0);
        BOOST_CHECK_THROW(db.set_C(C_bad), std::domain_error);
    }
    {
        Doubleton db(d, N0);
        VectorType r_bad(d+1);
        BOOST_CHECK_THROW(db.set_r(r_bad), std::domain_error);
    }

}

BOOST_AUTO_TEST_CASE(AssureOwnerSetXTest) {
    int d = 2;
    VectorType* x = new VectorType(d); (*x)[0] = 1.0;
    MatrixType* C = new MatrixType(d, d);
    VectorType* r0 = new VectorType(d);
    MatrixType* B = new MatrixType(d, d); B->setToIdentity();
    VectorType* r = new VectorType(d);

    Doubleton db(x, C, r0, B, r);
    BOOST_CHECK(db.common_x(x));

    VectorType new_x(d); new_x[0] = 2.0;
    db.set_x(new_x); // Should trigger assureOwner

    BOOST_CHECK(!db.common_x(x));
    BOOST_CHECK(db.get_x() == new_x);
    BOOST_CHECK((*x)[0] == 1.0);

    delete x; delete C; delete r0; delete B; delete r;
}

BOOST_AUTO_TEST_CASE(AddCoverageTest) {
    VectorType x(2);
    Doubleton db(x);
    VectorType v(3);
    BOOST_CHECK_THROW(db.add(v), std::logic_error);

    Doubleton db_bad(3, 3);
    BOOST_CHECK_THROW(db.add(db_bad), std::logic_error);

    Doubleton db_bad_N0(2, 5);
    BOOST_CHECK_THROW(db.add(db_bad_N0), std::logic_error);
}

BOOST_AUTO_TEST_SUITE_END()
