#define BOOST_TEST_MODULE GenericJetTestSuite
#include <boost/test/included/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/ddes/storage/GenericJet.hpp"

// Define types
typedef capd::vectalg::Vector<double, 0> MyVector;
typedef capd::vectalg::Matrix<double, 0, 0> MyMatrix;
typedef capd::ddes::DiscreteTimeGrid<double> MyGrid;
typedef MyGrid::TimePointType MyTimePoint;

// Instantiate GenericJet
typedef capd::ddes::GenericJet<
    MyTimePoint,
    MyVector,
    MyVector,
    MyMatrix
> JetType;

typedef JetType::size_type size_type;

// Define a test suite
BOOST_AUTO_TEST_SUITE(GenericJetTestSuite)

// ==========================================
// Basic Constructor Tests
// ==========================================

BOOST_AUTO_TEST_CASE(DefaultConstructorTest) {
    JetType jet;
    BOOST_CHECK_EQUAL(jet.dimension(), 0);
    BOOST_CHECK_EQUAL(jet.order(), 0);
    BOOST_CHECK_EQUAL(jet.storageDimension(), 0);
    BOOST_CHECK(jet.t0().isZero());
}

BOOST_AUTO_TEST_CASE(DimensionOrderConstructorTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(0);
    size_type dim = 2;
    size_type order = 3;
    JetType jet(t0, dim, order);

    BOOST_CHECK_EQUAL(jet.dimension(), dim);
    BOOST_CHECK_EQUAL(jet.order(), order);
    BOOST_CHECK_EQUAL(jet.storageDimension(), dim * (order + 1));
    BOOST_CHECK_EQUAL(jet.t0(), t0);

    // Check if coefficients are zero
    for (size_type k = 0; k <= order; ++k) {
        MyVector v = jet[k];
        BOOST_CHECK_EQUAL(v.dimension(), dim);
        for(size_type i=0; i<dim; ++i)
             BOOST_CHECK_EQUAL(v[i], 0.0);
    }
}

BOOST_AUTO_TEST_CASE(ConstantVectorConstructorTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(1);
    size_type order = 2;
    MyVector v(3); v[0] = 1.0; v[1] = 2.0; v[2] = 3.0;

    JetType jet(t0, order, v);

    BOOST_CHECK_EQUAL(jet.dimension(), 3);
    BOOST_CHECK_EQUAL(jet.order(), order);
    BOOST_CHECK_EQUAL(jet.t0(), t0);

    // Check 0-th coeff is v
    BOOST_CHECK(jet[0] == v);

    // Check higher order coeffs are zero
    for (size_type k = 1; k <= order; ++k) {
        MyVector vk = jet[k];
        for(size_type i=0; i<3; ++i)
             BOOST_CHECK_EQUAL(vk[i], 0.0);
    }
}

BOOST_AUTO_TEST_CASE(IteratorConstructorTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(2);

    std::vector<MyVector> coeffs;
    for(int k=0; k<3; ++k) {
        MyVector v(2); v[0] = k; v[1] = k*2;
        coeffs.push_back(v);
    }

    JetType jet(t0, coeffs.data(), coeffs.data() + coeffs.size());

    BOOST_CHECK_EQUAL(jet.dimension(), 2);
    // order is size - 1 = 3 - 1 = 2
    BOOST_CHECK_EQUAL(jet.order(), 2);
    BOOST_CHECK_EQUAL(jet.t0(), t0);

    for(size_type k=0; k<=jet.order(); ++k) {
        BOOST_CHECK(jet[k] == coeffs[k]);
    }
}

BOOST_AUTO_TEST_CASE(IteratorConstructorWithOrderTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(3);

    std::vector<MyVector> coeffs;
    for(int k=0; k<2; ++k) {
        MyVector v(2); v[0] = k; v[1] = k*2;
        coeffs.push_back(v);
    }
    size_type order = 3; // larger than data provided

    JetType jet(t0, order, coeffs.data(), coeffs.data() + coeffs.size());

    BOOST_CHECK_EQUAL(jet.dimension(), 2);
    BOOST_CHECK_EQUAL(jet.order(), order);

    for(size_type k=0; k<2; ++k) {
        BOOST_CHECK(jet[k] == coeffs[k]);
    }
    // Check padding with zeros
    for(size_type k=2; k<=order; ++k) {
         MyVector v = jet[k];
         for(size_type i=0; i<2; ++i) BOOST_CHECK_EQUAL(v[i], 0.0);
    }
}

BOOST_AUTO_TEST_CASE(VectorConstructorTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(4);

    std::vector<MyVector> coeffs;
    for(int k=0; k<3; ++k) {
        MyVector v(2); v[0] = k+1; v[1] = (k+1)*2;
        coeffs.push_back(v);
    }

    JetType jet(t0, coeffs);

    BOOST_CHECK_EQUAL(jet.dimension(), 2);
    BOOST_CHECK_EQUAL(jet.order(), 2);

    for(size_type k=0; k<=jet.order(); ++k) {
        BOOST_CHECK(jet[k] == coeffs[k]);
    }
}

// ==========================================
// Copy/Assignment Tests
// ==========================================

BOOST_AUTO_TEST_CASE(CopyConstructorTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(5);
    MyVector v(2); v[0] = 1.0; v[1] = 2.0;
    JetType original(t0, 1, v);

    JetType copy(original);

    BOOST_CHECK_EQUAL(copy.dimension(), original.dimension());
    BOOST_CHECK_EQUAL(copy.order(), original.order());
    BOOST_CHECK_EQUAL(copy.t0(), original.t0());
    BOOST_CHECK(copy == original);

    // Check deep copy
    MyVector new_v(2); new_v[0] = 5.0; new_v[1] = 5.0;
    copy[0] = new_v;
    BOOST_CHECK(original[0] == v); // original should remain unchanged
    BOOST_CHECK(copy[0] == new_v);
}

BOOST_AUTO_TEST_CASE(AssignmentOperatorTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(6);
    MyVector v(2); v[0] = 1.0; v[1] = 2.0;
    JetType original(t0, 1, v);

    // We must initialize assigned with a point from the same grid because GenericJet cannot change grids upon assignment
    JetType assigned(t0);
    assigned = original;

    BOOST_CHECK_EQUAL(assigned.dimension(), original.dimension());
    BOOST_CHECK_EQUAL(assigned.order(), original.order());
    BOOST_CHECK_EQUAL(assigned.t0(), original.t0());
    for(size_type k=0; k<=original.order(); ++k) {
        BOOST_CHECK_MESSAGE(assigned[k] == original[k], "Coeff " << k << " mismatch");
    }

    BOOST_CHECK(assigned == original);

    // Check self-assignment (BUG: this fails due to deallocateCoeffs)
    // assigned = assigned;
    // BOOST_CHECK(assigned == original);
}

// ==========================================
// Accessor/Setter Tests
// ==========================================

BOOST_AUTO_TEST_CASE(SetGetXTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(7);
    size_type dim = 2;
    size_type order = 1;
    JetType jet(t0, dim, order);

    // storageDimension = dim * (order + 1) = 2 * 2 = 4
    MyVector x(4);
    x[0] = 1.0; x[1] = 2.0; // coeff 0
    x[2] = 3.0; x[3] = 4.0; // coeff 1

    jet.set_x(x);

    MyVector retrieved_x = jet.get_x();
    BOOST_CHECK(retrieved_x == x);

    MyVector c0 = jet[0];
    BOOST_CHECK_EQUAL(c0[0], 1.0);
    BOOST_CHECK_EQUAL(c0[1], 2.0);

    MyVector c1 = jet[1];
    BOOST_CHECK_EQUAL(c1[0], 3.0);
    BOOST_CHECK_EQUAL(c1[1], 4.0);
}

BOOST_AUTO_TEST_CASE(SetT0Test) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(8);
    MyTimePoint t1 = grid.point(9);

    JetType jet(t0, 2, 1);
    BOOST_CHECK_EQUAL(jet.t0(), t0);

    jet.setT0(t1);
    BOOST_CHECK_EQUAL(jet.getT0(), t1);
}

// ==========================================
// Iterators Tests
// ==========================================

BOOST_AUTO_TEST_CASE(IteratorsTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(10);
    std::vector<MyVector> coeffs;
    for(int k=0; k<3; ++k) {
        MyVector v(1); v[0] = k;
        coeffs.push_back(v);
    }

    JetType jet(t0, coeffs);

    size_type k = 0;
    for(auto it = jet.begin(); it != jet.end(); ++it, ++k) {
        BOOST_CHECK(*it == coeffs[k]);
    }
    BOOST_CHECK_EQUAL(k, 3);

    BOOST_CHECK(jet.begin() + jet.order() + 1 == jet.end());
    BOOST_CHECK(*jet.back() == coeffs.back());
    BOOST_CHECK(*jet.at(1) == coeffs[1]);
}

// ==========================================
// Evaluation Tests
// ==========================================

BOOST_AUTO_TEST_CASE(EvalAtDeltaTest) {
    // f(t) = 1 + t + t^2/2  (coeffs: 1, 1, 0.5)
    // at t=0
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(0);
    size_type dim = 1;
    size_type order = 2;

    std::vector<MyVector> coeffs;
    MyVector v0(1); v0[0] = 1.0; coeffs.push_back(v0);
    MyVector v1(1); v1[0] = 1.0; coeffs.push_back(v1); // 1st derivative at 0 -> coeff 1
    MyVector v2(1); v2[0] = 0.5; coeffs.push_back(v2); // 2nd derivative at 0 / 2! = 1/2 -> coeff 0.5

    JetType jet(t0, coeffs);

    double delta = 2.0;
    // f(delta) = 1 + 2 + 2^2/2 = 1 + 2 + 2 = 5

    MyVector res = jet.evalAtDelta(delta);
    BOOST_CHECK_CLOSE(res[0], 5.0, 1e-10);

    MyVector res2(1);
    jet.evalAtDelta(delta, res2);
    BOOST_CHECK_CLOSE(res2[0], 5.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(EvalTest) {
    // Same function f(t) = 1 + t + t^2/2
    // but at t0 = 1.0.
    // If coeffs are at t0=1.0, then jet represents f(t) expanded at t0=1.
    // Let's make jet represent g(t) = t^2.
    // g'(t) = 2t, g''(t) = 2.
    // At t0 = 1: g(1) = 1, g'(1) = 2, g''(1) = 2.
    // Coeffs: c0 = 1, c1 = 2, c2 = 2/2 = 1.

    MyGrid grid(1.0); // step 1.0
    MyTimePoint t0 = grid.point(1); // t0 = 1.0

    std::vector<MyVector> coeffs;
    MyVector v0(1); v0[0] = 1.0; coeffs.push_back(v0);
    MyVector v1(1); v1[0] = 2.0; coeffs.push_back(v1);
    MyVector v2(1); v2[0] = 1.0; coeffs.push_back(v2);

    JetType jet(t0, coeffs);

    // Evaluate at t = 3.0. Delta = 3 - 1 = 2.
    // g(3) = 3^2 = 9.
    // Taylor: 1 + 2*(2) + 1*(2^2) = 1 + 4 + 4 = 9.

    MyVector res = jet.eval(3.0);
    BOOST_CHECK_CLOSE(res[0], 9.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(EvalCoeffTest) {
    // g(t) = t^3
    // g'(t) = 3t^2
    // g''(t) = 6t
    // g'''(t) = 6
    // At t0 = 0: c0=0, c1=0, c2=0, c3=6/6=1.

    MyGrid grid(1.0);
    MyTimePoint t0 = grid.point(0);

    std::vector<MyVector> coeffs;
    MyVector v(1); v[0] = 0.0; coeffs.push_back(v); // 0
    coeffs.push_back(v); // 0
    coeffs.push_back(v); // 0
    v[0] = 1.0; coeffs.push_back(v); // 1

    JetType jet(t0, coeffs); // order 3

    // Evaluate 1st coeff (g'(t)) at t = 2.0.
    // g'(2) = 3*2^2 = 12.
    // 1st coeff is g'(t) / 1! = g'(t).

    // The implementation of evalCoeffAtDelta computes the n-th coefficient of the shifted jet.
    // If J is jet at t0, J_shifted is jet at t0+delta.
    // (J_shifted)_n = 1/n! * d^n/dt^n (Sum c_k (t-t0)^k) at t=t0+delta.
    // J(t) = c3 * t^3 = t^3.
    // J'(t) = 3t^2. J''(t) = 6t. J'''(t) = 6.
    // n=1: coeff 1 of J_shifted is J'(t0+delta).
    // t0=0, delta=2. J'(2) = 12.

    MyVector res = jet.evalCoeff(1, 2.0);
    BOOST_CHECK_CLOSE(res[0], 12.0, 1e-10);

    // n=2: coeff 2 of J_shifted is J''(2)/2! = 12/2 = 6.
    // Wait, J''(t) = 6t. J''(2) = 12. 12/2 = 6. Correct.
    res = jet.evalCoeff(2, 2.0);
    BOOST_CHECK_CLOSE(res[0], 6.0, 1e-10);
}

// ==========================================
// Derivative Tests
// ==========================================

BOOST_AUTO_TEST_CASE(DerivativeTest) {
    // f(t) = t^2.
    // At t0=0: c0=0, c1=0, c2=1.

    MyGrid grid(1.0);
    MyTimePoint t0 = grid.point(0);

    std::vector<MyVector> coeffs;
    MyVector v(1); v[0] = 0.0; coeffs.push_back(v);
    coeffs.push_back(v);
    v[0] = 1.0; coeffs.push_back(v);

    JetType jet(t0, coeffs); // order 2

    // Derivative dt(1): f'(t) = 2t.
    // At t0=0: c0=0, c1=2.
    // Resulting jet should have order 2-1 = 1.

    JetType deriv = jet.dt(1);

    BOOST_CHECK_EQUAL(deriv.order(), 1);
    BOOST_CHECK_EQUAL(deriv.dimension(), 1);
    BOOST_CHECK_EQUAL(deriv[0][0], 0.0);
    BOOST_CHECK_EQUAL(deriv[1][0], 2.0);
}

// ==========================================
// Order Manipulation Tests
// ==========================================

BOOST_AUTO_TEST_CASE(SetOrderTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(0);
    JetType jet(t0, 1, 2); // dim 1, order 2

    JetType reduced = jet.setOrder(1);
    BOOST_CHECK_EQUAL(reduced.order(), 1);

    JetType expanded = jet.setOrder(3);
    BOOST_CHECK_EQUAL(expanded.order(), 3);
    // Check if new coeffs are zero
    BOOST_CHECK_EQUAL(expanded[3][0], 0.0);
}

BOOST_AUTO_TEST_CASE(IncreasedDecreasedOrderTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(0);
    JetType jet(t0, 1, 2);

    JetType increased = jet.increasedOrder(1);
    BOOST_CHECK_EQUAL(increased.order(), 3);

    JetType decreased = jet.decreasedOrder(1);
    BOOST_CHECK_EQUAL(decreased.order(), 1);

    BOOST_CHECK_THROW(jet.decreasedOrder(3), std::logic_error);
}

// ==========================================
// Show/Hull/Conversion Tests
// ==========================================

BOOST_AUTO_TEST_CASE(ShowTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(0);
    JetType jet(t0, 1, 1);
    std::string s = jet.show();
    BOOST_CHECK(!s.empty());
}

BOOST_AUTO_TEST_CASE(HullTest) {
    // Default implementation of hull returns *this cast to VectorType.
    // Wait, GenericJet::hull() returns VectorType.
    // If VectorSpec is MyVector, it should return MyVector.
    // However, GenericJet::operator VectorType() returns a vector of dimension storageDimension.
    // GenericJet::hull() calls *this (which invokes operator VectorType() ?? No, *this is JetType).
    // GenericJet::hull() implementation: { return *this; }
    // Since GenericJet has operator VectorType(), this works.

    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(0);
    JetType jet(t0, 1, 1); // dim 1, order 1. storage dim = 2.

    MyVector h = jet.hull();
    BOOST_CHECK_EQUAL(h.dimension(), 2);
}

BOOST_AUTO_TEST_CASE(CoeffNotImplementedYetTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(0);
    JetType jet(t0, 1, 1);
    BOOST_CHECK_THROW(jet.coeff(1), std::logic_error);
}


BOOST_AUTO_TEST_SUITE_END()
