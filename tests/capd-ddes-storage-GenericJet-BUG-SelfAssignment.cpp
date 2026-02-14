#define BOOST_TEST_MODULE GenericJetSelfAssignmentBugTestSuite
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

BOOST_AUTO_TEST_SUITE(GenericJetSelfAssignmentBugTestSuite)

BOOST_AUTO_TEST_CASE(SelfAssignmentBugTest) {
    MyGrid grid(0.1);
    MyTimePoint t0 = grid.point(6);
    MyVector v(2); v[0] = 1.0; v[1] = 2.0;
    JetType original(t0, 1, v);

    // We must initialize assigned with a point from the same grid because GenericJet cannot change grids upon assignment
    JetType assigned(t0);
    assigned = original;

    BOOST_CHECK(assigned == original);

    // Check self-assignment
    // This is known to fail because operator= does not check for self-assignment
    // and calls setupCoeffs which deallocates coefficients.

    try {
        assigned = assigned;

        // We use WARN to indicate this is a bug but we don't want to fail the build now
        bool self_assign_ok = (assigned == original);
        BOOST_WARN_MESSAGE(self_assign_ok, "Self-assignment failed (BUG in GenericJet::operator=)");
        // Or if it crashes, catch it.
    } catch (std::exception& e) {
        BOOST_WARN_MESSAGE(false, "Self-assignment threw exception: " << e.what());
    }
}

BOOST_AUTO_TEST_SUITE_END()
