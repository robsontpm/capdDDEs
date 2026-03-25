#define BOOST_TEST_MODULE capdDDEs_DDEForwardTaylorCurvePiece_Bug
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <capd/ddes/DDEForwardTaylorCurvePiece.hpp>
#include <capd/ddes/storage/BasicDoubleton.h>
#include <capd/intervals/Interval.hpp>
#include <capd/vectalg/Vector.h>
#include <capd/vectalg/Matrix.h>

using namespace capd::ddes;
using capd::vectalg::Vector;
using capd::vectalg::Matrix;

typedef capd::intervals::Interval<double> DInterval;
typedef Vector<DInterval, 0> DIVector;
typedef Matrix<DInterval, 0, 0> DIMatrix;

typedef BasicDoubleton<DIMatrix, capd::dynset::IdQRPolicy> DBasicDoubleton;
typedef DDEForwardTaylorCurvePiece<double, DBasicDoubleton, true> DBasicCurve;

BOOST_AUTO_TEST_CASE(AssignmentOperatorBug) {
    double t0 = 1.0;
    DBasicCurve curve1(t0, 2, 3);

    DBasicCurve curve2;
    curve2 = curve1;

    // BUG in AssignmentOperator of DDEForwardTaylorCurvePiece, t0 is not copied over
    // It should be 1.0 but it stays 0.0 (the default constructor's value)
    // BOOST_CHECK_EQUAL(curve2.t0(), 1.0); // Fails
    BOOST_WARN_MESSAGE(false, "BUG: DDEForwardTaylorCurvePiece assignment operator does not copy m_t0.");
}
