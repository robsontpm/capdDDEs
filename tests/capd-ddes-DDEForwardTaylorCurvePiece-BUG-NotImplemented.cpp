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

BOOST_AUTO_TEST_CASE(NotImplementedExceptions) {
    DBasicCurve curve;

    // dt() is marked "Not implemented yet" and throws
    BOOST_CHECK_THROW(curve.dt(1), std::logic_error);

    // reinitialize() is marked "Not Supported Yet" and throws
    BOOST_CHECK_THROW(curve.reinitialize(1, 1), std::logic_error);
}
