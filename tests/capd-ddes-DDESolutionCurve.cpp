#define BOOST_TEST_MODULE DDESolutionCurveTests
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <capd/ddes/DDESolutionCurve.h>
#include <capd/ddes/DDESolutionCurve.hpp>
#include <capd/ddes/DDECommon.h>
#include <capd/ddes/storage/SharedDoubleton.h>
#include <capd/ddes/storage/GenericJet.h>
#include <capd/ddes/storage/GenericJet.hpp>
#include <capd/capdlib.h>

using namespace capd::ddes;

typedef capd::DInterval Interval;
typedef capd::IMatrix IMatrix;
typedef capd::ddes::SharedDoubleton<IMatrix> DDESet;
typedef DiscreteTimeGrid<Interval> GridType;
typedef GridType::TimePointType TimePointType;

BOOST_AUTO_TEST_SUITE(DDESolutionCurve_Suite)

BOOST_AUTO_TEST_CASE(Constructors_Test) {
    GridType grid(Interval(0.1));

    // 1. DDESolutionCurve(grid, d, N0)
    DDESolutionCurve<DDESet> curve1(grid, 2, 3);
    BOOST_CHECK_EQUAL(curve1.dimension(), 2);
    BOOST_CHECK_EQUAL(curve1.storageN0(), 3);
    BOOST_CHECK_EQUAL(curve1.length(), 0);
    BOOST_CHECK(curve1.t0() == grid.point(0));

    // 2. DDESolutionCurve(t0, d, N0)
    DDESolutionCurve<DDESet> curve2(grid.point(5), 2, 3);
    BOOST_CHECK_EQUAL(curve2.dimension(), 2);
    BOOST_CHECK_EQUAL(curve2.storageN0(), 3);
    BOOST_CHECK(curve2.t0() == grid.point(5));
}

BOOST_AUTO_TEST_CASE(Constructors_Vector_Set_Test) {
    GridType grid(Interval(0.1));
    TimePointType t0 = grid.point(0);
    TimePointType t1 = grid.point(5); // t1 = 0.5

    // 3. DDESolutionCurve(t0, const SetType& value)
    DDESet set1(2, 3);
    capd::IVector set1x(2);
    set1x[0]=1.0; set1x[1]=2.0;
    set1.set_x(set1x);
    DDESolutionCurve<DDESet> curve3(t0, set1);
    BOOST_CHECK_EQUAL(curve3.dimension(), 2);
    BOOST_CHECK_EQUAL(curve3.storageN0(), 3);
    BOOST_CHECK_EQUAL(curve3.length(), 0);

    // 4. DDESolutionCurve(t0, const VectorType& value)
    capd::IVector vec1(3);
    vec1[0]=3.0; vec1[1]=4.0; vec1[2]=5.0;
    DDESolutionCurve<DDESet> curve4(t0, vec1);
    BOOST_CHECK_EQUAL(curve4.dimension(), 3);
    BOOST_CHECK_EQUAL(curve4.length(), 0);

    // 5. DDESolutionCurve(t0, t1, order, const SetType& value)
    DDESolutionCurve<DDESet> curve5(t0, t1, 2, set1);
    BOOST_CHECK_EQUAL(curve5.length(), 5);
    BOOST_CHECK_EQUAL(curve5.dimension(), 2);

    // 6. DDESolutionCurve(t0, t1, order, SetType& value) -> takes r0
    DDESet set2(2, 4);
    capd::IVector set2x(2);
    set2x[0]=5.0; set2x[1]=6.0;
    set2.set_x(set2x);
    DDESolutionCurve<DDESet> curve6(t0, t1, 3, set2);
    BOOST_CHECK_EQUAL(curve6.length(), 5);
    BOOST_CHECK_EQUAL(curve6.storageN0(), 4);

    // 7. DDESolutionCurve(t0, t1, order, const VectorType& value, N0)
    DDESolutionCurve<DDESet> curve7(t0, t1, 4, vec1, 5);
    BOOST_CHECK_EQUAL(curve7.length(), 5);
    BOOST_CHECK_EQUAL(curve7.storageN0(), 5);
}

BOOST_AUTO_TEST_CASE(Assignment_Copy_Test) {
    GridType grid(Interval(0.1));
    TimePointType t0 = grid.point(0);
    TimePointType t1 = grid.point(5); // t1 = 0.5

    capd::IVector vec1(3);
    vec1[0]=3.0; vec1[1]=4.0; vec1[2]=5.0;
    DDESolutionCurve<DDESet> orig(t0, t1, 4, vec1, 5);

    // Copy constructor
    DDESolutionCurve<DDESet> copy_curve(orig);
    BOOST_CHECK_EQUAL(copy_curve.length(), 5);
    BOOST_CHECK_EQUAL(copy_curve.dimension(), 3);
    BOOST_CHECK_EQUAL(copy_curve.storageN0(), 5);

    // Assignment operator
    DDESolutionCurve<DDESet> assign_curve(t0, 1, 1);
    assign_curve = orig;
    BOOST_CHECK_EQUAL(assign_curve.length(), 5);
    BOOST_CHECK_EQUAL(assign_curve.dimension(), 3);
    BOOST_CHECK_EQUAL(assign_curve.storageN0(), 5);

    // Assignment over different grid
    GridType grid2(Interval(0.2));
    DDESolutionCurve<DDESet> bad_grid_curve(grid2.point(0), 1, 1);
    BOOST_CHECK_THROW(bad_grid_curve = orig, std::logic_error);
}

BOOST_AUTO_TEST_CASE(Evaluation_Test) {
    GridType grid(Interval(0.1));
    TimePointType t0 = grid.point(0);
    TimePointType t1 = grid.point(5); // t1 = 0.5

    capd::IVector vec1(3);
    vec1[0]=3.0; vec1[1]=4.0; vec1[2]=5.0;
    DDESolutionCurve<DDESet> curve(t0, t1, 4, vec1, 5);

    // Test j() accessors
    BOOST_CHECK_EQUAL(curve.j(grid.point(0)).order(), 4);
    BOOST_CHECK_EQUAL(curve.j(grid.point(1)).order(), 4);
    BOOST_CHECK_EQUAL(curve.jetOrderAt(grid.point(2)), 4);

    BOOST_CHECK_THROW(curve.j(grid.point(10)), std::domain_error); // Out of bounds
}

BOOST_AUTO_TEST_CASE(Doubleton_Interface_Impl_Test) {
    GridType grid(Interval(0.1));
    TimePointType t0 = grid.point(0);
    TimePointType t1 = grid.point(2);

    capd::IVector vec1(2);
    vec1[0]=3.0; vec1[1]=4.0; // d=2
    DDESolutionCurve<DDESet> curve(t0, t1, 1, vec1, 3); // N0=3, length=2, order=1

    // storageDimension = current_value.dim + sum(pieces.storageDim)
    // piece.storageDim = d * (order + 1) = 2 * 2 = 4
    // total = 2 + 4 + 4 = 10
    BOOST_CHECK_EQUAL(curve.storageDimension(), 10);

    capd::IVector x = curve.get_x();
    BOOST_CHECK_EQUAL(x.dimension(), 10);

    IMatrix C = curve.get_C();
    BOOST_CHECK_EQUAL(C.numberOfRows(), 10);
    BOOST_CHECK_EQUAL(C.numberOfColumns(), 3);

    capd::IVector r0 = curve.get_r0();
    BOOST_CHECK_EQUAL(r0.dimension(), 3);

    IMatrix B = curve.get_B();
    BOOST_CHECK_EQUAL(B.numberOfRows(), 10);
    BOOST_CHECK_EQUAL(B.numberOfColumns(), 10);

    capd::IVector r = curve.get_r();
    BOOST_CHECK_EQUAL(r.dimension(), 10);

    IMatrix Binv = curve.get_Binv();
    BOOST_CHECK_EQUAL(Binv.numberOfRows(), 10);

    // Test set_x
    capd::IVector new_x(10);
    for(int i=0; i<10; ++i) new_x[i] = i;
    curve.set_x(new_x);
    BOOST_CHECK_EQUAL(curve.get_x()[5].leftBound() <= 5 && curve.get_x()[5].rightBound() >= 5, true);

    // Test set_C
    IMatrix new_C(10, 3);
    curve.set_C(new_C);

    // Test set_r0
    capd::IVector new_r0(3);
    curve.set_r0(new_r0);

    // Test set_Cr0
    curve.set_Cr0(new_C, new_r0);

    // Test affine transform
    IMatrix M_aff(10, 10);
    capd::IVector v_aff(10);
    for(int i=0; i<10; ++i) M_aff[i][i] = 2.0;
    BOOST_CHECK_THROW(curve.affineTransform(M_aff, v_aff), std::logic_error);

    BOOST_CHECK_THROW(curve.translate(v_aff), std::logic_error);

    // Test take_* semantics
    capd::IVector* p_x = curve.take_x();
    BOOST_CHECK(p_x != nullptr);
    delete p_x;

    IMatrix* p_C = curve.take_C();
    BOOST_CHECK(p_C != nullptr);
    delete p_C;
}

BOOST_AUTO_TEST_CASE(Operations_Test) {
    GridType grid(Interval(0.1));
    TimePointType t0 = grid.point(0);
    TimePointType t1 = grid.point(2);

    capd::IVector vec1(2);
    vec1[0]=3.0; vec1[1]=4.0;

    DDESolutionCurve<DDESet> curve(t0, t1, 1, vec1, 3); // N0=3, length=2, order=1

    // Add vector
    capd::IVector v_add(10);
    for(int i=0; i<10; ++i) v_add[i] = 1.0;
    BOOST_CHECK_THROW(curve.add(v_add), std::logic_error);

    // Add set
    DDESolutionCurve<DDESet> curve2(t0, t1, 1, vec1, 3); // same structure
    BOOST_CHECK_THROW(curve.add(curve2), std::logic_error);

    // Mid point
    capd::IVector mid = curve.midPoint();
    BOOST_CHECK_EQUAL(mid.dimension(), 10);
}

BOOST_AUTO_TEST_CASE(IO_Show_Test) {
    GridType grid(Interval(0.1));
    TimePointType t0 = grid.point(0);
    TimePointType t1 = grid.point(2);

    capd::IVector vec1(2);
    vec1[0]=3.0; vec1[1]=4.0;

    DDESolutionCurve<DDESet> curve(t0, t1, 1, vec1, 3);

    std::string s = curve.show();
    BOOST_CHECK(!s.empty());

    std::stringstream ss;
    ss << curve;
    BOOST_CHECK(!ss.str().empty());
}

BOOST_AUTO_TEST_CASE(Hull_Test) {
    GridType grid(Interval(0.1));
    TimePointType t0 = grid.point(0);
    TimePointType t1 = grid.point(2);

    capd::IVector vec1(2);
    vec1[0]=3.0; vec1[1]=4.0;

    DDESolutionCurve<DDESet> curve1(t0, t1, 1, vec1, 3);

    capd::IVector vec2(2);
    vec2[0]=5.0; vec2[1]=6.0;
    DDESolutionCurve<DDESet> curve2(t0, t1, 1, vec2, 3);

    DDESolutionCurve<DDESet> out(t0, t1, 1, vec2, 3);

    capd::ddes::hull(curve1, curve2, out);
    BOOST_CHECK_EQUAL(out.length(), 2);
}

BOOST_AUTO_TEST_CASE(Subcurve_Test) {
    GridType grid(Interval(0.1));
    TimePointType t0 = grid.point(0);
    TimePointType t1 = grid.point(5);

    capd::IVector vec1(2);
    vec1[0]=3.0; vec1[1]=4.0;

    DDESolutionCurve<DDESet> curve(t0, t1, 1, vec1, 3);

    // Subcurve via size_type index
    DDESolutionCurve<DDESet> sub1 = curve.subcurve(2, 4);
    BOOST_CHECK_EQUAL(sub1.length(), 2);

    DDESolutionCurve<DDESet> sub2 = curve.subcurve(3);
    BOOST_CHECK_EQUAL(sub2.length(), 2);

    // Subcurve via TimePointType
    DDESolutionCurve<DDESet> sub3 = curve.subcurve(grid.point(2), grid.point(4));
    BOOST_CHECK_EQUAL(sub3.length(), 2);

    DDESolutionCurve<DDESet> sub4 = curve.subcurve(grid.point(3));
    BOOST_CHECK_EQUAL(sub4.length(), 2);

    // domain errors
    BOOST_CHECK_THROW(curve.subcurve(grid.point(-1), grid.point(4)), std::domain_error);
    BOOST_CHECK_THROW(curve.subcurve(grid.point(2), grid.point(10)), std::domain_error);
    BOOST_CHECK_THROW(curve.subcurve(grid.point(10)), std::domain_error);

    // Reverse iterators
    auto rit = curve.rbegin();
    BOOST_CHECK(rit != curve.rend());

    // at() via time point
    auto it = curve.at(grid.point(1));
    BOOST_CHECK_EQUAL((*it)->order(), 1);

    BOOST_CHECK_THROW(curve.at(grid.point(10)), std::range_error);
}

BOOST_AUTO_TEST_CASE(Operations_And_MakeStorage_2) {
    GridType grid(Interval(0.1));
    TimePointType t0 = grid.point(0);
    TimePointType t1 = grid.point(2);

    capd::IVector vec1(2);
    vec1[0]=3.0; vec1[1]=4.0;

    DDESolutionCurve<DDESet> curve(t0, t1, 1, vec1, 3);

    curve.makeStorage_x();
    curve.makeStorage_r0();
    curve.makeStorage_B();
    curve.makeStorage_r();
    curve.makeStorage_C();

    BOOST_CHECK(curve.get_x().dimension() == 10);

    // mulThenAdd throws logic_error
    DDESolutionCurve<DDESet> curve2(t0, t1, 1, vec1, 3);
    BOOST_CHECK_THROW(curve.mulThenAdd(2.0, curve2), std::logic_error);

    // test the reinitialize stub
    BOOST_CHECK_THROW(curve.reinitialize(2, 3), std::logic_error);

    capd::IVector h = curve.hull();
    BOOST_CHECK_EQUAL(h.dimension(), 10);
}

BOOST_AUTO_TEST_CASE(Pointer_Setters) {
    GridType grid(Interval(0.1));
    TimePointType t0 = grid.point(0);
    TimePointType t1 = grid.point(2);

    capd::IVector vec1(2);
    vec1[0]=3.0; vec1[1]=4.0;

    DDESolutionCurve<DDESet> curve(t0, t1, 1, vec1, 3);

    capd::IVector* px = new capd::IVector(10);
    curve.set_x(px, true);

    IMatrix* pC = new IMatrix(10, 3);
    curve.set_C(pC, true);

    capd::IVector* pr0 = new capd::IVector(3);
    BOOST_CHECK_THROW(curve.set_r0(pr0, true), std::logic_error);

    IMatrix* pC2 = new IMatrix(10, 3);
    capd::IVector* pr02 = new capd::IVector(3);
    BOOST_CHECK_THROW(curve.set_Cr0(pC2, pr02, true, true), std::logic_error);

    IMatrix* pB = new IMatrix(10, 10);
    curve.set_B(pB, true);

    capd::IVector* pr = new capd::IVector(10);
    curve.set_r(pr, true);

    IMatrix* pBinv = new IMatrix(10, 10);
    curve.set_Binv(pBinv, true);

    // Non owning variants
    capd::IVector v(10);
    BOOST_CHECK_THROW(curve.set_x(&v, false), std::domain_error);
}

BOOST_AUTO_TEST_CASE(Misc_Test) {
    GridType grid(Interval(0.1));
    TimePointType t0 = grid.point(0);
    TimePointType t1 = grid.point(2);

    capd::IVector vec1(2);
    vec1[0]=3.0; vec1[1]=4.0;

    DDESolutionCurve<DDESet> curve_high(t0, t1, 3, vec1, 3);
    DDESolutionCurve<DDESet> curve_low(t0, t1, 1, vec1, 3);

    // operator=(orig) throws logic_error internally when grid mismatch? Yes, tested.
    // get_Binv is tested.

    // hull() with diff length
    DDESolutionCurve<DDESet> curve(t0, t1, 1, vec1, 3);
    DDESolutionCurve<DDESet> curve_diff(t0, grid.point(10), 1, vec1, 3);

    DDESolutionCurve<DDESet> out(t0, t1, 1, vec1, 3);
    BOOST_CHECK_THROW(capd::ddes::hull(curve, curve_diff, out), std::logic_error);
}

BOOST_AUTO_TEST_CASE(Dot_Impl_Coverage) {
    GridType grid(Interval(0.1));
    TimePointType t0 = grid.point(0);
    TimePointType t1 = grid.point(2);

    capd::IVector vec1(2);
    vec1[0]=3.0; vec1[1]=4.0;

    DDESolutionCurve<DDESet> curve(t0, t1, 1, vec1, 3);

    capd::IVector dot_v(curve.dimension());
    for(int i=0; i<curve.dimension(); ++i) dot_v[i] = 1.0;

    // Test the standard vector dot product
    Interval res = curve.dot(dot_v);
    BOOST_CHECK_EQUAL(res.leftBound() <= 7.0 && res.rightBound() >= 7.0, true);
}

BOOST_AUTO_TEST_SUITE_END()
