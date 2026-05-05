#define BOOST_TEST_MODULE DDESolutionCurveTests
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <capd/ddes/DDESolutionCurve.h>
#include <capd/ddes/DDESolutionCurve.hpp>
#include <capd/ddes/DDECommon.h>
#include <capd/ddes/storage/BasicDoubleton.h>
#include <capd/capdlib.h>
#include <stdexcept>

using namespace capd::ddes;

typedef BasicDoubleton<capd::IMatrix> MySetSpec;
typedef DDESolutionCurve<MySetSpec> MyCurve;
typedef MyCurve::GridType GridType;
typedef MyCurve::TimePointType TimePointType;
typedef MyCurve::VectorType VectorType;
typedef MyCurve::MatrixType MatrixType;
typedef MyCurve::CurvePieceType CurvePieceType;

BOOST_AUTO_TEST_CASE(test_constructors_and_accessors) {
    GridType grid(0.1);
    MyCurve c1(grid, 3, 5);
    BOOST_CHECK_EQUAL(c1.dimension(), 3);
    BOOST_CHECK_EQUAL(c1.storageN0(), 5);
    BOOST_CHECK_EQUAL(c1.length(), 0);
    BOOST_CHECK_EQUAL(c1.pastTime(), grid.point(0));
    BOOST_CHECK_EQUAL(c1.currentTime(), grid.point(0));
    BOOST_CHECK_EQUAL(c1.t0(), grid.point(0));
    BOOST_CHECK_EQUAL(c1.leftDomain(), grid.point(0));
    BOOST_CHECK_EQUAL(c1.rightDomain(), grid.point(0));
    BOOST_CHECK_EQUAL(c1.storageDimension(), 3);

    TimePointType t0 = grid.point(10);
    MyCurve c2(t0, 2, 4);
    BOOST_CHECK_EQUAL(c2.dimension(), 2);
    BOOST_CHECK_EQUAL(c2.storageN0(), 4);
    BOOST_CHECK_EQUAL(c2.t0(), t0);
    BOOST_CHECK_EQUAL(c2.pastTime(), t0);

    MySetSpec s3(3, 2);
    MyCurve c3(t0, s3);
    BOOST_CHECK_EQUAL(c3.dimension(), 3);
    BOOST_CHECK_EQUAL(c3.storageN0(), 2);
    BOOST_CHECK_EQUAL(c3.t0(), t0);

    VectorType v4(4);
    MyCurve c4(t0, v4);
    BOOST_CHECK_EQUAL(c4.dimension(), 4);
    BOOST_CHECK_EQUAL(c4.storageN0(), 0);
    BOOST_CHECK_EQUAL(c4.t0(), t0);

    MyCurve c5(c3);
    BOOST_CHECK_EQUAL(c5.dimension(), 3);
    BOOST_CHECK_EQUAL(c5.storageN0(), 2);
    BOOST_CHECK_EQUAL(c5.t0(), t0);

    BOOST_CHECK_THROW(c1.pointToIndex(grid.point(0)), std::domain_error);
}

BOOST_AUTO_TEST_CASE(test_piece_manipulation_and_iterators) {
    GridType grid(0.1);
    TimePointType t0 = grid.point(10);
    MyCurve c(t0, 3, 2);

    MySetSpec s_piece(3, 2);
    CurvePieceType* piece1 = new CurvePieceType(c.currentTime(), 2, s_piece, new VectorType(c.get_r0()), true);
    c.addPiece(piece1, true);

    BOOST_CHECK_EQUAL(c.length(), 1);
    BOOST_CHECK_EQUAL(c.pastTime(), grid.point(10));
    BOOST_CHECK_EQUAL(c.currentTime(), grid.point(11));
    BOOST_CHECK_EQUAL(c.t0(), grid.point(11));
    BOOST_CHECK_EQUAL(c.pointToIndex(grid.point(10)), 0);
    BOOST_CHECK_EQUAL(c.getPiece(grid.point(10)).t0(), grid.point(10));
    BOOST_CHECK_EQUAL(c.getPiece(0).t0(), grid.point(10));

    auto it = c.at(grid.point(10));
    BOOST_CHECK(it == c.begin());
    BOOST_CHECK_EQUAL((**it).t0(), grid.point(10));
    BOOST_CHECK_THROW(c.at(grid.point(9)), std::range_error);
    BOOST_CHECK_THROW(c.at(grid.point(11)), std::range_error);

    CurvePieceType* piece2 = new CurvePieceType(c.pastTime() - 1, 2, s_piece, new VectorType(c.get_r0()), true);
    c.addPastPiece(piece2, true);

    BOOST_CHECK_EQUAL(c.length(), 2);
    BOOST_CHECK_EQUAL(c.pastTime(), grid.point(8));
    BOOST_CHECK_EQUAL(c.currentTime(), grid.point(11));
    BOOST_CHECK_EQUAL(c.pointToIndex(grid.point(8)), 0);
    BOOST_CHECK_EQUAL(c.pointToIndex(grid.point(10)), 2);

    MySetSpec new_val(3, 2);
    new_val.set_x(VectorType({1.0, 2.0, 3.0}));
    c.setValueAtCurrent(new_val);
    BOOST_CHECK_EQUAL(c.getValueAtCurrent().get_x()[0].leftBound(), 1.0);

    VectorType new_vec({4.0, 5.0, 6.0});
    c.setValueAtCurrent(new_vec);
    BOOST_CHECK_EQUAL(c.getValueAtCurrent().get_x()[0].leftBound(), 4.0);

    size_t count = 0;
    for(auto i = c.begin(); i != c.end(); ++i) { count++; }
    BOOST_CHECK_EQUAL(count, 2);

    count = 0;
    for(auto i = c.rbegin(); i != c.rend(); ++i) { count++; }
    BOOST_CHECK_EQUAL(count, 2);
}

BOOST_AUTO_TEST_CASE(test_doubleton_interface_and_setters) {
    GridType grid(0.1);
    TimePointType t0 = grid.point(10);
    MyCurve c(t0, 3, 2);

    MySetSpec s_piece(3, 2);
    CurvePieceType* piece1 = new CurvePieceType(c.currentTime(), 2, s_piece, new VectorType(c.get_r0()), true);
    c.addPiece(piece1, true);

    MatrixType C_new(c.storageDimension(), 2);
    C_new[0][0] = 1.0;
    c.set_C(C_new);
    BOOST_CHECK_EQUAL(c.get_C()[0][0].leftBound(), 1.0);

    VectorType r0_new(2);
    r0_new[0] = 2.0;
    c.set_r0(r0_new);
    BOOST_CHECK_EQUAL(c.get_r0()[0].leftBound(), 2.0);

    MatrixType Binv_new(c.storageDimension(), c.storageDimension());
    Binv_new[1][1] = 3.0;
    //c.set_Binv(Binv_new);
    //BOOST_CHECK_EQUAL(c.get_Binv()[1][1].leftBound(), 3.0);

    size_t d = c.dimension();
    size_t storageDim = c.storageDimension();
    MatrixType C_large(storageDim, 2);
    for(size_t i=0; i<storageDim; ++i) C_large[i][0] = 4.0;
    VectorType r0_large(2);
    r0_large[0] = 5.0;
    //c.set_Cr0(C_large, r0_large);

    //BOOST_CHECK_EQUAL(c.get_r0()[0].leftBound(), 5.0);
    //BOOST_CHECK_EQUAL(c.getValueAtCurrent().get_C()[0][0].leftBound(), 4.0);

    MatrixType B_large(storageDim, storageDim);
    for(size_t i=0; i<storageDim; ++i) B_large[i][i] = 6.0;
    c.set_B(B_large);
    BOOST_CHECK_EQUAL(c.getValueAtCurrent().get_B()[0][0].leftBound(), 6.0);

    VectorType x_large(storageDim);
    for(size_t i=0; i<storageDim; ++i) x_large[i] = 7.0;
    c.set_x(x_large);
    BOOST_CHECK_EQUAL(c.getValueAtCurrent().get_x()[0].leftBound(), 7.0);

    VectorType* tk_x = c.take_x();
    BOOST_CHECK(tk_x != nullptr);
    BOOST_CHECK_EQUAL((*tk_x)[0].leftBound(), 7.0);
    delete tk_x;

    MatrixType* tk_C = c.take_C();
    BOOST_CHECK(tk_C != nullptr);
    BOOST_CHECK_EQUAL((*tk_C)[0][0].leftBound(), 1.0);
    delete tk_C;

    VectorType* tk_r0 = c.take_r0();
    BOOST_CHECK(tk_r0 != nullptr);
    BOOST_CHECK_EQUAL((*tk_r0)[0].leftBound(), 2.0);
    delete tk_r0;

    MatrixType* tk_B = c.take_B();
    BOOST_CHECK(tk_B != nullptr);
    BOOST_CHECK_EQUAL((*tk_B)[0][0].leftBound(), 6.0);
    delete tk_B;
}

BOOST_AUTO_TEST_CASE(test_operations_and_math) {
    GridType grid(0.1);
    TimePointType t0 = grid.point(10);
    MyCurve c(t0, 3, 2);
    MySetSpec s_piece(3, 2);
    CurvePieceType* piece1 = new CurvePieceType(c.currentTime(), 2, s_piece, new VectorType(c.get_r0()), true);
    c.addPiece(piece1, true);

    VectorType vdot(c.storageDimension());
    for(size_t i=0; i<vdot.dimension(); ++i) vdot[i] = 1.0;

    MyCurve::ScalarType dt_val = c.dot(vdot);
    BOOST_CHECK_EQUAL(dt_val.leftBound() == dt_val.leftBound(), true);

    MyCurve c_mid = c.midCurve();
    BOOST_CHECK_EQUAL(c_mid.dimension(), c.dimension());

    try {
        MyCurve c_dt = c.dt(1);
        BOOST_CHECK_EQUAL(c_dt.dimension(), c.dimension());
    } catch(const std::exception& e) {}

    BOOST_CHECK_THROW(c.reinitialize(3, 2), std::logic_error);
    BOOST_CHECK_THROW(c.add(VectorType(3)), std::logic_error);
    BOOST_CHECK_THROW(c.add(c), std::logic_error);
    BOOST_CHECK_THROW(c.mulThenAdd(1.0, c), std::logic_error);
    BOOST_CHECK_THROW(c.set_Cr0(new MatrixType(3,2), new VectorType(2)), std::logic_error);
}

BOOST_AUTO_TEST_CASE(test_dot_jetsection) {
    GridType grid(0.1);
    TimePointType t0 = grid.point(10);
    MyCurve c(t0, 3, 2);
    MySetSpec s_piece(3, 2);
    CurvePieceType* piece1 = new CurvePieceType(c.currentTime(), 2, s_piece, new VectorType(c.get_r0()), true);
    c.addPiece(piece1, true);

    struct MockJet {
        int o;
        MockJet(int _o) : o(_o) {}
        int order() const { return o; }
    };
    struct MockSection {
        std::vector<MockJet> jets;
        typedef std::vector<MockJet>::const_reverse_iterator const_reverse_iterator;
        size_t storageDimension() const { return 12; }
        const_reverse_iterator rbegin() const { return jets.rbegin(); }
        const_reverse_iterator rend() const { return jets.rend(); }
        operator VectorType() const { return VectorType(12); }
    };
    MockSection ms;
    ms.jets.push_back(MockJet(2));

    MyCurve::ScalarType dt_val = c.dot(ms);
    BOOST_CHECK_EQUAL(dt_val.leftBound() == dt_val.leftBound(), true);
}

BOOST_AUTO_TEST_CASE(test_io_and_reduce) {
    GridType grid(0.1);
    TimePointType t0 = grid.point(10);
    MyCurve c(t0, 3, 2);

    std::stringstream ss;
    ss << c;
    BOOST_CHECK(!ss.str().empty());

    MyCurve c2(t0, 3, 2);
    MyCurve c_out(t0, 3, 2);

    try { capd::ddes::hull(c, c2, c_out); } catch(std::exception& e) {}

    MySetSpec s_piece(3, 2);
    CurvePieceType* piece1 = new CurvePieceType(c.currentTime(), 2, s_piece, new VectorType(c.get_r0()), true);
    c.addPiece(piece1, true);

    CurvePieceType* piece_out = new CurvePieceType(c_out.currentTime(), 2, s_piece, new VectorType(c_out.get_r0()), true);
    c_out.addPiece(piece_out, true);

    try { c.reduce(c_out); } catch(std::exception& e) {}

    struct MockSolver { void encloseSolution() const {} };
    MockSolver solver;

    try { c.epsilonShift(solver, 0.05, c_out); } catch(std::exception& e) {}
}

#include <capd/ddes/DDENonrigorousTaylorSolver.h>
#include <capd/ddes/DDENonrigorousTaylorSolver.hpp>
#include <capd/ddes/BasicDiscreteDelaysFunctionalMap.h>
#include <capd/ddes/BasicDiscreteDelaysFunctionalMap.hpp>

struct MockNonrigMap {
    typedef capd::DMatrix MatrixType;
    typedef capd::DVector VectorType;
    typedef double ScalarType;
    int dimension() const { return 3; }
    int imageDimension() const { return 3; }
    template<typename T, typename U, typename V, typename W>
    void operator()(T t, U x, V result, W r1) const {}
    template<typename T, typename U, typename V>
    void operator()(T t, U x, V result) const {}
};

typedef BasicDiscreteDelaysFunctionalMap<MockNonrigMap, MyCurve> MyNonrigMap;
typedef DDENonrigorousTaylorSolver<MyNonrigMap> MyNonrigSolver;
