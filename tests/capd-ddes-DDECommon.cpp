#define BOOST_TEST_MODULE DDECommonTestSuite
#include <boost/test/included/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/ddes/DDECommon.h"
#include <sstream>
#include <vector>
#include <stdexcept>

BOOST_AUTO_TEST_SUITE(DDECommonTestSuite)

// ==========================================
// Helper Functions Tests
// ==========================================

BOOST_AUTO_TEST_CASE(HelperSafeDeleteTest) {
    int* ptr = new int(5);
    capd::ddes::helper_safe_delete(ptr, true);
    BOOST_CHECK(ptr == nullptr);

    int* ptr2 = new int(10);
    capd::ddes::helper_safe_delete(ptr2, false); // Should not delete
    BOOST_CHECK(ptr2 != nullptr);
    BOOST_CHECK_EQUAL(*ptr2, 10);
    delete ptr2;

    int* nullPtr = nullptr;
    capd::ddes::helper_safe_delete(nullPtr, true); // Should be safe
    BOOST_CHECK(nullPtr == nullptr);
}

BOOST_AUTO_TEST_CASE(HelperSafeArrayDeleteTest) {
    int* arr = new int[5];
    capd::ddes::helper_safe_array_delete(arr, true);
    BOOST_CHECK(arr == nullptr);

    int* arr2 = new int[5];
    capd::ddes::helper_safe_array_delete(arr2, false); // Should not delete
    BOOST_CHECK(arr2 != nullptr);
    delete[] arr2;

    int* nullArr = nullptr;
    capd::ddes::helper_safe_array_delete(nullArr, true); // Should be safe
    BOOST_CHECK(nullArr == nullptr);
}

BOOST_AUTO_TEST_CASE(HelperDumpTest) {
    std::stringstream ss;
    ss << "line1\nline2";
    capd::ddes::helper_dump_line(ss);
    std::string remaining;
    std::getline(ss, remaining);
    BOOST_CHECK_EQUAL(remaining, "line2");

    std::stringstream ss2;
    ss2 << "word1 word2";
    capd::ddes::helper_dump_badge(ss2);
    std::string remaining2;
    ss2 >> remaining2;
    BOOST_CHECK_EQUAL(remaining2, "word2");
}

BOOST_AUTO_TEST_CASE(EcloseStepTest) {
    double h_d = 0.1;
    BOOST_CHECK_EQUAL(capd::ddes::ecloseStep(h_d), 0.1);

    capd::interval h_i(0.1);
    capd::interval res = capd::ddes::ecloseStep(h_i);
    BOOST_CHECK(res.leftBound() <= 0.0);
    BOOST_CHECK(res.rightBound() >= 0.1);
}

BOOST_AUTO_TEST_CASE(ShowEnclosedIntervalTest) {
    double a = 1.0, b = 2.0;
    std::string s = capd::ddes::showEnclosedInterval(a, b);
    BOOST_CHECK(s.find("[1, 2)") != std::string::npos); // Rough check

    capd::interval i1(1.0), i2(2.0);
    std::string s2 = capd::ddes::showEnclosedInterval(i1, i2);
    BOOST_CHECK(s2.find("[1") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(RethrowTest) {
    try {
        try {
            throw std::runtime_error("Original error");
        } catch (const std::exception& e) {
            // Must specify both types because deduction fails for mixed types
            throw capd::ddes::rethrow<std::exception, std::runtime_error>("Wrapper error", e);
        }
    } catch (const std::runtime_error& e) {
        std::string msg = e.what();
        BOOST_CHECK(msg.find("Wrapper error") != std::string::npos);
        BOOST_CHECK(msg.find("Original error") != std::string::npos);
    }
}

BOOST_AUTO_TEST_CASE(ClosestIntTest) {
    BOOST_CHECK_EQUAL(capd::ddes::closestInt(1.1), 1);
    BOOST_CHECK_EQUAL(capd::ddes::closestInt(1.9), 1); // Cast to int truncates? No, wait.
    // implementation is return int(value); so 1.9 -> 1.
    
    // capd::interval i1(1.1);
    // BOOST_CHECK_EQUAL(capd::ddes::closestInt(i1), 1);
    // BUG: The above fails to compile because closestInt overload for Interval expects
    // capd::intervals::Interval but capd::interval is capd::filib::Interval.

    BOOST_CHECK_EQUAL(capd::ddes::closestSmallerInt(1.1), 1);
    BOOST_CHECK_EQUAL(capd::ddes::closestSmallerInt(-1.1), -2); // -1 - 1 = -2
    
    capd::interval i2(-1.1);
    // This should work because there is a specific overload for capd::interval declared (and likely defined in cpp)
    BOOST_CHECK_EQUAL(capd::ddes::closestSmallerInt(i2), -2);
}

// ==========================================
// DiscreteTimeGrid Tests
// ==========================================

BOOST_AUTO_TEST_CASE(DiscreteTimeGridTest) {
    typedef capd::ddes::DiscreteTimeGrid<double> Grid;
    double h = 0.1;
    Grid grid(h);
    BOOST_CHECK_EQUAL(grid.h(), h);

    Grid::TimePointType t0 = grid.point(0);
    BOOST_CHECK(t0.isZero());
    BOOST_CHECK_EQUAL(t0.toInt(), 0);

    Grid::TimePointType t1 = grid.point(1);
    BOOST_CHECK_EQUAL(t1.toInt(), 1);
    BOOST_CHECK_CLOSE((double)t1, 0.1, 1e-15);

    Grid::TimePointType t2 = t1 + 1;
    BOOST_CHECK_EQUAL(t2.toInt(), 2);

    Grid::TimePointType t3 = t1 + t2;
    BOOST_CHECK_EQUAL(t3.toInt(), 3);

    // Test different grids
    Grid grid2(0.2);
    Grid::TimePointType t_other = grid2.point(1);
    BOOST_CHECK_THROW(t1 + t_other, std::logic_error);

    // Test grid equality
    Grid grid3(h);
    // Warning in header says: two grids created with same physical constant might be !=
    // Implementation uses shared_ptr.
    BOOST_CHECK(grid != grid3);

    // Copy constructor should share the pointer
    Grid grid4(grid);
    BOOST_CHECK(grid == grid4);
    BOOST_CHECK(grid.point(1) + grid4.point(1) == grid.point(2));
}

BOOST_AUTO_TEST_CASE(DiscreteTimeGridSplitTest) {
    typedef capd::ddes::DiscreteTimeGrid<double> Grid;
    double h = 0.5;
    Grid grid(h);

    double t = 1.2;
    // TimePointType must be initialized with the grid because it holds a reference to it.
    // Default constructor creates a point on trivialGrid which cannot be assigned to.
    Grid::TimePointType ti = grid.point(0);
    double epsi;

    grid.split(t, ti, epsi);
    // 1.2 / 0.5 = 2.4 -> floor is 2.
    // ti should be 2.
    // epsi = 1.2 - (2 * 0.5) = 0.2

    BOOST_CHECK_EQUAL(ti.toInt(), 2);
    BOOST_CHECK_CLOSE(epsi, 0.2, 1e-12);
}

BOOST_AUTO_TEST_CASE(DiscreteTimePointOperatorsTest) {
    typedef capd::ddes::DiscreteTimeGrid<double> Grid;
    Grid grid(0.1);
    Grid::TimePointType t1 = grid(1);
    Grid::TimePointType t2 = grid(2);

    BOOST_CHECK(t1 < t2);
    BOOST_CHECK(t2 > t1);
    BOOST_CHECK(t1 <= t2);
    BOOST_CHECK(t1 != t2);
    BOOST_CHECK(t1 == grid(1));

    BOOST_CHECK_EQUAL((++t1).toInt(), 2);
    BOOST_CHECK_EQUAL(t1.toInt(), 2);

    BOOST_CHECK_EQUAL((t1--).toInt(), 2);
    BOOST_CHECK_EQUAL(t1.toInt(), 1);

    t1 += 2;
    BOOST_CHECK_EQUAL(t1.toInt(), 3);

    t1 -= 1;
    BOOST_CHECK_EQUAL(t1.toInt(), 2);
}

// ==========================================
// Taylor Sum Tests
// ==========================================

BOOST_AUTO_TEST_CASE(SumTaylorForwardTest) {
    std::vector<double> coeffs = {1.0, 2.0, 3.0}; // 1 + 2x + 3x^2
    double step = 2.0;
    double out = 0.0;
    
    capd::ddes::sumTaylorForward(coeffs.begin(), 2, step, out);
    // 1 + 2*2 + 3*4 = 1 + 4 + 12 = 17
    BOOST_CHECK_EQUAL(out, 17.0);
}

BOOST_AUTO_TEST_CASE(SumTaylorBackwardTest) {
    std::vector<double> coeffs = {1.0, 2.0, 3.0}; // 1 + 2x + 3x^2
    double step = 2.0;
    double out = 0.0;

    // Note: backward iterator usually goes from end to begin.
    // The implementation does --a inside the loop.
    // It starts at index n (highest power).
    // Loop j=0 to n.
    // out = *a + (step * out)
    
    // Iteration 0: a points to 3.0. out = 3 + 0 = 3. --a -> points to 2.0
    // Iteration 1: a points to 2.0. out = 2 + 2*3 = 8. --a -> points to 1.0
    // Iteration 2: a points to 1.0. out = 1 + 2*8 = 17. --a -> invalid (before begin)

    // So we need to pass iterator to the LAST element.
    std::vector<double>::iterator it = coeffs.end();
    --it; // Points to 3.0

    capd::ddes::sumTaylorBackward(it, 2, step, out);
    BOOST_CHECK_EQUAL(out, 17.0);
}

// ==========================================
// Matrix Block Tests
// ==========================================

BOOST_AUTO_TEST_CASE(ExtractDiagonalBlocksTest) {
    int d = 2;
    int rows = 4;
    int cols = 4;
    capd::DMatrix M(rows, cols);
    // Block 1 (top-left)
    M[0][0] = 1.0; M[0][1] = 2.0;
    M[1][0] = 3.0; M[1][1] = 4.0;

    // Block 2 (bottom-right)
    M[2][2] = 5.0; M[2][3] = 6.0;
    M[3][2] = 7.0; M[3][3] = 8.0;

    // Off-diagonal element
    M[0][3] = 9.0;

    int offCount = 0;
    std::vector<capd::DMatrix> blocks = capd::ddes::extractDiagonalBlocks(M, d, offCount);
    
    BOOST_CHECK_EQUAL(blocks.size(), 2);
    BOOST_CHECK_EQUAL(offCount, 1);

    BOOST_CHECK_EQUAL(blocks[0][0][0], 1.0);
    BOOST_CHECK_EQUAL(blocks[0][1][1], 4.0);

    BOOST_CHECK_EQUAL(blocks[1][0][0], 5.0); // Local indices
    BOOST_CHECK_EQUAL(blocks[1][1][1], 8.0);
    
    // Test invalid dimensions
    capd::DMatrix M_bad(3, 3);
    int offCountBad = 0;
    BOOST_CHECK_THROW(capd::ddes::extractDiagonalBlocks(M_bad, 2, offCountBad), std::range_error);
}

BOOST_AUTO_TEST_SUITE_END()
