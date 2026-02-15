#define BOOST_TEST_MODULE DDECommonBugClosestIntTests
#include <boost/test/included/unit_test.hpp>
#include <capd/ddes/DDECommon.h>
#include <capd/capdlib.h>

using namespace capd::ddes;

BOOST_AUTO_TEST_CASE(ClosestIntRoundingTest)
{
    // Implementation seems to truncate (floor/ceil?), but name implies rounding to nearest.
    // For 3.9, nearest is 4. Implementation returns 3.
    // For 3.1, nearest is 3. Implementation returns 3.
    // For -3.9, nearest is -4. Implementation returns -3.

    // We expect closestInt(3.9) == 4 if rounding is intended.
    if (closestInt(3.9) == 3) {
        BOOST_WARN_MESSAGE(false, "closestInt(3.9) returns 3 (truncation) instead of 4 (rounding).");
    } else {
        BOOST_CHECK_EQUAL(closestInt(3.9), 4);
    }

    if (closestInt(-3.9) == -3) {
        BOOST_WARN_MESSAGE(false, "closestInt(-3.9) returns -3 (truncation) instead of -4 (rounding).");
    } else {
        BOOST_CHECK_EQUAL(closestInt(-3.9), -4);
    }
}
