# Notes on Tests

## BasicDoubleton.h (DONE)
- Tests implemented in `tests/capd-ddes-storage-BasicDoubleton.cpp`. Coverage > 94%.
- **Known Bug 1:** Constructor `BasicDoubleton(size_type d, size_type N0 = -1)` fails due to unsigned underflow of `N0 = -1`.
- **Known Bug 2:** Constructor `BasicDoubleton(VectorType* x, ...)` fails to compile due to passing pointers instead of references to `setupFromData` and ambiguous `VectorType(0)`.
- These bugs are marked in the header file and corresponding tests are disabled or avoided.

## SharedDoubleton.h (DONE)
- **Status:** Tests Consolidated and Passing. Coverage Confirmed High.
- **Tests Implemented:** `tests/capd-ddes-storage-SharedDoubleton.cpp`
  - Consolidated all previous isolated tests (`Add`, `Mul`, `EdgeCases`, `Coverage`, `Constructors`).
  - All tests PASS.
  - **Fixed Crash:** `SanityCheckTest` was causing a crash because it reused a `SharedDoubleton` object after it was corrupted by exceptions (e.g. `set_x(nullptr)`). Fixed by using separate scopes for each check.
  - **Known Bug 3:** Constructor `SharedDoubleton(size_type d, size_type N0 = -1)` causes `std::bad_alloc`. This is reproduced in `tests/capd-ddes-storage-SharedDoubleton-BUG-Constructor.cpp` (catches exception).
  - **Confirmed Bug 4:** `SharedDoubleton` constructor with explicit data (e.g. `Doubleton(x, C, r0...)`) does not verify dimensions of `C` against `x` and `r0` before assignment. If dimensions mismatch, `IMatrix` assignment crashes (segfault). This is reproduced in `tests/capd-ddes-storage-SharedDoubleton-BUG-DataConstructor.cpp` (disabled by default to avoid crash).
  - **Coverage Tooling Note:**
  - **lcov 2.0+ Issue:** `lcov` 2.0+ reports artificially low coverage (~23%) due to issues mapping template header lines (`mismatch` errors).
  - **Resolution:** Downgrading to `lcov` 1.16 resolved the issue, reporting **98.0%** coverage, consistent with the user's `lcov` 1.14 report (98.2%).

## GenericJet.h (DONE)
- **Status:** Tests Implemented and Passing. Coverage logically high (report low due to lcov 2.0+).
- **Tests Implemented:** `tests/capd-ddes-storage-GenericJet.cpp`
  - Covers Constructors, Accessors, Modifiers, Evaluation, Derivative, Iterators.
- **Found Bug 5:** Self-assignment `jet = jet` fails because `operator=` deallocates coefficients before copying them (classic self-assignment issue).
  - Reproduced in `tests/capd-ddes-storage-GenericJet-BUG-SelfAssignment.cpp`.
  - The main test file `tests/capd-ddes-storage-GenericJet.cpp` avoids self-assignment to pass.

## DDECommon.h (DONE)
- **Status:** Tests Implemented and Passing. Coverage 99.3% (header), 100% (source).
- **Tests Implemented:** `tests/capd-ddes-DDECommon.cpp`
  - Covers `helper_dump_*`, `ecloseStep`, `showEnclosedInterval`, `rethrow`, `closestInt`, `closestSmallerInt`, `DiscreteTimeGrid`, `extractDiagonalBlocks`.
- **Known Issue:** `closestInt` (and related) implements truncation (casting to int) instead of rounding to nearest integer, despite the name suggesting otherwise.
  - Verified in `tests/capd-ddes-DDECommon-BUG-ClosestInt.cpp` which emits a warning if truncation is observed instead of rounding.
  - The implementation uses `int(value)` which truncates towards zero.
- **Coverage Tooling:** `lcov` 1.16 used.

## BasicDiscreteDelaysFunctionalMap.h (DONE)
- **Status:** Tests Implemented and Passing. Coverage 100% (header), 89% (hpp).
- **Tests Implemented:** `tests/capd-ddes-BasicDiscreteDelaysFunctionalMap.cpp`
  - Covers Constructors, Accessors, `operator()`, `collectComputationData`, `computeDDECoefficients` (with and without Jacobian).
  - Mocks `MapType`, `SolutionCurveSpec`, and `JetSpec` to isolate tests.
- **Found Bug 6:** `getMaxDelay` in `BasicDiscreteDelaysFunctionalMap.hpp` fails to compile because it compares an iterator (`tau`) with a value (`max_tau`). It should dereference the iterator (`*tau`).
  - The bug is marked in the header file.
  - The failing test case is extracted to `tests/capd-ddes-BasicDiscreteDelaysFunctionalMap-BUG-GetMaxDelay.cpp` (disabled by default to allow compilation).
- **Boost Test Framework:** Tests now use the precompiled `Boost::unit_test_framework` library (installed in `${HOME}/deps/boost`) instead of the header-only variant, to improve compilation speed and match project conventions. `BOOST_TEST_DYN_LINK` is defined in the new test files.

## DiscreteDelaysFunctionalMap.h (DONE)
- **Status:** Tests Implemented and Passing. Coverage 93.1%.
- **Tests Implemented:** `tests/capd-ddes-DiscreteDelaysFunctionalMap.cpp`
  - Covers Constructors, `operator()`, `collectComputationData`, `computeDDECoefficients`, `findRoughEnclosure`.
  - Uses `capd::Interval` and `capd::IVector` to verify rigorous logic.
- **Found Issue:** Template `checkDimension(AnyVector const&)` shadows `checkDimension(size_type)` if `AnyVector` type matches closer (or exactly) while `size_type` (usually `std::size_t`) requires conversion from `Vector::dimension()` return type (often `unsigned int` or `int`).
  - Workaround: In `MockSolutionCurve`, `size_type` was explicitly set to `unsigned int` to match `IVector::dimension()` and avoid ambiguity/shadowing.

## DDEForwardTaylorCurvePiece.h (DONE)
- **Status:** Tests Implemented and Passing. Coverage 89.8% (header), 92.0% (hpp).
- **Tests Implemented:** `tests/capd-ddes-DDEForwardTaylorCurvePiece.cpp`
  - Covers Constructors (Default, TimePoint, Copy, DimOrder, Vector, Set, iterators), Assignment, Evaluation methods (`taylor`, `summa`, `eval`, `evalCoeffs` and their Delta variants), operations (`mul`, `midCurve`), accessors/setters, iterators, exceptions, and `GenericJet` integration dot product.
- **Found Bug 7:** `AssignmentOperator` does not copy the base time `m_t0` from the source object.
  - Reproduced and verified in `tests/capd-ddes-DDEForwardTaylorCurvePiece-BUG-Assignment.cpp` (disabled by default, emits warning).
- **Found Issues:**
  - `dt(n)` derivative function is unimplemented and correctly throws `std::logic_error("Not implemented yet")`.
  - `reinitialize(d, N0)` is unimplemented and correctly throws `std::logic_error("Not Supported Yet")`.
  - Both these cases are isolated in `tests/capd-ddes-DDEForwardTaylorCurvePiece-BUG-NotImplemented.cpp`.
- **Note on `jetAt` Template Ambiguity:** In previous review it was identified that `jetAt(TimePointType)` and `jetAt(RealType)` cause compilation errors if `TimePointType` and `RealType` match. The user clarified this is intentional by design, as these types represent different domains in the application. Tests implement `MockTimePoint` to accurately reflect this.
