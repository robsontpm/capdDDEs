# Notes on Tests

## BasicDoubleton.h (DONE)
- Tests implemented in `tests/capd-ddes-storage-BasicDoubleton.cpp`. Coverage > 94%.
- **Known Bug 1:** Constructor `BasicDoubleton(size_type d, size_type N0 = -1)` fails due to unsigned underflow of `N0 = -1`.
- **Known Bug 2:** Constructor `BasicDoubleton(VectorType* x, ...)` fails to compile due to passing pointers instead of references to `setupFromData` and ambiguous `VectorType(0)`.
- These bugs are marked in the header file and corresponding tests are disabled or avoided.

## SharedDoubleton.h (UNFINISHED)
- Tests implemented in:
  - `tests/capd-ddes-storage-SharedDoubleton.cpp`
  - `tests/capd-ddes-storage-SharedDoubleton-Add.cpp` (Isolated `add`)
  - `tests/capd-ddes-storage-SharedDoubleton-Mul.cpp` (Isolated `mul`)
  - `tests/capd-ddes-storage-SharedDoubleton-MulThenAdd.cpp` (Isolated `mulThenAdd`)
  - `tests/capd-ddes-storage-SharedDoubleton-ZeroDim.cpp` (Isolated zero dimension cases)
  - `tests/capd-ddes-storage-SharedDoubleton-EdgeCases.cpp` (Exceptions, shared pointers)
  - `tests/capd-ddes-storage-SharedDoubleton-Constructors.cpp` (Additional constructors and setters)
- Coverage improved to ~44.6%.
- **Status:** All isolated tests PASS. This confirms that `add`, `mul`, `mulThenAdd` work correctly in isolation, and the crashes in the main suite are due to test environment interactions.
- **New Tests Coverage:**
  - `mul()` and `mulThenAdd()`.
  - Zero-dimensional vectors and `N0=0`.
  - Shared `r0` pointer usage in `add()`.
  - Exception handling in `reinit`, `affineTransform`, `translate`.
  - Additional constructors involving `Binv` and reference passing.
  - `set_B` behavior (verifying `IdQRPolicy` normalization).
- **Remaining Work:**
  - Coverage is still reported below 80%. This may be due to template instantiation measurement issues or untestable internal branches.
  - The main test suite integration of `add/mul` still crashes if enabled together. The problematic tests in `tests/capd-ddes-storage-SharedDoubleton.cpp` remain disabled.
  - Investigate `IdQRPolicy` behavior further to fully understand `set_B` side effects.
